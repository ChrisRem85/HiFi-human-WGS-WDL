# bash_pipeline: HiFi human WGS variant calling without a WDL engine

This is a bash re-implementation of the WDL workflows in
[workflows/](../workflows). It runs the exact same tools, via the exact same
pinned container images and (by default) the exact same parameters as the
WDL tasks, in a single sequential script per sample/family — no `miniwdl`,
`Cromwell`, or `sprocket` required. Docker is still used to run each
bioinformatics tool, since that's what makes the WDL tasks reproducible in
the first place; only the *workflow engine* is removed.

## Requirements

- `bash`, `docker` (with permission to pull `quay.io/pacbio/*`,
  `google/deepvariant`, `nvcr.io/nvidia/clara/clara-parabricks`)
- `python3` is **not** required on the host — all Python (plotting/stat
  extraction) runs inside the same containers the WDL tasks use.
- Optional: an NVIDIA GPU + `nvidia-container-toolkit` if you set
  `USE_GPU=true`.

## How it works

1. **`collect_manifest.sh`** — collects and validates every file the run
   needs *before* any analysis starts:
   - unpacks the pinned reference-data container (same one
     `unpack_container_manifest.wdl` uses) to resolve the reference FASTA,
     TRGT catalog, Sawfish exclude/expected-CN BEDs, MethBat region TSV,
     Paraphase genome build, etc.
   - validates every sample's HiFi reads BAM(s) exist.
   - writes `manifest.tsv` (every file used, with its sha256 + size) and
     `ref.env` (resolved reference paths, sourced by later steps).
2. **`run_pipeline.sh`** — runs every step in sequence, per sample:
   `pbmm2 align` → `mosdepth` → `DeepVariant` → `sawfish discover` →
   `paraphase` → `mitorsaw` → `kivvi` (KIV2, D4Z4), then, for `--mode family`,
   a joint-calling pass (`sawfish joint-call` + `GLnexus`, split back out per
   sample), then per-sample downstream: `hiphase` → `trgt genotype` →
   `pbjam bam-stats` → `bcftools stats/roh` → SV stats → `methbat`
   pileup/profile → `pbstarphase` (PGx).

Every tool call, image digest, and default parameter lives in **`config.sh`**
and **`steps.sh`** — a single source of truth used identically for every
sample in a run, which is what guarantees identical output files and
identical parameters across samples.

## Usage

```bash
# 1. Everything referenced below (BAMs, and the --out directory) must live
#    under a common --data-root, since it's bind-mounted into every
#    container 1:1 (so absolute host paths resolve unchanged inside them).
DATA_ROOT=/data

# 2. Collect + validate the manifest (downloads/unpacks reference assets).
bash_pipeline/collect_manifest.sh \
  --samples bash_pipeline/samples.example.tsv \
  --data-root "$DATA_ROOT" \
  --out "$DATA_ROOT/manifest"

# 3. Run the pipeline.
bash_pipeline/run_pipeline.sh \
  --manifest-dir "$DATA_ROOT/manifest" \
  --data-root "$DATA_ROOT" \
  --out "$DATA_ROOT/results" \
  --mode family   # or: --mode singleton
```

Outputs land in `$DATA_ROOT/results/<sample_id>/{upstream,downstream}/...`
and, for family runs, `$DATA_ROOT/results/<family_id>.joint/...` — the same
filenames (`<sample_id>.<ref_name>.<suffix>`) that the WDL tasks produce.

### Overriding parameters

Every knob in `config.sh` can be overridden via environment variable, e.g.:

```bash
# DeepVariant on GPU, without Parabricks (uses docker's --runtime=nvidia)
THREADS=64 USE_GPU=true \
  bash_pipeline/run_pipeline.sh --manifest-dir ... --data-root ... --out ... --mode family

# DeepVariant via Parabricks on GPU instead
THREADS=64 USE_GPU=true USE_PARABRICKS_DEEPVARIANT=true \
  bash_pipeline/run_pipeline.sh --manifest-dir ... --data-root ... --out ... --mode family
```

`USE_GPU=true` alone runs DeepVariant's own `-gpu` image with
`DEEPVARIANT_GPU_DOCKER_ARGS` (default: `--runtime=nvidia -e NVIDIA_VISIBLE_DEVICES=0`).
Adding `USE_PARABRICKS_DEEPVARIANT=true` instead runs Parabricks with
`PARABRICKS_GPU_DOCKER_ARGS` (default: `--gpus all`). Override either
variable if your host's NVIDIA container setup differs (e.g. to target a
different GPU index or use Docker's native `--gpus` flag for both).

Because `config.sh` is sourced once per script invocation and every sample in
a run shares the same environment, changing a variable changes it for every
sample in that run — there's no per-sample parameter drift.

## Known simplifications vs. the WDL workflows

- **DeepVariant** is invoked via the standard single-node `run_deepvariant`
  wrapper (`--num_shards=$THREADS`) instead of WDL's manual 64-way
  make_examples/call_variants/postprocess_variants sharding, which exists
  purely to parallelize across many small cloud VMs. Same model
  (`--model_type=PACBIO`), same output files.
- Per-task CPU/memory tuned for specific cloud VM shapes (e.g. sawfish's
  128 GiB) are not reproduced; use the global `THREADS` (and, where a tool
  takes an explicit memory flag, `MERGE_MEM_GB`/`GLNEXUS_MEM_GB`) knobs for
  your own hardware instead.
- `fail_reads` (optional low-quality/failed-read rescue for TRGT) is not
  implemented; it's an optional input in the WDL workflows too.
- Alignment chunking (`--chunk N/M`) exists in the WDL purely to scatter
  `pbmm2 align` across many separate cloud/HPC nodes in parallel. This
  pipeline runs everything on one node, and its chunk loop is sequential, so
  chunking would only add repeated startup/index-load overhead with no
  speedup — it therefore defaults to `USE_ALIGNMENT_CHUNKING=false` (single
  `pbmm2 align` pass per input BAM). Set it to `true` only if you also
  parallelize the chunk loop yourself.
