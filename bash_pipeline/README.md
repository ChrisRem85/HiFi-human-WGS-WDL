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

## Tools

Every tool below runs inside its own pinned Docker image (see `IMG_*` in
`config.sh`); this is the same list of tools the WDL tasks in
[workflows/](../workflows) invoke. See [docs/tools.md](../docs/tools.md) and
[docs/tools_containers.md](../docs/tools_containers.md) for full version and
container details.

- **pbindex** (`pbtk`) — indexes an unaligned HiFi reads BAM so `pbmm2` can
  scatter it into chunks; only used when `USE_ALIGNMENT_CHUNKING=true`.
- **pbmm2** — PacBio's `minimap2`-based aligner tuned for HiFi reads; indexes
  the reference and aligns each input BAM to it.
- **pbsamoa** — merges multiple per-movie aligned BAMs for a sample into a
  single sorted, indexed BAM (only run when a sample has >1 input BAM).
- **mosdepth** — computes per-region depth from the aligned BAM, used for
  mean coverage and chrY-based sex inference.
- **DeepVariant** — deep-learning small variant (SNV/indel) caller; run
  either on CPU, on GPU (`USE_GPU=true`), or via NVIDIA **Parabricks**
  (`USE_PARABRICKS_DEEPVARIANT=true`) for GPU-accelerated calling.
- **sawfish** — structural variant and copy-number discovery/calling;
  `sawfish discover` runs per sample, then `sawfish joint-call` merges one or
  more samples' discover output into a genotyped SV VCF (a single sample for
  `--mode singleton`, the whole family for `--mode family`).
- **sawshark** — post-processes `sawfish joint-call` output into the final
  genotyped SV VCF.
- **paraphase** — paralog-aware variant calling for genes with high sequence
  similarity (e.g. _SMN1_/_SMN2_, _PMS2_/_PMS2CL_); uses `minimap2`
  internally for realignment.
- **mitorsaw** — mitochondrial haplotype variant caller and haplotype
  statistics.
- **kivvi** — targeted genotyping for the KIV2 (_LPA_) and D4Z4 (FSHD)
  repeat regions.
- **GLnexus** — joint genotyper that merges per-sample GVCFs into a
  family-level small-variant callset (`--mode family` only).
- **HiPhase** — haplotype-phases small variants and SVs together and
  haplotags the aligned BAM.
- **TRGT** — targeted tandem repeat genotyper; also reports repeat
  "dropouts" (regions with insufficient haplotype coverage) via its bundled
  `find_trgt_dropouts.py`.
- **pbjam** — computes HiFi read/alignment QC statistics (read length,
  quality, mapped %, etc.) and length/quality distribution plots.
- **bcftools** — used throughout for VCF manipulation: small-variant
  stats/plots, runs-of-homozygosity (`bcftools roh`) detection, and SV
  summary stats (counts/length distributions by SV type).
- **samtools** — used throughout for BAM indexing/sorting (e.g. TRGT's
  spanning-reads BAM) and header inspection.
- **MethBat** — 5mC/5hmC methylation pileup (`methbat pileup`) from the
  haplotagged BAM, and region-level methylation profiling (`methbat
  profile`).
- **StarPhase** (`pbstarphase`) — HLA typing and pharmacogenomic (PGx)
  diplotype calling from the phased small-variant/SV VCFs and aligned BAM.

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

### Resuming an interrupted or already-finished run

By default, re-running `run_pipeline.sh` unconditionally re-runs and
overwrites every step's output. Pass `--resume` (or set `RESUME=true`) to
skip any step whose expected output file(s) already exist from a prior run:

```bash
bash_pipeline/run_pipeline.sh \
  --manifest-dir "$DATA_ROOT/manifest" \
  --data-root "$DATA_ROOT" \
  --out "$DATA_ROOT/results" \
  --mode family --resume
```

Each step checks its own final output(s) (e.g. the phased VCF + index for
`hiphase`, the VCF + index for `sawfish call`/`GLnexus`) and is skipped only
if every one of them already exists and is non-empty; otherwise it re-runs
and overwrites as usual. This is a per-step check, not a whole-run
checkpoint, so a run interrupted partway through a step will simply re-run
that one step in full on resume.

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
- `collect_manifest.sh` is not covered by `--resume`: it always re-unpacks
  the reference bundle and re-validates/re-hashes every sample's BAM(s). This
  is cheap relative to the analysis steps, so it isn't worth skipping.
