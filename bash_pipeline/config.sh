#!/usr/bin/env bash
# Single source of truth for every tool parameter, container image and
# reference asset used by the pipeline. Every step script sources this file,
# so every sample/family run with the same config.sh is guaranteed to use
# identical parameters and identical container image digests.
#
# Values below mirror the defaults baked into the upstream WDL tasks
# (workflows/wdl-common/wdl/tasks/*.wdl and workflows/wdl-common/wdl/workflows/*).
# Do not change an image digest without also updating the matching WDL task,
# and vice versa, or the bash and WDL pipelines will silently diverge.

set -euo pipefail

# ---------------------------------------------------------------------------
# Container registry / images (pinned by digest, same as the WDL `runtime` blocks)
# ---------------------------------------------------------------------------
: "${CONTAINER_REGISTRY:=quay.io/pacbio}"

IMG_PBMM2="${CONTAINER_REGISTRY}/pbmm2@sha256:86e39aa67fa5d385769d5f119739e8811ce550163e3b9dfc42bd58d1fecdf3a8"
IMG_PBTK="${CONTAINER_REGISTRY}/pbtk@sha256:f27bafa0ae6ffff6170d45a86a5677406320c7866760391f89ebe900a74d2039"
IMG_PBSAMOA="${CONTAINER_REGISTRY}/pbsamoa@sha256:cf70c89422d63c3a4e4b2ffd912160c3e7481303a211cd6134236fe6e9f9461f"
IMG_BASE="${CONTAINER_REGISTRY}/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae"
IMG_MOSDEPTH="${CONTAINER_REGISTRY}/mosdepth@sha256:f4edf52e6a31eb18f1755e8cb90224fec6aa88e4a4353a673b16a89f59cbe822"
IMG_SAWFISH="${CONTAINER_REGISTRY}/sawfish@sha256:5cf8f02790ecb89e652885c57f103fc9597c41f75ee71be2c05510cbbdf68f59"
IMG_PARAPHASE="${CONTAINER_REGISTRY}/paraphase@sha256:665ac4fefcef92e0023395eae9e0ca3f6e65d1228c49afba3ad3f8e7b3f3d1cb"
IMG_MITORSAW="${CONTAINER_REGISTRY}/mitorsaw@sha256:9b610f343ae018b21c86977259d91dba587e1f22d7a73e115937e428d3c7adeb"
IMG_KIVVI="${CONTAINER_REGISTRY}/kivvi@sha256:9e4a390821ea999af9b3faa483c5f7dce58041c439a35e4af3676f25edecd204"
IMG_HIPHASE="${CONTAINER_REGISTRY}/hiphase@sha256:41ebe22b55c66e2e78da2013f7fffaecc02a8b4e980400c3ea8d03c87330522e"
IMG_TRGT="${CONTAINER_REGISTRY}/trgt@sha256:648aee4a2c9d7371a48e454a7143861a242b853d81ff5453924cc0095d207824"
IMG_PBJAM="${CONTAINER_REGISTRY}/pbjam@sha256:1b90537ed6683ac7b89002c0f59e940b4ea46c600d292904f549af1ae2760bba"
IMG_METHBAT="${CONTAINER_REGISTRY}/methbat@sha256:281569947c0a6154f9a6bdaa1bb5f1e67dd160ccb49028d3218a49a37a0488fb"
IMG_PBSTARPHASE="${CONTAINER_REGISTRY}/pbstarphase@sha256:86aeeb3a9663c22c135d08353de700050e1d268406c1f63832cad5f139cdd390"
IMG_GLNEXUS="${CONTAINER_REGISTRY}/glnexus@sha256:ce6fecf59dddc6089a8100b31c29c1e6ed50a0cf123da9f2bc589ee4b0c69c8e"

IMG_DEEPVARIANT_CPU="google/deepvariant:1.10.0"
IMG_DEEPVARIANT_GPU="google/deepvariant:1.10.0-gpu"
IMG_PARABRICKS="nvcr.io/nvidia/clara/clara-parabricks:4.7.0-1"

# Reference data bundle container. Choose the variant matching REF_NAME.
: "${REF_NAME:=GRCh38_GIABv3}"
IMG_REF_DATA_GRCH38="${CONTAINER_REGISTRY}/workflow-data-container-hifi-human-wgs-wdl-grch38@sha256:5e3f23e44c1c09762838e81677f7d136367be2a4063efc608129841cef745b42"
IMG_REF_DATA_GRCH38_GIABV3="${CONTAINER_REGISTRY}/workflow-data-container-hifi-human-wgs-wdl-grch38_giabv3@sha256:f3799c1eff1b816a95ea06ce567d3324d06e945428567077d96eca07c76ad0aa"
case "${REF_NAME}" in
  GRCh38_GIABv3) IMG_REF_DATA="${IMG_REF_DATA_GRCH38_GIABV3}" ;;
  GRCh38)        IMG_REF_DATA="${IMG_REF_DATA_GRCH38}" ;;
  *) echo "Unknown REF_NAME='${REF_NAME}'. Set IMG_REF_DATA explicitly." >&2 ;;
esac

# ---------------------------------------------------------------------------
# Compute / resource parameters (uniform across every sample in the run)
# ---------------------------------------------------------------------------
: "${THREADS:=$(nproc 2>/dev/null || echo 8)}"
: "${USE_GPU:=false}"                    # true to use DeepVariant/Parabricks GPU path
: "${USE_PARABRICKS_DEEPVARIANT:=false}" # true to use Parabricks instead of DeepVariant
: "${USE_ALIGNMENT_CHUNKING:=false}"     # WDL scatters chunks across cloud nodes; this pipeline runs them sequentially on one node, so chunking only adds overhead here
: "${ALIGNMENT_CHUNKS:=16}"              # used only if USE_ALIGNMENT_CHUNKING=true (matches WDL's fixed 16-way chunking)
: "${GLNEXUS_MEM_GB:=60}"

# `docker run` GPU-exposure flags, split by GPU consumer since they expect
# different container GPU runtimes. Override either if your host uses a
# different NVIDIA setup (e.g. Docker's native `--gpus` support everywhere).
: "${DEEPVARIANT_GPU_DOCKER_ARGS:=--runtime=nvidia -e NVIDIA_VISIBLE_DEVICES=0}"
: "${PARABRICKS_GPU_DOCKER_ARGS:=--gpus all}"
: "${MERGE_MEM_GB:=16}"
: "${DEEPVARIANT_GVCF_OUTPUT:=true}"

# ---------------------------------------------------------------------------
# Tool parameters (verbatim defaults from the WDL task definitions)
# ---------------------------------------------------------------------------
# pbmm2 align
PBMM2_MIN_LENGTH=50
PBMM2_STRIP_KINETICS=true
PBMM2_KEEP_UNMAPPED=true

# mosdepth
MOSDEPTH_MAX_NORM_FEMALE_CHRY_DEPTH_DEFAULT=0.1   # overridden by reference bundle value if present

# trgt genotype
TRGT_MIN_READ_QUALITY=-1.0
TRGT_MAX_DEPTH=150
TRGT_HAPLOTYPE_COVERAGE_THRESHOLD=2

# bcftools_stats_roh_small_variants
BCFTOOLS_ROH_MIN_LENGTH=100000
BCFTOOLS_ROH_MIN_QUAL=20

# sv_stats
SV_STATS_MIN_LENGTH=50
SV_STATS_MAX_SCAR_LENGTH=10

# methbat pileup
METHBAT_MIN_MAPQ=1
METHBAT_MIN_COVERAGE=4
METHBAT_EDGE_TRIMMING_SIZE=20
METHBAT_PHASE_SET_MIN_FRACTION=0.7
METHBAT_SKIP_5MC=false
METHBAT_SKIP_5HMC=false
METHBAT_SKIP_6MA=true

# hiphase (all unset => tool defaults, matching WDL's optional inputs left undefined)
HIPHASE_PRESET=""
HIPHASE_MIN_GQ=""
HIPHASE_MIN_MAPQ=""
HIPHASE_DISABLE_GLOBAL_REALIGNMENT=false
HIPHASE_NO_SUPPLEMENTAL_JOINS=false
HIPHASE_PHASE_SINGLETONS=false

export CONTAINER_REGISTRY REF_NAME IMG_REF_DATA
export IMG_PBMM2 IMG_PBTK IMG_PBSAMOA IMG_BASE IMG_MOSDEPTH IMG_SAWFISH IMG_PARAPHASE \
  IMG_MITORSAW IMG_KIVVI IMG_HIPHASE IMG_TRGT IMG_PBJAM IMG_METHBAT IMG_PBSTARPHASE \
  IMG_GLNEXUS IMG_DEEPVARIANT_CPU IMG_DEEPVARIANT_GPU IMG_PARABRICKS
export THREADS USE_GPU USE_PARABRICKS_DEEPVARIANT USE_ALIGNMENT_CHUNKING ALIGNMENT_CHUNKS \
  GLNEXUS_MEM_GB MERGE_MEM_GB DEEPVARIANT_GVCF_OUTPUT \
  DEEPVARIANT_GPU_DOCKER_ARGS PARABRICKS_GPU_DOCKER_ARGS
export PBMM2_MIN_LENGTH PBMM2_STRIP_KINETICS PBMM2_KEEP_UNMAPPED \
  MOSDEPTH_MAX_NORM_FEMALE_CHRY_DEPTH_DEFAULT \
  TRGT_MIN_READ_QUALITY TRGT_MAX_DEPTH TRGT_HAPLOTYPE_COVERAGE_THRESHOLD \
  BCFTOOLS_ROH_MIN_LENGTH BCFTOOLS_ROH_MIN_QUAL \
  SV_STATS_MIN_LENGTH SV_STATS_MAX_SCAR_LENGTH \
  METHBAT_MIN_MAPQ METHBAT_MIN_COVERAGE METHBAT_EDGE_TRIMMING_SIZE METHBAT_PHASE_SET_MIN_FRACTION \
  METHBAT_SKIP_5MC METHBAT_SKIP_5HMC METHBAT_SKIP_6MA \
  HIPHASE_PRESET HIPHASE_MIN_GQ HIPHASE_MIN_MAPQ HIPHASE_DISABLE_GLOBAL_REALIGNMENT \
  HIPHASE_NO_SUPPLEMENTAL_JOINS HIPHASE_PHASE_SINGLETONS
