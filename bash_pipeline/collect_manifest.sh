#!/usr/bin/env bash
# Step 0: collect and validate a manifest of every file the pipeline needs
# before any analysis runs: the reference asset bundle (unpacked from its
# container, exactly as unpack_container_manifest.wdl does) plus every
# sample's HiFi reads BAM(s).
#
# Writes:
#   <out>/ref.env        - resolved reference file paths + parameters (sourced by later steps)
#   <out>/manifest.tsv    - every file used by the run, with checksum + size
#
# Usage:
#   collect_manifest.sh --samples samples.tsv --data-root /data --out /data/manifest
#
# samples.tsv columns (tab-separated, no header):
#   family_id  sample_id  sex(MALE|FEMALE|.)  hifi_reads(comma-separated absolute paths)
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=config.sh
source "${SCRIPT_DIR}/config.sh"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"

SAMPLES_TSV=""
OUT_DIR=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --samples) SAMPLES_TSV="$2"; shift 2 ;;
    --data-root) DATA_ROOT="$2"; shift 2 ;;
    --out) OUT_DIR="$2"; shift 2 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

[[ -n "${SAMPLES_TSV}" ]] || die "--samples <samples.tsv> is required"
[[ -n "${DATA_ROOT:-}" ]] || die "--data-root <dir> is required"
[[ -n "${OUT_DIR}" ]] || die "--out <dir> is required"
require_cmd docker
require_cmd sha256sum
require_file "${SAMPLES_TSV}"

DATA_ROOT="$(cd "${DATA_ROOT}" && pwd)"
mkdir -p "${OUT_DIR}"
OUT_DIR="$(cd "${OUT_DIR}" && pwd)"
assert_under_data_root "${OUT_DIR}"

MANIFEST="${OUT_DIR}/manifest.tsv"
echo -e "category\tsample_id\tpath\tsha256\tsize_bytes" > "${MANIFEST}"

manifest_add() {
  local category="$1" sample_id="$2" path="$3"
  require_file "${path}"
  assert_under_data_root "${path}"
  local sha size
  sha="$(sha256sum "${path}" | cut -d' ' -f1)"
  size="$(stat -c%s "${path}" 2>/dev/null || stat -f%z "${path}")"
  echo -e "${category}\t${sample_id}\t${path}\t${sha}\t${size}" >> "${MANIFEST}"
}

# ---------------------------------------------------------------------------
# 1. Unpack the reference asset bundle (same container image + script the
#    WDL's unpack_container_manifest task uses).
# ---------------------------------------------------------------------------
REF_DIR="${OUT_DIR}/reference"
mkdir -p "${REF_DIR}"
log "Unpacking reference bundle (${IMG_REF_DATA}) into ${REF_DIR}"
run_docker "${IMG_REF_DATA}" "${REF_DIR}" \
  "python3 /opt/scripts/unpack_container.py --manifest /opt/manifests/manifest.json --output-dir ."

REF_FASTA="$(find "${REF_DIR}/ref_fasta" -type f | head -n1)"
REF_INDEX="$(find "${REF_DIR}/ref_index" -type f | head -n1)"
TRGT_BED="$(find "${REF_DIR}/trgt_tandem_repeat_bed" -type f | head -n1)"
SAWFISH_EXCLUDE_BED="$(find "${REF_DIR}/sawfish_exclude_bed" -type f | head -n1)"
SAWFISH_EXCLUDE_BED_INDEX="$(find "${REF_DIR}/sawfish_exclude_bed_index" -type f | head -n1)"
SAWFISH_EXPECTED_BED_MALE="$(find "${REF_DIR}/sawfish_expected_bed_male" -type f | head -n1)"
SAWFISH_EXPECTED_BED_FEMALE="$(find "${REF_DIR}/sawfish_expected_bed_female" -type f | head -n1)"
METHBAT_REGION_TSV="$(find "${REF_DIR}/methbat_region_tsv" -type f | head -n1)"
REF_NAME_RESOLVED="$(cat "${REF_DIR}/ref_name")"
MAX_NORM_FEMALE_CHRY_DEPTH="$(cat "${REF_DIR}/max_norm_female_chrY_depth")"
PARAPHASE_GENOME_BUILD="$(cat "${REF_DIR}/paraphase_genome_build")"
RUN_STARPHASE="$(cat "${REF_DIR}/run_starphase")"

for f in "${REF_FASTA}" "${REF_INDEX}" "${TRGT_BED}" "${SAWFISH_EXCLUDE_BED}" \
  "${SAWFISH_EXCLUDE_BED_INDEX}" "${SAWFISH_EXPECTED_BED_MALE}" "${SAWFISH_EXPECTED_BED_FEMALE}" \
  "${METHBAT_REGION_TSV}"; do
  [[ -n "${f}" ]] || die "Reference unpack did not produce an expected asset; check ${REF_DIR}"
  manifest_add "reference" "-" "${f}"
done

cat > "${OUT_DIR}/ref.env" <<EOF
REF_NAME="${REF_NAME_RESOLVED}"
REF_FASTA="${REF_FASTA}"
REF_INDEX="${REF_INDEX}"
TRGT_BED="${TRGT_BED}"
SAWFISH_EXCLUDE_BED="${SAWFISH_EXCLUDE_BED}"
SAWFISH_EXCLUDE_BED_INDEX="${SAWFISH_EXCLUDE_BED_INDEX}"
SAWFISH_EXPECTED_BED_MALE="${SAWFISH_EXPECTED_BED_MALE}"
SAWFISH_EXPECTED_BED_FEMALE="${SAWFISH_EXPECTED_BED_FEMALE}"
METHBAT_REGION_TSV="${METHBAT_REGION_TSV}"
MAX_NORM_FEMALE_CHRY_DEPTH="${MAX_NORM_FEMALE_CHRY_DEPTH}"
PARAPHASE_GENOME_BUILD="${PARAPHASE_GENOME_BUILD}"
RUN_STARPHASE="${RUN_STARPHASE}"
EOF
log "Reference assets resolved -> ${OUT_DIR}/ref.env"

# ---------------------------------------------------------------------------
# 2. Validate every sample's HiFi reads BAM(s) and record them.
# ---------------------------------------------------------------------------
while IFS=$'\t' read -r family_id sample_id sex hifi_reads_csv; do
  [[ -z "${family_id}" || "${family_id}" == \#* ]] && continue
  IFS=',' read -r -a bams <<< "${hifi_reads_csv}"
  [[ "${#bams[@]}" -gt 0 ]] || die "Sample ${sample_id}: no hifi_reads listed"
  for bam in "${bams[@]}"; do
    manifest_add "hifi_reads" "${sample_id}" "${bam}"
  done
  log "Sample ${sample_id} (family ${family_id}, sex ${sex}): ${#bams[@]} hifi_reads BAM(s) validated"
done < "${SAMPLES_TSV}"

cp "${SAMPLES_TSV}" "${OUT_DIR}/samples.tsv"

log "Manifest complete: ${MANIFEST}"
log "Run 'run_pipeline.sh --manifest-dir ${OUT_DIR} --data-root ${DATA_ROOT} --mode <singleton|family>' next."
