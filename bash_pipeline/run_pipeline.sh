#!/usr/bin/env bash
# Main entrypoint: runs the full PacBio HiFi human WGS pipeline, in sequence,
# for every sample listed in the manifest produced by collect_manifest.sh.
#
# Usage:
#   run_pipeline.sh --manifest-dir <dir> --data-root <dir> --out <dir> --mode singleton|family
#
# --manifest-dir must be the --out directory previously passed to
# collect_manifest.sh (it contains ref.env and samples.tsv).
#
# Every sample is run through the identical sequence of steps, using the
# identical config.sh parameters and identical output filenames (see
# steps.sh), which is what guarantees "same output files" / "same parameters"
# across every sample in a run.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=config.sh
source "${SCRIPT_DIR}/config.sh"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"
# shellcheck source=steps.sh
source "${SCRIPT_DIR}/steps.sh"

MANIFEST_DIR=""
OUT_ROOT=""
MODE=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --manifest-dir) MANIFEST_DIR="$2"; shift 2 ;;
    --data-root) DATA_ROOT="$2"; shift 2 ;;
    --out) OUT_ROOT="$2"; shift 2 ;;
    --mode) MODE="$2"; shift 2 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

[[ -n "${MANIFEST_DIR}" ]] || die "--manifest-dir <dir> is required"
[[ -n "${DATA_ROOT:-}" ]] || die "--data-root <dir> is required"
[[ -n "${OUT_ROOT}" ]] || die "--out <dir> is required"
[[ "${MODE}" == "singleton" || "${MODE}" == "family" ]] || die "--mode must be 'singleton' or 'family'"
require_cmd docker

DATA_ROOT="$(cd "${DATA_ROOT}" && pwd)"
mkdir -p "${OUT_ROOT}"
OUT_ROOT="$(cd "${OUT_ROOT}" && pwd)"
assert_under_data_root "${OUT_ROOT}"

require_file "${MANIFEST_DIR}/ref.env"
require_file "${MANIFEST_DIR}/samples.tsv"
# shellcheck source=/dev/null
source "${MANIFEST_DIR}/ref.env"

MAX_NORM_FEMALE_CHRY_DEPTH="${MAX_NORM_FEMALE_CHRY_DEPTH:-${MOSDEPTH_MAX_NORM_FEMALE_CHRY_DEPTH_DEFAULT}}"

log "=== HiFi human WGS bash pipeline ==="
log "mode=${MODE} ref_name=${REF_NAME} data_root=${DATA_ROOT} out=${OUT_ROOT} threads=${THREADS}"

declare -a ALL_SAMPLE_IDS=()
declare -A SAMPLE_SEX=()
declare -A SAMPLE_ALIGNED_BAM=()
declare -A SAMPLE_ALIGNED_BAI=()
declare -A SAMPLE_GVCF=()
declare -A SAMPLE_GVCF_TBI=()
declare -A SAMPLE_DISCOVER_DIR=()
declare -A SAMPLE_SV_VCF=()
declare -A SAMPLE_SV_VCF_TBI=()
declare -A SAMPLE_SMALL_VCF=()
declare -A SAMPLE_SMALL_VCF_TBI=()
declare -A SAMPLE_FAMILY=()
declare -a FAMILY_IDS=()

# ---------------------------------------------------------------------------
# Per-sample upstream: align, coverage, small variants, SV discovery,
# paraphase, mitorsaw, kivvi.
# ---------------------------------------------------------------------------
run_upstream_for_sample() {
  local family_id="$1" sample_id="$2" sex="$3" hifi_reads_csv="$4"
  local sample_out="${OUT_ROOT}/${sample_id}/upstream"
  mkdir -p "${sample_out}"

  IFS=',' read -r -a bams <<< "${hifi_reads_csv}"

  log "[${sample_id}] pbmm2 align (${#bams[@]} input BAM(s))"
  local align_out
  align_out="$(step_pbmm2_align "${sample_id}" "${REF_FASTA}" "${REF_NAME}" "${sample_out}/align" "${bams[@]}")"
  local aligned_bam aligned_bai
  IFS=$'\t' read -r aligned_bam aligned_bai <<< "${align_out}"

  log "[${sample_id}] mosdepth"
  local mosdepth_out
  mosdepth_out="$(step_mosdepth "${sample_id}" "${REF_NAME}" "${aligned_bam}" "${aligned_bai}" "${sample_out}/mosdepth" "${MAX_NORM_FEMALE_CHRY_DEPTH}")"
  local inferred_sex
  inferred_sex="$(echo "${mosdepth_out}" | cut -f6)"
  if [[ -z "${sex}" || "${sex}" == "." ]]; then
    sex="${inferred_sex}"
    log "[${sample_id}] sex not provided; inferred sex=${sex:-UNKNOWN}"
  fi

  log "[${sample_id}] DeepVariant"
  local dv_out
  dv_out="$(step_deepvariant "${sample_id}" "${REF_NAME}" "${REF_FASTA}" "${REF_INDEX}" "${aligned_bam}" "${aligned_bai}" "${sample_out}/deepvariant")"
  local dv_vcf dv_vcf_tbi dv_gvcf dv_gvcf_tbi
  IFS=$'\t' read -r dv_vcf dv_vcf_tbi dv_gvcf dv_gvcf_tbi <<< "${dv_out}"

  log "[${sample_id}] sawfish discover"
  local discover_dir
  discover_dir="$(step_sawfish_discover "${sample_id}" "${sex}" "${aligned_bam}" "${aligned_bai}" "${REF_FASTA}" "${REF_INDEX}" "${SAWFISH_EXCLUDE_BED}" "${SAWFISH_EXCLUDE_BED_INDEX}" "${SAWFISH_EXPECTED_BED_MALE}" "${SAWFISH_EXPECTED_BED_FEMALE}" "${sample_out}/sawfish_discover")"

  log "[${sample_id}] paraphase"
  step_paraphase "${sample_id}" "${aligned_bam}" "${aligned_bai}" "${REF_FASTA}" "${REF_INDEX}" "${PARAPHASE_GENOME_BUILD}" "${sample_out}/paraphase" >/dev/null

  log "[${sample_id}] mitorsaw"
  step_mitorsaw "${aligned_bam}" "${aligned_bai}" "${REF_FASTA}" "${REF_INDEX}" "${sample_id}.${REF_NAME}" "${sample_out}/mitorsaw" >/dev/null

  log "[${sample_id}] kivvi (kiv2, d4z4)"
  step_kivvi kiv2 "${aligned_bam}" "${aligned_bai}" "${sample_id}.${REF_NAME}" "${sample_out}/kivvi_kiv2" >/dev/null
  step_kivvi d4z4 "${aligned_bam}" "${aligned_bai}" "${sample_id}.${REF_NAME}" "${sample_out}/kivvi_d4z4" >/dev/null

  ALL_SAMPLE_IDS+=("${sample_id}")
  SAMPLE_FAMILY["${sample_id}"]="${family_id}"
  SAMPLE_SEX["${sample_id}"]="${sex}"
  SAMPLE_ALIGNED_BAM["${sample_id}"]="${aligned_bam}"
  SAMPLE_ALIGNED_BAI["${sample_id}"]="${aligned_bai}"
  SAMPLE_GVCF["${sample_id}"]="${dv_gvcf}"
  SAMPLE_GVCF_TBI["${sample_id}"]="${dv_gvcf_tbi}"
  SAMPLE_DISCOVER_DIR["${sample_id}"]="${discover_dir}"

  if [[ "${MODE}" == "singleton" ]]; then
    log "[${sample_id}] sawfish call (single-sample)"
    local sc_out
    sc_out="$(step_sawfish_call "${sample_id}.${REF_NAME}.structural_variants" "${sample_out}/sawfish_call" "${REF_FASTA}" "${REF_INDEX}" "${sample_id}" "${discover_dir}" "${aligned_bam}" "${aligned_bai}")"
    local sv_vcf sv_vcf_tbi
    IFS=$'\t' read -r sv_vcf sv_vcf_tbi _ <<< "${sc_out}"
    SAMPLE_SV_VCF["${sample_id}"]="${sv_vcf}"
    SAMPLE_SV_VCF_TBI["${sample_id}"]="${sv_vcf_tbi}"
    SAMPLE_SMALL_VCF["${sample_id}"]="${dv_vcf}"
    SAMPLE_SMALL_VCF_TBI["${sample_id}"]="${dv_vcf_tbi}"
  fi
}

# ---------------------------------------------------------------------------
# Per-family joint calling (family mode only): sawfish joint-call + GLnexus,
# each split back out by sample.
# ---------------------------------------------------------------------------
run_joint_for_family() {
  local family_id="$1"
  shift
  local sample_ids=("$@")
  local joint_out="${OUT_ROOT}/${family_id}.joint"
  mkdir -p "${joint_out}"

  local discover_dirs="" aligned_bams="" aligned_bam_indices="" sample_ids_csv="" gvcfs=()
  local sid first=1
  for sid in "${sample_ids[@]}"; do
    if [[ ${first} -eq 1 ]]; then
      sample_ids_csv="${sid}"; discover_dirs="${SAMPLE_DISCOVER_DIR[$sid]}"
      aligned_bams="${SAMPLE_ALIGNED_BAM[$sid]}"; aligned_bam_indices="${SAMPLE_ALIGNED_BAI[$sid]}"
      first=0
    else
      sample_ids_csv+=",${sid}"; discover_dirs+=",${SAMPLE_DISCOVER_DIR[$sid]}"
      aligned_bams+=",${SAMPLE_ALIGNED_BAM[$sid]}"; aligned_bam_indices+=",${SAMPLE_ALIGNED_BAI[$sid]}"
    fi
    gvcfs+=("${SAMPLE_GVCF[$sid]}")
  done

  log "[${family_id}] sawfish joint-call (${#sample_ids[@]} samples)"
  local sc_out
  sc_out="$(step_sawfish_call "${family_id}.joint.${REF_NAME}.structural_variants" "${joint_out}/sawfish_call" "${REF_FASTA}" "${REF_INDEX}" "${sample_ids_csv}" "${discover_dirs}" "${aligned_bams}" "${aligned_bam_indices}")"
  local joint_sv_vcf joint_sv_vcf_tbi
  IFS=$'\t' read -r joint_sv_vcf joint_sv_vcf_tbi _ <<< "${sc_out}"

  log "[${family_id}] splitting joint SV VCF by sample"
  step_split_vcf_by_sample "${joint_sv_vcf}" "${joint_sv_vcf_tbi}" false "${joint_out}/split_sv" "${sample_ids[@]}"

  log "[${family_id}] GLnexus joint genotyping"
  local gl_out
  gl_out="$(step_glnexus "${family_id}.joint" "${REF_NAME}" "${joint_out}/glnexus" "${gvcfs[@]}")"
  local joint_small_vcf joint_small_vcf_tbi
  IFS=$'\t' read -r joint_small_vcf joint_small_vcf_tbi <<< "${gl_out}"

  log "[${family_id}] splitting joint small-variant VCF by sample"
  step_split_vcf_by_sample "${joint_small_vcf}" "${joint_small_vcf_tbi}" true "${joint_out}/split_small" "${sample_ids[@]}"

  local sv_base small_base
  sv_base="$(basename "${joint_sv_vcf}")"
  small_base="$(basename "${joint_small_vcf}")"
  for sid in "${sample_ids[@]}"; do
    SAMPLE_SV_VCF["${sid}"]="${joint_out}/split_sv/${sid}.${sv_base}"
    SAMPLE_SV_VCF_TBI["${sid}"]="${joint_out}/split_sv/${sid}.${sv_base}.tbi"
    SAMPLE_SMALL_VCF["${sid}"]="${joint_out}/split_small/${sid}.${small_base}"
    SAMPLE_SMALL_VCF_TBI["${sid}"]="${joint_out}/split_small/${sid}.${small_base}.tbi"
  done
}

# ---------------------------------------------------------------------------
# Per-sample downstream: phasing, TRGT, stats, methylation, PGx.
# ---------------------------------------------------------------------------
run_downstream_for_sample() {
  local sample_id="$1"
  local sex="${SAMPLE_SEX[$sample_id]}"
  local out="${OUT_ROOT}/${sample_id}/downstream"
  mkdir -p "${out}"

  log "[${sample_id}] hiphase"
  local hp_out
  hp_out="$(step_hiphase "${sample_id}" "${REF_NAME}" "${REF_FASTA}" "${REF_INDEX}" \
    "${SAMPLE_SMALL_VCF[$sample_id]}" "${SAMPLE_SMALL_VCF_TBI[$sample_id]}" \
    "${SAMPLE_SV_VCF[$sample_id]}" "${SAMPLE_SV_VCF_TBI[$sample_id]}" \
    "${SAMPLE_ALIGNED_BAM[$sample_id]}" "${SAMPLE_ALIGNED_BAI[$sample_id]}" "${out}/hiphase")"
  local phased_small_vcf phased_small_tbi phased_sv_vcf phased_sv_tbi haplotagged_bam haplotagged_bai
  IFS=$'\t' read -r phased_small_vcf phased_small_tbi phased_sv_vcf phased_sv_tbi haplotagged_bam haplotagged_bai _ <<< "${hp_out}"

  log "[${sample_id}] trgt genotype"
  step_trgt "${sample_id}" "${REF_NAME}" "${sex}" "${haplotagged_bam}" "${haplotagged_bai}" \
    "${REF_FASTA}" "${REF_INDEX}" "${TRGT_BED}" "${SAWFISH_EXPECTED_BED_MALE}" "${SAWFISH_EXPECTED_BED_FEMALE}" \
    "${out}/trgt" >/dev/null

  log "[${sample_id}] pbjam bam-stats"
  step_pbjam_bam_stats "${sample_id}" "${REF_NAME}" "${haplotagged_bam}" "${haplotagged_bai}" "${out}/pbjam" >/dev/null

  log "[${sample_id}] bcftools stats/roh"
  step_bcftools_stats_roh "${sample_id}" "${REF_NAME}" "${phased_small_vcf}" "${REF_FASTA}" "${out}/bcftools_small_variants" >/dev/null

  log "[${sample_id}] sv_stats"
  step_sv_stats "${sample_id}" "${REF_NAME}" "${phased_sv_vcf}" "${out}/sv_stats" >/dev/null

  log "[${sample_id}] methbat pileup"
  local mb_out
  mb_out="$(step_methbat_pileup "${haplotagged_bam}" "${haplotagged_bai}" "${sample_id}.${REF_NAME}" "${out}/methbat_pileup")"
  local cpg_bed
  cpg_bed="$(printf '%s\n' "${mb_out}" | head -n1 | cut -f1)"

  log "[${sample_id}] methbat profile"
  step_methbat_profile "${cpg_bed}" "${METHBAT_REGION_TSV}" "${sample_id}.${REF_NAME}" "${out}/methbat_profile" >/dev/null

  if [[ "${RUN_STARPHASE}" == "true" ]]; then
    log "[${sample_id}] pbstarphase"
    step_pbstarphase "${sample_id}" "${phased_small_vcf}" "${phased_small_tbi}" "${phased_sv_vcf}" "${phased_sv_tbi}" \
      "${haplotagged_bam}" "${haplotagged_bai}" "${REF_FASTA}" "${REF_INDEX}" "${out}/pbstarphase" >/dev/null
  fi

  log "[${sample_id}] done -> ${out}"
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
while IFS=$'\t' read -r family_id sample_id sex hifi_reads_csv; do
  [[ -z "${family_id}" || "${family_id}" == \#* ]] && continue
  run_upstream_for_sample "${family_id}" "${sample_id}" "${sex}" "${hifi_reads_csv}"
  if [[ "${MODE}" == "family" ]]; then
    if [[ ! " ${FAMILY_IDS[*]-} " == *" ${family_id} "* ]]; then
      FAMILY_IDS+=("${family_id}")
    fi
  fi
done < "${MANIFEST_DIR}/samples.tsv"

if [[ "${MODE}" == "family" ]]; then
  for fam in "${FAMILY_IDS[@]}"; do
    fam_samples=()
    for sid in "${ALL_SAMPLE_IDS[@]}"; do
      [[ "${SAMPLE_FAMILY[$sid]}" == "${fam}" ]] && fam_samples+=("${sid}")
    done
    run_joint_for_family "${fam}" "${fam_samples[@]}"
  done
fi

for sid in "${ALL_SAMPLE_IDS[@]}"; do
  run_downstream_for_sample "${sid}"
done

log "Pipeline complete. Outputs under ${OUT_ROOT}/<sample_id>/{upstream,downstream}$( [[ "${MODE}" == "family" ]] && echo " and ${OUT_ROOT}/<family_id>.joint" )."
