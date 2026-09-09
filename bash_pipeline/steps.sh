#!/usr/bin/env bash
# One bash function per WDL task/workflow step. Every function:
#   - takes plain positional arguments (no WDL structs)
#   - writes outputs into a caller-supplied output directory using the same
#     naming convention the WDL tasks use, so every sample produces the same
#     set of output files under the same names
#   - reads its tool parameters from config.sh, so every sample/run sourcing
#     the same config.sh uses identical parameters
#
# All functions assume config.sh and lib/common.sh have already been sourced,
# and that DATA_ROOT is set.
set -euo pipefail

# ===========================================================================
# Alignment: pbmm2 index + align (+ optional chunking) + pbsamoa merge
# ===========================================================================
# step_pbmm2_align <sample_id> <ref_fasta> <ref_name> <outdir> <bam...>
# Prints "aligned_bam<TAB>aligned_bam_index" on success.
step_pbmm2_align() {
  local sample_id="$1" ref_fasta="$2" ref_name="$3" outdir="$4"
  shift 4
  local bams=("$@")

  mkdir -p "${outdir}"
  local idx_dir="${outdir}/pbmm2_index"
  local ref_base mmi
  ref_base="$(basename "${ref_fasta}")"
  mmi="${ref_base%.*}.mmi"
  if [[ ! -f "${idx_dir}/${mmi}" ]]; then
    link_input "${ref_fasta}" "${idx_dir}"
    run_docker "${IMG_PBMM2}" "${idx_dir}" \
      "pbmm2 --version >&2; pbmm2 index --num-threads ${THREADS} --log-level INFO --preset HIFI '${ref_base}' '${mmi}'"
  fi
  local pbmm2_index="${idx_dir}/${mmi}"

  local aligned_bams=()
  local movie_idx=0
  local bam movie work bam_base is_aligned num_chunks
  for bam in "${bams[@]}"; do
    movie_idx=$((movie_idx + 1))
    movie="$(basename "${bam}" .bam)"
    work="${outdir}/align_${movie_idx}_${movie}"
    mkdir -p "${work}"
    link_input "${bam}" "${work}"
    link_input "${pbmm2_index}" "${work}"
    bam_base="$(basename "${bam}")"

    is_aligned="$(run_docker "${IMG_BASE}" "${work}" \
      "samtools view -H '${bam_base}' | grep -c '^@SQ' || true" | tr -d '[:space:]')"

    num_chunks=0
    if [[ "${USE_ALIGNMENT_CHUNKING}" == "true" && "${is_aligned}" == "0" ]]; then
      if ! resume_available "${work}/${bam_base}.pbi"; then
        run_docker "${IMG_PBTK}" "${work}" \
          "pbindex --version >&2; pbindex --num-threads ${THREADS} '${bam_base}'"
      fi
      num_chunks="${ALIGNMENT_CHUNKS}"
    fi

    local strip_flag="" unmapped_flag=""
    [[ "${PBMM2_STRIP_KINETICS}" == "true" ]] && strip_flag="--strip"
    [[ "${PBMM2_KEEP_UNMAPPED}" == "true" ]] && unmapped_flag="--unmapped"

    if [[ "${num_chunks}" -eq 0 ]]; then
      local out_bam="${sample_id}.${movie}.chunk_0.${ref_name}.aligned.bam"
      if resume_available "${work}/${out_bam}" "${work}/${out_bam}.bai"; then
        log "[resume] pbmm2 align: ${out_bam} already present, skipping"
      else
        run_docker "${IMG_PBMM2}" "${work}" \
          "pbmm2 --version >&2; pbmm2 align --num-threads ${THREADS} --sort-memory 4G --preset HIFI --sample '${sample_id}' --log-level INFO --sort --strip-tags HP,PS,PC ${strip_flag} ${unmapped_flag} --min-length ${PBMM2_MIN_LENGTH} '${mmi}' '${bam_base}' '${out_bam}'"
      fi
      aligned_bams+=("${work}/${out_bam}")
    else
      local c
      for c in $(seq 1 "${num_chunks}"); do
        local out_bam="${sample_id}.${movie}.chunk_${c}.${ref_name}.aligned.bam"
        if resume_available "${work}/${out_bam}" "${work}/${out_bam}.bai"; then
          log "[resume] pbmm2 align: ${out_bam} already present, skipping"
        else
          run_docker "${IMG_PBMM2}" "${work}" \
            "pbmm2 --version >&2; pbmm2 align --num-threads ${THREADS} --sort-memory 4G --preset HIFI --sample '${sample_id}' --log-level INFO --sort --strip-tags HP,PS,PC ${strip_flag} ${unmapped_flag} --min-length ${PBMM2_MIN_LENGTH} --chunk '${c}/${num_chunks}' --chunk-mode scatter '${mmi}' '${bam_base}' '${out_bam}'"
        fi
        aligned_bams+=("${work}/${out_bam}")
      done
    fi
  done

  local final_bam final_bai
  if [[ "${#aligned_bams[@]}" -gt 1 ]]; then
    local merge_dir="${outdir}/merge"
    mkdir -p "${merge_dir}"
    local out_prefix="${sample_id}.${ref_name}.hifi_reads"
    final_bam="${merge_dir}/${out_prefix}.bam"
    final_bai="${merge_dir}/${out_prefix}.bam.bai"
    if resume_available "${final_bam}" "${final_bai}"; then
      log "[resume] pbsamoa merge: ${out_prefix}.bam already present, skipping"
    else
      local merge_list="" b
      for b in "${aligned_bams[@]}"; do
        link_input "${b}" "${merge_dir}"
        link_input "${b}.bai" "${merge_dir}"
        merge_list+=" $(basename "${b}")"
      done
      run_docker "${IMG_PBSAMOA}" "${merge_dir}" \
        "pbsamoa merge --compress-threads $(( (THREADS * 3 + 3) / 4 )) --decode-threads $(( THREADS / 4 > 0 ? THREADS / 4 : 1 )) --memory '$(( MERGE_MEM_GB / 2 > 0 ? MERGE_MEM_GB / 2 : 1 ))G' --compression 6 --bai '${out_prefix}.bam'${merge_list}"
    fi
  else
    final_bam="${aligned_bams[0]}"
    final_bai="${aligned_bams[0]}.bai"
  fi

  printf '%s\t%s\n' "${final_bam}" "${final_bai}"
}

# ===========================================================================
# mosdepth: coverage + sex inference
# ===========================================================================
# step_mosdepth <sample_id> <ref_name> <aligned_bam> <aligned_bam_index> <outdir> <max_norm_female_chrY_depth>
step_mosdepth() {
  local sample_id="$1" ref_name="$2" aligned_bam="$3" aligned_bam_index="$4" outdir="$5" max_norm_female_chrY_depth="$6"
  mkdir -p "${outdir}"
  link_input "${aligned_bam}" "${outdir}"
  link_input "${aligned_bam_index}" "${outdir}"
  local bam_base out_prefix
  bam_base="$(basename "${aligned_bam}")"
  out_prefix="${bam_base%.bam}"

  local summary="${sample_id}.${ref_name}.mosdepth.summary.txt"
  local region_bed="${sample_id}.${ref_name}.mosdepth.regions.bed.gz"
  local region_bed_index="${sample_id}.${ref_name}.mosdepth.regions.bed.gz.csi"
  if resume_available "${outdir}/${summary}" "${outdir}/${region_bed}" "${outdir}/${region_bed_index}"; then
    log "[resume] mosdepth: ${summary} already present, skipping"
  else
    run_docker "${IMG_MOSDEPTH}" "${outdir}" \
      "mosdepth --version >&2; mosdepth --threads $(( THREADS > 1 ? THREADS - 1 : 0 )) --by 500 --no-per-base --use-median '${out_prefix}' '${bam_base}'"

    if [[ ! -f "${outdir}/${summary}" ]]; then
      mv --verbose "${outdir}/${out_prefix}.mosdepth.summary.txt" "${outdir}/${summary}"
      mv --verbose "${outdir}/${out_prefix}.regions.bed.gz" "${outdir}/${region_bed}"
      mv --verbose "${outdir}/${out_prefix}.regions.bed.gz.csi" "${outdir}/${region_bed_index}"
    fi
  fi

  local plot="${sample_id}.${ref_name}.mosdepth.depth_distribution.png"
  if ! resume_available "${outdir}/${plot}"; then
    run_docker "${IMG_BASE}" "${outdir}" \
      "python3 /pipeline-scripts/plot_depth_distribution.py '${region_bed}' '${sample_id}.${ref_name}' '${plot}'"
  fi

  if ! resume_available "${outdir}/mean_depth.txt" "${outdir}/inferred_sex.txt"; then
    run_docker "${IMG_BASE}" "${outdir}" \
      "python3 /pipeline-scripts/mosdepth_stats.py '${summary}' '${max_norm_female_chrY_depth}' mean_depth.txt inferred_sex.txt"
  fi

  local mean_depth inferred_sex
  mean_depth="$(cat "${outdir}/mean_depth.txt")"
  inferred_sex="$(cat "${outdir}/inferred_sex.txt")"

  printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${outdir}/${summary}" "${outdir}/${region_bed}" "${outdir}/${region_bed_index}" \
    "${outdir}/${plot}" "${mean_depth}" "${inferred_sex}"
}

# ===========================================================================
# DeepVariant (CPU/GPU) or Parabricks DeepVariant
# ===========================================================================
# step_deepvariant <sample_id> <ref_name> <ref_fasta> <ref_index> <aligned_bam> <aligned_bam_index> <outdir>
# Prints "vcf<TAB>vcf_index<TAB>gvcf<TAB>gvcf_index"
step_deepvariant() {
  local sample_id="$1" ref_name="$2" ref_fasta="$3" ref_index="$4" aligned_bam="$5" aligned_bam_index="$6" outdir="$7"
  mkdir -p "${outdir}"
  link_input "${ref_fasta}" "${outdir}"
  link_input "${ref_index}" "${outdir}"
  link_input "${aligned_bam}" "${outdir}"
  link_input "${aligned_bam_index}" "${outdir}"

  local ref_base bam_base
  ref_base="$(basename "${ref_fasta}")"
  bam_base="$(basename "${aligned_bam}")"
  local vcf="${sample_id}.${ref_name}.small_variants.vcf.gz"
  local gvcf="${sample_id}.${ref_name}.small_variants.g.vcf.gz"
  local gvcf_flags=""
  [[ "${DEEPVARIANT_GVCF_OUTPUT}" == "true" ]] && gvcf_flags="--output_gvcf=${gvcf}"

  local -a dv_expected=("${outdir}/${vcf}" "${outdir}/${vcf}.tbi")
  [[ "${DEEPVARIANT_GVCF_OUTPUT}" == "true" ]] && dv_expected+=("${outdir}/${gvcf}" "${outdir}/${gvcf}.tbi")

  if resume_available "${dv_expected[@]}"; then
    log "[resume] DeepVariant: ${vcf} already present, skipping"
  elif [[ "${USE_GPU}" == "true" && "${USE_PARABRICKS_DEEPVARIANT}" == "true" ]]; then
    # Parabricks path (requires an NVIDIA GPU + nvidia-container-toolkit)
    local raw_vcf raw_gvcf
    raw_vcf="${sample_id}.${ref_name}.small_variants.vcf"
    raw_gvcf="${sample_id}.${ref_name}.small_variants.g.vcf"
    local pb_gvcf_flag="" pb_out="${raw_vcf}"
    if [[ "${DEEPVARIANT_GVCF_OUTPUT}" == "true" ]]; then
      pb_gvcf_flag="--gvcf"
      pb_out="${raw_gvcf}"
    fi
    run_docker_gpu "${IMG_PARABRICKS}" "${outdir}" \
      "export TCMALLOC_MAX_TOTAL_THREAD_CACHE_BYTES=268435456; /usr/local/parabricks/pbrun deepvariant --num-gpus 1 --preserve-file-symlinks --run-partition --mode pacbio --ref '${ref_base}' ${pb_gvcf_flag} --in-bam '${bam_base}' --out-variants '${pb_out}'" \
      "${PARABRICKS_GPU_DOCKER_ARGS}"
    if [[ "${DEEPVARIANT_GVCF_OUTPUT}" == "true" ]]; then
      run_docker "${IMG_BASE}" "${outdir}" \
        "bcftools view --output-type z --output-file '${gvcf}' '${raw_gvcf}'; bcftools index --tbi --force '${gvcf}'; rm '${raw_gvcf}'"
    fi
    run_docker "${IMG_BASE}" "${outdir}" \
      "bcftools view --exclude-uncalled --output-type z --output-file '${vcf}' '${raw_vcf}'; bcftools index --tbi --force '${vcf}'; rm '${raw_vcf}'"
  elif [[ "${USE_GPU}" == "true" ]]; then
    # DeepVariant GPU path, without Parabricks
    run_docker_gpu "${IMG_DEEPVARIANT_GPU}" "${outdir}" \
      "/opt/deepvariant/bin/run_deepvariant --model_type=PACBIO --ref='${ref_base}' --reads='${bam_base}' --sample_name='${sample_id}' --output_vcf='${vcf}' ${gvcf_flags} --num_shards=${THREADS}" \
      "${DEEPVARIANT_GPU_DOCKER_ARGS}"
  else
    run_docker "${IMG_DEEPVARIANT_CPU}" "${outdir}" \
      "/opt/deepvariant/bin/run_deepvariant --model_type=PACBIO --ref='${ref_base}' --reads='${bam_base}' --sample_name='${sample_id}' --output_vcf='${vcf}' ${gvcf_flags} --num_shards=${THREADS}"
  fi

  local gvcf_index=""
  [[ "${DEEPVARIANT_GVCF_OUTPUT}" == "true" ]] && gvcf_index="${outdir}/${gvcf}.tbi"
  printf '%s\t%s\t%s\t%s\n' "${outdir}/${vcf}" "${outdir}/${vcf}.tbi" \
    "${outdir}/${gvcf}" "${gvcf_index}"
}

# ===========================================================================
# Sawfish: SV discovery + (single-sample or joint) calling
# ===========================================================================
# step_sawfish_discover <sample_id> <sex> <aligned_bam> <aligned_bam_index> <ref_fasta> <ref_index>
#   <exclude_bed> <exclude_bed_index> <expected_male_bed> <expected_female_bed> <outdir>
# Prints "discover_dir" (the pipeline keeps discover output as a directory
# instead of the WDL's tar, since we never leave the local filesystem)
step_sawfish_discover() {
  local sample_id="$1" sex="$2" aligned_bam="$3" aligned_bam_index="$4" ref_fasta="$5" ref_index="$6" \
    exclude_bed="$7" exclude_bed_index="$8" expected_male_bed="$9" expected_female_bed="${10}" outdir="${11}"
  mkdir -p "${outdir}"
  link_input "${aligned_bam}" "${outdir}"
  link_input "${aligned_bam_index}" "${outdir}"
  link_input "${ref_fasta}" "${outdir}"
  link_input "${ref_index}" "${outdir}"
  link_input "${exclude_bed}" "${outdir}"
  link_input "${exclude_bed_index}" "${outdir}"

  local expected_bed="${expected_female_bed}"
  [[ "${sex}" == "MALE" ]] && expected_bed="${expected_male_bed}"
  link_input "${expected_bed}" "${outdir}"

  local out_prefix="${sample_id}"
  if resume_available "${outdir}/${out_prefix}"; then
    log "[resume] sawfish discover: ${out_prefix} already present, skipping"
  else
    run_docker "${IMG_SAWFISH}" "${outdir}" \
      "sawfish --version >&2; sawfish discover --threads ${THREADS} --disable-path-canonicalization --ref '$(basename "${ref_fasta}")' --bam '$(basename "${aligned_bam}")' --expected-cn '$(basename "${expected_bed}")' --cnv-excluded-regions '$(basename "${exclude_bed}")' --output-dir '${out_prefix}'"
  fi

  printf '%s\n' "${outdir}/${out_prefix}"
}

# step_sawfish_call <out_prefix> <outdir> <ref_fasta> <ref_index> <sample_ids_csv> <discover_dirs_csv> <aligned_bams_csv> <aligned_bam_indices_csv>
# Used both for the singleton per-sample call and the family joint call.
# Prints "vcf<TAB>vcf_index<TAB>supporting_reads"
step_sawfish_call() {
  local out_prefix="$1" outdir="$2" ref_fasta="$3" ref_index="$4" \
    sample_ids_csv="$5" discover_dirs_csv="$6" aligned_bams_csv="$7" aligned_bam_indices_csv="$8"
  mkdir -p "${outdir}"

  local vcf="${outdir}/${out_prefix}.vcf.gz"
  local vcf_index="${outdir}/${out_prefix}.vcf.gz.tbi"

  if resume_available "${vcf}" "${vcf_index}"; then
    log "[resume] sawfish call: ${out_prefix}.vcf.gz already present, skipping"
  else
    link_input "${ref_fasta}" "${outdir}"
    link_input "${ref_index}" "${outdir}"

    local sample_ids discover_dirs aligned_bams aligned_bam_indices
    IFS=',' read -r -a sample_ids <<< "${sample_ids_csv}"
    IFS=',' read -r -a discover_dirs <<< "${discover_dirs_csv}"
    IFS=',' read -r -a aligned_bams <<< "${aligned_bams_csv}"
    IFS=',' read -r -a aligned_bam_indices <<< "${aligned_bam_indices_csv}"

    local b i
    for b in "${aligned_bams[@]}" "${aligned_bam_indices[@]}"; do
      link_input "${b}" "${outdir}"
    done
    # sawfish reads discover output from a directory named after the sample; copy
    # (not symlink) each discover directory into the call working dir under its
    # own name, matching the WDL's `tar --extract` of one dir per sample.
    local samples_flags=""
    for i in "${!sample_ids[@]}"; do
      local sid="${sample_ids[$i]}" ddir="${discover_dirs[$i]}"
      rm --recursive --force "${outdir}/${sid}"
      cp --recursive "${ddir}" "${outdir}/${sid}"
      samples_flags+=" --sample '${sid}'"
    done

    run_docker "${IMG_SAWFISH}" "${outdir}" \
      "sawfish --version >&2; sawfish joint-call --threads ${THREADS} --report-supporting-reads${samples_flags} --output-dir '${out_prefix}'"

    run_docker "${IMG_SAWFISH}" "${outdir}" \
      "sawshark --threads $(( THREADS / 2 > 0 ? THREADS / 2 : 1 )) --vcf '${out_prefix}/genotyped.sv.vcf.gz' | bcftools view - --output-type z --write-index=tbi --output '${out_prefix}.vcf.gz'"

    mv --verbose "${outdir}/${out_prefix}/supporting_reads.json.gz" "${outdir}/${out_prefix}.supporting_reads.json.gz" 2>/dev/null || true

    for i in "${!sample_ids[@]}"; do
      local sid="${sample_ids[$i]}" prefix
      if [[ "${#sample_ids[@]}" -gt 1 ]]; then
        prefix="${sid}.${out_prefix}"
      else
        prefix="${out_prefix}"
      fi
      local sampledir
      sampledir="$(find "${outdir}/${out_prefix}/samples" -maxdepth 1 -type d -name "sample????_${sid}" | head -n1)"
      if [[ -n "${sampledir}" ]]; then
        mv --verbose "${sampledir}/copynum.bedgraph" "${outdir}/${prefix}.copynum.bedgraph"
        mv --verbose "${sampledir}/depth.bw" "${outdir}/${prefix}.depth.bw"
        mv --verbose "${sampledir}/gc_bias_corrected_depth.bw" "${outdir}/${prefix}.gc_bias_corrected_depth.bw"
        mv --verbose "${sampledir}/copynum.summary.json" "${outdir}/${prefix}.copynum.summary.json"
      fi
    done
    for sid in "${sample_ids[@]}"; do rm --recursive --force "${outdir}/${sid}"; done
  fi

  printf '%s\t%s\t%s\n' "${vcf}" "${vcf_index}" \
    "${outdir}/${out_prefix}.supporting_reads.json.gz"
}

# ===========================================================================
# Paraphase / Mitorsaw / Kivvi (segmental-dup genes, mtDNA, KIV2 & D4Z4 repeats)
# ===========================================================================
# step_paraphase <sample_id> <aligned_bam> <aligned_bam_index> <ref_fasta> <ref_index> <genome_build> <outdir>
step_paraphase() {
  local sample_id="$1" aligned_bam="$2" aligned_bam_index="$3" ref_fasta="$4" ref_index="$5" genome_build="$6" outdir="$7"
  mkdir -p "${outdir}"
  link_input "${aligned_bam}" "${outdir}"
  link_input "${aligned_bam_index}" "${outdir}"
  link_input "${ref_fasta}" "${outdir}"
  link_input "${ref_index}" "${outdir}"

  if resume_available "${outdir}/${sample_id}.paraphase.json"; then
    log "[resume] paraphase: ${sample_id}.paraphase.json already present, skipping"
  else
    run_docker "${IMG_PARAPHASE}" "${outdir}" \
      "paraphase --version >&2; paraphase --threads ${THREADS} --bam '$(basename "${aligned_bam}")' --reference '$(basename "${ref_fasta}")' --genome '${genome_build}' --out ./ 2>&1 | tee paraphase.log || echo 'Paraphase failed for sample ${sample_id}' >> messages.txt; if ls '${sample_id}_paraphase_vcfs'/*.vcf &> /dev/null; then tar --gzip --create --file '${sample_id}.paraphase_vcfs.tar.gz' '${sample_id}_paraphase_vcfs'/*.vcf; fi"
  fi

  printf '%s\t%s\t%s\t%s\n' \
    "${outdir}/${sample_id}.paraphase.json" "${outdir}/${sample_id}.paraphase.bam" \
    "${outdir}/${sample_id}.paraphase.bam.bai" "${outdir}/${sample_id}.paraphase_vcfs.tar.gz"
}

# step_mitorsaw <aligned_bam> <aligned_bam_index> <ref_fasta> <ref_index> <out_prefix> <outdir>
step_mitorsaw() {
  local aligned_bam="$1" aligned_bam_index="$2" ref_fasta="$3" ref_index="$4" out_prefix="$5" outdir="$6"
  mkdir -p "${outdir}"
  link_input "${aligned_bam}" "${outdir}"
  link_input "${aligned_bam_index}" "${outdir}"
  link_input "${ref_fasta}" "${outdir}"
  link_input "${ref_index}" "${outdir}"

  if resume_available "${outdir}/${out_prefix}.mitorsaw.vcf.gz" "${outdir}/${out_prefix}.mitorsaw.vcf.gz.tbi" "${outdir}/${out_prefix}.mitorsaw.json"; then
    log "[resume] mitorsaw: ${out_prefix}.mitorsaw.vcf.gz already present, skipping"
  else
    run_docker "${IMG_MITORSAW}" "${outdir}" \
      "mitorsaw --version >&2; mitorsaw haplotype --reference '$(basename "${ref_fasta}")' --bam '$(basename "${aligned_bam}")' --output-vcf '${out_prefix}.mitorsaw.vcf.gz' --output-hap-stats '${out_prefix}.mitorsaw.json'"
  fi

  printf '%s\t%s\t%s\n' "${outdir}/${out_prefix}.mitorsaw.vcf.gz" "${outdir}/${out_prefix}.mitorsaw.vcf.gz.tbi" \
    "${outdir}/${out_prefix}.mitorsaw.json"
}

# step_kivvi <mode: kiv2|d4z4> <aligned_bam> <aligned_bam_index> <out_prefix> <outdir>
step_kivvi() {
  local mode="$1" aligned_bam="$2" aligned_bam_index="$3" out_prefix="$4" outdir="$5"
  mkdir -p "${outdir}"
  link_input "${aligned_bam}" "${outdir}"
  link_input "${aligned_bam_index}" "${outdir}"

  if resume_available "${outdir}/messages.txt"; then
    log "[resume] kivvi ${mode}: messages.txt already present, skipping"
  else
    run_docker "${IMG_KIVVI}" "${outdir}" \
      "touch messages.txt; kivvi --version >&2; kivvi --bam '$(basename "${aligned_bam}")' --out . --prefix '${out_prefix}' ${mode} || echo 'kivvi ${mode} failed, presumably due to low coverage.' >> messages.txt; if [ -f '${out_prefix}.kivvi.${mode}.vcf' ]; then bgzip '${out_prefix}.kivvi.${mode}.vcf'; tabix --preset vcf '${out_prefix}.kivvi.${mode}.vcf.gz'; fi"
  fi

  printf '%s\t%s\n' "${outdir}/${out_prefix}.kivvi.${mode}.vcf.gz" "${outdir}/${out_prefix}.kivvi.${mode}.vcf.gz.tbi"
}

# ===========================================================================
# HiPhase: phase small variants + SVs, haplotag alignments
# ===========================================================================
# step_hiphase <sample_id> <ref_name> <ref_fasta> <ref_index> <small_variant_vcf> <small_variant_vcf_index>
#   <sv_vcf> <sv_vcf_index> <aligned_bam> <aligned_bam_index> <outdir>
step_hiphase() {
  local sample_id="$1" ref_name="$2" ref_fasta="$3" ref_index="$4" \
    small_variant_vcf="$5" small_variant_vcf_index="$6" sv_vcf="$7" sv_vcf_index="$8" \
    aligned_bam="$9" aligned_bam_index="${10}" outdir="${11}"
  mkdir -p "${outdir}"
  link_input "${ref_fasta}" "${outdir}"
  link_input "${ref_index}" "${outdir}"
  link_input "${small_variant_vcf}" "${outdir}"
  link_input "${small_variant_vcf_index}" "${outdir}"
  link_input "${sv_vcf}" "${outdir}"
  link_input "${sv_vcf_index}" "${outdir}"
  link_input "${aligned_bam}" "${outdir}"
  link_input "${aligned_bam_index}" "${outdir}"

  local phased_small_vcf phased_sv_vcf haplotagged_bam
  phased_small_vcf="$(basename "${small_variant_vcf}" .vcf.gz).phased.vcf.gz"
  phased_sv_vcf="$(basename "${sv_vcf}" .vcf.gz).phased.vcf.gz"
  haplotagged_bam="${sample_id}.${ref_name}.haplotagged.bam"

  local opt_flags=""
  [[ -n "${HIPHASE_PRESET}" ]] && opt_flags+=" --preset '${HIPHASE_PRESET}'"
  [[ -n "${HIPHASE_MIN_MAPQ}" ]] && opt_flags+=" --min-mapq '${HIPHASE_MIN_MAPQ}'"
  [[ -n "${HIPHASE_MIN_GQ}" ]] && opt_flags+=" --min-vcf-qual '${HIPHASE_MIN_GQ}'"
  [[ "${HIPHASE_NO_SUPPLEMENTAL_JOINS}" == "true" ]] && opt_flags+=" --no-supplemental-joins"
  [[ "${HIPHASE_PHASE_SINGLETONS}" == "true" ]] && opt_flags+=" --phase-singletons"
  [[ "${HIPHASE_DISABLE_GLOBAL_REALIGNMENT}" == "true" ]] && opt_flags+=" --disable-global-realignment"

  local stats="${sample_id}.${ref_name}.hiphase.stats.tsv"
  local blocks="${sample_id}.${ref_name}.hiphase.blocks.tsv"
  if resume_available "${outdir}/${phased_small_vcf}" "${outdir}/${phased_sv_vcf}" "${outdir}/${haplotagged_bam}" "${outdir}/${stats}" "${outdir}/${blocks}"; then
    log "[resume] hiphase: ${haplotagged_bam} already present, skipping"
  else
    run_docker "${IMG_HIPHASE}" "${outdir}" \
      "hiphase --version >&2; hiphase --threads ${THREADS}${opt_flags} --sample-name '${sample_id}' --vcf '$(basename "${small_variant_vcf}")' --vcf '$(basename "${sv_vcf}")' --output-vcf '${phased_small_vcf}' --output-vcf '${phased_sv_vcf}' --bam '$(basename "${aligned_bam}")' --output-bam '${haplotagged_bam}' --reference '$(basename "${ref_fasta}")' --summary-file '${stats}' --blocks-file '${blocks}' --haplotag-file '${sample_id}.${ref_name}.hiphase.haplotags.tsv'; gzip '${sample_id}.${ref_name}.hiphase.haplotags.tsv'"
  fi

  if ! resume_available "${outdir}/${haplotagged_bam}.bai"; then
    run_docker "${IMG_BASE}" "${outdir}" \
      "samtools index '${haplotagged_bam}'"
  fi

  local phased_bp phase_ng50
  phased_bp="$(run_docker "${IMG_BASE}" "${outdir}" "python3 /pipeline-scripts/hiphase_get_stat.py '${stats}' basepairs_per_block_sum")"
  phase_ng50="$(run_docker "${IMG_BASE}" "${outdir}" "python3 /pipeline-scripts/hiphase_get_stat.py '${stats}' block_ng50")"

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${outdir}/${phased_small_vcf}" "${outdir}/${phased_small_vcf}.tbi" \
    "${outdir}/${phased_sv_vcf}" "${outdir}/${phased_sv_vcf}.tbi" \
    "${outdir}/${haplotagged_bam}" "${outdir}/${haplotagged_bam}.bai" \
    "${outdir}/${stats}" "${outdir}/${sample_id}.${ref_name}.hiphase.blocks.tsv" \
    "${phased_bp}" "${phase_ng50}"
}

# ===========================================================================
# TRGT: tandem repeat genotyping
# ===========================================================================
# step_trgt <sample_id> <ref_name> <sex> <aligned_bam> <aligned_bam_index> <ref_fasta> <ref_index>
#   <trgt_bed> <expected_male_bed> <expected_female_bed> <outdir>
step_trgt() {
  local sample_id="$1" ref_name="$2" sex="$3" aligned_bam="$4" aligned_bam_index="$5" \
    ref_fasta="$6" ref_index="$7" trgt_bed="$8" expected_male_bed="$9" expected_female_bed="${10}" outdir="${11}"
  mkdir -p "${outdir}"
  link_input "${aligned_bam}" "${outdir}"
  link_input "${aligned_bam_index}" "${outdir}"
  link_input "${ref_fasta}" "${outdir}"
  link_input "${ref_index}" "${outdir}"
  link_input "${trgt_bed}" "${outdir}"

  local karyotype="XX" expected_bed="${expected_female_bed}"
  if [[ "${sex}" == "MALE" ]]; then
    karyotype="XY"
    expected_bed="${expected_male_bed}"
  fi
  link_input "${expected_bed}" "${outdir}"

  local out_prefix="${sample_id}.${ref_name}"
  if resume_available "${outdir}/${out_prefix}.trgt.sorted.vcf.gz" "${outdir}/${out_prefix}.trgt.sorted.vcf.gz.tbi" \
    "${outdir}/${out_prefix}.trgt.spanning.sorted.bam" "${outdir}/${out_prefix}.trgt.spanning.sorted.bam.bai" \
    "${outdir}/${out_prefix}.trgt.dropouts.txt"; then
    log "[resume] trgt genotype: ${out_prefix}.trgt.sorted.vcf.gz already present, skipping"
  else
    run_docker "${IMG_TRGT}" "${outdir}" \
      "trgt --version >&2; trgt genotype --threads ${THREADS} --karyotype '${karyotype}' --genome '$(basename "${ref_fasta}")' --repeats '$(basename "${trgt_bed}")' --reads '$(basename "${aligned_bam}")' --max-depth ${TRGT_MAX_DEPTH} --min-read-quality=${TRGT_MIN_READ_QUALITY} --output-prefix '${out_prefix}.trgt'"

    run_docker "${IMG_TRGT}" "${outdir}" \
      "bcftools sort --output-type z --output '${out_prefix}.trgt.sorted.vcf.gz' --write-index=tbi '${out_prefix}.trgt.vcf.gz'"

    run_docker "${IMG_TRGT}" "${outdir}" \
      "samtools sort --threads ${THREADS} -m 800M --write-index -o '${out_prefix}.trgt.spanning.sorted.bam##idx##${out_prefix}.trgt.spanning.sorted.bam.bai' '${out_prefix}.trgt.spanning.bam'"

    run_docker "${IMG_TRGT}" "${outdir}" \
      "find_trgt_dropouts.py --ploidybed '$(basename "${expected_bed}")' --coverage ${TRGT_HAPLOTYPE_COVERAGE_THRESHOLD} '$(basename "${trgt_bed}")' '${out_prefix}.trgt.spanning.sorted.bam' > '${out_prefix}.trgt.dropouts.txt'"
  fi

  local genotyped uncalled
  genotyped="$(run_docker "${IMG_TRGT}" "${outdir}" \
    "bcftools view --no-header --exclude-uncalled '${out_prefix}.trgt.sorted.vcf.gz' | wc --lines" | tr -d '[:space:]')"
  uncalled="$(run_docker "${IMG_TRGT}" "${outdir}" \
    "bcftools view --no-header --uncalled '${out_prefix}.trgt.sorted.vcf.gz' | wc --lines" | tr -d '[:space:]')"

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${outdir}/${out_prefix}.trgt.sorted.vcf.gz" "${outdir}/${out_prefix}.trgt.sorted.vcf.gz.tbi" \
    "${outdir}/${out_prefix}.trgt.spanning.sorted.bam" "${outdir}/${out_prefix}.trgt.spanning.sorted.bam.bai" \
    "${outdir}/${out_prefix}.trgt.dropouts.txt" "${genotyped}" "${uncalled}"
}

# ===========================================================================
# pbjam bam-stats, bcftools stats/roh, sv_stats, methbat, pbstarphase
# ===========================================================================
# step_pbjam_bam_stats <sample_id> <ref_name> <bam> <bam_index> <outdir>
step_pbjam_bam_stats() {
  local sample_id="$1" ref_name="$2" bam="$3" bam_index="$4" outdir="$5"
  mkdir -p "${outdir}"
  link_input "${bam}" "${outdir}"
  link_input "${bam_index}" "${outdir}"

  local json="${sample_id}.${ref_name}.pbjam.json"
  if resume_available "${outdir}/${json}" "${outdir}/${sample_id}.read_length_distribution.png" "${outdir}/${sample_id}.read_quality_distribution.png"; then
    log "[resume] pbjam bam-stats: ${json} already present, skipping"
  else
    run_docker "${IMG_PBJAM}" "${outdir}" \
      "pbjam bam-stats --threads ${THREADS} --include-unmapped --input-bam '$(basename "${bam}")' --output-json '${json}' --plot-label '${sample_id}.${ref_name}' --image-prefix '${sample_id}.${ref_name}'; mv --verbose '${sample_id}.${ref_name}.read_length_distribution.png' '${sample_id}.read_length_distribution.png'; mv --verbose '${sample_id}.${ref_name}.read_quality_distribution.png' '${sample_id}.read_quality_distribution.png'"
  fi

  local key val results=()
  for key in read_count read_length_mean read_length_median read_length_n50 \
    read_quality_mean read_quality_median mapped_read_count mapped_read_percentage \
    gap_compress_identity_mean gap_compress_identity_median; do
    val="$(run_docker "${IMG_PBJAM}" "${outdir}" "python3 /pipeline-scripts/pbjam_get_value.py '${json}' '${key}'")"
    results+=("${val}")
  done

  printf '%s\t%s\t%s\t%s\t' "${outdir}/${sample_id}.read_length_distribution.png" \
    "${outdir}/${sample_id}.read_quality_distribution.png" \
    "${outdir}/${sample_id}.${ref_name}.mapq_distribution.png" \
    "${outdir}/${sample_id}.${ref_name}.mg_distribution.png"
  ( IFS=$'\t'; printf '%s\n' "${results[*]}" )
}

# step_bcftools_stats_roh <sample_id> <ref_name> <vcf> <ref_fasta> <outdir>
step_bcftools_stats_roh() {
  local sample_id="$1" ref_name="$2" vcf="$3" ref_fasta="$4" outdir="$5"
  mkdir -p "${outdir}"
  link_input "${vcf}" "${outdir}"
  link_input "${ref_fasta}" "${outdir}"

  local vcf_base ref_base stats_txt
  vcf_base="$(basename "${vcf}")"
  ref_base="$(basename "${ref_fasta}")"
  stats_txt="${sample_id}.${ref_name}.small_variants.vcf.stats.txt"

  if resume_available "${outdir}/${stats_txt}" "${outdir}/snv_count.txt" "${outdir}/indel_count.txt" "${outdir}/tstv_ratio.txt" "${outdir}/hethom_ratio.txt"; then
    log "[resume] bcftools stats: ${stats_txt} already present, skipping"
  else
    run_docker "${IMG_BASE}" "${outdir}" \
      "bcftools --version >&2; bcftools norm --fasta-ref '${ref_base}' --multiallelics - '${vcf_base}' 2>/dev/null | bcftools view --apply-filters .,PASS --exclude 'GQ<20.0 || GT=\"ref\" || GT=\"mis\" || ALT=\".\"' --trim-alt-alleles - | bcftools stats --samples '${sample_id}' --fasta-ref '${ref_base}' - > '${stats_txt}'"

    run_docker "${IMG_BASE}" "${outdir}" \
      "grep -w '^SN' '${stats_txt}' | grep 'number of SNPs:' | cut -f4 > snv_count.txt; grep -w '^SN' '${stats_txt}' | grep 'number of indels:' | cut -f4 > indel_count.txt; grep -w '^TSTV' '${stats_txt}' | cut -f5 > tstv_ratio.txt; nHets=\$(grep -w '^PSC' '${stats_txt}' | cut -f6); nNonRefHom=\$(grep -w '^PSC' '${stats_txt}' | cut -f5); printf %.2f \"\$((10**2 * nHets / nNonRefHom))e-2\" > hethom_ratio.txt"
  fi

  local snv_plot="${sample_id}.${ref_name}.small_variants.snv_distribution.png"
  local indel_plot="${sample_id}.${ref_name}.small_variants.indel_distribution.png"
  if ! resume_available "${outdir}/${snv_plot}"; then
    run_docker "${IMG_BASE}" "${outdir}" \
      "grep -w '^ST' '${stats_txt}' | cut -f3,4 | awk -v OFS='\t' 'BEGIN {print \"type\", \"count\"} {print \$1, \$2}' > st.tsv; python3 /pipeline-scripts/plot_snvs.py st.tsv '${sample_id}.${ref_name}' '${snv_plot}'"
  fi
  if ! resume_available "${outdir}/${indel_plot}"; then
    run_docker "${IMG_BASE}" "${outdir}" \
      "grep -w '^IDD' '${stats_txt}' | cut -f3,4 | awk -v OFS='\t' 'BEGIN {print \"length\", \"count\"} {print \$1, \$2}' > idd.tsv; python3 /pipeline-scripts/plot_indels.py idd.tsv '${sample_id}.${ref_name}' '${indel_plot}'"
  fi

  local roh_out="${sample_id}.${ref_name}.bcftools_roh.out"
  local roh_bed="${sample_id}.${ref_name}.bcftools_roh.bed"
  if resume_available "${outdir}/${roh_out}.gz" "${outdir}/${roh_bed}"; then
    log "[resume] bcftools roh: ${roh_bed} already present, skipping"
  else
    run_docker "${IMG_BASE}" "${outdir}" \
      "bcftools roh --threads $(( THREADS > 1 ? THREADS - 1 : 0 )) --AF-dflt 0.4 '${vcf_base}' > '${roh_out}'"
    run_docker "${IMG_BASE}" "${outdir}" \
      "python3 /pipeline-scripts/roh_bed.py '${roh_out}' '${BCFTOOLS_ROH_MIN_LENGTH}' '${BCFTOOLS_ROH_MIN_QUAL}' '${roh_bed}'; gzip -f '${roh_out}'"
  fi

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${outdir}/${stats_txt}" "${outdir}/${roh_out}.gz" "${outdir}/${roh_bed}" \
    "$(cat "${outdir}/snv_count.txt")" "$(cat "${outdir}/indel_count.txt")" \
    "$(cat "${outdir}/tstv_ratio.txt")" "$(cat "${outdir}/hethom_ratio.txt")" \
    "${outdir}/${snv_plot}"
  printf '%s\n' "${outdir}/${indel_plot}"
}

# step_sv_stats <sample_id> <ref_name> <vcf> <outdir>
step_sv_stats() {
  local sample_id="$1" ref_name="$2" vcf="$3" outdir="$4"
  mkdir -p "${outdir}"
  link_input "${vcf}" "${outdir}"
  local vcf_base
  vcf_base="$(basename "${vcf}")"

  if resume_available "${outdir}/stat_DUP.txt" "${outdir}/stat_DEL.txt" "${outdir}/stat_INS.txt" \
    "${outdir}/stat_INV.txt" "${outdir}/stat_BND.txt" "${outdir}/stat_SWAP.txt"; then
    log "[resume] sv_stats: stats already present, skipping"
  else
    run_docker "${IMG_BASE}" "${outdir}" "bcftools --version
bcftools query -f '%INFO/SVLEN\n' --include '(GT!=\"ref\" & GT!=\"./.\" & GT!=\".\") & FILTER=\"PASS\" & ABS(SVLEN)>=${SV_STATS_MIN_LENGTH} & SVTYPE=\"DUP\"' '${vcf_base}' > length_DUP.txt || : > length_DUP.txt
wc --lines < length_DUP.txt > stat_DUP.txt
bcftools query -f '%INFO/SVLEN\n' --include 'GT=\"alt\" & FILTER=\"PASS\" & ABS(SVLEN)>=${SV_STATS_MIN_LENGTH} & SVTYPE=\"DEL\" & STRLEN(ALT[0])<=${SV_STATS_MAX_SCAR_LENGTH}' '${vcf_base}' > length_DEL.txt || : > length_DEL.txt
wc --lines < length_DEL.txt > stat_DEL.txt
bcftools query -f '%INFO/SVLEN\n' --include 'GT=\"alt\" & FILTER=\"PASS\" & ABS(SVLEN)>=${SV_STATS_MIN_LENGTH} & SVTYPE=\"INS\" & STRLEN(REF)<=${SV_STATS_MAX_SCAR_LENGTH}' '${vcf_base}' > length_INS.txt || : > length_INS.txt
wc --lines < length_INS.txt > stat_INS.txt
bcftools view --no-header --include 'GT=\"alt\" & FILTER=\"PASS\" & ABS(SVLEN)>=${SV_STATS_MIN_LENGTH} & (SVTYPE=\"INS\" | SVTYPE=\"DEL\") & STRLEN(REF)>${SV_STATS_MAX_SCAR_LENGTH} & STRLEN(ALT[0])>${SV_STATS_MAX_SCAR_LENGTH}' '${vcf_base}' | wc --lines > stat_SWAP.txt || echo 0 > stat_SWAP.txt
bcftools query -f '%INFO/SVLEN\n' --include 'GT=\"alt\" & FILTER=\"PASS\" & SVTYPE=\"INV\"' '${vcf_base}' > length_INV.txt || : > length_INV.txt
wc --lines < length_INV.txt > stat_INV.txt
bcftools view --no-header --include 'GT=\"alt\" & FILTER=\"PASS\" & SVTYPE=\"BND\"' '${vcf_base}' | wc --lines > stat_BND.txt || echo 0 > stat_BND.txt"
  fi

  local plot="${sample_id}.${ref_name}.sv_stats.png"
  if ! resume_available "${outdir}/${plot}"; then
    run_docker "${IMG_BASE}" "${outdir}" \
      "python3 /pipeline-scripts/plot_sv_stats.py length_INS.txt length_DEL.txt length_DUP.txt length_INV.txt ${SV_STATS_MIN_LENGTH} '${sample_id}.${ref_name}' '${plot}'"
  fi

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$(cat "${outdir}/stat_DUP.txt")" "$(cat "${outdir}/stat_DEL.txt")" "$(cat "${outdir}/stat_INS.txt")" \
    "$(cat "${outdir}/stat_INV.txt")" "$(cat "${outdir}/stat_BND.txt")" "$(cat "${outdir}/stat_SWAP.txt")" \
    "${outdir}/${plot}"
}

# step_methbat_pileup <haplotagged_bam> <haplotagged_bam_index> <out_prefix> <outdir>
step_methbat_pileup() {
  local haplotagged_bam="$1" haplotagged_bam_index="$2" out_prefix="$3" outdir="$4"
  mkdir -p "${outdir}"
  link_input "${haplotagged_bam}" "${outdir}"
  link_input "${haplotagged_bam_index}" "${outdir}"

  local skip5mc="" skip5hmc="" skip6ma=""
  [[ "${METHBAT_SKIP_5MC}" == "true" ]] && skip5mc="--skip-5mC"
  [[ "${METHBAT_SKIP_5HMC}" == "true" ]] && skip5hmc="--skip-5hmC"
  [[ "${METHBAT_SKIP_6MA}" == "true" ]] && skip6ma="--skip-6mA"

  if resume_available "${outdir}/messages.txt" "${outdir}/${out_prefix}.combined.bed.count" "${outdir}/${out_prefix}.hap1.bed.count" "${outdir}/${out_prefix}.hap2.bed.count"; then
    log "[resume] methbat pileup: ${out_prefix} already present, skipping"
  else
    run_docker "${IMG_METHBAT}" "${outdir}" \
      "touch messages.txt; methbat --version >&2; methbat pileup --threads ${THREADS} --input-bam '$(basename "${haplotagged_bam}")' --min-mapq ${METHBAT_MIN_MAPQ} --min-coverage ${METHBAT_MIN_COVERAGE} --edge-trimming-size ${METHBAT_EDGE_TRIMMING_SIZE} --phase-set-min-fraction ${METHBAT_PHASE_SET_MIN_FRACTION} ${skip5mc} ${skip5hmc} ${skip6ma} --output-prefix '${out_prefix}' || echo 'MethBat pileup failed' >> messages.txt; echo 0 > '${out_prefix}.combined.bed.count'; echo 0 > '${out_prefix}.hap1.bed.count'; echo 0 > '${out_prefix}.hap2.bed.count'; if [ -f '${out_prefix}.5mC.bed.gz' ]; then zgrep -v '^#' '${out_prefix}.5mC.bed.gz' | grep -c Total > '${out_prefix}.combined.bed.count' || true; zgrep -v '^#' '${out_prefix}.5mC.bed.gz' | grep -c hap1 > '${out_prefix}.hap1.bed.count' || true; zgrep -v '^#' '${out_prefix}.5mC.bed.gz' | grep -c hap2 > '${out_prefix}.hap2.bed.count' || true; fi"
  fi

  printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${outdir}/${out_prefix}.5mC.bed.gz" "${outdir}/${out_prefix}.5mC.bed.gz.tbi" \
    "${outdir}/${out_prefix}.5hmC.bed.gz" "${outdir}/${out_prefix}.5hmC.bed.gz.tbi" \
    "$(cat "${outdir}/${out_prefix}.hap1.bed.count")" "$(cat "${outdir}/${out_prefix}.hap2.bed.count")"
  printf '%s\n' "$(cat "${outdir}/${out_prefix}.combined.bed.count")"
}

# step_methbat_profile <cpg_pileup_bed> <region_tsv> <out_prefix> <outdir>
step_methbat_profile() {
  local cpg_pileup_bed="$1" region_tsv="$2" out_prefix="$3" outdir="$4"
  mkdir -p "${outdir}"
  [[ -f "${cpg_pileup_bed}" ]] || { printf '\t0\t0\t0\n'; return 0; }
  link_input "${cpg_pileup_bed}" "${outdir}"
  link_input "${region_tsv}" "${outdir}"

  local profile="${out_prefix}.methbat.profile.tsv"
  if resume_available "${outdir}/${profile}" "${outdir}/methylated_count.txt" "${outdir}/unmethylated_count.txt" "${outdir}/asm_count.txt"; then
    log "[resume] methbat profile: ${profile} already present, skipping"
  else
    run_docker "${IMG_METHBAT}" "${outdir}" \
      "methbat --version >&2; methbat profile --input-pileup '$(basename "${cpg_pileup_bed}")' --input-regions '$(basename "${region_tsv}")' --output-region-profile '${profile}'; awk '\$5==\"Methylated\" {print}' '${profile}' | wc -l > methylated_count.txt; awk '\$5==\"Unmethylated\" {print}' '${profile}' | wc -l > unmethylated_count.txt; awk '\$5==\"AlleleSpecificMethylation\" {print}' '${profile}' | wc -l > asm_count.txt"
  fi

  printf '%s\t%s\t%s\t%s\n' "${outdir}/${profile}" \
    "$(cat "${outdir}/methylated_count.txt")" "$(cat "${outdir}/unmethylated_count.txt")" \
    "$(cat "${outdir}/asm_count.txt")"
}

# step_pbstarphase <out_prefix> <phased_small_vcf> <phased_small_vcf_index> <phased_sv_vcf> <phased_sv_vcf_index>
#   <aligned_bam> <aligned_bam_index> <ref_fasta> <ref_index> <outdir>
step_pbstarphase() {
  local out_prefix="$1" phased_small_vcf="$2" phased_small_vcf_index="$3" phased_sv_vcf="$4" phased_sv_vcf_index="$5" \
    aligned_bam="$6" aligned_bam_index="$7" ref_fasta="$8" ref_index="$9" outdir="${10}"
  mkdir -p "${outdir}"
  link_input "${phased_small_vcf}" "${outdir}"
  link_input "${phased_small_vcf_index}" "${outdir}"
  link_input "${phased_sv_vcf}" "${outdir}"
  link_input "${phased_sv_vcf_index}" "${outdir}"
  link_input "${aligned_bam}" "${outdir}"
  link_input "${aligned_bam_index}" "${outdir}"
  link_input "${ref_fasta}" "${outdir}"
  link_input "${ref_index}" "${outdir}"

  if resume_available "${outdir}/${out_prefix}.pbstarphase.json" "${outdir}/${out_prefix}.pbstarphase.tsv"; then
    log "[resume] pbstarphase: ${out_prefix}.pbstarphase.json already present, skipping"
  else
    run_docker "${IMG_PBSTARPHASE}" "${outdir}" \
      "pbstarphase --version >&2; pbstarphase diplotype --database /opt/pbstarphase_db.json.gz --reference '$(basename "${ref_fasta}")' --vcf '$(basename "${phased_small_vcf}")' --sv-vcf '$(basename "${phased_sv_vcf}")' --bam '$(basename "${aligned_bam}")' --output-calls '${out_prefix}.pbstarphase.json' --pharmcat-tsv '${out_prefix}.pbstarphase.tsv'"
  fi

  printf '%s\t%s\n' "${outdir}/${out_prefix}.pbstarphase.json" "${outdir}/${out_prefix}.pbstarphase.tsv"
}

# ===========================================================================
# GLnexus + split-by-sample (family/joint mode only)
# ===========================================================================
# step_glnexus <cohort_id> <ref_name> <outdir> <gvcf...>
step_glnexus() {
  local cohort_id="$1" ref_name="$2" outdir="$3"
  shift 3
  local gvcfs=("$@")
  mkdir -p "${outdir}"

  local vcf="${cohort_id}.${ref_name}.small_variants.vcf.gz"

  if resume_available "${outdir}/${vcf}" "${outdir}/${vcf}.tbi"; then
    log "[resume] GLnexus: ${vcf} already present, skipping"
  else
    local g gvcf_list=""
    for g in "${gvcfs[@]}"; do
      link_input "${g}" "${outdir}"
      link_input "${g}.tbi" "${outdir}"
      gvcf_list+=" $(basename "${g}")"
    done

    cat > "${outdir}/config.yml" <<'YAML'
unifier_config:
  min_AQ1: 0
  min_AQ2: 0
  min_GQ: 0
  monoallelic_sites_for_lost_alleles: true
  max_alleles_per_site: 32
genotyper_config:
  required_dp: 1
  revise_genotypes: false
  allow_partial_data: true
  more_PL: true
  trim_uncalled_alleles: true
  liftover_fields:
    - orig_names: [MIN_DP, DP]
      name: DP
      description: '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">'
      type: int
      combi_method: min
      number: basic
      count: 1
      ignore_non_variants: true
    - orig_names: [AD]
      name: AD
      description: '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">'
      type: int
      number: alleles
      combi_method: min
      default_type: zero
      count: 0
    - orig_names: [GQ]
      name: GQ
      description: '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">'
      type: int
      number: basic
      combi_method: min
      count: 1
      ignore_non_variants: true
    - orig_names: [PL]
      name: PL
      description: '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype Likelihoods">'
      type: int
      number: genotype
      combi_method: missing
      count: 0
      ignore_non_variants: true
YAML

    local bcf="${cohort_id}.${ref_name}.small_variants.bcf"
    run_docker "${IMG_GLNEXUS}" "${outdir}" \
      "glnexus_cli --help 2>&1 | grep -Eo 'glnexus_cli release v[0-9a-f.-]+' >&2; bcftools --version >&2; glnexus_cli --threads ${THREADS} --mem-gbytes ${GLNEXUS_MEM_GB} --dir '${cohort_id}.${ref_name}.GLnexus.DB' --config ./config.yml${gvcf_list} > '${bcf}'; bcftools view --threads $(( THREADS > 1 ? THREADS - 1 : 0 )) --output-type z --output-file '${vcf}' '${bcf}'; bcftools index --threads $(( THREADS > 1 ? THREADS - 1 : 0 )) --tbi '${vcf}'; rm --recursive --force '${cohort_id}.${ref_name}.GLnexus.DB' '${bcf}'"
  fi

  printf '%s\t%s\n' "${outdir}/${vcf}" "${outdir}/${vcf}.tbi"
}

# step_split_vcf_by_sample <vcf> <vcf_index> <exclude_uncalled: true|false> <outdir> <sample_id...>
# Writes "<sample_id>.<vcf_basename>" per sample into outdir.
step_split_vcf_by_sample() {
  local vcf="$1" vcf_index="$2" exclude_uncalled="$3" outdir="$4"
  shift 4
  local sample_ids=("$@")
  mkdir -p "${outdir}"

  local vcf_base
  vcf_base="$(basename "${vcf}")"
  local exclude_flag=""
  [[ "${exclude_uncalled}" == "true" ]] && exclude_flag="--exclude-uncalled"

  local sid all_resumed=true
  for sid in "${sample_ids[@]}"; do
    resume_available "${outdir}/${sid}.${vcf_base}" "${outdir}/${sid}.${vcf_base}.tbi" || all_resumed=false
  done
  if [[ "${all_resumed}" != "true" ]]; then
    link_input "${vcf}" "${outdir}"
    link_input "${vcf_index}" "${outdir}"
  fi

  for sid in "${sample_ids[@]}"; do
    if resume_available "${outdir}/${sid}.${vcf_base}" "${outdir}/${sid}.${vcf_base}.tbi"; then
      log "[resume] split_vcf_by_sample: ${sid}.${vcf_base} already present, skipping"
      continue
    fi
    run_docker "${IMG_BASE}" "${outdir}" \
      "bcftools --version >&2; bcftools view --threads $(( THREADS > 1 ? THREADS - 1 : 0 )) --samples '${sid}' ${exclude_flag} --output-type z --output '${sid}.${vcf_base}' --write-index=tbi '${vcf_base}'"
  done
}
