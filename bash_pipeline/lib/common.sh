#!/usr/bin/env bash
# Shared helper functions used by every step script.
set -euo pipefail

log()  { printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2; }
die()  { log "ERROR: $*"; exit 1; }

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "Required command '$1' not found on PATH."
}

require_file() {
  [[ -f "$1" ]] || die "Required file not found: $1"
}

require_dir() {
  [[ -d "$1" ]] || die "Required directory not found: $1"
}

# Every input/output file used by the pipeline must live under DATA_ROOT.
# We bind-mount DATA_ROOT to the *same* absolute path inside every container,
# so tools can be invoked with plain host-absolute paths (mirroring the way
# the WDL tasks symlink localized files into their own cwd and reference them
# by path) without juggling a different mount per file.
assert_under_data_root() {
  local path
  path="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
  case "${path}" in
    "${DATA_ROOT}"/*|"${DATA_ROOT}") ;;
    *) die "Path '$1' is not under DATA_ROOT='${DATA_ROOT}'. Re-run with --data-root set to a common ancestor of all inputs/outputs." ;;
  esac
}

# run_docker <image> <workdir> <command>
# Runs a single command string inside <image>, with DATA_ROOT mounted 1:1 and
# <workdir> (which must be under DATA_ROOT) as the container's cwd. The
# pipeline's py/ helper scripts (used to reproduce the WDL tasks' inline
# plotting/parsing logic) are mounted read-only at /pipeline-scripts.
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

run_docker() {
  local image="$1" workdir="$2" cmd="$3"
  assert_under_data_root "${workdir}"
  mkdir -p "${workdir}"
  docker run --rm \
    --user "$(id -u):$(id -g)" \
    -v "${DATA_ROOT}:${DATA_ROOT}" \
    -v "${PIPELINE_DIR}/py:/pipeline-scripts:ro" \
    -w "${workdir}" \
    "${image}" \
    bash -o pipefail -c "${cmd}"
}

run_docker_gpu() {
  local image="$1" workdir="$2" cmd="$3" gpu_docker_args="${4:---gpus all}"
  assert_under_data_root "${workdir}"
  mkdir -p "${workdir}"
  # shellcheck disable=SC2086 # gpu_docker_args is intentionally word-split (e.g. "-e FOO=bar")
  docker run --rm ${gpu_docker_args} \
    --user "$(id -u):$(id -g)" \
    -v "${DATA_ROOT}:${DATA_ROOT}" \
    -v "${PIPELINE_DIR}/py:/pipeline-scripts:ro" \
    -w "${workdir}" \
    "${image}" \
    bash -o pipefail -c "${cmd}"
}

# link_input <target-file> <link-dir>
# Symlinks a file into a step's working directory using its basename, same
# idiom the WDL tasks use ("ln --symbolic --verbose ... ."), so every step's
# command line can refer to inputs by basename regardless of their original
# location on DATA_ROOT.
link_input() {
  local target="$1" dir="$2"
  mkdir -p "${dir}"
  # No --verbose: this function's stdout must stay silent because callers
  # like step_pbmm2_align/step_mosdepth/etc. are invoked as "$(step_xxx ...)",
  # and any extra stdout here would corrupt their tab-separated return value.
  ln --symbolic --force "$(cd "$(dirname "${target}")" && pwd)/$(basename "${target}")" "${dir}/"
}
