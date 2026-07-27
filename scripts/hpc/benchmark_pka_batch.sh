#!/usr/bin/env bash
# Benchmark pKa batch array submission vs individual submissions on HPC.
#
# Uses inputs from $PROJECT/chemsmart/test_pka/gaussian (pka_scale.csv + *.xyz).
# Run from sandbox (login node kills chemsmart).
#
# Usage:
#   ./benchmark_pka_batch.sh                     # all modes, wait for completion
#   MODE=array_parallel ./benchmark_pka_batch.sh
#   DRY_RUN=1 ./benchmark_pka_batch.sh           # --test only, no sbatch
#   WAIT=0 ./benchmark_pka_batch.sh                # submit timing only
#
# Modes (MODE):
#   individual     one chemsmart sub per CSV row (legacy-style)
#   array_serial   one SLURM array, -M 1 (one task at a time)
#   array_parallel one SLURM array, -M NUM_JOBS (all rows concurrently)
#   all            run the three modes above sequentially (default)
#
# Environment:
#   PROJECT          default: /project/xlzhang/jingyi
#   PKA_DIR          default: $PROJECT/chemsmart/test_pka/gaussian
#   RUN_ROOT         default: $PROJECT/chemsmart/test_batch/pka_benchmark
#   SERVER           default: cu_dev
#   GAUSSIAN_PROJECT default: pka_wat
#   PKA_TABLE        default: pka_scale.csv
#   SCHEME           default: direct
#   ARRAY_PARALLEL_M default: 0 (= number of CSV data rows)
#   POLL_SEC         default: 30
#   DRY_RUN          default: 0  (set 1 to pass --test)
#   WAIT             default: 1   (set 0 to skip queue wait)
#   FORCE_RERUN      default: 1   (pass -R / --no-skip-completed)

set -euo pipefail

PROJECT="${PROJECT:-/project/xlzhang/jingyi}"
PKA_DIR="${PKA_DIR:-${PROJECT}/chemsmart/test_pka/gaussian}"
RUN_ROOT="${RUN_ROOT:-${PROJECT}/chemsmart/test_batch/pka_benchmark}"
SERVER="${SERVER:-cu_dev}"
GAUSSIAN_PROJECT="${GAUSSIAN_PROJECT:-pka_wat}"
PKA_TABLE="${PKA_TABLE:-pka_scale.csv}"
SCHEME="${SCHEME:-direct}"
MODE="${MODE:-all}"
POLL_SEC="${POLL_SEC:-30}"
DRY_RUN="${DRY_RUN:-0}"
WAIT="${WAIT:-1}"
FORCE_RERUN="${FORCE_RERUN:-1}"
ARRAY_PARALLEL_M="${ARRAY_PARALLEL_M:-0}"

TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
SUMMARY_FILE="${RUN_ROOT}/benchmark_summary_${TIMESTAMP}.tsv"
LOG_DIR="${RUN_ROOT}/logs_${TIMESTAMP}"

mkdir -p "${RUN_ROOT}" "${LOG_DIR}"

if ! command -v chemsmart >/dev/null 2>&1; then
  echo "error: chemsmart not found on PATH" >&2
  exit 1
fi

if [[ ! -f "${PKA_DIR}/${PKA_TABLE}" ]]; then
  echo "error: missing ${PKA_DIR}/${PKA_TABLE}" >&2
  exit 1
fi

export CHEMSMART_CONFIG_DIR="${CHEMSMART_CONFIG_DIR:-${HOME}/.chemsmart}"

NUM_JOBS="$(
  python3 - "${PKA_DIR}/${PKA_TABLE}" <<'PY'
import csv, sys
with open(sys.argv[1], newline="") as fh:
    print(len(list(csv.DictReader(fh))))
PY
)"

if [[ "${ARRAY_PARALLEL_M}" -eq 0 ]]; then
  ARRAY_PARALLEL_M="${NUM_JOBS}"
fi

SUB_TEST_FLAG=()
if [[ "${DRY_RUN}" == "1" ]]; then
  SUB_TEST_FLAG=(--test)
fi

RERUN_FLAG=()
if [[ "${FORCE_RERUN}" == "1" ]]; then
  RERUN_FLAG=(-R)
fi

log() {
  printf '[%s] %s\n' "$(date '+%F %T')" "$*" | tee -a "${LOG_DIR}/benchmark.log"
}

elapsed_since() {
  local start="$1"
  date +%s | awk -v s="${start}" '{printf "%.1f", $1 - s}'
}

extract_submitted_job_ids() {
  local logfile="$1"
  grep -Eo 'Submitted batch job [0-9]+' "${logfile}" | awk '{print $NF}' | sort -u
}

wait_for_jobs() {
  local -a ids=("$@")
  if [[ "${#ids[@]}" -eq 0 ]]; then
    log "no job ids to wait for"
    return 0
  fi
  log "waiting for job(s): ${ids[*]}"
  while true; do
    local pending=0
    for id in "${ids[@]}"; do
      if squeue -j "${id}" -h 2>/dev/null | grep -q .; then
        pending=1
        break
      fi
    done
    if [[ "${pending}" -eq 0 ]]; then
      break
    fi
    sleep "${POLL_SEC}"
  done
}

max_sacct_elapsed_seconds() {
  local -a ids=("$@")
  local max_elapsed=0
  local id elapsed
  for id in "${ids[@]}"; do
    elapsed="$(
      sacct -j "${id}" --parsable2 --noheader \
        --format=JobID,Elapsed \
        2>/dev/null \
        | awk -F'|' '
            $1 ~ /^[0-9]+(_[0-9]+)?$/ && $2 != "" {
              split($2, p, ":")
              if (length(p) == 2) sec = p[1]*60 + p[2]
              else if (length(p) == 3) sec = p[1]*3600 + p[2]*60 + p[3]
              else sec = 0
              if (sec > max) max = sec
            }
            END { print max+0 }'
    )"
    if [[ "${elapsed}" -gt "${max_elapsed}" ]]; then
      max_elapsed="${elapsed}"
    fi
  done
  echo "${max_elapsed}"
}

prepare_run_dir() {
  local mode="$1"
  local run_dir="${RUN_ROOT}/${mode}_${TIMESTAMP}"
  mkdir -p "${run_dir}"
  ln -sf "${PKA_DIR}/${PKA_TABLE}" "${run_dir}/${PKA_TABLE}"
  for xyz in "${PKA_DIR}"/*.xyz; do
    [[ -f "${xyz}" ]] || continue
    ln -sf "${xyz}" "${run_dir}/$(basename "${xyz}")"
  done
  printf '%s\n' "${run_dir}"
}

write_summary_row() {
  printf '%s\n' "$1" >> "${SUMMARY_FILE}"
}

run_individual() {
  local run_dir submit_start submit_secs wall_secs
  local logfile="${LOG_DIR}/individual_${TIMESTAMP}.log"
  run_dir="$(prepare_run_dir individual)"
  log "mode=individual run_dir=${run_dir}"

  cd "${run_dir}"
  submit_start="$(date +%s)"
  : > "${logfile}"

  while IFS=$'\t' read -r filepath proton_index charge multiplicity; do
    log "individual submit: ${filepath} pi=${proton_index}"
    chemsmart sub \
      "${SUB_TEST_FLAG[@]}" \
      -s "${SERVER}" \
      gaussian \
      -p "${GAUSSIAN_PROJECT}" \
      -f "${filepath}" \
      --charge "${charge}" \
      --multiplicity "${multiplicity}" \
      pka \
      --scheme "${SCHEME}" \
      submit \
      -pi "${proton_index}" \
      "${RERUN_FLAG[@]}" \
      2>&1 | tee -a "${logfile}"
  done < <(
    python3 - "${PKA_TABLE}" <<'PY'
import csv, sys

with open(sys.argv[1], newline="") as fh:
    for row in csv.DictReader(fh):
        print(
            "\t".join(
                [
                    row["filepath"].strip(),
                    row["proton_index"].strip(),
                    row["charge"].strip(),
                    row["multiplicity"].strip(),
                ]
            )
        )
PY
  )

  submit_secs="$(elapsed_since "${submit_start}")"
  mapfile -t job_ids < <(extract_submitted_job_ids "${logfile}")

  wall_secs="${submit_secs}"
  if [[ "${WAIT}" == "1" && "${DRY_RUN}" != "1" && "${#job_ids[@]}" -gt 0 ]]; then
    wait_for_jobs "${job_ids[@]}"
    wall_secs="$(elapsed_since "${submit_start}")"
  fi

  compute_secs="NA"
  if [[ "${DRY_RUN}" != "1" && "${#job_ids[@]}" -gt 0 ]]; then
    compute_secs="$(max_sacct_elapsed_seconds "${job_ids[@]}")"
  fi

  write_summary_row "$(printf '%s\t%s\t%s\t%s\t%s\t%s\t%s' \
    individual NA "${#job_ids[@]}" "${submit_secs}" "${wall_secs}" \
    "${compute_secs}" "${run_dir}")"
  log "individual done: jobs=${#job_ids[@]} submit_s=${submit_secs} wall_s=${wall_secs}"
}

run_array() {
  local throttle="$1"
  local mode="$2"
  local run_dir submit_start submit_secs wall_secs
  local logfile="${LOG_DIR}/${mode}_${TIMESTAMP}.log"
  run_dir="$(prepare_run_dir "${mode}")"
  log "mode=${mode} throttle=-M ${throttle} run_dir=${run_dir}"

  cd "${run_dir}"
  submit_start="$(date +%s)"

  chemsmart sub \
    "${SUB_TEST_FLAG[@]}" \
    -s "${SERVER}" \
    --no-run-in-parallel \
    -M "${throttle}" \
    gaussian \
    -p "${GAUSSIAN_PROJECT}" \
    -f "${PKA_TABLE}" \
    pka \
    --scheme "${SCHEME}" \
    batch \
    "${RERUN_FLAG[@]}" \
    2>&1 | tee "${logfile}"

  submit_secs="$(elapsed_since "${submit_start}")"
  mapfile -t job_ids < <(extract_submitted_job_ids "${logfile}")

  wall_secs="${submit_secs}"
  if [[ "${WAIT}" == "1" && "${DRY_RUN}" != "1" && "${#job_ids[@]}" -gt 0 ]]; then
    wait_for_jobs "${job_ids[@]}"
    wall_secs="$(elapsed_since "${submit_start}")"
  fi

  compute_secs="NA"
  if [[ "${DRY_RUN}" != "1" && "${#job_ids[@]}" -gt 0 ]]; then
    compute_secs="$(max_sacct_elapsed_seconds "${job_ids[@]}")"
  fi

  write_summary_row "$(printf '%s\t%s\t%s\t%s\t%s\t%s\t%s' \
    "${mode}" "${throttle}" "${#job_ids[@]}" "${submit_secs}" "${wall_secs}" \
    "${compute_secs}" "${run_dir}")"
  log "${mode} done: array_jobs=${#job_ids[@]} submit_s=${submit_secs} wall_s=${wall_secs}"
}

init_summary() {
  cat > "${SUMMARY_FILE}" <<EOF
# pKa batch benchmark ${TIMESTAMP}
# PKA_DIR=${PKA_DIR}  table=${PKA_TABLE}  rows=${NUM_JOBS}
# SERVER=${SERVER}  project=${GAUSSIAN_PROJECT}  scheme=${SCHEME}
# DRY_RUN=${DRY_RUN}  WAIT=${WAIT}
mode	throttle	num_sbatch_jobs	submit_s	wall_s	compute_s	run_dir
EOF
}

print_footer() {
  cat <<EOF

Speedup vs array_serial (wall clock):
  speedup = wall_s(array_serial) / wall_s(mode)

Quick commands (run inside ${PKA_DIR}):

  # Batch array (parallel, ${NUM_JOBS} tasks at once)
  chemsmart sub -s ${SERVER} --no-run-in-parallel -M ${ARRAY_PARALLEL_M} \\
    gaussian -p ${GAUSSIAN_PROJECT} -f ${PKA_TABLE} \\
    pka --scheme ${SCHEME} batch -R

  # Batch array (serial throttle)
  chemsmart sub -s ${SERVER} --no-run-in-parallel -M 1 \\
    gaussian -p ${GAUSSIAN_PROJECT} -f ${PKA_TABLE} \\
    pka --scheme ${SCHEME} batch -R

  # Single row (example: 1a.xyz)
  chemsmart sub -s ${SERVER} gaussian -p ${GAUSSIAN_PROJECT} -f 1a.xyz \\
    --charge 0 --multiplicity 1 pka --scheme ${SCHEME} submit -pi 10 -R
EOF
}

main() {
  log "starting pKa batch benchmark (MODE=${MODE}, rows=${NUM_JOBS})"
  init_summary
  case "${MODE}" in
    all)
      run_individual
      run_array 1 array_serial
      run_array "${ARRAY_PARALLEL_M}" array_parallel
      ;;
    individual)
      run_individual
      ;;
    array_serial)
      run_array 1 array_serial
      ;;
    array_parallel)
      run_array "${ARRAY_PARALLEL_M}" array_parallel
      ;;
    *)
      echo "MODE must be individual, array_serial, array_parallel, or all" >&2
      exit 2
      ;;
  esac
  print_footer | tee -a "${SUMMARY_FILE}"
  log "summary written to ${SUMMARY_FILE}"
}

main "$@"
