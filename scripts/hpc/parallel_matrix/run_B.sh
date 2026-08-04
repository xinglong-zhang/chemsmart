#!/usr/bin/env bash
# BatchJob parallel matrix: B1–B7 (pKa CSV, 5 rows).
#
# Usage (from sandbox):
#   DRY_RUN=1 ./run_B.sh
#   CASES=B1,B6 DRY_RUN=1 ./run_B.sh
#   DRY_RUN=0 WAIT=1 ./run_B.sh

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=common.sh
source "${SCRIPT_DIR}/common.sh"

ensure_chemsmart
sub_flags
rerun_flags
mkdir -p "${MATRIX_ROOT}/logs"
SUMMARY="${MATRIX_ROOT}/logs/B_summary_${TIMESTAMP}.tsv"
CASES="${CASES:-B1,B2,B3,B4,B5,B6,B7}"

NUM_ROWS="$(
  python3 - "${PKA_DIR}/${PKA_TABLE}" <<'PY'
import csv, sys
with open(sys.argv[1], newline="") as fh:
    print(len(list(csv.DictReader(fh))))
PY
)"

run_pka_batch() {
  local case_id="$1"
  shift
  local run_dir logfile
  run_dir="$(prepare_case_dir "${case_id}")"
  logfile="${MATRIX_ROOT}/logs/${case_id}_${TIMESTAMP}.log"
  log "${case_id}: run_dir=${run_dir}"
  (
    cd "${run_dir}"
    chemsmart sub \
      "${SUB_FLAGS[@]}" \
      -s "${SERVER}" \
      "$@" \
      gaussian \
      -p "${GAUSSIAN_PROJECT_BATCH}" \
      -f "${PKA_TABLE}" \
      pka \
      --scheme "${SCHEME}" \
      batch \
      "${RERUN_FLAGS[@]}" \
      >"${logfile}" 2>&1
  )
  # Keep a short console breadcrumb without polluting command substitution.
  tail -n 5 "${logfile}" >&2 || true
  printf '%s\n' "${run_dir}"
}

check_b_case() {
  local case_id="$1"
  local expected_re="$2"
  local run_dir="$3"
  local _old
  _old="$(pwd)"
  cd "${run_dir}"
  expect_array "${case_id}" "${expected_re}" || true
  expect_shared_runscript "${case_id}.runscript" || true
  cd "${_old}"
}

run_case() {
  local case_id="$1"
  local run_dir line _old
  case "${case_id}" in
    B1)
      run_dir="$(run_pka_batch B1 --no-run-in-parallel)"
      check_b_case B1 '^#SBATCH --array=1-'"${NUM_ROWS}"'%1$' "${run_dir}"
      ;;
    B2)
      # -M should be ignored under serial default
      run_dir="$(run_pka_batch B2 --no-run-in-parallel -M "${NUM_ROWS}")"
      check_b_case B2 '^#SBATCH --array=1-'"${NUM_ROWS}"'%1$' "${run_dir}"
      ;;
    B3)
      run_dir="$(run_pka_batch B3 --run-in-parallel)"
      _old="$(pwd)"
      cd "${run_dir}"
      line="$(extract_array_line "$(find_newest 'chemsmart_sub_array_*.sh')")"
      if [[ "${line}" == "#SBATCH --array=1-${NUM_ROWS}" ]] \
        || [[ "${line}" =~ ^#SBATCH\ --array=1-${NUM_ROWS}%[0-9]+$ ]]; then
        record PASS B3 "${line}"
      else
        record FAIL B3 "unexpected array line: ${line}"
      fi
      expect_shared_runscript B3.runscript || true
      cd "${_old}"
      ;;
    B4)
      run_dir="$(run_pka_batch B4 --run-in-parallel -M 1)"
      check_b_case B4 '^#SBATCH --array=1-'"${NUM_ROWS}"'%1$' "${run_dir}"
      ;;
    B5)
      run_dir="$(run_pka_batch B5 --run-in-parallel -M 2)"
      check_b_case B5 '^#SBATCH --array=1-'"${NUM_ROWS}"'%2$' "${run_dir}"
      ;;
    B6)
      run_dir="$(run_pka_batch B6 --run-in-parallel -M "${NUM_ROWS}")"
      check_b_case B6 '^#SBATCH --array=1-'"${NUM_ROWS}"'%'"${NUM_ROWS}"'$' "${run_dir}"
      ;;
    B7)
      run_dir="$(run_pka_batch B7 --run-in-parallel -M $((NUM_ROWS * 2)))"
      check_b_case B7 '^#SBATCH --array=1-'"${NUM_ROWS}"'%'"${NUM_ROWS}"'$' "${run_dir}"
      ;;
    *)
      record SKIP "${case_id}" "unknown B-case"
      ;;
  esac
}

log "B-series: rows=${NUM_ROWS} cases=${CASES} DRY_RUN=${DRY_RUN} WAIT=${WAIT}"
IFS=',' read -r -a CASE_LIST <<< "${CASES}"
for case_id in "${CASE_LIST[@]}"; do
  case_id="$(echo "${case_id}" | tr -d '[:space:]')"
  [[ -n "${case_id}" ]] || continue
  run_case "${case_id}"
  if [[ "${DRY_RUN}" == "0" && "${WAIT}" == "1" ]]; then
    wait_for_job_ids_from_log \
      "${MATRIX_ROOT}/logs/${case_id}_${TIMESTAMP}.log"
  fi
done
print_summary "${SUMMARY}"
[[ "${FAIL}" -eq 0 ]]
