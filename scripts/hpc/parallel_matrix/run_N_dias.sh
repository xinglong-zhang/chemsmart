#!/usr/bin/env bash
# Nestable DIAS parallel matrix: N1–N4.
#
# Uses TS-mode DIAS on an IRC/opt log (3 children: mol + f1 + f2).
#
# Usage (from sandbox):
#   DRY_RUN=1 ./run_N_dias.sh
#   CASES=N1,N2 DRY_RUN=1 ./run_N_dias.sh

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=common.sh
source "${SCRIPT_DIR}/common.sh"

ensure_chemsmart
sub_flags
rerun_flags
mkdir -p "${MATRIX_ROOT}/logs"
SUMMARY="${MATRIX_ROOT}/logs/N_dias_summary_${TIMESTAMP}.tsv"
CASES="${CASES:-N1,N2,N3,N4}"
DIAS_N_CHILDREN="${DIAS_N_CHILDREN:-3}"

run_dias() {
  local case_id="$1"
  shift
  local run_dir logfile
  run_dir="$(prepare_case_dir "${case_id}")"
  logfile="${MATRIX_ROOT}/logs/${case_id}_${TIMESTAMP}.log"
  log "${case_id}: run_dir=${run_dir}"
  if [[ ! -f "${run_dir}/${DIAS_LOG}" ]]; then
    record SKIP "${case_id}" "missing ${DIAS_LOG} in ${TEST_BATCH_DIR}"
    printf '%s\n' "${run_dir}"
    return 0
  fi
  (
    cd "${run_dir}"
    chemsmart sub \
      "${SUB_FLAGS[@]}" \
      -s "${SERVER}" \
      "$@" \
      gaussian \
      -p "${GAUSSIAN_PROJECT_DIAS}" \
      -f "${DIAS_LOG}" \
      --charge "${DIAS_CHARGE}" \
      --multiplicity "${DIAS_MULTIPLICITY}" \
      -l ts \
      dias \
      -i "${DIAS_FRAGMENTS}" \
      -n "${DIAS_EVERY_N}" \
      --mode "${DIAS_MODE}" \
      -c1 "${DIAS_FRAG1_CHARGE}" \
      -m1 "${DIAS_FRAG1_MULTIPLICITY}" \
      -c2 "${DIAS_FRAG2_CHARGE}" \
      -m2 "${DIAS_FRAG2_MULTIPLICITY}" \
      "${RERUN_FLAGS[@]}" \
      >"${logfile}" 2>&1
  )
  tail -n 5 "${logfile}" >&2 || true
  printf '%s\n' "${run_dir}"
}

run_case() {
  local case_id="$1"
  local run_dir line _old
  case "${case_id}" in
    N1)
      run_dir="$(run_dias N1 --no-run-in-parallel)"
      _old="$(pwd)"
      cd "${run_dir}"
      if [[ -f "${DIAS_LOG}" ]]; then
        if ls chemsmart_sub_array_*.sh >/dev/null 2>&1; then
          record FAIL N1 "unexpected array submit under --no-run-in-parallel"
        else
          # Submit script is named from -l ts, not the dias job token.
          expect_single_submit N1 'chemsmart_sub_ts*.sh' || \
            expect_single_submit N1 'chemsmart_sub_*.sh' || true
        fi
      fi
      cd "${_old}"
      ;;
    N2)
      run_dir="$(run_dias N2 --run-in-parallel)"
      _old="$(pwd)"
      cd "${run_dir}"
      if [[ -f "${DIAS_LOG}" ]]; then
        line="$(extract_array_line "$(find_newest 'chemsmart_sub_array_*.sh')")"
        if [[ "${line}" == "#SBATCH --array=1-${DIAS_N_CHILDREN}" ]] \
          || [[ "${line}" =~ ^#SBATCH\ --array=1-${DIAS_N_CHILDREN}%[0-9]+$ ]]; then
          record PASS N2 "${line}"
        else
          record FAIL N2 "expected array 1-${DIAS_N_CHILDREN}[%M], got '${line}'"
        fi
        expect_shared_runscript N2.runscript || true
        expect_child_index N2.child_index || true
      fi
      cd "${_old}"
      ;;
    N3)
      run_dir="$(run_dias N3 --run-in-parallel -M 1)"
      _old="$(pwd)"
      cd "${run_dir}"
      if [[ -f "${DIAS_LOG}" ]]; then
        expect_array N3 '^#SBATCH --array=1-'"${DIAS_N_CHILDREN}"'%1$' || true
        expect_shared_runscript N3.runscript || true
        expect_child_index N3.child_index || true
      fi
      cd "${_old}"
      ;;
    N4)
      run_dir="$(run_dias N4 --run-in-parallel -M 2)"
      _old="$(pwd)"
      cd "${run_dir}"
      if [[ -f "${DIAS_LOG}" ]]; then
        expect_array N4 '^#SBATCH --array=1-'"${DIAS_N_CHILDREN}"'%2$' || true
        expect_shared_runscript N4.runscript || true
        expect_child_index N4.child_index || true
      fi
      cd "${_old}"
      ;;
    *)
      record SKIP "${case_id}" "unknown N-case"
      ;;
  esac
}

log "N-dias: children=${DIAS_N_CHILDREN} mode=${DIAS_MODE} log=${DIAS_LOG} cases=${CASES} DRY_RUN=${DRY_RUN} WAIT=${WAIT}"
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
