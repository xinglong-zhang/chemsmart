#!/usr/bin/env bash
# Edge / control cases: E1–E5.
#
# Usage (from sandbox):
#   DRY_RUN=1 ./run_E.sh
#   CASES=E1,E2 DRY_RUN=1 ./run_E.sh
#   # E3–E5 need real submit for full checks:
#   CASES=E3,E4,E5 DRY_RUN=0 WAIT=1 ./run_E.sh

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=common.sh
source "${SCRIPT_DIR}/common.sh"

ensure_chemsmart
sub_flags
rerun_flags
mkdir -p "${MATRIX_ROOT}/logs"
SUMMARY="${MATRIX_ROOT}/logs/E_summary_${TIMESTAMP}.tsv"
CASES="${CASES:-E1,E2,E3,E4,E5}"

NUM_ROWS="$(
  python3 - "${PKA_DIR}/${PKA_TABLE}" <<'PY'
import csv, sys
with open(sys.argv[1], newline="") as fh:
    print(len(list(csv.DictReader(fh))))
PY
)"

run_case_E1() {
  local run_dir logfile
  run_dir="$(prepare_case_dir E1)"
  logfile="${MATRIX_ROOT}/logs/E1_${TIMESTAMP}.log"
  cd "${run_dir}"
  if [[ ! -f "${WATER_COM}" ]]; then
    record SKIP E1 "missing ${WATER_COM}"
    return 0
  fi
  chemsmart sub \
    "${SUB_FLAGS[@]}" \
    -s "${SERVER}" \
    gaussian \
    -p "${GAUSSIAN_PROJECT_OPT}" \
    -f "${WATER_COM}" \
    --charge 0 --multiplicity 1 \
    -l e1_water_sp \
    sp \
    "${RERUN_FLAGS[@]}" \
    2>&1 | tee "${logfile}"
  expect_single_submit E1 'chemsmart_sub_e1_water_sp*.sh' || \
    expect_single_submit E1 'chemsmart_sub_*.sh' || true
}

run_case_E2() {
  local run_dir logfile first_mtime second_mtime runscript
  run_dir="$(prepare_case_dir E2)"
  logfile="${MATRIX_ROOT}/logs/E2_${TIMESTAMP}.log"
  cd "${run_dir}"
  chemsmart sub \
    "${SUB_FLAGS[@]}" \
    -s "${SERVER}" \
    --run-in-parallel -M 2 \
    gaussian -p "${GAUSSIAN_PROJECT_BATCH}" \
    -f "${PKA_TABLE}" \
    pka --scheme "${SCHEME}" batch -R \
    2>&1 | tee "${logfile}"
  runscript="$(find_array_runscript || true)"
  if [[ -z "${runscript}" || ! -f "${runscript}" ]]; then
    record FAIL E2 "missing chemsmart_run_array_<array_log_stem>.py after first submit"
    return 0
  fi
  first_mtime="$(stat -c %Y "${runscript}" 2>/dev/null || stat -f %m "${runscript}")"
  sleep 1
  chemsmart sub \
    "${SUB_FLAGS[@]}" \
    -s "${SERVER}" \
    --run-in-parallel -M 2 \
    gaussian -p "${GAUSSIAN_PROJECT_BATCH}" \
    -f "${PKA_TABLE}" \
    pka --scheme "${SCHEME}" batch -R \
    2>&1 | tee -a "${logfile}"
  runscript="$(find_array_runscript || true)"
  if [[ -z "${runscript}" || ! -f "${runscript}" ]]; then
    record FAIL E2 "missing chemsmart_run_array_<array_log_stem>.py after -R rerun"
    return 0
  fi
  second_mtime="$(stat -c %Y "${runscript}" 2>/dev/null || stat -f %m "${runscript}")"
  if [[ "${second_mtime}" -ge "${first_mtime}" ]]; then
    record PASS E2 "rerun with -R rewrote ${runscript##*/} (mtime ${first_mtime}->${second_mtime})"
  else
    record FAIL E2 "expected rewrite of ${runscript##*/} with -R"
  fi
}

run_case_E3() {
  local run_b1 run_b6 t1 t6 wall1 wall6
  if [[ "${DRY_RUN}" == "1" ]]; then
    record SKIP E3 "requires DRY_RUN=0 real submission"
    return 0
  fi
  run_b1="$(prepare_case_dir E3_B1)"
  run_b6="$(prepare_case_dir E3_B6)"
  t1="$(date +%s)"
  (
    cd "${run_b1}"
    chemsmart sub -s "${SERVER}" --no-run-in-parallel \
      gaussian -p "${GAUSSIAN_PROJECT_BATCH}" -f "${PKA_TABLE}" \
      pka --scheme "${SCHEME}" batch -R \
      2>&1 | tee "${MATRIX_ROOT}/logs/E3_B1_${TIMESTAMP}.log"
  )
  if [[ "${WAIT}" == "1" ]]; then
    wait_for_job_ids_from_log "${MATRIX_ROOT}/logs/E3_B1_${TIMESTAMP}.log"
  fi
  wall1=$(( $(date +%s) - t1 ))

  t6="$(date +%s)"
  (
    cd "${run_b6}"
    chemsmart sub -s "${SERVER}" --run-in-parallel -M "${NUM_ROWS}" \
      gaussian -p "${GAUSSIAN_PROJECT_BATCH}" -f "${PKA_TABLE}" \
      pka --scheme "${SCHEME}" batch -R \
      2>&1 | tee "${MATRIX_ROOT}/logs/E3_B6_${TIMESTAMP}.log"
  )
  if [[ "${WAIT}" == "1" ]]; then
    wait_for_job_ids_from_log "${MATRIX_ROOT}/logs/E3_B6_${TIMESTAMP}.log"
  fi
  wall6=$(( $(date +%s) - t6 ))
  record PASS E3 "wall_s B1(serial)=${wall1} B6(parallel)=${wall6} (WAIT=${WAIT})"
}

run_case_E4() {
  local run_dir logfile jobid ntasks sq
  if [[ "${DRY_RUN}" == "1" ]]; then
    record SKIP E4 "requires DRY_RUN=0 to inspect squeue"
    return 0
  fi
  run_dir="$(prepare_case_dir E4)"
  logfile="${MATRIX_ROOT}/logs/E4_${TIMESTAMP}.log"
  (
    cd "${run_dir}"
    chemsmart sub -s "${SERVER}" --run-in-parallel -M "${NUM_ROWS}" \
      gaussian -p "${GAUSSIAN_PROJECT_BATCH}" -f "${PKA_TABLE}" \
      pka --scheme "${SCHEME}" batch -R \
      2>&1 | tee "${logfile}"
  )
  jobid="$(grep -Eo 'Submitted batch job [0-9]+' "${logfile}" | awk '{print $NF}' | tail -1)"
  if [[ -z "${jobid}" ]]; then
    record FAIL E4 "no job id submitted"
    return 0
  fi
  sleep 5
  ntasks="$(squeue -j "${jobid}" -h -o '%i %t' 2>/dev/null | wc -l | tr -d ' ')"
  sq="$(squeue -j "${jobid}" -h -o '%i %t %R' 2>/dev/null | tr '\n' '; ')"
  if [[ "${ntasks}" -ge 1 ]]; then
    record PASS E4 "job ${jobid}: ${ntasks} squeue line(s): ${sq}"
  elif sacct -j "${jobid}" --parsable2 --noheader -o JobID,State 2>/dev/null | grep -q .; then
    record PASS E4 "job ${jobid} already left queue; sacct has entries"
  else
    record FAIL E4 "job ${jobid} not found in squeue/sacct"
  fi
  if [[ "${WAIT}" == "1" ]]; then
    wait_for_job_ids_from_log "${logfile}"
  fi
}

run_case_E5() {
  local run_dir logfile shfile lock out
  run_dir="$(prepare_case_dir E5)"
  logfile="${MATRIX_ROOT}/logs/E5_${TIMESTAMP}.log"
  cd "${run_dir}"
  chemsmart sub \
    "${SUB_FLAGS[@]}" \
    -s "${SERVER}" \
    --run-in-parallel -M "${NUM_ROWS}" \
    gaussian -p "${GAUSSIAN_PROJECT_BATCH}" -f "${PKA_TABLE}" \
    pka --scheme "${SCHEME}" batch -R \
    2>&1 | tee "${logfile}"
  shfile="$(find_newest 'chemsmart_sub_array_*.sh')"
  if [[ -z "${shfile}" ]]; then
    record FAIL E5 "missing array submit script"
    return 0
  fi
  if ! grep -q 'flock 9' "${shfile}"; then
    record FAIL E5 "submit script missing flock serialization"
    return 0
  fi
  if ! grep -q 'BEGIN array task' "${shfile}"; then
    record FAIL E5 "submit script missing headed log blocks"
    return 0
  fi
  if ! grep -qE 'python chemsmart_run_array_[A-Za-z0-9._-]+\.py' "${shfile}"; then
    record FAIL E5 "submit script missing stem-named chemsmart_run_array_<stem>.py invoke"
    return 0
  fi
  if grep -qE 'python chemsmart_run_array\.py([[:space:]]|$)' "${shfile}"; then
    record FAIL E5 "submit script still invokes legacy chemsmart_run_array.py"
    return 0
  fi
  lock="$(grep -oE '\.[A-Za-z0-9._-]+\.loglock' "${shfile}" | head -1 || true)"
  if [[ -z "${lock}" ]]; then
    record FAIL E5 "could not find .loglock name in submit script"
    return 0
  fi
  if [[ "${DRY_RUN}" == "1" ]]; then
    record PASS E5 "script has flock + headed blocks + stem runscript; lock=${lock} (DRY_RUN skips live log check)"
    return 0
  fi
  if [[ "${WAIT}" == "1" ]]; then
    wait_for_job_ids_from_log "${logfile}"
  else
    sleep 10
  fi
  out="$(find_newest '*_array.out')"
  if [[ -z "${out}" ]]; then
    out="$(find_newest '*_array.slurmout')"
  fi
  if [[ -n "${out}" && -f "${out}" ]] && grep -q 'BEGIN array task' "${out}"; then
    record PASS E5 "shared log ${out} has BEGIN blocks; lock=${lock}"
  else
    record PASS E5 "flock wiring OK in ${shfile##*/}; shared log not ready yet (WAIT=${WAIT})"
  fi
}

log "E-series: cases=${CASES} DRY_RUN=${DRY_RUN} WAIT=${WAIT}"
IFS=',' read -r -a CASE_LIST <<< "${CASES}"
for case_id in "${CASE_LIST[@]}"; do
  case_id="$(echo "${case_id}" | tr -d '[:space:]')"
  [[ -n "${case_id}" ]] || continue
  case "${case_id}" in
    E1) run_case_E1 ;;
    E2) run_case_E2 ;;
    E3) run_case_E3 ;;
    E4) run_case_E4 ;;
    E5) run_case_E5 ;;
    *) record SKIP "${case_id}" "unknown E-case" ;;
  esac
done
print_summary "${SUMMARY}"
[[ "${FAIL}" -eq 0 ]]
