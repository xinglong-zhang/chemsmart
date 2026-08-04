#!/usr/bin/env bash
# Minimum parallel matrix: 7 cases × 28 top-level job types (covering map).
#
# Usage (from sandbox):
#   DRY_RUN=1 ./run_minimum.sh
#   CASES=B1,B6 DRY_RUN=1 ./run_minimum.sh
#   JOBS=gaussian/sp,orca/sp DRY_RUN=1 ./run_minimum.sh
#   DRY_RUN=0 WAIT=1 CASES=E3 ./run_minimum.sh   # live SP wall B1 vs B6 shapes

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../common.sh
source "${SCRIPT_DIR}/../common.sh"

export MINIMUM_INPUTS="${MINIMUM_INPUTS:-${SCRIPT_DIR}/inputs}"
MAP_FILE="${MAP_FILE:-${SCRIPT_DIR}/JOB_CASE_MAP.tsv}"
GAUSSIAN_PROJECT_MIN="${GAUSSIAN_PROJECT_MIN:-${GAUSSIAN_PROJECT_OPT}}"
ORCA_PROJECT="${ORCA_PROJECT:-test}"
PKA_TABLE_MIN="${PKA_TABLE_MIN:-pka_scale.csv}"
NUM_FIVE="${NUM_FIVE:-5}"
NUM_TWO="${NUM_TWO:-2}"
NUM_PKA_MIN="${NUM_PKA_MIN:-2}"

ensure_chemsmart
sub_flags
rerun_flags
mkdir -p "${MATRIX_ROOT}/logs"
SUMMARY="${MATRIX_ROOT}/logs/minimum_summary_${TIMESTAMP}.tsv"
CASES="${CASES:-B1,B2,B6,O2,N1,N2,E3}"
JOBS="${JOBS:-}"

case_selected() {
  local case_id="$1"
  local c
  IFS=',' read -r -a _cases <<< "${CASES}"
  for c in "${_cases[@]}"; do
    c="$(echo "${c}" | tr -d '[:space:]')"
    [[ "${c}" == "${case_id}" ]] && return 0
  done
  return 1
}

job_selected() {
  local prog="$1"
  local job="$2"
  local key="${prog}/${job}"
  [[ -z "${JOBS}" ]] && return 0
  local j
  IFS=',' read -r -a _jobs <<< "${JOBS}"
  for j in "${_jobs[@]}"; do
    j="$(echo "${j}" | tr -d '[:space:]')"
    [[ "${j}" == "${key}" ]] && return 0
  done
  return 1
}

case_sub_flags() {
  # Populate CASE_FLAGS for the matrix case policy.
  local case_id="$1"
  CASE_FLAGS=()
  case "${case_id}" in
    B1) CASE_FLAGS+=(--no-run-in-parallel) ;;
    B2) CASE_FLAGS+=(--no-run-in-parallel -M "${NUM_FIVE}") ;;
    B6) CASE_FLAGS+=(--run-in-parallel -M "${NUM_FIVE}") ;;
    O2) CASE_FLAGS+=(--run-in-parallel -M "${NUM_TWO}") ;;
    N1) CASE_FLAGS+=(--no-run-in-parallel) ;;
    N2) CASE_FLAGS+=(--run-in-parallel -M 2) ;;
    E3) ;;
    *)
      echo "error: unknown case '${case_id}'" >&2
      return 1
      ;;
  esac
}

resolve_filename() {
  local file="$1"
  if [[ "${file}" == "DIAS_LOG" ]]; then
    printf '%s\n' "${DIAS_LOG}"
  else
    printf '%s\n' "${file}"
  fi
}

build_program_cli() {
  # Sets PROGRAM_CLI array: program -p ... -f ... [ -i ... ] [charge/mult] [dias opts]
  local case_id="$1"
  local program="$2"
  local job="$3"
  local filename="$4"
  local indices="$5"
  local row_label="$6"
  PROGRAM_CLI=("${program}")
  if [[ "${program}" == "gaussian" ]]; then
    if [[ "${job}" == "dias" ]]; then
      PROGRAM_CLI+=(-p "${GAUSSIAN_PROJECT_DIAS}")
    elif [[ "${job}" == "pka" ]]; then
      PROGRAM_CLI+=(-p "${GAUSSIAN_PROJECT_MIN}")
    else
      PROGRAM_CLI+=(-p "${GAUSSIAN_PROJECT_MIN}")
    fi
  else
    PROGRAM_CLI+=(-p "${ORCA_PROJECT}")
  fi
  PROGRAM_CLI+=(-f "${filename}")
  if [[ "${indices}" != "-" && -n "${indices}" ]]; then
    PROGRAM_CLI+=(-i "${indices}")
  fi
  if [[ "${job}" == "dias" ]]; then
    PROGRAM_CLI+=(
      --charge "${DIAS_CHARGE}"
      --multiplicity "${DIAS_MULTIPLICITY}"
    )
  else
    PROGRAM_CLI+=(--charge 0 --multiplicity 1)
  fi
  # Colliding labels (job-token names) regression-test CLI rewrite anchoring.
  case "${job}" in
    qrc) PROGRAM_CLI+=(-l ts) ;;
    dias) PROGRAM_CLI+=(-l irc) ;;
    crest|traj|opt) PROGRAM_CLI+=(-l opt) ;;
    *) PROGRAM_CLI+=(-l "${row_label}") ;;
  esac
}

build_job_cli() {
  # Sets JOB_CLI: job token + optional args (+ dias fragment opts).
  local job="$1"
  local job_args="$2"
  JOB_CLI=("${job}")
  if [[ "${job}" == "dias" ]]; then
    JOB_CLI+=(
      -i "${DIAS_FRAGMENTS}"
      -n "${DIAS_EVERY_N}"
      --mode "${DIAS_MODE}"
      -c1 "${DIAS_FRAG1_CHARGE}"
      -m1 "${DIAS_FRAG1_MULTIPLICITY}"
      -c2 "${DIAS_FRAG2_CHARGE}"
      -m2 "${DIAS_FRAG2_MULTIPLICITY}"
    )
  elif [[ "${job_args}" != "-" && -n "${job_args}" ]]; then
    local token
    while IFS= read -r token; do
      [[ -n "${token}" ]] && JOB_CLI+=("${token}")
    done < <(python3 -c 'import shlex,sys; print("\n".join(shlex.split(sys.argv[1])))' "${job_args}")
  fi
}

assert_expect() {
  local tag="$1"
  local expect="$2"
  local n_tasks="${3:-${NUM_FIVE}}"
  case "${expect}" in
    array_serial)
      expect_array "${tag}" "^#SBATCH --array=1-${n_tasks}%1$" || true
      expect_shared_runscript "${tag}.runscript" || true
      ;;
    array_serial_M)
      expect_array "${tag}" "^#SBATCH --array=1-${n_tasks}%1$" || true
      expect_shared_runscript "${tag}.runscript" || true
      ;;
    array_parallel)
      expect_array "${tag}" "^#SBATCH --array=1-${n_tasks}%${n_tasks}$" || true
      expect_shared_runscript "${tag}.runscript" || true
      ;;
    array_serial_pka)
      expect_array "${tag}" "^#SBATCH --array=1-${NUM_PKA_MIN}%1$" || true
      expect_shared_runscript "${tag}.runscript" || true
      ;;
    o2_labels)
      expect_array "${tag}" "^#SBATCH --array=1-${NUM_TWO}%${NUM_TWO}$" || true
      expect_shared_runscript "${tag}.runscript" || true
      expect_o2_labels "${tag}.labels" || true
      ;;
    n1_nestable)
      if ls chemsmart_sub_array_*.sh >/dev/null 2>&1; then
        record FAIL "${tag}" "unexpected array submit under N1 / --no-run-in-parallel"
      else
        expect_single_submit "${tag}" 'chemsmart_sub_*.sh' || true
      fi
      ;;
    n2_nestable)
      expect_array "${tag}" '^#SBATCH --array=1-[0-9]+(%[0-9]+)?$' || true
      expect_shared_runscript "${tag}.runscript" || true
      expect_child_index "${tag}.child" || true
      ;;
    single)
      expect_single_submit "${tag}" 'chemsmart_sub_*.sh' || true
      ;;
    *)
      record FAIL "${tag}" "unknown expect kind '${expect}'"
      ;;
  esac
}

run_map_row() {
  local case_id="$1"
  local program="$2"
  local job="$3"
  local input_kind="$4"
  local file="$5"
  local indices="$6"
  local job_args="$7"
  local expect="$8"

  local tag="${case_id}/${program}/${job}"
  local filename row_label run_dir logfile rc _old
  filename="$(resolve_filename "${file}")"
  row_label="min_${case_id}_${program}_${job}"
  # Unique workdir per job type (cases are reused across rows).
  run_dir="$(prepare_case_dir "min_${case_id}_${program}_${job}" "${input_kind}")"
  logfile="${MATRIX_ROOT}/logs/min_${case_id}_${program}_${job}_${TIMESTAMP}.log"

  case_sub_flags "${case_id}" || {
    record FAIL "${tag}" "bad case flags"
    return 0
  }
  build_program_cli "${case_id}" "${program}" "${job}" "${filename}" "${indices}" "${row_label}"
  build_job_cli "${job}" "${job_args}"

  if [[ ! -e "${run_dir}/${filename}" ]]; then
    record SKIP "${tag}" "missing input ${filename} (kind=${input_kind})"
    return 0
  fi

  log "${tag}: run_dir=${run_dir}"
  rc=0
  # Singles are script-gen checks only; live compute is reserved for E3 timing.
  local row_sub_flags=("${SUB_FLAGS[@]}")
  if [[ "${expect}" == "single" && "${DRY_RUN}" != "1" ]]; then
    row_sub_flags=(--test)
  fi
  (
    cd "${run_dir}"
    chemsmart sub \
      "${row_sub_flags[@]}" \
      -s "${SERVER}" \
      "${CASE_FLAGS[@]}" \
      "${PROGRAM_CLI[@]}" \
      "${JOB_CLI[@]}" \
      "${RERUN_FLAGS[@]}" \
      >"${logfile}" 2>&1
  ) || rc=$?

  if [[ "${rc}" -ne 0 ]]; then
    record FAIL "${tag}" "chemsmart sub exit ${rc}; see ${logfile}"
    tail -n 20 "${logfile}" >&2 || true
    return 0
  fi

  _old="$(pwd)"
  cd "${run_dir}"
  assert_expect "${tag}" "${expect}"
  cd "${_old}"
  tail -n 3 "${logfile}" >&2 || true
}

run_e3_timing() {
  # Live wall time: five_mols gaussian SP serial (%1) vs parallel (%5).
  local run_b1 run_b6 t1 t6 wall1 wall6 logfile1 logfile6
  if [[ "${DRY_RUN}" == "1" ]]; then
    record SKIP "E3/timing" "requires DRY_RUN=0 WAIT=1"
    return 0
  fi
  if [[ "${WAIT}" != "1" ]]; then
    record SKIP "E3/timing" "requires WAIT=1"
    return 0
  fi

  run_b1="$(prepare_case_dir "min_E3_timing_B1" five_mols)"
  run_b6="$(prepare_case_dir "min_E3_timing_B6" five_mols)"
  logfile1="${MATRIX_ROOT}/logs/min_E3_timing_B1_${TIMESTAMP}.log"
  logfile6="${MATRIX_ROOT}/logs/min_E3_timing_B6_${TIMESTAMP}.log"

  if [[ ! -f "${run_b1}/five_mols.xyz" ]]; then
    record SKIP "E3/timing" "missing five_mols.xyz"
    return 0
  fi

  t1="$(date +%s)"
  (
    cd "${run_b1}"
    chemsmart sub -s "${SERVER}" --no-run-in-parallel \
      gaussian -p "${GAUSSIAN_PROJECT_MIN}" -f five_mols.xyz -i 1,2,3,4,5 \
      --charge 0 --multiplicity 1 -l min_E3_sp_serial \
      sp -R \
      >"${logfile1}" 2>&1
  )
  wait_for_job_ids_from_log "${logfile1}"
  wall1=$(( $(date +%s) - t1 ))

  t6="$(date +%s)"
  (
    cd "${run_b6}"
    chemsmart sub -s "${SERVER}" --run-in-parallel -M "${NUM_FIVE}" \
      gaussian -p "${GAUSSIAN_PROJECT_MIN}" -f five_mols.xyz -i 1,2,3,4,5 \
      --charge 0 --multiplicity 1 -l min_E3_sp_parallel \
      sp -R \
      >"${logfile6}" 2>&1
  )
  wait_for_job_ids_from_log "${logfile6}"
  wall6=$(( $(date +%s) - t6 ))

  record PASS "E3/timing" "wall_s serial(B1)=${wall1} parallel(B6)=${wall6}"
}

# ---- main ----
if [[ ! -f "${MAP_FILE}" ]]; then
  echo "error: map not found: ${MAP_FILE}" >&2
  exit 1
fi

log "minimum matrix: map=${MAP_FILE} cases=${CASES} DRY_RUN=${DRY_RUN} WAIT=${WAIT}"
log "MINIMUM_INPUTS=${MINIMUM_INPUTS}"

while IFS=$'\t' read -r case_id program job input_kind file indices job_args expect || [[ -n "${case_id:-}" ]]; do
  [[ -z "${case_id}" || "${case_id}" =~ ^# || "${case_id}" == "case" ]] && continue
  case_selected "${case_id}" || continue
  job_selected "${program}" "${job}" || continue
  run_map_row "${case_id}" "${program}" "${job}" "${input_kind}" "${file}" "${indices}" "${job_args}" "${expect}"
done < "${MAP_FILE}"

if case_selected E3; then
  run_e3_timing
fi

print_summary "${SUMMARY}"
[[ "${FAIL}" -eq 0 ]]
