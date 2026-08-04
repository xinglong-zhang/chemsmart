#!/usr/bin/env bash
# Shared helpers for parallel-matrix HPC tests.
# shellcheck disable=SC2034

set -euo pipefail

PROJECT="${PROJECT:-/project/xlzhang/jingyi}"
PROJECT="${PROJECT%/}"
TEST_BATCH_DIR="${TEST_BATCH_DIR:-${PROJECT}/chemsmart/test_batch}"
PKA_DIR="${PKA_DIR:-${PROJECT}/chemsmart/test_pka/gaussian}"
MATRIX_ROOT="${MATRIX_ROOT:-${TEST_BATCH_DIR}/parallel_matrix}"
SERVER="${SERVER:-cu_batch}"
GAUSSIAN_PROJECT_BATCH="${GAUSSIAN_PROJECT_BATCH:-pka_wat}"
GAUSSIAN_PROJECT_OPT="${GAUSSIAN_PROJECT_OPT:-test}"
GAUSSIAN_PROJECT_DIAS="${GAUSSIAN_PROJECT_DIAS:-test}"
PKA_TABLE="${PKA_TABLE:-pka_scale.csv}"
SCHEME="${SCHEME:-direct}"
WATER_COM="${WATER_COM:-water.com}"
DIAS_LOG="${DIAS_LOG:-model_R_c5_ts_irc_flatr_flat_opt.log}"
DIAS_FRAGMENTS="${DIAS_FRAGMENTS:-1-10}"
DIAS_EVERY_N="${DIAS_EVERY_N:-20}"
DIAS_MODE="${DIAS_MODE:-ts}"
DIAS_CHARGE="${DIAS_CHARGE:-2}"
DIAS_MULTIPLICITY="${DIAS_MULTIPLICITY:-2}"
DIAS_FRAG1_CHARGE="${DIAS_FRAG1_CHARGE:-0}"
DIAS_FRAG1_MULTIPLICITY="${DIAS_FRAG1_MULTIPLICITY:-1}"
DIAS_FRAG2_CHARGE="${DIAS_FRAG2_CHARGE:-2}"
DIAS_FRAG2_MULTIPLICITY="${DIAS_FRAG2_MULTIPLICITY:-2}"
# Minimum suite fixtures (override when inputs live under MATRIX_ROOT on HPC)
MINIMUM_INPUTS="${MINIMUM_INPUTS:-}"
ORCA_PROJECT="${ORCA_PROJECT:-test}"
DRY_RUN="${DRY_RUN:-1}"
WAIT="${WAIT:-0}"
FORCE_RERUN="${FORCE_RERUN:-1}"
POLL_SEC="${POLL_SEC:-30}"
TIMESTAMP="${TIMESTAMP:-$(date +%Y%m%d_%H%M%S)}"

export CHEMSMART_CONFIG_DIR="${CHEMSMART_CONFIG_DIR:-${HOME}/.chemsmart}"

PASS=0
FAIL=0
SKIP=0
RESULTS=()

ensure_chemsmart() {
  if ! command -v chemsmart >/dev/null 2>&1; then
    echo "error: chemsmart not found on PATH (use sandbox + conda activate chemsmart)" >&2
    exit 1
  fi
}

log() {
  # stderr so command substitutions only capture explicit stdout (e.g. run_dir)
  printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2
}

record() {
  local status="$1"
  local case_id="$2"
  local detail="$3"
  RESULTS+=("${status}"$'\t'"${case_id}"$'\t'"${detail}")
  case "${status}" in
    PASS) PASS=$((PASS + 1)) ;;
    FAIL) FAIL=$((FAIL + 1)) ;;
    SKIP) SKIP=$((SKIP + 1)) ;;
  esac
  log "${status}: ${case_id} — ${detail}"
}

link_pka_inputs() {
  local run_dir="$1"
  if [[ ! -f "${PKA_DIR}/${PKA_TABLE}" ]]; then
    return 0
  fi
  ln -sfn "${PKA_DIR}/${PKA_TABLE}" "${run_dir}/${PKA_TABLE}"
  local xyz
  for xyz in "${PKA_DIR}"/*.xyz; do
    [[ -f "${xyz}" ]] || continue
    ln -sfn "${xyz}" "${run_dir}/$(basename "${xyz}")"
  done
}

link_water_inputs() {
  local run_dir="$1"
  local f
  for f in "${WATER_COM}" water.xyz; do
    if [[ -f "${TEST_BATCH_DIR}/${f}" ]]; then
      ln -sfn "${TEST_BATCH_DIR}/${f}" "${run_dir}/${f}"
    fi
  done
}

link_dias_inputs() {
  local run_dir="$1"
  if [[ -f "${TEST_BATCH_DIR}/${DIAS_LOG}" ]]; then
    ln -sfn "${TEST_BATCH_DIR}/${DIAS_LOG}" "${run_dir}/${DIAS_LOG}"
  fi
}

link_minimum_file() {
  local run_dir="$1"
  local name="$2"
  local src=""
  if [[ -n "${MINIMUM_INPUTS}" && -f "${MINIMUM_INPUTS}/${name}" ]]; then
    src="${MINIMUM_INPUTS}/${name}"
  elif [[ -f "${MATRIX_ROOT}/minimum/inputs/${name}" ]]; then
    src="${MATRIX_ROOT}/minimum/inputs/${name}"
  elif [[ -f "${TEST_BATCH_DIR}/parallel_matrix/minimum/inputs/${name}" ]]; then
    src="${TEST_BATCH_DIR}/parallel_matrix/minimum/inputs/${name}"
  fi
  if [[ -n "${src}" ]]; then
    ln -sfn "${src}" "${run_dir}/${name}"
  fi
}

link_five_mols_inputs() {
  local run_dir="$1"
  link_minimum_file "${run_dir}" five_mols.xyz
}

link_two_mols_inputs() {
  local run_dir="$1"
  link_minimum_file "${run_dir}" two_mols.xyz
}

link_single_inputs() {
  local run_dir="$1"
  local f
  for f in h2.xyz water.com water.xyz water.inp; do
    link_minimum_file "${run_dir}" "${f}"
  done
}

link_minimum_pka_inputs() {
  local run_dir="$1"
  local f
  for f in pka_scale.csv acid1.xyz acid2.xyz; do
    link_minimum_file "${run_dir}" "${f}"
  done
}

link_qrc_inputs() {
  local run_dir="$1"
  link_minimum_file "${run_dir}" ts_freq.log
}

link_orca_qrc_inputs() {
  local run_dir="$1"
  link_minimum_file "${run_dir}" orca_freq.out
}

link_neb_inputs() {
  local run_dir="$1"
  link_minimum_file "${run_dir}" h2.xyz
  link_minimum_file "${run_dir}" h2_product.xyz
}

case_input_kind() {
  # Map matrix case id → input set (override with prepare_case_dir kind arg).
  local case_id="$1"
  case "${case_id}" in
    B*|E2*|E3*|E4*|E5*) printf '%s\n' pka ;;
    E1) printf '%s\n' water ;;
    N*) printf '%s\n' dias ;;
    O2) printf '%s\n' two_mols ;;
    *)
      echo "error: unknown case id '${case_id}' for input staging" >&2
      return 1
      ;;
  esac
}

prepare_case_dir() {
  local case_id="$1"
  local kind="${2:-}"
  local run_dir="${MATRIX_ROOT}/runs/${case_id}_${TIMESTAMP}"
  mkdir -p "${run_dir}"
  if [[ -z "${kind}" ]]; then
    kind="$(case_input_kind "${case_id}")"
  fi
  case "${kind}" in
    pka) link_pka_inputs "${run_dir}" ;;
    water) link_water_inputs "${run_dir}" ;;
    dias) link_dias_inputs "${run_dir}" ;;
    five_mols) link_five_mols_inputs "${run_dir}" ;;
    two_mols) link_two_mols_inputs "${run_dir}" ;;
    single) link_single_inputs "${run_dir}" ;;
    pka_min|minimum_pka) link_minimum_pka_inputs "${run_dir}" ;;
    qrc) link_qrc_inputs "${run_dir}" ;;
    orca_qrc) link_orca_qrc_inputs "${run_dir}" ;;
    neb) link_neb_inputs "${run_dir}" ;;
    *)
      echo "error: unknown input kind '${kind}'" >&2
      return 1
      ;;
  esac
  printf '%s\n' "${run_dir}"
}

sub_flags() {
  SUB_FLAGS=()
  if [[ "${DRY_RUN}" == "1" ]]; then
    SUB_FLAGS+=(--test)
  fi
}

rerun_flags() {
  RERUN_FLAGS=()
  if [[ "${FORCE_RERUN}" == "1" ]]; then
    RERUN_FLAGS+=(-R)
  fi
}

extract_array_line() {
  local shfile="$1"
  if [[ ! -f "${shfile}" ]]; then
    echo ""
    return 0
  fi
  grep -E '^#SBATCH --array=' "${shfile}" | head -1 | tr -d '\r' || true
}

find_newest() {
  local pattern="$1"
  # shellcheck disable=SC2086
  ls -t ${pattern} 2>/dev/null | head -1 || true
}

find_array_runscript() {
  # Shared per-array run script: chemsmart_run_array_<array_log_stem>.py
  # (not legacy chemsmart_run_array.py or per-task chemsmart_run_array_N.py).
  local f base
  local -a shared=()
  local -a leftover=()
  for f in chemsmart_run_array_*.py; do
    [[ -f "${f}" ]] || continue
    base="${f#chemsmart_run_array_}"
    base="${base%.py}"
    if [[ "${base}" =~ ^[0-9]+$ ]]; then
      leftover+=("${f}")
    else
      shared+=("${f}")
    fi
  done
  if [[ ${#leftover[@]} -gt 0 ]]; then
    printf 'LEFTOVER:%s\n' "${leftover[*]}"
    return 1
  fi
  if [[ ${#shared[@]} -eq 0 ]]; then
    printf '\n'
    return 1
  fi
  ls -t "${shared[@]}" 2>/dev/null | head -1
}

expect_array() {
  # expect_array CASE_ID EXPECTED_REGEX [SUBMIT_GLOB]
  local case_id="$1"
  local expected_re="$2"
  local submit_glob="${3:-chemsmart_sub_array_*.sh}"
  local shfile
  shfile="$(find_newest "${submit_glob}")"
  local line
  line="$(extract_array_line "${shfile}")"
  if [[ -z "${shfile}" ]]; then
    record FAIL "${case_id}" "missing submit script matching ${submit_glob}"
    return 1
  fi
  if [[ "${line}" =~ ${expected_re} ]]; then
    record PASS "${case_id}" "${shfile##*/}: ${line}"
    return 0
  fi
  record FAIL "${case_id}" "expected /${expected_re}/, got '${line}' in ${shfile##*/}"
  return 1
}

expect_single_submit() {
  local case_id="$1"
  local submit_glob="${2:-chemsmart_sub_*.sh}"
  local shfile
  shfile="$(find_newest "${submit_glob}")"
  if [[ -z "${shfile}" ]]; then
    record FAIL "${case_id}" "missing single submit script ${submit_glob}"
    return 1
  fi
  if grep -qE '^#SBATCH --array=' "${shfile}"; then
    record FAIL "${case_id}" "expected non-array submit, found array in ${shfile##*/}"
    return 1
  fi
  record PASS "${case_id}" "non-array submit: ${shfile##*/}"
  return 0
}

expect_shared_runscript() {
  local case_id="$1"
  local runscript status
  status=0
  runscript="$(find_array_runscript)" || status=$?
  if [[ "${runscript}" == LEFTOVER:* ]]; then
    record FAIL "${case_id}" "found leftover per-task ${runscript#LEFTOVER:}"
    return 1
  fi
  if [[ "${status}" -ne 0 || -z "${runscript}" || ! -f "${runscript}" ]]; then
    record FAIL "${case_id}" "missing chemsmart_run_array_<array_log_stem>.py"
    return 1
  fi
  if [[ -f chemsmart_run_array.py ]]; then
    record FAIL "${case_id}" "legacy chemsmart_run_array.py still present alongside ${runscript##*/}"
    return 1
  fi
  record PASS "${case_id}" "shared ${runscript##*/} present"
  return 0
}

expect_o2_labels() {
  # Distinct --label values in TASK_CLI for a 2-mol O2 batch.
  local case_id="$1"
  local runscript status
  status=0
  runscript="$(find_array_runscript)" || status=$?
  if [[ "${status}" -ne 0 || -z "${runscript}" || ! -f "${runscript}" ]]; then
    record FAIL "${case_id}" "missing chemsmart_run_array_<stem>.py for O2 label check"
    return 1
  fi
  local nlabels
  nlabels="$(
    python3 - "${runscript}" <<'PY'
import ast, pathlib, re, sys
text = pathlib.Path(sys.argv[1]).read_text()
m = re.search(r"TASK_CLI\s*=\s*(\{.*?\})", text, re.S)
if not m:
    print(0)
    raise SystemExit
task_cli = ast.literal_eval(m.group(1))
labels = set()
for args in task_cli.values():
    for i, tok in enumerate(args):
        if tok in ("-l", "--label") and i + 1 < len(args):
            labels.add(args[i + 1])
print(len(labels))
PY
  )"
  if [[ "${nlabels}" -ge 2 ]]; then
    record PASS "${case_id}" "TASK_CLI has ${nlabels} distinct --label values (${runscript##*/})"
    return 0
  fi
  record FAIL "${case_id}" "expected ≥2 distinct --label in TASK_CLI, got ${nlabels}"
  return 1
}

expect_child_index() {
  # Require --child-index after the nestable job token in every TASK_CLI entry.
  local case_id="$1"
  local runscript status
  status=0
  runscript="$(find_array_runscript)" || status=$?
  if [[ "${status}" -ne 0 || -z "${runscript}" || ! -f "${runscript}" ]]; then
    record FAIL "${case_id}" "missing chemsmart_run_array_<stem>.py for --child-index check"
    return 1
  fi
  local check_status
  check_status="$(
    python3 - "${runscript}" <<'PY'
import ast, pathlib, re, sys
text = pathlib.Path(sys.argv[1]).read_text()
m = re.search(r"TASK_CLI\s*=\s*(\{.*?\})", text, re.S)
if not m:
    print("missing_task_cli")
    raise SystemExit
task_cli = ast.literal_eval(m.group(1))
nestable = {"qrc", "dias", "crest", "traj"}
for args in task_cli.values():
    try:
        ci = args.index("--child-index")
    except ValueError:
        print("missing_child_index")
        raise SystemExit
    job_pos = [i for i, tok in enumerate(args) if tok in nestable]
    if not job_pos or max(job_pos) > ci:
        print("child_index_before_job")
        raise SystemExit
print("ok")
PY
  )"
  if [[ "${check_status}" == "ok" ]]; then
    record PASS "${case_id}" "TASK_CLI has --child-index after nestable job token (${runscript##*/})"
    return 0
  fi
  record FAIL "${case_id}" "bad --child-index placement (${check_status})"
  return 1
}

print_summary() {
  local out="${1:-}"
  {
    echo
    echo "==== parallel matrix summary (${TIMESTAMP}) ===="
    echo -e "status\tcase\tdetail"
    local i=0
    while [[ $i -lt ${#RESULTS[@]} ]]; do
      echo -e "${RESULTS[$i]}"
      i=$((i + 1))
    done
    echo "PASS=${PASS} FAIL=${FAIL} SKIP=${SKIP}"
    echo "DRY_RUN=${DRY_RUN} WAIT=${WAIT} SERVER=${SERVER}"
  } | tee ${out:+"${out}"}
}

wait_for_job_ids_from_log() {
  local logfile="$1"
  mapfile -t ids < <(grep -Eo 'Submitted batch job [0-9]+' "${logfile}" | awk '{print $NF}' | sort -u)
  if [[ "${#ids[@]}" -eq 0 ]]; then
    log "no submitted job ids in ${logfile}"
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
    [[ "${pending}" -eq 0 ]] && break
    sleep "${POLL_SEC}"
  done
}
