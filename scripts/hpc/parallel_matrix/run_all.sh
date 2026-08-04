#!/usr/bin/env bash
# Run B + E + N-dias parallel matrix suites.
#
# Usage (sandbox):
#   DRY_RUN=1 ./run_all.sh
#   DRY_RUN=1 CASES_B=B1,B6 CASES_E=E1,E5 CASES_N=N1,N2 ./run_all.sh
#   DRY_RUN=1 RUN_MINIMUM=1 ./run_all.sh

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

DRY_RUN="${DRY_RUN:-1}"
WAIT="${WAIT:-0}"
CASES_B="${CASES_B:-B1,B2,B3,B4,B5,B6,B7}"
CASES_E="${CASES_E:-E1,E2,E3,E4,E5}"
CASES_N="${CASES_N:-N1,N2,N3,N4}"
RUN_MINIMUM="${RUN_MINIMUM:-0}"

echo "=== parallel matrix ALL (DRY_RUN=${DRY_RUN} WAIT=${WAIT}) ==="

status=0
CASES="${CASES_B}" DRY_RUN="${DRY_RUN}" WAIT="${WAIT}" \
  bash "${SCRIPT_DIR}/run_B.sh" || status=1

CASES="${CASES_E}" DRY_RUN="${DRY_RUN}" WAIT="${WAIT}" \
  bash "${SCRIPT_DIR}/run_E.sh" || status=1

CASES="${CASES_N}" DRY_RUN="${DRY_RUN}" WAIT="${WAIT}" \
  bash "${SCRIPT_DIR}/run_N_dias.sh" || status=1

if [[ "${RUN_MINIMUM}" == "1" ]]; then
  DRY_RUN="${DRY_RUN}" WAIT="${WAIT}" \
    bash "${SCRIPT_DIR}/minimum/run_minimum.sh" || status=1
fi

echo "=== ALL finished (aggregate exit=${status}) ==="
exit "${status}"
