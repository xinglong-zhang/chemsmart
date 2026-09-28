#!/bin/bash
# R11 truth-4 item 3 records check (read-only, no python): for every goal
# ledger under the roots given, each analysis_evidence_recorded stream and
# whether a session_stream_recorded row of the same ledger names it.
# One line per ledger that holds analysis evidence:
#   <first row date> ev=<evidence streams> named=<of them session-named>
#   session_rows=<session streams named> word=<last settled state> <ledger>
for root in "$@"; do
  find "$root" -path '*/.chemsmart-agent/goals/*/ledger.jsonl' \
    -not -path '*sealed*' -not -path '*/claude/*' -not -path '*private*' \
    2>/dev/null
done | sort | while read -r L; do
  ev=$(grep '"kind": "analysis_evidence_recorded"' "$L" \
    | grep -o '"evidence": "[^"]*"' | sed 's/"evidence": "//; s/"$//; s#^runs/##' | sort -u)
  [ -z "$ev" ] && continue
  ss=$(grep '"kind": "session_stream_recorded"' "$L" \
    | grep -o '"run_id": "[^"]*"' | sed 's/"run_id": "//; s/"$//' | sort -u)
  n_ev=$(printf '%s\n' "$ev" | grep -c .)
  n_ss=$(printf '%s\n' "$ss" | grep -c .)
  n_named=$(comm -12 <(printf '%s\n' "$ev") <(printf '%s\n' "$ss") | grep -c .)
  first=$(head -1 "$L" | grep -o '"at": "[^"]*"' | head -1 | cut -d'"' -f4 | cut -c1-10)
  word=$(grep '"kind": "goal_settled"' "$L" | tail -1 | grep -o '"state": "[^"]*"' | tail -1 | cut -d'"' -f4)
  echo "$first ev=$n_ev named=$n_named session_rows=$n_ss word=${word:-unsettled} $L"
done
