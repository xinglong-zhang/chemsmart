"""The registered dissent replication's blind packet (R11 behav, rule R-DISSENT).

Collects, for every D2 sample of the baseline (d2-head-*) and of arm M
(d2-armM-*), the declared observables that carry an expectation (a non-empty
expected_sign, or both band ends), from the successful
declare_requested_observable calls after the cut -- host records only. Each
sample gets a random code; the packet prints each code's observables
(meaning, sign, band, basis, role) and nothing about the model, the arm or
the sample; the key (code -> cell, sample) is written to a file and never
printed. A sample with no such observable still gets a code.

    python dissent_packet.py OUT_DIR

writes OUT_DIR/dissent-packet.txt and OUT_DIR/dissent-key.json.
"""

import json
import random
import sys
from pathlib import Path

S = Path(
    "/private/tmp/claude-501/-Users-hongjiseung-developer-chemsmart/"
    "c946a1a6-a959-400c-ba03-fd5f4544cb43/scratchpad/behav"
)
CELLS = (
    ("baseline", "d2-head-deepseek"),
    ("baseline", "d2-head-qwen"),
    ("arm-m", "d2-armM-deepseek"),
    ("arm-m", "d2-armM-qwen"),
)
SEED = 20260928


def expectations(row):
    found = []
    for act in row.get("acts_after_cut") or ():
        if act.get("name") != "declare_requested_observable":
            continue
        if (act.get("status") or ["", ""])[0] != "ok":
            continue
        try:
            args = json.loads(act.get("arguments") or "{}")
        except ValueError:
            continue
        for item in args.get("observables") or ():
            if not isinstance(item, dict):
                continue
            banded = (
                item.get("expected_low") is not None
                and item.get("expected_high") is not None
            )
            if not (str(item.get("expected_sign") or "").strip() or banded):
                continue
            found.append(item)
    return found


def main(out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    items = []
    for root, cell in CELLS:
        path = S / root / cell / "samples.jsonl"
        for line in path.read_text().splitlines():
            row = json.loads(line)
            if row.get("infra"):
                continue
            items.append((cell, row["sample"], expectations(row)))
    random.Random(SEED).shuffle(items)
    key = {}
    lines = []
    for number, (cell, sample, found) in enumerate(items):
        code = f"D{number:02d}"
        key[code] = [cell, sample]
        lines.append(
            f"===== {code} ({len(found)} observables with an expectation)"
        )
        for item in found:
            lines.append(
                f"  - id={item.get('observable_id')} role={item.get('role') or 'requested'}"
                f" sign={item.get('expected_sign')}"
                f" band=[{item.get('expected_low')}, {item.get('expected_high')}]"
                f" unit={item.get('unit')}"
            )
            lines.append(f"    meaning: {item.get('meaning')}")
            lines.append(f"    basis: {item.get('expectation_basis')}")
            if item.get("failure_update_rule"):
                lines.append(
                    f"    if falsified: {item.get('failure_update_rule')}"
                )
    (out / "dissent-packet.txt").write_text("\n".join(lines) + "\n")
    (out / "dissent-key.json").write_text(json.dumps(key, indent=1))
    print("codes:", len(key))


if __name__ == "__main__":
    main(sys.argv[1])
