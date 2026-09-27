"""Typed hypothesis and noticing acts per archived session, by model (scratch).

Unit: a live-session stream (runs/live-*/events.jsonl) with at least one
provider turn whose observed model is known; the session's model is the
most common observed model in it. Acts, from the host's own events:

  DECLARED  any requested_observable_declared
  HYP       a declared observable with role "diagnostic" and a non-empty
            failure_update_rule
  DECISION  any scientific_decision_recorded
  FINDING   a decision carrying non-empty "findings" anywhere in its payload
  UNREACH   a decision carrying non-empty "unreachable_observables" or
            "unreachable_observable_ids" anywhere in its payload

Only sessions on or after the first day any session declared a diagnostic
are compared (the affordance existed); rows split by model and by the
campaign directory two levels under the root.

    python hyp_census.py ROOT [ROOT ...]
"""

import json
import os
import sys
from collections import Counter, defaultdict


def streams(root):
    for base, _dirs, files in os.walk(root):
        if "events.jsonl" in files and "/runs/live-" in base + "/":
            yield os.path.join(base, "events.jsonl")


def found(value, keys):
    if isinstance(value, dict):
        for key, item in value.items():
            if key in keys and item not in (None, "", [], {}, ()):
                return True
            if found(item, keys):
                return True
    elif isinstance(value, list):
        return any(found(item, keys) for item in value)
    return False


def census(path):
    models = Counter()
    first = ""
    acts = dict.fromkeys(("DECLARED", "HYP", "DECISION", "FINDING", "UNREACH"), False)
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip():
                continue
            try:
                event = json.loads(line)
            except ValueError:
                continue
            first = first or str(event.get("timestamp") or "")[:10]
            kind = event.get("kind")
            payload = event.get("payload") or {}
            if kind == "provider_turn_observed":
                model = str(payload.get("observed_model") or "")
                if model and model != "not_observed":
                    models[model] += 1
            elif kind == "requested_observable_declared":
                acts["DECLARED"] = True
                for item in payload.get("observables") or ():
                    if (
                        isinstance(item, dict)
                        and item.get("role") == "diagnostic"
                        and str(item.get("failure_update_rule") or "").strip()
                    ):
                        acts["HYP"] = True
            elif kind == "scientific_decision_recorded":
                acts["DECISION"] = True
                acts["FINDING"] |= found(payload, {"findings"})
                acts["UNREACH"] |= found(
                    payload, {"unreachable_observables", "unreachable_observable_ids"}
                )
    if not models:
        return None
    return models.most_common(1)[0][0], first, acts


def main(roots):
    rows = []
    seen = set()
    for root in roots:
        for path in streams(root):
            real = os.path.realpath(path)
            if real in seen:
                continue
            seen.add(real)
            result = census(path)
            if result is None:
                continue
            model, day, acts = result
            rel = os.path.relpath(path, root).split(os.sep)
            family = "/".join(rel[:2]) if len(rel) > 2 else rel[0]
            rows.append((model, day, os.path.basename(root.rstrip("/")) + ":" + family, acts))
    start = min((day for _m, day, _f, acts in rows if acts["HYP"]), default="")
    print("first day any session declared a diagnostic:", start)
    names = ("sessions", "DECLARED", "HYP", "DECISION", "FINDING", "UNREACH")
    by_family = defaultdict(Counter)
    totals = defaultdict(Counter)
    for model, day, family, acts in rows:
        if day < start:
            continue
        for key, bucket in (((model, family), by_family), ((model,), totals)):
            bucket[key]["sessions"] += 1
            for name, present in acts.items():
                bucket[key][name] += int(present)
    print("model family | " + " ".join(names))
    for key, counts in sorted(by_family.items()):
        print(" ".join(key) + " | " + " ".join(str(counts[n]) for n in names))
    for key, counts in sorted(totals.items()):
        print("TOTAL " + key[0] + " | " + " ".join(str(counts[n]) for n in names))


if __name__ == "__main__":
    main(sys.argv[1:])
