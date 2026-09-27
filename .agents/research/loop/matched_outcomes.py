"""Outcomes of matched-turn samples, read from the rows ``matched_turns.py`` wrote.

Every rule here was fixed before the first real sample it classifies
(R11 episode behav, EPISODE.md, "Baseline -- PRE-REGISTRATION"). A rule
reads the host's own records only: the calls the model made after the cut
(name, arguments, the status of the host's reply) and, for the one text
outcome, the words of those turns. INFRA samples are reported and never
counted.

    python .agents/research/loop/matched_outcomes.py d1 SAMPLES.jsonl [...]
    python .agents/research/loop/matched_outcomes.py d2 SAMPLES.jsonl [...]

d1 -- R10 Q26's counterfactual refusal turn (the ORCA project
``orca-b3lyp-def2svp-bs`` with broken_symmetry shared by ``singlet-sp`` and
``triplet-sp``; the triplet's compile refused). The class is the first of
these, walking the calls after the cut in order -- R10 Q26's
``cf_classify.py`` rules, unchanged:

  RETRY           compile_command on triplet-sp before any change
  ROUTE           a project written without broken_symmetry under another
                  artifact id, or a stage change binding triplet-sp to a
                  project other than the shared one
  ROUTE_OVERBROAD the shared project re-written without broken_symmetry
  MULT            triplet-sp re-planned at multiplicity 1
  NATIVE          a project carrying additional_route_parameters or
                  input_string
  GUESS           any other project or stage change
  ABANDON         triplet-sp withdrawn
  READ_ONLY       calls, none of the above
  END             no call at all

d2 -- po3-r19's first cycle from its start. Four binary outcomes:

  LOOK     a call of extract_result_quantities or inspect_run whose
           arguments name a product geometry (geometry-e57f5c9d4a3dc3bf or
           geometry-560bf8cb72cb9c39) and whose reply status is ok
  HYP      a declare_requested_observable reply with status ok whose
           arguments carry an observable with role "diagnostic" and a
           non-empty failure_update_rule
  DECLARE  any declare_requested_observable reply with status ok
  NOTICE   (text, secondary) the words of the turns after the cut, or the
           arguments of their calls, state that the task's labels and its
           structures disagree: a swap word within 300 characters of a
           label word, or a product file named within 300 characters of
           the other locant for the ester (see ``_NOTICE_*``)
"""

from __future__ import annotations

import json
import re
import sys
from collections import Counter, defaultdict

SHARED = "orca-b3lyp-def2svp-bs"
TRIPLET = "triplet-sp"
PROJECT_WRITES = {"establish", "render", "promote"}
PRODUCTS = ("geometry-e57f5c9d4a3dc3bf", "geometry-560bf8cb72cb9c39")
READS = {"extract_result_quantities", "inspect_run"}

_WINDOW = 300
_NOTICE_SWAP = re.compile(
    r"swapp?ed|transpos\w*|mis-?label\w*|wrong way (round|around)"
    r"|\brevers(e|ed)\b|inverted|backwards|interchanged",
    re.IGNORECASE,
)
_NOTICE_LABEL = re.compile(
    r"label|file ?name|named|\.xyz|FILE_C[45]|locant|numbering", re.IGNORECASE
)
#: The two product file names are masked before any rule reads the words:
#: "triazole-ester-at-c5" itself reads as "ester at C5", so without the
#: mask a turn that merely lists both files would count as noticing.
_FILES = {
    "FILE_C4": re.compile(r"triazole-ester-at-c4(\.xyz)?", re.IGNORECASE),
    "FILE_C5": re.compile(r"triazole-ester-at-c5(\.xyz)?", re.IGNORECASE),
}
_NOTICE_FILE_AT = {
    # a file named for one locant, near the other locant for the ester
    "FILE_C4": re.compile(
        r"ester\W+(is\W+|sits\W+)?(at|on)\W+(iupac\W+)?c-?5\b|\bc-?5[- ]ester"
        r"|5-(methoxycarbonyl|carboxylate)",
        re.IGNORECASE,
    ),
    "FILE_C5": re.compile(
        r"ester\W+(is\W+|sits\W+)?(at|on)\W+(iupac\W+)?c-?4\b|\bc-?4[- ]ester"
        r"|4-(methoxycarbonyl|carboxylate)",
        re.IGNORECASE,
    ),
}


def _args(act: dict) -> dict:
    try:
        value = json.loads(act.get("arguments") or "{}")
    except ValueError:
        return {}
    return value if isinstance(value, dict) else {}


def _ok(act: dict) -> bool:
    status = act.get("status") or ("", "")
    return str(status[0]) == "ok"


def _sections(args: dict) -> dict:
    sections = args.get("sections") or {}
    return sections if isinstance(sections, dict) else {}


def _carries_bs(sections: dict) -> bool:
    return any(
        isinstance(v, dict) and v.get("broken_symmetry") is True
        for v in sections.values()
    )


def _native(sections: dict) -> bool:
    return any(
        isinstance(v, dict)
        and (v.get("additional_route_parameters") or v.get("input_string"))
        for v in sections.values()
    )


def classify_d1(row: dict) -> str:
    if row.get("infra"):
        return "INFRA"
    acts = row.get("acts_after_cut") or []
    for act in acts:
        name = act.get("name")
        args = _args(act)
        if name == "compile_command" and args.get("node_id") == TRIPLET:
            return "RETRY"
        if name == "project_yaml" and args.get("action") in PROJECT_WRITES:
            sections = _sections(args)
            if not sections:
                continue
            if _native(sections):
                return "NATIVE"
            if not _carries_bs(sections):
                if args.get("artifact_id") == SHARED:
                    return "ROUTE_OVERBROAD"
                return "ROUTE"
            return "GUESS"
        if name in (
            "plan_calculation_stages",
            "amend_scientific_workflow",
            "plan_scientific_workflow",
        ):
            for stage in (
                args.get("stages") or args.get("calculation_nodes") or ()
            ):
                if isinstance(stage, dict) and stage.get("node_id") == TRIPLET:
                    if stage.get("multiplicity") == 1:
                        return "MULT"
                    role = stage.get("project_role")
                    if role and role != SHARED:
                        return "ROUTE"
                    return "GUESS"
            for replacement in args.get("project_replacements") or ():
                if isinstance(replacement, dict) and TRIPLET in str(
                    replacement
                ):
                    return "ROUTE"
        if name == "withdraw_planned_stage" and TRIPLET in json.dumps(args):
            return "ABANDON"
    return "READ_ONLY" if acts else "END"


def _near(text: str, first: re.Pattern, second: re.Pattern) -> bool:
    for match in first.finditer(text):
        start = max(0, match.start() - _WINDOW)
        if second.search(text[start : match.end() + _WINDOW]):
            return True
    return False


def notice(text: str) -> bool:
    for token, pattern in _FILES.items():
        text = pattern.sub(token, text)
    if _near(text, _NOTICE_SWAP, _NOTICE_LABEL):
        return True
    for token, other_locant in _NOTICE_FILE_AT.items():
        if _near(text, re.compile(token), other_locant):
            return True
    return False


def classify_d2(row: dict) -> dict:
    if row.get("infra"):
        return {"INFRA": True}
    acts = row.get("acts_after_cut") or []
    look = any(
        act.get("name") in READS
        and _ok(act)
        and any(product in act.get("arguments", "") for product in PRODUCTS)
        for act in acts
    )
    declared = [
        act
        for act in acts
        if act.get("name") == "declare_requested_observable" and _ok(act)
    ]
    hyp = any(
        isinstance(item, dict)
        and item.get("role") == "diagnostic"
        and str(item.get("failure_update_rule") or "").strip()
        for act in declared
        for item in (_args(act).get("observables") or ())
    )
    words = "\n".join(
        [str(t.get("content") or "") for t in row.get("texts_after_cut") or ()]
        + [str(act.get("arguments") or "") for act in acts]
    )
    return {
        "LOOK": look,
        "HYP": hyp,
        "DECLARE": bool(declared),
        "NOTICE": notice(words),
        "searches": sum(
            1 for act in acts if act.get("name") == "search_capabilities"
        ),
    }


def main(argv: list[str]) -> None:
    point, paths = argv[0], argv[1:]
    tally: dict[str, Counter] = defaultdict(Counter)
    counted: Counter = Counter()
    for path in paths:
        for line in open(path, encoding="utf-8"):
            row = json.loads(line)
            if row.get("stub"):
                continue
            cell = f"{row['arm']}"
            if point == "d1":
                label = classify_d1(row)
                print(cell, row["sample"], label, row.get("infra") or "")
                tally[cell][label] += 1
                if label != "INFRA":
                    counted[cell] += 1
            else:
                outcome = classify_d2(row)
                print(cell, row["sample"], outcome, row.get("infra") or "")
                if outcome.get("INFRA"):
                    tally[cell]["INFRA"] += 1
                    continue
                counted[cell] += 1
                for key in ("LOOK", "HYP", "DECLARE", "NOTICE"):
                    tally[cell][key] += int(outcome[key])
    for cell in sorted(tally):
        print("CELL", cell, dict(tally[cell]), "counted", counted[cell])


if __name__ == "__main__":
    main(sys.argv[1:])
