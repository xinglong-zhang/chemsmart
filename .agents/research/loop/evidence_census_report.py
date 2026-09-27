"""Census E, read: the tallies EPISODE.md quotes, recomputed from outputs.

``evidence_census.py`` writes rows beside the records; this reads them back
(and, for the two questions rows do not carry, the records themselves). It
grades nothing. Every figure is a count with its denominator, and each
subcommand says which surface it reads.

    python evidence_census_report.py predictions OUT_DIR [OUT_DIR ...]
    python evidence_census_report.py gaps OUT_DIR [OUT_DIR ...]
    python evidence_census_report.py routes OUT_DIR [OUT_DIR ...]
    python evidence_census_report.py claim-constants ROOT [ROOT ...]
    python evidence_census_report.py carried-literals ROOT [ROOT ...]

``predictions`` evaluates the pre-registered P1-P6 per output directory and
pooled. ``gaps`` groups the Agent's own statements of a gap by their reason.
``routes`` summarises accepted plans and executed goal routes. The last two
read event streams and public transcripts directly: claims whose value
stands on a model-authored constant, and literal nodes whose value re-types
a host number (in-session) or carries one only the task or wake showed.
"""

from __future__ import annotations

import collections
import importlib.util
import json
import os
import re
import sys
from pathlib import Path

_CENSUS = Path(__file__).with_name("evidence_census.py")
_spec = importlib.util.spec_from_file_location("evidence_census", _CENSUS)
census = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(census)

DELIVERED = ("final_message", "finding", "unreachable")
PLAN_TOOLS = (
    "tool:plan_scientific_workflow",
    "tool:plan_result_extraction",
    "tool:plan_quantity_expression",
    "tool:plan_thermochemistry",
    "tool:plan_claim_rendering",
    "tool:plan_scientific_validation",
    "tool:plan_unsupported_external",
)
SELECTOR_CLASSES = (
    "selector_undeclared_for_jobtype",
    "selector_not_in_enum",
    "reader_does_not_provide",
    "value_absent_in_result",
)


def infrastructure(record):
    return census.infrastructure(record)


def load(dirs):
    sessions, rows, goals = {}, [], {}
    for d in dirs:
        for line in open(os.path.join(d, "sessions.jsonl"), encoding="utf-8"):
            record = json.loads(line)
            if not infrastructure(record):
                sessions[record["session"]] = record
        for line in open(os.path.join(d, "rows.jsonl"), encoding="utf-8"):
            row = json.loads(line)
            if (
                row.get("session") in sessions
                or row.get("detector") == "project_file"
            ):
                rows.append(row)
        path = os.path.join(d, "goals.json")
        if os.path.exists(path):
            goals.update(json.load(open(path, encoding="utf-8")))
    return sessions, rows, goals


def pct(count, total):
    return f"{count} ({100.0 * count / (total or 1):.1f}%)"


def prose_leaf(row):
    """Prose-ness from the leaf's own key, with the census's current table."""

    leaf = [
        p for p in str(row.get("leaf") or "").split(".") if p and p != "[]"
    ]
    return bool(leaf) and leaf[-1] in census.PROSE_KEYS


def predictions(dirs):
    sessions, rows, _goals = load(dirs)
    n = len(sessions)
    calls = sum(
        (s.get("counts") or {}).get("tool_calls", 0) for s in sessions.values()
    )
    hatch = [r for r in rows if r["detector"] == "hatch"]
    hatch_sessions = {r["session"] for r in hatch}
    hatch_calls = {
        (r["session"], r["message_index"], r["call_index"]) for r in hatch
    }
    operative = [
        r
        for r in rows
        if r["detector"] == "string_exit"
        and not prose_leaf(r)
        and not str(r.get("leaf", "")).startswith("sections.")
    ]
    operative_calls = {
        (r["session"], r["message_index"], r["call_index"]) for r in operative
    }
    selector = [
        r
        for r in rows
        if r["detector"] == "refusal" and r.get("class") in SELECTOR_CLASSES
    ]
    places = collections.defaultdict(set)
    for r in selector:
        for name, where in r.get("placed") or ():
            places[name].add(str(where))
    nowhere = sorted(
        name for name, where in places.items() if where == {"declared_nowhere"}
    )
    planners = {
        s
        for s, record in sessions.items()
        if any(tool in (record.get("counts") or {}) for tool in PLAN_TOOLS)
    }
    gaps = {r["session"] for r in rows if r["detector"] == "declared_gap"}
    reextract = {
        r["session"]
        for r in rows
        if r["detector"] == "citation"
        and r["origin"] == "context"
        and r["tool"] == "extract_result_quantities"
        and r["reply_status"] == "ok"
    }
    cited = collections.defaultdict(set)
    for r in rows:
        if r["detector"] == "citation":
            cited[(r["session"], r["message_index"], r.get("call_index"))].add(
                r["origin"]
            )
    unknown_inputs = [
        r
        for r in rows
        if r["detector"] == "refusal"
        and r.get("class") == "receipt_unknown"
        and r["tool"]
        in ("evaluate_quantity_expression", "record_analysis_claims")
    ]
    cross_run = [
        r
        for r in unknown_inputs
        if cited.get((r["session"], r["message_index"], r.get("call_index")))
    ]
    real, control = collections.Counter(), collections.Counter()
    for r in rows:
        if r["detector"] == "number" and r["surface"] in DELIVERED:
            real[str(r.get("class")).split(":")[0]] += 1
            control[str(r.get("control_class")).split(":")[0]] += 1
    total, ctotal = sum(real.values()), sum(control.values())
    host = ("bound", "bound_other", "host_text")
    print(f"behavioural sessions {n}; tool calls {calls}")
    print(
        f"P1 sessions with a hatch or free-word key {pct(len(hatch_sessions), n)};"
        f" hatch calls {len(hatch_calls)}"
    )
    print(
        f"P2 calls with an operative native/path/shell string outside project"
        f" sections {pct(len(operative_calls), calls)}:"
        f" {dict(collections.Counter((r['tool'], r['pattern'], r['reply_status']) for r in operative))}"
    )
    print(
        f"P3 selector refusals {len(selector)}; distinct names {len(places)};"
        f" declared nowhere {nowhere}"
    )
    print(
        f"P4 analysis-planning sessions {len(planners)}; with an Agent-declared gap"
        f" {pct(len(gaps & planners), len(planners))}"
    )
    print(
        f"P5 sessions re-extracting a context-shown artifact {len(reextract)};"
        f" expression/claim inputs refused as unknown receipts {len(unknown_inputs)},"
        f" of which citing a digest no earlier reply carried {len(cross_run)}"
    )
    print(
        f"P6 delivered numbers {total}: typed record {pct(real['bound'], total)};"
        f" any host value {pct(sum(real[k] for k in host), total)};"
        f" prose-computed {pct(real['prose_computed'], total)};"
        f" context {pct(real['context'], total)};"
        f" own + unmatched {pct(real['own_argument'] + real['unmatched'], total)}"
    )
    print(
        f"   control: typed record {pct(control['bound'], ctotal)};"
        f" any host value {pct(sum(control[k] for k in host), ctotal)};"
        f" prose-computed {pct(control['prose_computed'], ctotal)};"
        f" context {pct(control['context'], ctotal)};"
        f" unmatched {pct(control['unmatched'], ctotal)}"
    )
    by_band = collections.defaultdict(collections.Counter)
    for r in rows:
        if r["detector"] == "number" and r["surface"] in DELIVERED:
            by_band[(r["surface"], r.get("precision"))][
                str(r.get("class")).split(":")[0]
            ] += 1
    for (surface, band), counter in sorted(by_band.items()):
        size = sum(counter.values())
        print(
            f"   {surface} {band}: typed record {pct(counter['bound'], size)} of {size}"
        )


#: A declared gap's reason, typed by the first pattern it matches.
GAP_TYPES = (
    # First: a reason about a held, non-stationary structure also names the
    # job types whose thermochemistry axis is unsupported (r10/q21 g1-hooh),
    # and would otherwise read as the ORCA sp frequency gap.
    ("gate: thermochemistry at a non-stationary structure", r"non-stationary"),
    (
        "vocabulary: adaptive mode selection",
        r"mode selection|adaptive assignment",
    ),
    (
        "vocabulary: ORCA sp frequencies",
        r"vibrational_frequencies.*\bsp\b|\bsp\b.*vibrational_frequencies|freq-enabled single point",
    ),
    (
        "vocabulary: frontier orbitals",
        r"\bhomo\b|\blumo\b|frontier[- ]orbital|orbital-eigenvalue",
    ),
    ("vocabulary: atomic spin populations", r"spin[- ]population|spin densit"),
    (
        "vocabulary: IRC path geometry",
        r"irc.*(positions|trajectory|endpoint geometry)|trajectory selectors",
    ),
    ("vocabulary: NEB energies", r"\bneb\b"),
    (
        "vocabulary: energy to wavenumber",
        r"energy[- >]*to[- ]*(cm|wavenumber)|hartree->cm|cm\^-1.*different dimensions",
    ),
    (
        "vocabulary: categorical claim",
        r"categorical|text label|numerical claims only",
    ),
    ("vocabulary: spatial extent", r"spatial extent|<r\*\*2>|<r\^2>"),
    ("vocabulary: excited-state spin", r"excited[- _]state[- _]spin|per-root"),
    (
        "vocabulary: scan coverage",
        r"scan.*parser coverage|torsional-map extraction",
    ),
    ("settings: counterpoise ghost atoms", r"ghost|counterpoise"),
    ("settings: functional not written", r"mn15|not materiali[sz]able"),
    (
        "constant: literature value",
        r"g_aq\(h\+\)|calibrat|ferrocen|fc\+|fc0|reference constant|absolute potential",
    ),
    (
        "gate: thermochemistry at a non-minimum",
        r"imaginary (harmonic )?(mode|frequenc)|not a (true|validated) minimum|non-stationary",
    ),
    (
        "execution: engine, envelope or preview",
        r"engine binding|execution-unsupported|preview.*cannot pass|mpirun|no executable"
        r"|unsupported_jobtype|capability red|binary is missing|missing from the bound",
    ),
    (
        "evidence: missing producer or input",
        r"no (approved|admissible|semi-empirical|gas-phase)|not registered|absent"
        r"|upstream|consumes the blocked|no .* result|required producer|grant|budget"
        r"|never produced|cannot exist|not materialised|cannot be admitted",
    ),
)


def gap_type(text):
    low = str(text or "").lower()
    for name, pattern in GAP_TYPES:
        if re.search(pattern, low):
            return name
    return "other"


def workspace_label(path):
    for marker in ("ax41-refine-100/", "experiments-public/", "jiseung/"):
        if marker in path:
            path = path.split(marker, 1)[1]
            break
    return path.split("/.chemsmart-agent")[0].split("/workspace/")[0]


def gaps(dirs):
    _sessions, rows, _goals = load(dirs)
    per_type = collections.defaultdict(set)
    for r in rows:
        if r["detector"] != "declared_gap":
            continue
        texts = (
            [r["blocked_reason"]]
            if r.get("blocked_reason")
            else re.findall(
                r'"blocked_reason": "((?:[^"\\]|\\.)*)',
                r.get("arguments_excerpt") or "",
            )
        )
        for text in texts or [""]:
            per_type[gap_type(text)].add(workspace_label(r["path"]))
    for name, where in sorted(
        per_type.items(), key=lambda kv: (-len(kv[1]), kv[0])
    ):
        print(f"{len(where):4d} workspaces  {name}")
        if not name.startswith(("evidence", "gate", "execution", "other")):
            for label in sorted(where):
                print("         ", label)


def routes(dirs):
    sessions, _rows, goals = load(dirs)
    plans = collections.Counter()
    sizes = collections.Counter()
    for record in sessions.values():
        planned = record.get("planned_routes") or []
        if not planned or not (planned[-1].get("nodes") or {}):
            continue
        nodes = planned[-1]["nodes"]
        plans["+".join(sorted({k.split(":")[0] for k in nodes}))] += 1
        sizes[min(sum(nodes.values()), 10)] += 1
    executed = collections.Counter()
    edges = collections.Counter()
    for goal in goals.values():
        route = goal.get("route") or {}
        if not route.get("executed"):
            continue
        executed[
            "+".join(sorted({k.split(":")[0] for k in route["executed"]}))
        ] += 1
        edges.update(
            {
                k: v
                for k, v in (route.get("edges") or {}).items()
                if not k.startswith("analysis:")
            }
        )
    print(
        f"accepted plans with calculation nodes: {sum(plans.values())}; programs {dict(plans.most_common())}"
    )
    print(
        f"  calculation nodes per plan (10 = 10 or more): {dict(sorted(sizes.items()))}"
    )
    print(
        f"goals with executed nodes: {sum(executed.values())}; programs {dict(executed.most_common())}"
    )
    print(f"  recorded edges: {dict(edges)}")


PRUNE_PREFIX = census.PRUNE_PREFIX
PRUNE_EXACT = census.PRUNE_EXACT


def files(roots, predicate):
    for root in roots:
        for dirpath, names in census.walk(root, ()):
            for name in names:
                if predicate(name):
                    yield os.path.join(dirpath, name)


def claim_constants(roots):
    """Claims whose value stands on an expression output that names a
    model-authored constant, joined within each event stream."""

    seen, total, resting, streams = set(), 0, 0, set()
    values = collections.Counter()
    for path in files(roots, lambda name: name == "events.jsonl"):
        constants, claims = {}, []
        for row in census.load_jsonl(path):
            payload = row.get("payload") or {}
            if row.get("kind") == "quantity_expression_evaluated":
                for dep in (payload.get("record") or {}).get(
                    "output_dependencies"
                ) or ():
                    constants[
                        (payload.get("receipt_sha256"), dep.get("output_id"))
                    ] = (dep.get("model_authored_constants") or ())
            elif row.get("kind") == "analysis_claims_recorded":
                if payload.get("receipt_sha256") in seen:
                    continue
                seen.add(payload.get("receipt_sha256"))
                claims.extend(
                    (payload.get("record") or {}).get("claims") or ()
                )
        for claim in claims:
            total += 1
            named = constants.get(
                (claim.get("source_receipt_sha256"), claim.get("quantity_id"))
            )
            if named:
                resting += 1
                streams.add(path)
                for item in named:
                    values[json.dumps(item, sort_keys=True)[:100]] += 1
    print(
        f"claims {total}; whose value names a model-authored constant {pct(resting, total)}; streams {len(streams)}"
    )
    for text, count in values.most_common(25):
        print(f"  {count:4d} {text}")


#: Conditions and physical constants: a literal carrying one of these is the
#: constants registry's business, not a carried result.
COMMON = (
    298.15,
    101325.0,
    1.01325,
    8.314462618,
    8.31446261815324,
    0.008314462618,
    0.00831446261815324,
    0.00198720425864083,
    1.98720425864083,
    627.509474,
    627.5095,
    2625.49964,
    2625.5,
    27.211386,
    27.2114,
    219474.63,
    0.0433641,
    23.060548,
    96.485332,
    4.184,
    1.772453850905516,
    3.141592653589793,
    6.283185307179586,
    0.69503476,
    1.380649e-23,
    6.62607015e-34,
    6.02214076e23,
    1.6605390666e-27,
    2.99792458e10,
    1.007825,
    15.994915,
    14.003074,
    0.529177,
    1.889726,
    1239.84198,
    0.0119627,
    83.593,
    349.755,
    1.4812,
    2.4789,
    0.592,
    0.025609,
    1.9872,
    0.0019872,
    2.302585,
)


def carried_literals(roots):
    """Literal nodes (>= 5 significant digits, not a condition or constant)
    whose value equals a number an earlier typed reply of the session
    carried, or one only the task or wake message carried."""

    def numbers(value):
        return census.numeric_values(value)

    def near(x, pool):
        return any(abs(x - v) <= max(1e-9, 1e-7 * abs(x)) for v in pool)

    digests = set()
    kept = retyped = 0
    carried = []
    for path in files(
        roots,
        lambda n: n.startswith("public-transcript") and n.endswith(".json"),
    ):
        try:
            data = json.load(open(path, encoding="utf-8"))
        except Exception:  # noqa: BLE001
            continue
        if data.get("transcript_sha256") in digests:
            continue
        digests.add(data.get("transcript_sha256"))
        context, typed = [], []
        for message in data.get("transcript") or []:
            role = message.get("role")
            parsed = census.parse_json(message.get("content"))
            if role == "user":
                context.extend(
                    numbers(parsed)
                    if parsed is not None
                    else census.numbers_in_text(str(message.get("content")))
                )
            elif role == "tool" and parsed is not None:
                typed.extend(numbers(parsed))
            elif role == "assistant":
                for call in message.get("tool_calls") or []:
                    arguments = census.parse_json(
                        (call.get("function") or {}).get("arguments") or "{}"
                    )
                    for path_, leaf in census.leaves(
                        arguments if isinstance(arguments, dict) else {}
                    ):
                        if not path_ or path_[-1] != "literal_value":
                            continue
                        if isinstance(leaf, bool) or not isinstance(
                            leaf, (int, float)
                        ):
                            continue
                        x = float(leaf)
                        if x == 0 or any(
                            abs(abs(x) - c) <= 1e-6 * c for c in COMMON
                        ):
                            continue
                        if census.significant_digits(repr(abs(x))) < 5:
                            continue
                        kept += 1
                        if near(x, typed):
                            retyped += 1
                        elif near(x, context):
                            carried.append((workspace_label(path), x))
    print(
        f"literal values (>= 5 significant digits, not a condition or constant): {kept}"
    )
    print(
        f"  equal to a number an earlier typed reply of the session carried: {retyped}"
    )
    print(
        f"  equal only to a number the task or wake message carried: {len(carried)}"
    )
    for label, x in sorted(set(carried)):
        print(f"    {label} = {x}")


def main(argv):
    if len(argv) < 3:
        print(__doc__)
        return 2
    command, targets = argv[1], argv[2:]
    handlers = {
        "predictions": predictions,
        "gaps": gaps,
        "routes": routes,
        "claim-constants": claim_constants,
        "carried-literals": carried_literals,
    }
    if command not in handlers:
        print(__doc__)
        return 2
    if command == "predictions":
        for target in targets:
            print("==", target.rstrip("/").split("/")[-1])
            predictions([target])
        if len(targets) > 1:
            print("== pooled")
            predictions(targets)
        return 0
    handlers[command](targets)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
