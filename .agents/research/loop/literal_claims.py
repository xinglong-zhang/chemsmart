"""Which claims in a host-rendered report rest on a model-authored literal?

    python .agents/research/loop/literal_claims.py OUT_DIR ROOT [ROOT ...]

Every completed or partial analysis report under the roots (sealed,
private and ``claude/`` paths excluded) renders a table of "host-rendered
numerical claims": value, unit, source receipt. For each row the run's own
stream is read: the claim record gives the quantity the claim carries, and
an expression receipt's ``output_dependencies`` say, per output, which
model-authored constants (``literal`` nodes) it stands on; the literal's
own node value gives its dimension. Each row is classed:

- host: the source is an extraction or a thermochemistry receipt, or an
  expression output that stands on no model-authored constant;
- count: it stands only on dimensionless literals (a coefficient, an
  electron count);
- condition: it also stands on a temperature or a pressure literal;
- physical: it stands on a literal with any other dimension -- an energy,
  a potential, a length: a number the model supplied. Such a literal is
  "carried" when its value equals, to 1e-9 relative, a quantity an
  extraction or a thermochemistry receipt of the same workspace minted
  (a host number typed again), and "foreign" otherwise (a literature
  value, an estimate).

A row is also "pure" when its output does no arithmetic on any receipt
(``arithmetic_node_count`` 0 and no source receipts): the value the table
shows is the model's own number. Imports nothing from ``chemsmart``.
Written for R11 episode `truth` (Phase II item 5).
"""

from __future__ import annotations

import collections
import json
import re
import sys
from pathlib import Path

EXCLUDED = ("sealed", "/private", "/claude/", "/grade-", "/grading/")
_ROW = re.compile(
    r"^\|\s*([^|]+?)\s*\|\s*`([^`]*)`\s*\|\s*([^|]*?)\s*\|\s*`([0-9a-f]{64})`\s*\|"
)
_TEMPERATURE = (0, 0, 1, 0, 0, 0)


def rows(path: Path) -> list[dict]:
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return []
    out = []
    for line in text.splitlines():
        if line.strip():
            try:
                out.append(json.loads(line))
            except json.JSONDecodeError:
                continue
    return out


def host_values(agent: Path | None, run_stream: Path) -> list[float]:
    """Every number an extraction or thermochemistry receipt minted in the
    run's stream and in the workspace's recorded streams."""

    streams = {run_stream}
    if agent is not None:
        streams.update(agent.glob("runs/*/events.jsonl"))
        streams.update(agent.glob("goals/*/runs/*/events.jsonl"))
    values = []
    for stream in streams:
        for event in rows(stream):
            if event.get("kind") not in {
                "result_quantities_extracted",
                "thermochemistry_derived",
            }:
                continue
            record = (event.get("payload") or {}).get("record") or {}
            for quantity in record.get("quantities") or ():
                value = quantity.get("value")
                if isinstance(value, (int, float)) and not isinstance(
                    value, bool
                ):
                    values.append(float(value))
    return values


def carried(value: float, pool: list[float]) -> bool:
    return any(
        abs(value - other) <= max(1e-12, 1e-9 * abs(other)) for other in pool
    )


def main() -> None:
    out_dir = Path(sys.argv[1])
    out_dir.mkdir(parents=True, exist_ok=True)
    reports = []
    for root in sys.argv[2:]:
        for pattern in (
            "completed-analysis-report.md",
            "partial-analysis-report.md",
        ):
            for path in Path(root).rglob(pattern):
                if any(part in str(path) for part in EXCLUDED):
                    continue
                reports.append(path)
    results = []
    counts: collections.Counter = collections.Counter()
    for report in sorted(set(reports)):
        run = report.parent.parent
        stream = run / "events.jsonl"
        events = rows(stream)
        agent = next(
            (p for p in report.parents if p.name == ".chemsmart-agent"), None
        )
        expressions = {}
        claims = {}
        for event in events:
            payload = event.get("payload") or {}
            if event.get("kind") == "quantity_expression_evaluated":
                expressions[str(payload.get("receipt_sha256") or "")] = (
                    payload.get("record") or {}
                )
            elif event.get("kind") == "analysis_claims_recorded":
                for claim in (payload.get("record") or {}).get("claims") or ():
                    claims[str(claim.get("claim_id") or "")] = claim
        pool = None
        counts["reports"] += 1
        for line in report.read_text(errors="replace").splitlines():
            match = _ROW.match(line)
            if not match:
                continue
            claim_id, value, unit, source = match.groups()
            counts["rows"] += 1
            row = {
                "report": str(report),
                "claim_id": claim_id,
                "value": value,
                "unit": unit,
                "source": source[:12],
            }
            record = expressions.get(source)
            claim = claims.get(claim_id) or {}
            if record is None:
                row["class"] = "host"
                counts["host"] += 1
                results.append(row)
                continue
            output_id = str(claim.get("quantity_id") or "")
            dependency = next(
                (
                    d
                    for d in record.get("output_dependencies") or ()
                    if str(d.get("output_id") or "") == output_id
                ),
                None,
            )
            if dependency is None:
                row["class"] = "unread: no dependency row for the output"
                counts["unread"] += 1
                results.append(row)
                continue
            constants = dependency.get("model_authored_constants") or ()
            nodes = {
                str(n.get("quantity_id") or ""): n
                for n in record.get("node_values") or ()
            }
            pure = not dependency.get("arithmetic_node_count") and not (
                dependency.get("source_receipt_sha256s")
            )
            kinds = []
            literals = []
            for constant in constants:
                node = nodes.get(str(constant.get("node_id") or "")) or {}
                dimension = tuple(int(x) for x in node.get("dimension") or ())
                padded = dimension + (0,) * (6 - len(dimension))
                if not any(padded):
                    kinds.append("count")
                elif padded[:6] == _TEMPERATURE:
                    kinds.append("condition")
                else:
                    if pool is None:
                        pool = host_values(agent, stream)
                    value_c = node.get("value")
                    is_carried = isinstance(value_c, (int, float)) and carried(
                        float(value_c), pool
                    )
                    kinds.append("physical")
                    literals.append(
                        {
                            "node": constant.get("node_id"),
                            "value": node.get("source_value"),
                            "unit": node.get("source_unit"),
                            "carried": is_carried,
                        }
                    )
            if "physical" in kinds:
                klass = "physical"
            elif "condition" in kinds:
                klass = "condition"
            elif kinds:
                klass = "count"
            else:
                klass = "host"
            row.update(
                {
                    "class": klass,
                    "pure": pure,
                    "literals": literals,
                    "expression": record.get("expression_id"),
                    "output": output_id,
                }
            )
            counts[klass] += 1
            if pure:
                counts["pure"] += 1
            for literal in literals:
                counts[
                    "physical literal "
                    + ("carried" if literal["carried"] else "foreign")
                ] += 1
            results.append(row)
    (out_dir / "literal_claims.json").write_text(json.dumps(results, indent=1))
    for key, n in sorted(counts.items()):
        print(f"{n:6d} {key}")
    for row in results:
        if row.get("class") == "physical" or row.get("pure"):
            print(
                f"  {row['claim_id'][:34]:34} {row['value'][:14]:>14} "
                f"{row['unit'][:9]:9} pure={row.get('pure')} "
                + "; ".join(
                    f"{lit['node']}={lit['value']} {lit['unit']}"
                    f"{' (carried)' if lit['carried'] else ''}"
                    for lit in row.get("literals") or ()
                )[:140]
                + f"  <- {Path(row['report']).parts[-5]}/{Path(row['report']).parts[-4]}"
            )


if __name__ == "__main__":
    main()
