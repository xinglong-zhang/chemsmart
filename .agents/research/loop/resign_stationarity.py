"""Re-sign the stationarity words of archived goals on the imported tree.

    HOME=<fence> PYTHONPATH=<tree> python \
        .agents/research/loop/resign_stationarity.py OUT_DIR label=ROOT [...]

Two classes of word say what a structure is (``signed_words.py`` W10, W11):
a stationary-point characterisation (an order, on a stationary structure)
and a free energy (defined at a stationary point of a named surface). Both
are pure functions of the program's result file, so they can be signed
again on any tree: for every ``stationary_point_characterised`` and every
``thermochemistry_derived`` receipt a goal recorded, the result file it
names (by digest, through the goal's verified results; a recorded path is
mapped onto the archive's copy of the workspace and its digest checked) is
opened with the imported tree's reader and asked again --
``build_stationary_point_characterisation`` for an order,
``free_energy_surface`` for a free energy. A result that is absent or whose
bytes differ is reported, never read.

For each free energy the census also says whether a claim rendered it (a
claim whose source receipt is the thermochemistry receipt, or an
expression that bound it where the expression recorded its inputs), and
the goal's settlement word, so a free energy the imported tree would refuse
is counted where it was delivered. Written for R11 episode `truth`.
"""

from __future__ import annotations

import collections
import hashlib
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from word_reader import Goal, discover  # noqa: E402


def recorded_prefix(goal: Goal) -> str:
    """The workspace root the goal's records were written under."""

    for entry in goal.ledger:
        recorded = str((entry.get("payload") or {}).get("approval_file") or "")
        if "/.chemsmart-agent/" in recorded:
            return recorded.split("/.chemsmart-agent/", 1)[0]
    for _stamp, _role, _path, event in goal.events:
        text = json.dumps(event.get("payload") or {})
        marker = text.find("/.chemsmart-agent/")
        if marker > 0:
            start = text.rfind('"', 0, marker)
            return text[start + 1 : marker]
    return ""


def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    from chemsmart.agent._contracts import TrustedArtifactRefV1
    from chemsmart.agent.execution import build_stationary_point_characterisation
    from chemsmart.analysis.result_quantities import free_energy_surface
    from chemsmart.analysis.result_readers import reader_for
    import chemsmart

    print("chemsmart imported from", chemsmart.__file__, flush=True)
    out_dir = Path(sys.argv[1])
    out_dir.mkdir(parents=True, exist_ok=True)
    rows_out: list[dict] = []
    for item in discover(sys.argv[2:]):
        goal = Goal(Path(item["agent"]), item["goal_id"])
        prefix = recorded_prefix(goal)
        local_root = goal.agent.parent
        artifacts: dict[str, dict] = {}
        for _s, _r, _p, payload in goal.of("program_result_verified"):
            record = payload.get("record") or {}
            for artifact in record.get("output_artifacts") or ():
                artifacts.setdefault(str(artifact.get("sha256") or ""), {
                    "path": str(artifact.get("path") or ""),
                    "kind": str(artifact.get("kind") or ""),
                    "size": int(artifact.get("size_bytes") or 0),
                    "node": str(payload.get("node_id") or record.get("node_id") or ""),
                    "verified": str(payload.get("status") or ""),
                })
        claimed_from: set[str] = set()
        for _stamp, claim in goal.claims():
            claimed_from.add(str(claim.get("source_receipt_sha256") or ""))
        for _s, _r, _p, payload in goal.of("quantity_expression_evaluated"):
            if str(payload.get("receipt_sha256") or "") in claimed_from:
                for binding in (payload.get("record") or {}).get("input_bindings") or payload.get("input_bindings") or ():
                    if isinstance(binding, dict):
                        claimed_from.add(str(binding.get("source_receipt_sha256") or ""))

        def resolve(sha: str) -> tuple[Path | None, dict, str]:
            meta = artifacts.get(sha)
            if meta is None:
                return None, {}, "digest names no verified result"
            recorded = meta["path"]
            candidate = Path(recorded)
            if prefix and recorded.startswith(prefix):
                candidate = local_root / recorded[len(prefix):].lstrip("/")
            if not candidate.is_file():
                return None, meta, "file absent"
            if sha256_of(candidate) != sha:
                return None, meta, "bytes differ"
            return candidate, meta, ""

        seen: set[str] = set()
        for _s, _r, _p, payload in goal.of("thermochemistry_derived", "stationary_point_characterised"):
            receipt = str(payload.get("receipt_sha256") or "")
            if receipt in seen:
                continue
            seen.add(receipt)
            record = payload.get("record") or {}
            characterisation = "order_claimed" in record
            sha = str(record.get("result_artifact_sha256") if characterisation else (payload.get("artifact_sha256") or record.get("artifact_sha256") or ""))
            program = str(record.get("program") or "").lower()
            row = {
                "label": item["label"], "goal_id": item["goal_id"], "word": goal.word,
                "class": "W10" if characterisation else "W11",
                "receipt": receipt[:12], "program": program,
                "delivered": receipt in claimed_from,
            }
            path, meta, why = resolve(sha)
            row["node"] = meta.get("node", "")
            row["verified"] = meta.get("verified", "")
            if path is None:
                row["pin"] = "unread: " + why
                rows_out.append(row)
                continue
            reader = reader_for(program)
            try:
                if characterisation:
                    ref = TrustedArtifactRefV1(
                        artifact_id="archived-result", kind=reader.artifact_kind,
                        sha256=sha, size_bytes=path.stat().st_size,
                        path=str(path.resolve()), cli_value=str(path.resolve()),
                    )
                    build_stationary_point_characterisation(
                        result_artifact=ref, program=program,
                        order_claimed=int(record.get("order_claimed")),
                    )
                    row["pin"] = "certified"
                else:
                    output = reader.open_output(str(path))
                    surface = free_energy_surface(program, output)
                    row["pin"] = f"surface {surface.surface}"
                    row["stationarity"] = surface.stationarity.stationarity
                    row["basis"] = str(getattr(surface.stationarity, "basis", ""))
            except Exception as exc:  # noqa: BLE001 - the refusal is the finding
                gate = getattr(exc, "gate", "")
                row["pin"] = "refused" + (f" ({gate})" if gate else "")
                row["reason"] = str(exc)[:300]
            rows_out.append(row)
    (out_dir / "stationarity.json").write_text(json.dumps(rows_out, indent=1))
    counts = collections.Counter(
        (r["class"], r["pin"].split(":")[0] if r["pin"].startswith("unread") else r["pin"], r["delivered"])
        for r in rows_out
    )
    for (cls, pin, delivered), n in sorted(counts.items()):
        print(f"{cls:4} {pin:60} delivered={delivered!s:5} {n}")


if __name__ == "__main__":
    main()
