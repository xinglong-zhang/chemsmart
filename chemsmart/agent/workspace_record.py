"""What this workspace has learned: a host-written record that spans goals.

Nothing in the goal machinery spanned goals. Two sealed windows ran the
same six blind queries; the sulfone's gauche preference in acetonitrile
came back at -2.23 kJ/mol under def2-SVP and -0.91 under def2-TZVP, on
either side of the author's stated 2 kJ/mol design cutoff, and the only
reader who could see both numbers was the scientist reading the two
windows by hand (NOVEL-1/2 po2, 2026-09-04). Every goal started with no
memory of any other run.

This record is append-only, receipt-derived and host-written: per
recorded run, one row per verified result (which input, which bound
state, which program and job type, at which level of theory, what
energy, what ending) and one row per delivered claim (which id, what
value, in what unit, from which run). Rendered into every wake, it lets
a session see that the same input was computed before, at what level,
and what was claimed -- and a divergence line names two delivered values
of one claim that disagree. The host names the disagreement and never
its meaning; a model never writes here; nothing here grades anything.
The owner ruled it built (2026-09-05) and ruled a host-seeded prior run
admissible in a sealed window when it is written from receipts and named
by digest.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping

from chemsmart.agent.terminal_states import GEOMETRY_SEARCH_JOBTYPES

WORKSPACE_RECORD_FILE = "workspace-record.jsonl"
WORKSPACE_RECORD_SCHEMA = "chemsmart.workspace-record.v1"

#: The project settings that name a level of theory in one line.
_LEVEL_KEYS = (
    "functional",
    "ab_initio",
    "basis",
    "aux_basis",
    "dispersion",
    "solvent_model",
    "solvent",
    "freq",
    # A correlated method's level is not named without its frozen-core
    # convention (PySCF correlates every electron unless told otherwise;
    # ORCA and Gaussian freeze the core by default), and an excited-surface
    # optimisation's level is the response it ran on and the root it
    # followed: two results with one functional and one basis are at
    # different levels when one of these differs.
    "frozen_core",
    "response_method",
    "state_manifold",
    "nstates",
    "excited_state_root",
)

#: Two delivered values of one claim differ when they disagree by more
#: than this fraction of the larger magnitude. The number is displayed;
#: whether the disagreement matters is the scientist's.
DIVERGENCE_RELATIVE_TOLERANCE = 0.05


def workspace_record_path(workspace: str | Path) -> Path:
    return Path(workspace) / ".chemsmart-agent" / WORKSPACE_RECORD_FILE


#: The output a program's reader opens, by suffix (the executor files every
#: native output under one kind).
_RESULT_SUFFIXES: dict[str, tuple[str, ...]] = {
    "orca": (".out",),
    "gaussian": (".log", ".out"),
    "xtb": (".out",),
    "pyscf": (".h5",),
}


def _result_artifact_sha256(program: str, artifacts: Any) -> str:
    """The digest of the one output a reader opens for this result."""

    suffixes = _RESULT_SUFFIXES.get(str(program).lower(), ())
    for item in artifacts:
        if not isinstance(item, Mapping):
            continue
        kind = str(item.get("kind") or "")
        path = str(item.get("path") or "")
        if kind.endswith("_output") and kind != "program_output":
            return str(item.get("sha256") or "")
        if kind == "program_output" and any(
            path.endswith(s) for s in suffixes
        ):
            return str(item.get("sha256") or "")
    return ""


def _level_from_settings_text(text: str) -> dict[str, Any]:
    try:
        settings = json.loads(text) if text else {}
    except json.JSONDecodeError:
        return {}
    if not isinstance(settings, Mapping):
        return {}
    return {
        key: settings[key]
        for key in _LEVEL_KEYS
        if settings.get(key) not in (None, "", False)
    }


def _review_rows(review_file: str | Path | None) -> dict[str, dict]:
    if review_file is None:
        return {}
    try:
        review = json.loads(Path(review_file).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    packet = review.get("workflow_execution_review") or review
    rows: dict[str, dict] = {}
    for row in packet.get("node_reviews") or ():
        if isinstance(row, Mapping) and row.get("node_id"):
            rows[str(row["node_id"])] = dict(row)
    return rows


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _assessment_with_provenance(
    assessment: Mapping[str, Any] | None,
    provenance: Mapping[str, Any] | None,
) -> dict[str, Any] | None:
    """The assessment, carrying what the host resolved to reach it."""

    if assessment is None:
        return None
    row = dict(assessment)
    if provenance:
        reference = str(provenance.get("uncertainty_reference") or "")
        if reference:
            row["uncertainty_reference"] = reference
        components = tuple(provenance.get("uncertainty_components") or ())
        if components:
            row["uncertainty_components"] = components
        observations = tuple(provenance.get("uncertainty_observations") or ())
        combination = provenance.get("uncertainty_combination")
        if combination:
            row["uncertainty_combination"] = dict(combination)
        if observations:
            row["uncertainty_observations"] = observations
    return row


def record_run(
    workspace: str | Path,
    *,
    goal_id: str,
    cycle: int,
    run_events_path: str | Path,
    run: str = "",
    review_file: str | Path | None = None,
    recorded_at: str | None = None,
) -> int:
    """Append what one run's stream proves to the workspace record.

    Returns the number of rows appended. Reads only receipts the stream
    carries and the review the human decided on; a stream it cannot read
    appends nothing.
    """

    try:
        lines = Path(run_events_path).read_text(encoding="utf-8").splitlines()
    except OSError:
        return 0
    rows = _review_rows(review_file)
    stamp = recorded_at or _utc_now()
    entries: list[dict[str, Any]] = []
    levels_in_run: set[str] = set()
    # Which node's reached geometry each consumer was handed, so a later
    # reader can join a Hessian's characterisation onto the optimisation.
    handoffs: dict[str, str] = {}
    for line in lines:
        try:
            event = json.loads(line)
        except json.JSONDecodeError:
            continue
        if event.get("kind") != "optimized_geometry_handed_off":
            continue
        payload = event.get("payload") or {}
        if str(payload.get("status") or "") == "validated_handoff":
            handoffs[str(payload.get("consumer_node_id") or "")] = str(
                payload.get("producer_node_id") or ""
            )
    for line in lines:
        try:
            event = json.loads(line)
        except json.JSONDecodeError:
            continue
        kind = event.get("kind")
        payload = event.get("payload") or {}
        if kind == "program_result_verified":
            record = payload.get("record") or {}
            node_id = str(record.get("node_id") or "")
            observations = record.get("observations") or {}
            program = str(record.get("program") or "")
            block = observations.get(program) or {}
            if not isinstance(block, Mapping):
                block = {}
            row = rows.get(node_id, {})
            identity = row.get("molecular_identity") or {}
            coordinate = identity.get("coordinate_identity") or {}
            level_sha256 = str(row.get("project_settings_text_sha256") or "")
            if level_sha256:
                levels_in_run.add(level_sha256)
            entries.append(
                {
                    "kind": "result",
                    "goal_id": goal_id,
                    "cycle": int(cycle),
                    "run": run,
                    "node_id": node_id,
                    "program": program,
                    "jobtype": str(
                        observations.get("jobtype") or record.get("jobtype")
                    ),
                    "input_artifact_sha256": str(
                        record.get("input_artifact_sha256") or ""
                    ),
                    "geometry_artifact_sha256": str(
                        coordinate.get("geometry_artifact_sha256") or ""
                    ),
                    "formula": str(identity.get("formula") or ""),
                    "charge": block.get("charge", identity.get("charge")),
                    "multiplicity": block.get(
                        "multiplicity", identity.get("multiplicity")
                    ),
                    "level": _level_from_settings_text(
                        str(row.get("project_settings_text") or "")
                    ),
                    "level_sha256": level_sha256,
                    "energy_hartree": block.get("energy_hartree"),
                    "vibrational_mode_count": block.get(
                        "vibrational_mode_count"
                    ),
                    "printed_modes": printed_modes(record),
                    "geometry_producer_node_id": handoffs.get(node_id, ""),
                    "state": str(record.get("state") or ""),
                    "result_receipt_sha256": str(
                        record.get("receipt_sha256") or ""
                    ),
                    "output_artifact_sha256s": tuple(
                        str(item.get("sha256") or "")
                        for item in record.get("output_artifacts") or ()
                        if item.get("sha256")
                    ),
                    "result_artifact_sha256": _result_artifact_sha256(
                        program, record.get("output_artifacts") or ()
                    ),
                    "recorded_at": stamp,
                }
            )
        elif kind == "analysis_claims_recorded":
            record = payload.get("record") or {}
            # The assessment travels with the number it judges. Without
            # it a later cycle reading an earlier claim through this
            # record found no assessment at all, and a missing
            # assessment read as no open requirement -- the mirror of
            # the sticky union, and the more dangerous half.
            assessments = {
                str(row.get("observable_id") or ""): row
                for row in payload.get("sufficiency") or ()
                if isinstance(row, Mapping) and row.get("observable_id")
            }
            # And what the word rests on. The host resolves a citation
            # to grant `met` and the record kept a bare boolean, so
            # neither a human nor the next cycle could ask which
            # receipt backed the number that discharged the contract --
            # the same value-resolved-then-discarded defect one layer
            # out. No gate can decide whether a number is an
            # uncertainty; that is chemistry. What the host owes is
            # that its word is auditable, and that means carrying the
            # citation and the terms beside it (owner ruling,
            # 2026-09-09).
            provenance = {
                str(claim.get("claim_id") or ""): {
                    "uncertainty_reference": str(
                        claim.get("uncertainty_reference") or ""
                    ),
                    "uncertainty_components": tuple(
                        dict(item)
                        for item in claim.get("uncertainty_components") or ()
                        if isinstance(item, Mapping)
                    ),
                    "uncertainty_observations": tuple(
                        str(item)
                        for item in claim.get("uncertainty_observations") or ()
                    ),
                    "uncertainty_combination": (
                        dict(claim.get("uncertainty_combination") or {})
                        or None
                    ),
                }
                for claim in record.get("claims") or ()
                if claim.get("uncertainty_reference")
                or claim.get("uncertainty_components")
                or claim.get("uncertainty_observations")
            }
            for claim in record.get("claims") or ():
                value = claim.get("display_value")
                if isinstance(value, bool) or not isinstance(
                    value, (int, float)
                ):
                    continue
                entries.append(
                    {
                        "kind": "claim",
                        "goal_id": goal_id,
                        "cycle": int(cycle),
                        "run": run,
                        "claim_id": str(claim.get("claim_id") or ""),
                        "quantity_id": str(claim.get("quantity_id") or ""),
                        "value": float(value),
                        "unit": str(claim.get("display_unit") or ""),
                        # The dimension travels with the number: the
                        # settlement judges a declared observable in it,
                        # and a row that carried only a display unit let
                        # six 'e' declarations be answered by
                        # dimensionless counts (OPEN-1 ino3, 2026-09-07).
                        "dimension": tuple(
                            int(item) for item in claim.get("dimension") or ()
                        ),
                        "sufficiency": _assessment_with_provenance(
                            assessments.get(str(claim.get("claim_id") or ""))
                            or assessments.get(
                                str(claim.get("quantity_id") or "")
                            ),
                            provenance.get(str(claim.get("claim_id") or "")),
                        ),
                        "claim_receipt_sha256": str(
                            payload.get("receipt_sha256") or ""
                        ),
                        "source_receipt_sha256": str(
                            claim.get("source_receipt_sha256") or ""
                        ),
                        "recorded_at": stamp,
                    }
                )
    for entry in entries:
        if entry["kind"] == "claim":
            entry["level_sha256s"] = tuple(sorted(levels_in_run))
    if not entries:
        return 0
    path = workspace_record_path(workspace)
    # Durable, because the invariant is durable. The driver held a
    # per-(cycle, stream) marker in a set on the instance, so a second
    # process -- a restart, a wake, any new driver -- started empty and
    # appended the same stream's rows again: 17 rows became 34 with ten
    # claim ids doubled. Every row already carries what identifies it.
    already = {
        (
            str(row.get("goal_id") or ""),
            int(row.get("cycle") or 0),
            str(row.get("kind") or ""),
            str(row.get("claim_id") or row.get("node_id") or ""),
            str(
                row.get("claim_receipt_sha256")
                or row.get("result_artifact_sha256")
                or ""
            ),
        )
        for row in read_workspace_record(workspace)
    }
    entries = [
        entry
        for entry in entries
        if (
            str(entry.get("goal_id") or ""),
            int(entry.get("cycle") or 0),
            str(entry.get("kind") or ""),
            str(entry.get("claim_id") or entry.get("node_id") or ""),
            str(
                entry.get("claim_receipt_sha256")
                or entry.get("result_artifact_sha256")
                or ""
            ),
        )
        not in already
    ]
    if not entries:
        return 0
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        for entry in entries:
            handle.write(
                json.dumps(
                    {"schema_version": WORKSPACE_RECORD_SCHEMA, **entry},
                    sort_keys=True,
                )
                + "\n"
            )
    return len(entries)


def printed_modes(record: Mapping[str, Any]) -> bool:
    """Whether a verified result carries a printed frequency block.

    Read from the validator's own observations, program-neutrally: the
    block that records ``vibrational_mode_count`` answers; a block the
    neutral sensor step wrote carries ``consequential_imaginary_mode_count``
    only when the reader served frequencies, so its presence answers for
    a record written before the count travelled.
    """

    observations = record.get("observations") or {}
    for value in observations.values():
        if isinstance(value, Mapping) and "vibrational_mode_count" in value:
            try:
                return int(value.get("vibrational_mode_count") or 0) > 0
            except (TypeError, ValueError):
                return False
    return any(
        isinstance(value, Mapping)
        and "consequential_imaginary_mode_count" in value
        for value in observations.values()
    )


def _row_printed_modes(entry: Mapping[str, Any]) -> bool | None:
    """A recorded row's word on its modes; None when the row cannot say."""

    if isinstance(entry.get("printed_modes"), bool):
        return bool(entry["printed_modes"])
    try:
        return int(entry.get("vibrational_mode_count") or 0) != 0
    except (TypeError, ValueError):
        return None


def failed_artifacts(workspace: str | Path) -> tuple[str, ...]:
    """Output artifacts of recorded results the validator did not pass.

    A claim recorded in a later cycle's session on a result an earlier
    run typed failed carries no verified record in its own stream; the
    workspace record does, so the settlement can still say the number
    stands on a node that did not meet its promise, or on one the
    session had the host characterise (PySCF round g5, 2026-09-12: six
    claims on a characterised saddle, and the settlement said nothing).
    """

    digests: set[str] = set()
    for entry in read_workspace_record(workspace):
        if entry.get("kind") != "result":
            continue
        if str(entry.get("state") or "") == "valid":
            continue
        digests.update(
            str(item) for item in entry.get("output_artifact_sha256s") or ()
        )
    return tuple(sorted(digests))


def uncharacterised_artifacts(workspace: str | Path) -> tuple[str, ...]:
    """Output artifacts of recorded opt/ts results that printed no modes.

    A result verified in an earlier cycle is read as a registered result
    later, and that later stream carries no verified record for it; the
    workspace record does, so a number claimed on it still says its
    stationary point was never characterised.

    An optimisation whose reached geometry a validated Hessian consumed
    through the handoff edge, in the same run, is characterised by that
    Hessian: ORCA's opt+freq is one node and PySCF's or xTB's are two,
    and the same physics gets the same word (PySCF round, 2026-09-12).
    """

    entries = read_workspace_record(workspace)
    characterised: set[tuple[str, str, str]] = set()
    for entry in entries:
        if entry.get("kind") != "result":
            continue
        if str(entry.get("state") or "") != "valid":
            continue
        producer = str(entry.get("geometry_producer_node_id") or "")
        if producer and _row_printed_modes(entry):
            characterised.add(
                (
                    str(entry.get("goal_id") or ""),
                    str(entry.get("run") or ""),
                    producer,
                )
            )
    digests: set[str] = set()
    for entry in entries:
        if entry.get("kind") != "result":
            continue
        if str(entry.get("jobtype") or "") not in GEOMETRY_SEARCH_JOBTYPES:
            continue
        printed = _row_printed_modes(entry)
        if printed is None or printed:
            continue
        key = (
            str(entry.get("goal_id") or ""),
            str(entry.get("run") or ""),
            str(entry.get("node_id") or ""),
        )
        if key in characterised:
            continue
        digests.update(
            str(item) for item in entry.get("output_artifact_sha256s") or ()
        )
    return tuple(sorted(digests))


def read_workspace_record(workspace: str | Path) -> tuple[dict, ...]:
    path = workspace_record_path(workspace)
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError:
        return ()
    entries = []
    for line in lines:
        try:
            entries.append(json.loads(line))
        except json.JSONDecodeError:
            continue
    return tuple(entries)


def _normalised(identifier: str) -> str:
    """An id with case and separators removed.

    Two windows named one delivered number dg-sulfone-mecn and
    dg_sulfone_mecn, and a join on the raw id saw two claims and no
    disagreement. Case and separators are not identity; the raw ids
    are displayed beside the join.
    """

    return "".join(ch for ch in str(identifier).lower() if ch.isalnum())


def _claim_keys(entry: Mapping[str, Any]) -> tuple[str, ...]:
    keys: list[str] = []
    for field in ("claim_id", "quantity_id"):
        value = _normalised(str(entry.get(field) or ""))
        if value and value not in keys:
            keys.append(value)
    return tuple(keys)


def _where(entry: Mapping[str, Any], value: float) -> dict[str, Any]:
    return {
        "value": value,
        "claim_id": entry.get("claim_id"),
        "quantity_id": entry.get("quantity_id"),
        "goal_id": entry.get("goal_id"),
        "run": entry.get("run"),
        "level_sha256s": entry.get("level_sha256s"),
    }


def divergences(entries: tuple[dict, ...]) -> tuple[dict, ...]:
    """Pairs of delivered values of one claim that disagree.

    Same claim (by either id, case and separators ignored), same unit,
    from different runs or levels, differing by more than the stated
    fraction of the larger magnitude. Both numbers are shown with where
    they came from; the host says that they differ and never why.
    """

    claims = [entry for entry in entries if entry.get("kind") == "claim"]
    by_key: dict[tuple[str, str], list[dict]] = {}
    for entry in claims:
        for key in _claim_keys(entry):
            by_key.setdefault((key, str(entry.get("unit"))), []).append(entry)
    found: list[dict] = []
    seen: set[tuple[str, int, int]] = set()
    for (claim_id, unit), group in sorted(by_key.items()):
        for index, first in enumerate(group):
            for second in group[index + 1 :]:
                pair = (unit, id(first), id(second))
                if pair in seen:
                    continue
                seen.add(pair)
                same_run = first.get("goal_id") == second.get(
                    "goal_id"
                ) and first.get("run") == second.get("run")
                same_level = tuple(first.get("level_sha256s") or ()) == tuple(
                    second.get("level_sha256s") or ()
                )
                if same_run and same_level:
                    continue
                a, b = float(first["value"]), float(second["value"])
                scale = max(abs(a), abs(b))
                if scale == 0.0:
                    continue
                if abs(a - b) <= DIVERGENCE_RELATIVE_TOLERANCE * scale:
                    continue
                found.append(
                    {
                        "claim_id": claim_id,
                        "unit": unit,
                        "values": (_where(first, a), _where(second, b)),
                        "difference": a - b,
                        "relative_tolerance": DIVERGENCE_RELATIVE_TOLERANCE,
                    }
                )
    return tuple(found)


def render_workspace_record(workspace: str | Path) -> dict[str, Any]:
    """The record as a wake context reads it; empty when nothing is recorded."""

    entries = read_workspace_record(workspace)
    if not entries:
        return {}
    results = tuple(e for e in entries if e.get("kind") == "result")
    claims = tuple(e for e in entries if e.get("kind") == "claim")
    levels: dict[str, dict] = {}
    for entry in results:
        digest = str(entry.get("level_sha256") or "")
        if digest and digest not in levels and entry.get("level"):
            levels[digest] = dict(entry["level"])
    return {
        "schema_version": WORKSPACE_RECORD_SCHEMA,
        "meaning": (
            "host-written from receipts of this workspace's recorded runs, "
            "across goals; a level is named by the digest of the project "
            "settings that produced it; a divergence is two delivered "
            "values of one claim that disagree beyond the stated tolerance, "
            "stated and never explained"
        ),
        "levels": levels,
        "results": tuple(
            {
                **{
                    key: entry.get(key)
                    for key in (
                        "goal_id",
                        "cycle",
                        "run",
                        "node_id",
                        "program",
                        "jobtype",
                        "input_artifact_sha256",
                        "formula",
                        "charge",
                        "multiplicity",
                        "level_sha256",
                        "energy_hartree",
                        "vibrational_mode_count",
                        "state",
                    )
                },
                # The id a reading tool accepts for this result, when the
                # bytes are registered: the host's own rule, stated where
                # the number is shown (REACH-1 ino3 passed digests and was
                # refused ten times).
                "artifact_id": (
                    f"{entry.get('program')}-result-"
                    f"{str(entry.get('result_artifact_sha256'))[:16]}"
                    if entry.get("result_artifact_sha256")
                    else ""
                ),
            }
            for entry in results
        ),
        "claims": tuple(
            {
                key: entry.get(key)
                for key in (
                    "goal_id",
                    "cycle",
                    "run",
                    "claim_id",
                    "quantity_id",
                    "value",
                    "unit",
                    # The settlement and the completion gate both join
                    # a claim to its declaration by id *and* dimension,
                    # and the projection the wake hands the model
                    # carried the unit and not the dimension: a session
                    # told to read the record could not see why one of
                    # its own rows failed to join.
                    "dimension",
                    "level_sha256s",
                    # What this goal established about the number's
                    # precision, so a session reading the record sees
                    # the current assessment rather than only the value.
                    "sufficiency",
                )
            }
            for entry in claims
        ),
        "divergences": divergences(entries),
    }


__all__ = [
    "DIVERGENCE_RELATIVE_TOLERANCE",
    "WORKSPACE_RECORD_FILE",
    "WORKSPACE_RECORD_SCHEMA",
    "divergences",
    "read_workspace_record",
    "record_run",
    "render_workspace_record",
    "uncharacterised_artifacts",
    "workspace_record_path",
]
