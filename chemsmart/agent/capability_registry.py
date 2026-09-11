"""One schema for every capability of the agent, and a ladder computed
from what exists rather than hand-listed.

A capability is a program job type, a typed tool, a result selector, an
expression operation, a validation predicate, a literature constant, an
advisory skill, a guide, or a natural-language rule. Each climbs the same
ladder -- declared, wired, advertised, tested, qualified -- and each rung
is derived: declared from the registry that owns it, wired from the code
that joins it, advertised from the surface or prompt that shows it,
tested from the ``capability`` markers tests carry, qualified from the
release record and the host's own qualification store. Adding a
capability is therefore a local edit at its own registry; this module
notices it.
"""

from __future__ import annotations

import json
import re
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Mapping

CAPABILITY_KINDS = (
    "program_jobtype",
    "tool",
    "selector",
    "setting",
    "operation",
    "predicate",
    "signal",
    "gate",
    "constant",
    "skill",
    "guide",
    "rule",
)

LADDER = ("declared", "wired", "advertised", "tested", "qualified")

RELEASE_RECORD = Path(__file__).with_name("qualification") / "release.json"
HOST_QUALIFICATION_STORE = (
    Path.home() / ".chemsmart" / "agent" / "qualification.jsonl"
)

#: ``@pytest.mark.capability(...)`` and nothing else. The bare form
#: matched ``program_capability(...)`` and ``engine_capability(...)``
#: too, so the ladder's own input carried sixteen tokens no capability
#: could ever have -- ``cpu``, ``opt``, ``ts``, ``xtb``, ``2``, ``3`` --
#: scooped out of unrelated call sites. An instrument whose input is
#: polluted reports coverage it never measured.
_MARKER = re.compile(r"(?<!\w)capability\(([^)]*)\)")


@dataclass(frozen=True)
class CapabilityV1:
    kind: str
    id: str
    family: str
    tier: str
    declared_by: str
    wired_by: str = ""
    advertised_in: str = ""
    tested_by: tuple[str, ...] = ()
    #: Files whose coverage of this capability comes only from a family
    #: wildcard (``selector:*``). Deliberately does NOT lift `status` to
    #: `tested`: nine blanket markers were reporting tested for 557 of
    #: 712 capabilities, and `selector:*` alone certified all 393
    #: selectors from two files. A family test is real evidence about
    #: the family and no evidence about any one member.
    family_tested_by: tuple[str, ...] = ()
    qualified_by: tuple[str, ...] = ()

    @property
    def status(self) -> str:
        if self.qualified_by:
            return "qualified"
        if self.tested_by:
            return "tested"
        if self.advertised_in:
            return "advertised"
        if self.wired_by:
            return "wired"
        return "declared"

    @property
    def key(self) -> str:
        return f"{self.kind}:{self.id}"

    def public_record(self) -> dict:
        return {**asdict(self), "status": self.status}


def test_markers(tests_root: Path | None) -> dict[str, tuple[str, ...]]:
    """capability ids each test file pins, from ``capability(...)`` markers.
    A trailing ``*`` matches every id of that prefix."""

    found: dict[str, list[str]] = {}
    if tests_root is None or not Path(tests_root).is_dir():
        return {}
    for path in sorted(Path(tests_root).rglob("test_*.py")):
        text = path.read_text(encoding="utf-8")
        for match in _MARKER.finditer(text):
            for token in re.findall(r"['\"]([^'\"]+)['\"]", match.group(1)):
                found.setdefault(token, []).append(str(path.name))
    return {token: tuple(sorted(set(files))) for token, files in found.items()}


_SIGNAL_EMISSION = re.compile(r'"signal_id":\s*"([a-z_.0-9]+)"')


def _signal_emitters(package_root: Path) -> dict[str, tuple[str, ...]]:
    """Which module raises each anomaly signal, read from the source.

    Derived rather than declared, for the same reason ``tested_by`` is:
    a rung that is asserted cannot fail, and a sensor nobody raises must
    be able to report itself unwired.
    """

    found: dict[str, list[str]] = {}
    for path in sorted(Path(package_root).rglob("*.py")):
        try:
            text = path.read_text(encoding="utf-8")
        except OSError:
            continue
        for signal in set(_SIGNAL_EMISSION.findall(text)):
            found.setdefault(signal, []).append(path.name)
    return {
        signal: tuple(sorted(set(names))) for signal, names in found.items()
    }


def _gate_raisers(
    package_root: Path, gates: tuple[tuple[str, str], ...]
) -> dict[str, tuple[str, ...]]:
    """Which module raises each code gate, read from the source."""

    ids = {gate_id for gate_id, _ in gates}
    found: dict[str, list[str]] = {}
    for path in sorted(Path(package_root).rglob("*.py")):
        if path.name == "rules.py":
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except OSError:
            continue
        for gate_id in ids:
            if gate_id in text:
                found.setdefault(gate_id, []).append(path.name)
    return {
        gate_id: tuple(sorted(set(names))) for gate_id, names in found.items()
    }


def _matches(pattern: str, key: str) -> bool:
    if pattern.endswith("*"):
        return key.startswith(pattern[:-1])
    return pattern == key


def load_release_records(path: Path = RELEASE_RECORD) -> dict[str, dict]:
    try:
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    return {
        f"{item['kind']}:{item['id']}": dict(item)
        for item in payload.get("records", ())
        if item.get("kind") and item.get("id")
    }


def load_host_qualifications(
    path: Path = HOST_QUALIFICATION_STORE,
) -> dict[str, list[dict]]:
    found: dict[str, list[dict]] = {}
    try:
        lines = Path(path).read_text(encoding="utf-8").splitlines()
    except OSError:
        return {}
    for line in lines:
        try:
            item = json.loads(line)
        except json.JSONDecodeError:
            continue
        if item.get("kind") and item.get("id"):
            found.setdefault(f"{item['kind']}:{item['id']}", []).append(item)
    return found


def record_host_qualification(
    entries: Iterable[Mapping[str, object]],
    path: Path = HOST_QUALIFICATION_STORE,
) -> int:
    """Append qualification entries to the host store; returns how many."""

    rows = [json.dumps(dict(item), sort_keys=True) for item in entries]
    if not rows:
        return 0
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        for row in rows:
            handle.write(row + "\n")
    return len(rows)


def build_capability_registry(
    *,
    tests_root: Path | None = None,
    release_path: Path = RELEASE_RECORD,
    host_store: Path | None = HOST_QUALIFICATION_STORE,
) -> tuple[CapabilityV1, ...]:
    """Every capability, from the registries that already own each kind."""

    from chemsmart.agent.capabilities import load_program_capabilities
    from chemsmart.agent.execution import ANOMALY_SIGNALS
    from chemsmart.agent.guides import GUIDES, LEAF_OPERATIONS, LEAF_TOOLS
    from chemsmart.agent.rules import CODE_GATES, POLICY_RULES
    from chemsmart.agent.scientific_toolchain import (
        ANALYSIS_VALIDATION_PREDICATES,
    )
    from chemsmart.agent.skills import available_skill_ids
    from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1
    from chemsmart.agent.tool_specs import (
        build_approved_execution_tool_surface,
        build_command_compiled_tool_surface,
    )
    from chemsmart.analysis.literature_constants import LITERATURE_CONSTANTS
    from chemsmart.analysis.quantity_expressions import OPERATION_DESCRIPTIONS
    from chemsmart.analysis.result_readers import RESULT_READERS

    markers = test_markers(tests_root)
    release = load_release_records(release_path)
    host = load_host_qualifications(host_store) if host_store else {}

    def tested(key: str) -> tuple[str, ...]:
        """Files naming this capability exactly. Wildcards excluded."""

        files: list[str] = []
        for pattern, names in markers.items():
            if not pattern.endswith("*") and _matches(pattern, key):
                files.extend(names)
        return tuple(sorted(set(files)))

    def family_tested(key: str) -> tuple[str, ...]:
        """Files covering this capability only through a wildcard."""

        files: list[str] = []
        for pattern, names in markers.items():
            if pattern.endswith("*") and _matches(pattern, key):
                files.extend(names)
        return tuple(sorted(set(files)))

    def qualified(key: str) -> tuple[str, ...]:
        refs: list[str] = []
        record = release.get(key)
        if record:
            refs.append(
                f"release:{record.get('status', 'claimed')}:"
                f"{record.get('run') or record.get('source') or ''}"
            )
        for item in host.get(key, ()):
            refs.append(f"host:{item.get('run', '')}")
        return tuple(refs)

    records: list[CapabilityV1] = []

    registry = load_program_capabilities()
    for program in registry.programs:
        for engine, jobtype in sorted(set(program.preview_engine_job_pairs)):
            key_id = f"{program.program}:{engine}:{jobtype}"
            executable = (engine, jobtype) in set(
                program.execution_engine_job_pairs
            )
            key = f"program_jobtype:{key_id}"
            records.append(
                CapabilityV1(
                    kind="program_jobtype",
                    id=key_id,
                    family=program.program,
                    tier="T0",
                    declared_by="chemsmart.settings.capabilities",
                    wired_by="live Click schema + preview overlay",
                    advertised_in=(
                        "inspect_program" if executable else "preview only"
                    ),
                    tested_by=tested(key),
                    family_tested_by=family_tested(key),
                    qualified_by=qualified(key) if executable else (),
                )
            )

    handled = set(CommandCompiledToolHostV1.TOOL_HANDLERS)
    every_leaf = tuple(guide.guide_id for guide in GUIDES)
    planning = {
        item["function"]["name"]
        for item in build_command_compiled_tool_surface(
            guides=every_leaf
        ).tool_definitions
    }
    execution = {
        item["function"]["name"]
        for item in build_approved_execution_tool_surface().tool_definitions
    }
    for name in sorted(handled | planning | execution):
        key = f"tool:{name}"
        advertised = []
        if name in planning:
            advertised.append(
                f"planning (leaf {LEAF_TOOLS[name]})"
                if name in LEAF_TOOLS
                else "planning (stem)"
            )
        if name in execution:
            advertised.append("execution")
        records.append(
            CapabilityV1(
                kind="tool",
                id=name,
                family=LEAF_TOOLS.get(name, "stem"),
                tier="T3" if name in LEAF_TOOLS else "T0",
                declared_by="chemsmart.agent.tool_specs",
                wired_by=(
                    "CommandCompiledToolHostV1.TOOL_HANDLERS"
                    if name in handled
                    else ""
                ),
                advertised_in=", ".join(advertised),
                tested_by=tested(key),
                family_tested_by=family_tested(key),
                qualified_by=qualified(key),
            )
        )

    for program, reader in sorted(RESULT_READERS.items()):
        accessors = set(getattr(reader, "accessors", ()) or ())
        for jobtype, selectors in reader.jobtype_selectors:
            for selector in selectors:
                key_id = f"{program}:{jobtype}:{selector}"
                records.append(
                    CapabilityV1(
                        kind="selector",
                        id=key_id,
                        family=program,
                        tier="T0",
                        declared_by="chemsmart.analysis.result_readers",
                        # Computed, not asserted. This was the constant
                        # string "reader accessor" for every selector,
                        # so the ladder reported `wired` for a selector
                        # whose accessor had been deleted -- a rung that
                        # cannot fail measures nothing.
                        wired_by=(
                            f"{program} accessor {selector}"
                            if selector in accessors
                            else ""
                        ),
                        advertised_in="extract_result_quantities",
                        tested_by=tested(f"selector:{key_id}"),
                        family_tested_by=family_tested(f"selector:{key_id}"),
                        qualified_by=qualified(f"selector:{key_id}"),
                    )
                )

    # Project-owned settings are the load-bearing surface of the hub
    # thesis -- the model's whole method vocabulary travels through
    # project YAML -- and they were not a capability kind at all, so two
    # tests already carried `setting:` markers that matched nothing and
    # `opt_convergence` could be advertised, accepted, normalised and
    # read back while the reader had no property for it.
    for program in registry.programs:
        parameters = tuple(
            sorted(getattr(program, "project_owned_parameters", ()) or ())
        )
        domains = dict(getattr(program, "project_parameter_domains", ()) or ())
        for parameter in parameters:
            key_id = f"{program.program}:{parameter}"
            records.append(
                CapabilityV1(
                    kind="setting",
                    id=key_id,
                    family=program.program,
                    tier="T0",
                    declared_by="chemsmart.settings.capabilities",
                    # A declared domain is what makes the round trip
                    # generable without a human writing chemistry into
                    # a test, so it is the rung, not a decoration.
                    wired_by=(
                        f"domain {len(domains[parameter])} values"
                        if parameter in domains
                        else ""
                    ),
                    advertised_in="inspect_program",
                    tested_by=tested(f"setting:{key_id}")
                    + tested(f"setting:{parameter}"),
                    family_tested_by=family_tested(f"setting:{key_id}"),
                    qualified_by=qualified(f"setting:{key_id}"),
                )
            )

    # Anomaly sensors. `wired_by` is computed by looking for the id at
    # an emitting site, the same way `tested_by` is computed by looking
    # for a marker in a test: a signal nothing raises is declared and
    # not wired, and the ladder must be able to say so.
    emitters = _signal_emitters(Path(__file__).parent)
    for signal in ANOMALY_SIGNALS:
        records.append(
            CapabilityV1(
                kind="signal",
                id=signal,
                family=signal.split(".", 1)[0],
                tier="T0",
                declared_by="chemsmart.agent.execution.ANOMALY_SIGNALS",
                wired_by=", ".join(emitters.get(signal, ())),
                advertised_in="anomalies_observed / settlement word",
                tested_by=tested(f"signal:{signal}"),
                family_tested_by=family_tested(f"signal:{signal}"),
                qualified_by=qualified(f"signal:{signal}"),
            )
        )

    # Code gates. `wired_by` is the source that raises the gate id, so
    # a gate declared in the charter's list and raised by nothing
    # reports itself declared-and-unwired instead of passing silently.
    raisers = _gate_raisers(Path(__file__).parent, CODE_GATES)
    for gate_id, invariant in CODE_GATES:
        records.append(
            CapabilityV1(
                kind="gate",
                id=gate_id,
                family=gate_id.split(".", 1)[0],
                tier="T0",
                declared_by="chemsmart.agent.rules.CODE_GATES",
                wired_by=", ".join(raisers.get(gate_id, ())),
                advertised_in=invariant,
                tested_by=tested(f"gate:{gate_id}"),
                family_tested_by=family_tested(f"gate:{gate_id}"),
                qualified_by=qualified(f"gate:{gate_id}"),
            )
        )

    for name in sorted(OPERATION_DESCRIPTIONS):
        leaf = LEAF_OPERATIONS.get(name)
        records.append(
            CapabilityV1(
                kind="operation",
                id=name,
                family=leaf or "stem",
                tier="T3" if leaf else "T1",
                declared_by="chemsmart.analysis.quantity_expressions",
                wired_by="evaluate_quantity_expression",
                advertised_in=(
                    f"evaluate_quantity_expression (leaf {leaf})"
                    if leaf
                    else "evaluate_quantity_expression (stem)"
                ),
                tested_by=tested(f"operation:{name}"),
                family_tested_by=family_tested(f"operation:{name}"),
                qualified_by=qualified(f"operation:{name}"),
            )
        )
    for name in ANALYSIS_VALIDATION_PREDICATES:
        records.append(
            CapabilityV1(
                kind="predicate",
                id=name,
                family="validation",
                tier="T1",
                declared_by="chemsmart.agent.scientific_toolchain",
                wired_by="evaluate_scientific_validation",
                advertised_in="plan_scientific_workflow.validation_rules",
                tested_by=tested(f"predicate:{name}"),
                family_tested_by=family_tested(f"predicate:{name}"),
                qualified_by=qualified(f"predicate:{name}"),
            )
        )
    for name in sorted(LITERATURE_CONSTANTS):
        records.append(
            CapabilityV1(
                kind="constant",
                id=name,
                family="constants",
                tier="T3",
                declared_by="chemsmart.analysis.literature_constants",
                wired_by="constant operation",
                advertised_in="evaluate_quantity_expression (leaf constants)",
                tested_by=tested(f"constant:{name}"),
                family_tested_by=family_tested(f"constant:{name}"),
                qualified_by=qualified(f"constant:{name}"),
            )
        )
    for name in available_skill_ids():
        records.append(
            CapabilityV1(
                kind="skill",
                id=name,
                family="skills",
                tier="T1",
                declared_by="chemsmart/agent/skills",
                wired_by="open_guide",
                advertised_in="system prompt skill index",
                tested_by=tested(f"skill:{name}"),
                family_tested_by=family_tested(f"skill:{name}"),
                qualified_by=qualified(f"skill:{name}"),
            )
        )
    for guide in GUIDES:
        records.append(
            CapabilityV1(
                kind="guide",
                id=guide.guide_id,
                family=guide.guide_id,
                tier=guide.tier,
                declared_by="chemsmart.agent.guides",
                wired_by="activate_guides",
                advertised_in="system prompt guide index; open_guide",
                tested_by=tested(f"guide:{guide.guide_id}"),
                family_tested_by=family_tested(f"guide:{guide.guide_id}"),
                qualified_by=qualified(f"guide:{guide.guide_id}"),
            )
        )
    for rule in POLICY_RULES:
        records.append(
            CapabilityV1(
                kind="rule",
                id=rule.rule_id,
                family=rule.placement,
                tier=rule.tier,
                declared_by="chemsmart.agent.rules",
                wired_by="render_rules",
                advertised_in=rule.placement,
                tested_by=tested(f"rule:{rule.rule_id}"),
                family_tested_by=family_tested(f"rule:{rule.rule_id}"),
                qualified_by=qualified(f"rule:{rule.rule_id}"),
            )
        )
    return tuple(records)


def render_capability_matrix(records: Iterable[CapabilityV1]) -> str:
    """A human table: one row per capability, the ladder said out loud."""

    rows = list(records)
    width = max((len(item.key) for item in rows), default=10)
    lines = [f"{'capability':<{width}}  status      tested  qualified"]
    for item in rows:
        lines.append(
            f"{item.key:<{width}}  {item.status:<10}  "
            f"{len(item.tested_by):>6}  "
            f"{'; '.join(item.qualified_by) or 'unsupported'}"
        )
    counts: dict[str, int] = {}
    for item in rows:
        counts[item.status] = counts.get(item.status, 0) + 1
    lines.append("")
    lines.append(
        "by status: "
        + ", ".join(f"{stage} {counts.get(stage, 0)}" for stage in LADDER)
    )
    return "\n".join(lines)


__all__ = [
    "CAPABILITY_KINDS",
    "HOST_QUALIFICATION_STORE",
    "LADDER",
    "RELEASE_RECORD",
    "CapabilityV1",
    "build_capability_registry",
    "load_host_qualifications",
    "load_release_records",
    "record_host_qualification",
    "render_capability_matrix",
    "test_markers",
]
