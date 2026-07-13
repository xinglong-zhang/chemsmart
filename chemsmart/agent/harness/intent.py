"""Canonical intent assertions for ChemSmart commands.

The command semantic gate answers whether ChemSmart can execute a command.
This module answers the independent question: did the command preserve the
observable research intent?  Runtime requests use conservative extraction;
benchmark fixtures may provide a complete :class:`IntentSpec` directly.
"""

from __future__ import annotations

import re
import shlex
from dataclasses import asdict, dataclass, field
from pathlib import PurePath
from typing import Any, Literal

from chemsmart.agent.harness.workflow_state import project_name_from_request
from chemsmart.agent.model_command_parser import parse_model_command

IntentVerdict = Literal["ok", "reject"]


@dataclass(frozen=True)
class IntentSpec:
    action: str | None = None
    program: str | None = None
    kind: str | None = None
    project: str | None = None
    server: str | None = None
    input_path: str | None = None
    output_path: str | None = None
    charge: int | str | None = None
    multiplicity: int | str | None = None
    execution_mode: str | None = None
    chemistry: dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_dict(cls, value: dict[str, Any]) -> "IntentSpec":
        known = {
            key: value.get(key)
            for key in (
                "action",
                "program",
                "kind",
                "project",
                "server",
                "input_path",
                "output_path",
                "charge",
                "multiplicity",
                "execution_mode",
            )
        }
        chemistry = value.get("chemistry")
        known["chemistry"] = dict(chemistry) if isinstance(chemistry, dict) else {}
        return cls(**known)

    @classmethod
    def from_request(cls, request: str) -> "IntentSpec":
        text = str(request or "")
        lowered = text.lower()
        program = "gaussian" if "gaussian" in lowered else "orca" if "orca" in lowered else None
        kind = _kind_from_request(lowered, program)
        action = "sub" if re.search(r"\b(submit|submission|queue|cluster|hpc|chemsmart\s+sub)\b", lowered) else (
            "run" if re.search(r"\b(run|locally|chemsmart\s+run)\b", lowered) else None
        )
        input_match = re.search(
            r"(?<![\w./-])([\w./-]+\.(?:xyz|log|out|com|gjf|inp|db|hess))\b",
            text,
            flags=re.IGNORECASE,
        )
        charge = _number_after(text, ("charge", "-c"), integer=True)
        multiplicity = _number_after(text, ("multiplicity", "-m"), integer=True)
        chemistry: dict[str, Any] = {}
        for key, labels, integer in (
            ("nstates", ("nstates", "excited states"), True),
            ("root", ("root", "target state"), True),
            ("nimages", ("nimages", "images"), True),
            ("num_steps", ("steps",), True),
            ("step_size", ("step size", "increment"), False),
            ("recalc_hess", ("recalc hess", "recalculate hessian every"), True),
            ("trust_radius", ("trust radius",), False),
            ("maxstep", ("maxstep", "maximum step"), True),
            ("num_confs_to_run", ("lowest conformers", "conformers"), True),
        ):
            value = _number_after(text, labels, integer=integer)
            if value is not None:
                chemistry[key] = value
        states = re.search(r"\b(singlets?|triplets?|50-50)\b", lowered)
        if states and kind == "gaussian.tddft":
            chemistry["states"] = states.group(1).rstrip("s") + "s" if states.group(1) != "50-50" else "50-50"
        if "eqsolv" in lowered or "equilibrium solvation" in lowered:
            chemistry["eqsolv"] = True
        direction = re.search(r"\b(forward|reverse|both)\b(?:\s+directions?)?", lowered)
        if direction:
            chemistry["direction"] = direction.group(1)
        bond = re.search(r"\b(?:scan|stretch|vary)\s+(?:the\s+)?bond\s+(\d+)\s*[-–]\s*(\d+)", lowered)
        if bond:
            chemistry["coordinates"] = [int(bond.group(1)), int(bond.group(2))]
        frozen = re.search(r"\b(?:freeze|frozen|keep)\s+(?:the\s+)?bond\s+(\d+)\s*[-–]\s*(\d+)", lowered)
        if frozen:
            chemistry["constrained_coordinates"] = [[int(frozen.group(1)), int(frozen.group(2))]]
        qm_atoms = re.search(r"\b(?:qm|high[- ]level)\s+atoms?\s+([0-9,\-–\s]+)", lowered)
        if qm_atoms:
            chemistry["high_level_atoms"] = re.sub(r"\s+", "", qm_atoms.group(1).replace("–", "-"))
        endpoint_matches = re.findall(r"([\w./-]+\.xyz)\b", text, flags=re.IGNORECASE)
        if len(endpoint_matches) > 1:
            chemistry["ending_xyzfile"] = endpoint_matches[1]
        joboption = re.search(r"\b(neb(?:-ci|-ts)?|fast-neb-ts|tight-neb-ts)\b", lowered)
        if joboption:
            chemistry["joboption"] = joboption.group(1).upper()
        return cls(
            action=action,
            program=program,
            kind=kind,
            project=project_name_from_request(text) or None,
            server=_server_from_request(text),
            input_path=input_match.group(1) if input_match else None,
            charge=charge,
            multiplicity=multiplicity,
            execution_mode="submit" if action == "sub" else "local" if action == "run" else None,
            chemistry=chemistry,
        )

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class ObservedIntent:
    action: str | None
    program: str | None
    kind: str | None
    project: str | None
    server: str | None
    input_path: str | None
    output_path: str | None
    charge: str | None
    multiplicity: str | None
    execution_mode: str | None
    chemistry: dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_command(cls, command: str, *, cwd: str | None = None) -> "ObservedIntent":
        parsed = parse_model_command(command, cwd=cwd)
        try:
            tokens = shlex.split(command)
        except ValueError:
            tokens = []
        job = parsed.job
        kind = f"{parsed.program}.{job}" if parsed.program and job else None
        if parsed.program and "qmmm" in tokens:
            kind = f"{parsed.program}.qmmm"
        if kind == "gaussian.td":
            kind = "gaussian.tddft"
        chemistry = dict(parsed.structural_options)
        if "num_steps" not in chemistry and "num_steps_or_every_n_points" in chemistry:
            chemistry["num_steps"] = chemistry["num_steps_or_every_n_points"]
        if "step_size" not in chemistry and "step_size_or_solv" in chemistry:
            chemistry["step_size"] = chemistry["step_size_or_solv"]
        if parsed.route_parameters:
            chemistry["route_parameters"] = parsed.route_parameters
        if parsed.opt_options:
            chemistry["opt_options"] = parsed.opt_options
            maxstep = re.search(
                r"\bmaxstep\s*=\s*(\d+)", parsed.opt_options, re.IGNORECASE
            )
            if maxstep:
                chemistry["maxstep"] = maxstep.group(1)
        for key, aliases, flag_value in _CHEMISTRY_OPTIONS:
            value = _option_value(tokens, aliases, flag_value=flag_value)
            if value is not None and key not in chemistry:
                chemistry[key] = value
        return cls(
            action=parsed.action,
            program=parsed.program,
            kind=kind,
            project=parsed.project,
            server=parsed.server,
            input_path=parsed.filename,
            output_path=None,
            charge=parsed.charge,
            multiplicity=parsed.multiplicity,
            execution_mode="submit" if parsed.action == "sub" else "local" if parsed.action == "run" else None,
            chemistry=chemistry,
        )

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class IntentAssertion:
    id: str
    expected: Any
    observed: Any
    status: Literal["pass", "fail"]

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class IntentResult:
    verdict: IntentVerdict
    expected: IntentSpec
    observed: ObservedIntent
    assertions: tuple[IntentAssertion, ...]

    @property
    def failed_rule_ids(self) -> list[str]:
        return [row.id for row in self.assertions if row.status == "fail"]

    def to_dict(self) -> dict[str, Any]:
        return {
            "verdict": self.verdict,
            "expected": self.expected.to_dict(),
            "observed": self.observed.to_dict(),
            "assertions": [row.to_dict() for row in self.assertions],
            "failed_rule_ids": self.failed_rule_ids,
        }


def evaluate_intent(
    command: str,
    expected: IntentSpec | dict[str, Any],
    *,
    cwd: str | None = None,
) -> IntentResult:
    spec = expected if isinstance(expected, IntentSpec) else IntentSpec.from_dict(expected)
    observed = ObservedIntent.from_command(command, cwd=cwd)
    rows: list[IntentAssertion] = []
    for field_name in (
        "action",
        "program",
        "kind",
        "project",
        "server",
        "input_path",
        "output_path",
        "charge",
        "multiplicity",
        "execution_mode",
    ):
        expected_value = getattr(spec, field_name)
        if expected_value is None:
            continue
        observed_value = getattr(observed, field_name)
        _append_assertion(rows, f"intent.{field_name}", expected_value, observed_value, path=field_name.endswith("path"))
    for key, expected_value in spec.chemistry.items():
        if expected_value is None:
            continue
        if key == "route_contains":
            route = str(observed.chemistry.get("route_parameters") or "")
            rows.append(
                IntentAssertion(
                    id="intent.chemistry.route_contains",
                    expected=expected_value,
                    observed=route,
                    status="pass"
                    if str(expected_value).lower() in route.lower()
                    else "fail",
                )
            )
        else:
            _append_assertion(rows, f"intent.chemistry.{key}", expected_value, observed.chemistry.get(key), path=key.endswith("file"))
    verdict: IntentVerdict = "reject" if any(row.status == "fail" for row in rows) else "ok"
    return IntentResult(verdict=verdict, expected=spec, observed=observed, assertions=tuple(rows))


def _append_assertion(rows: list[IntentAssertion], rule_id: str, expected: Any, observed: Any, *, path: bool = False) -> None:
    same = _equivalent(expected, observed, path=path)
    rows.append(IntentAssertion(id=rule_id, expected=expected, observed=observed, status="pass" if same else "fail"))


def _equivalent(expected: Any, observed: Any, *, path: bool = False) -> bool:
    if path and expected is not None and observed is not None:
        return str(PurePath(str(expected))) == str(PurePath(str(observed)))
    if isinstance(expected, (list, tuple)):
        return _number_list(expected) == _number_list(observed)
    if isinstance(expected, bool):
        return expected is _as_bool(observed)
    try:
        return float(expected) == float(observed)
    except (TypeError, ValueError):
        left = re.sub(r"[\s\[\](){}]", "", str(expected).lower())
        right = re.sub(r"[\s\[\](){}]", "", str(observed).lower())
        return left == right


def _number_list(value: Any) -> list[str]:
    return re.findall(r"-?\d+(?:\.\d+)?", str(value))


def _as_bool(value: Any) -> bool | None:
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in {"true", "yes", "1", "eqsolv"}:
        return True
    if text in {"false", "no", "0", "noneqsolv"}:
        return False
    return None


def _number_after(text: str, labels: tuple[str, ...], *, integer: bool) -> int | float | None:
    number = r"[-+]?\d+" if integer else r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)"
    for label in labels:
        match = re.search(rf"\b{re.escape(label)}\b\s*(?:=|:|of|is)?\s*({number})", text, flags=re.IGNORECASE)
        if match:
            return int(match.group(1)) if integer else float(match.group(1))
    return None


# The literal marker table below misses common natural phrasings of a
# coordinate scan ("scan the bond", "relaxed PES scan of the C-C bond") and of
# a frozen constraint ("hold the 1-2 bond fixed"), so the derived intent kind
# came back None and a silent collapse to ``opt`` went uncaught. These regexes
# recover those phrasings. The TS guard prevents a "…with the bond frozen"
# transition-state request from being misread as a plain modred.
_TS_INTENT_RE = re.compile(
    r"transition[-\s]*state|saddle|\boptts\b|\bopt[-\s]*ts\b|\bts\b", re.I
)
_SCAN_VERB_RE = re.compile(
    r"\bscan\w*\b|\bpes\b|potential energy surface", re.I
)
_SCAN_INTENT_RE = re.compile(
    r"(?:\bpes\b|potential energy surface|relaxed\s+scan|coordinate\s+scan|\bscan\w*\b)"
    r"[^.]*\b(?:bond|angle|dihedral|coordinate|distance|length|torsion|[a-z][-–][a-z])\b"
    r"|\b(?:bond|angle|dihedral|coordinate|distance|torsion)\b[^.]*\bscan\w*\b",
    re.I,
)
_FREEZE_INTENT_RE = re.compile(
    r"\b(?:freez\w*|constrain\w*|fix\w*|hold\w*|keep\w*|lock\w*)\b"
    r"[^.]*\b(?:bond|angle|dihedral|coordinate|distance|length|torsion)\b"
    r"|\b(?:bond|angle|dihedral|distance|length|torsion)\b"
    r"[^.]*\b(?:fixed|frozen|constrained|held|locked)\b"
    r"|\bmodred\w*\b",
    re.I,
)


def _kind_from_request(lowered: str, program: str | None) -> str | None:
    if not program:
        return None
    if not _TS_INTENT_RE.search(lowered):
        if _SCAN_INTENT_RE.search(lowered):
            return f"{program}.scan"
        if _FREEZE_INTENT_RE.search(lowered) and not _SCAN_VERB_RE.search(
            lowered
        ):
            return f"{program}.modred"
    markers = (
        ("qmmm", ("qm/mm", "qmmm", "oniom")),
        ("tddft", ("td-dft", "tddft", "excited state", "absorption spectrum")),
        ("scan", ("relaxed scan", "coordinate scan", "scan bond", "scan angle", "scan dihedral")),
        ("modred", ("modred", "constrained optimization", "freeze bond")),
        ("neb", ("neb-ts", "neb-ci", "nudged elastic band", " neb ")),
        ("irc", ("irc", "intrinsic reaction coordinate")),
        ("qrc", ("qrc", "quasi-reaction coordinate")),
        ("dias", ("dias", "distortion interaction", "distortion-interaction")),
        ("crest", ("crest", "conformer refinement", "lowest conformers")),
        ("traj", ("trajectory frames", "trajectory refinement")),
        ("resp", ("resp",)),
        ("nci", ("nci", "non-covalent interaction")),
        ("wbi", ("wbi", "wiberg bond")),
        ("ts", ("transition state", "transition-state", "optts", "scants")),
        ("sp", ("single point", "single-point")),
        ("opt", ("optimization", "optimisation", "optimize", "optimise")),
    )
    padded = f" {lowered} "
    for suffix, words in markers:
        if any(word in padded for word in words):
            if suffix == "tddft" and program != "gaussian":
                return None
            if suffix in {"dias", "crest", "traj", "resp", "nci", "wbi"} and program != "gaussian":
                return None
            if suffix == "neb" and program != "orca":
                return None
            return f"{program}.{suffix}"
    return None


def _server_from_request(text: str) -> str | None:
    for pattern in (
        r"(?:^|\s)(?:-s|--server)(?:=|\s+)([A-Za-z0-9_.-]+)",
        r"\b(?:server|cluster)\s+([A-Za-z0-9_.-]+)\b",
        r"\b(?:submit|queue)\s+(?:it\s+)?(?:to|on)\s+([A-Za-z0-9_.-]+)\b",
    ):
        match = re.search(pattern, text, flags=re.IGNORECASE)
        if match:
            return match.group(1)
    return None


def _option_value(tokens: list[str], aliases: tuple[str, ...], *, flag_value: Any = None) -> Any:
    for index, token in enumerate(tokens):
        if token in aliases:
            if flag_value is not None:
                return flag_value
            if index + 1 < len(tokens):
                return tokens[index + 1]
        for alias in aliases:
            prefix = f"{alias}="
            if token.startswith(prefix):
                return token[len(prefix):]
    return None


_CHEMISTRY_OPTIONS = (
    ("states", ("--states",), None),
    ("root", ("--root",), None),
    ("nstates", ("--nstates",), None),
    ("eqsolv", ("--eqsolv",), True),
    ("direction", ("--direction",), None),
    ("recalc_hess", ("--recalc-hess",), None),
    ("trust_radius", ("-t", "--trust-radius"), None),
    ("tssearch_type", ("-ts", "--tssearch-type"), None),
    ("ending_xyzfile", ("-e", "--ending-xyzfile"), None),
    ("nimages", ("--nimages",), None),
    ("joboption", ("--joboption",), None),
    ("constrained_coordinates", ("-cc", "--constrained-coordinates"), None),
    ("fragment_indices", ("-i", "--fragment-indices"), None),
    ("every_n_points", ("-n", "--every-n-points"), None),
    ("mode", ("--mode",), None),
    ("high_level_atoms", ("-ha", "--high-level-atoms"), None),
    ("low_level_atoms", ("-la", "--low-level-atoms"), None),
    ("low_level_method", ("-lm", "--low-level-method", "--low-level-force-field"), None),
    ("num_confs_to_run", ("-nc", "--num-confs-to-run"), None),
    ("num_structures_to_run", ("-ns", "--num-structures-to-run"), None),
    ("grouping_strategy", ("-g", "--grouping-strategy"), None),
)


__all__ = [
    "IntentAssertion",
    "IntentResult",
    "IntentSpec",
    "ObservedIntent",
    "evaluate_intent",
]
