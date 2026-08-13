"""Explicit deterministic preview-verifier adapters for CLI programs."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from chemsmart.agent._contracts import (
    ContractError,
    TrustedArtifactRefV1,
    canonical_data,
    canonical_sha256,
    file_sha256,
)
from chemsmart.agent.projects import ProjectValidationReceiptV1


@dataclass(frozen=True)
class ProgramVerifierAdapterV1:
    program: str
    adapter_id: str
    adapter_version: str
    preview_artifact_contract: str
    adapter_sha256: str


@dataclass(frozen=True)
class ProgramSupportStatusV1:
    program: str
    status: str
    adapter_sha256: str
    rule_ids: tuple[str, ...]


@dataclass(frozen=True)
class PreviewValidationFindingV1:
    rule_id: str
    field: str
    expected: Any
    observed: Any
    evidence_ref: str


@dataclass(frozen=True)
class ProgramPreviewExpectationV1:
    schema_version: str
    program: str
    jobtype: str
    input_artifact: TrustedArtifactRefV1
    project_receipt_sha256: str
    project_sha256: str
    settings: tuple[tuple[str, Any], ...]
    charge: int
    multiplicity: int
    expectation_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.program-preview-expectation.v1":
            raise ContractError("unsupported preview expectation schema")
        body = {
            "schema_version": self.schema_version,
            "program": self.program,
            "jobtype": self.jobtype,
            "input_artifact": self.input_artifact,
            "project_receipt_sha256": self.project_receipt_sha256,
            "project_sha256": self.project_sha256,
            "settings": self.settings,
            "charge": self.charge,
            "multiplicity": self.multiplicity,
        }
        if self.expectation_sha256 != canonical_sha256(body):
            raise ContractError("preview expectation digest mismatch")


@dataclass(frozen=True)
class ProgramPreviewValidationReceiptV1:
    schema_version: str
    expectation_sha256: str
    adapter_sha256: str
    workspace_artifact_set_sha256: str
    findings: tuple[PreviewValidationFindingV1, ...]
    status: str
    receipt_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.program-preview-validation.v1":
            raise ContractError("unsupported preview validation schema")
        if self.status not in {"valid", "invalid"}:
            raise ContractError("invalid preview validation status")
        body = {
            "schema_version": self.schema_version,
            "expectation_sha256": self.expectation_sha256,
            "adapter_sha256": self.adapter_sha256,
            "workspace_artifact_set_sha256": (
                self.workspace_artifact_set_sha256
            ),
            "findings": self.findings,
            "status": self.status,
        }
        if self.receipt_sha256 != canonical_sha256(body):
            raise ContractError("preview validation receipt digest mismatch")


def _adapter(
    program: str, adapter_id: str, contract: str
) -> ProgramVerifierAdapterV1:
    body = {
        "program": program,
        "adapter_id": adapter_id,
        "adapter_version": "1",
        "preview_artifact_contract": contract,
        "module_sha256": file_sha256(__file__),
    }
    return ProgramVerifierAdapterV1(
        program=program,
        adapter_id=adapter_id,
        adapter_version="1",
        preview_artifact_contract=contract,
        adapter_sha256=canonical_sha256(body),
    )


PROGRAM_PREVIEW_VERIFIERS = {
    "gaussian": _adapter(
        "gaussian",
        "gaussian_native_input_v1",
        "parsed .com/.gjf plus geometry",
    ),
    "orca": _adapter(
        "orca", "orca_native_input_v1", "parsed .inp plus geometry"
    ),
    "xtb": _adapter(
        "xtb",
        "xtb_preview_receipt_v1",
        "hash-bound preview receipt plus generated XYZ and canonical argv",
    ),
    "pyscf": _adapter(
        "pyscf", "pyscf_receipt_hdf5_v1", "input/run receipts plus HDF5"
    ),
}


def program_support_status_matrix(
    registry,
) -> tuple[ProgramSupportStatusV1, ...]:
    rows = []
    for capability in registry.programs:
        adapter = PROGRAM_PREVIEW_VERIFIERS.get(capability.program)
        rows.append(
            ProgramSupportStatusV1(
                program=capability.program,
                status="preview_only" if adapter else "reference_only",
                adapter_sha256=adapter.adapter_sha256 if adapter else "",
                rule_ids=(
                    ("agent.verifier.explicit_adapter",)
                    if adapter
                    else ("agent.verifier.adapter_missing",)
                ),
            )
        )
    return tuple(rows)


def build_preview_expectation(
    *,
    program: str,
    jobtype: str,
    input_artifact: TrustedArtifactRefV1,
    project: ProjectValidationReceiptV1 | None,
    charge: int,
    multiplicity: int,
) -> ProgramPreviewExpectationV1:
    body = {
        "schema_version": "chemsmart.program-preview-expectation.v1",
        "program": program,
        "jobtype": jobtype,
        "input_artifact": input_artifact,
        "project_receipt_sha256": project.receipt_sha256 if project else "",
        "project_sha256": project.project_sha256 if project else "",
        "settings": project.settings if project else (),
        "charge": int(charge),
        "multiplicity": int(multiplicity),
    }
    return ProgramPreviewExpectationV1(
        **body, expectation_sha256=canonical_sha256(body)
    )


def validate_preview_workspace(
    expectation: ProgramPreviewExpectationV1,
    workspace: str | Path,
) -> ProgramPreviewValidationReceiptV1:
    adapter = PROGRAM_PREVIEW_VERIFIERS.get(expectation.program)
    if adapter is None:
        raise ContractError("program has no explicit preview verifier adapter")
    root = Path(workspace)
    paths = tuple(sorted(path for path in root.rglob("*") if path.is_file()))
    findings = []
    try:
        if expectation.program in {"gaussian", "orca"}:
            findings.extend(_validate_native_input(expectation, paths))
        elif expectation.program == "xtb":
            findings.extend(_validate_xtb_preview(expectation, paths))
        elif expectation.program == "pyscf":
            findings.extend(_validate_pyscf_preview(expectation, paths))
    except Exception as exc:
        findings.append(
            PreviewValidationFindingV1(
                rule_id="preview.verifier.exception",
                field="adapter",
                expected="deterministic parser completion",
                observed=type(exc).__name__,
                evidence_ref=f"adapter:{adapter.adapter_id}",
            )
        )
    artifact_rows = tuple(
        (
            path.relative_to(root).as_posix(),
            path.stat().st_size,
            file_sha256(path),
        )
        for path in paths
    )
    body = {
        "schema_version": "chemsmart.program-preview-validation.v1",
        "expectation_sha256": expectation.expectation_sha256,
        "adapter_sha256": adapter.adapter_sha256,
        "workspace_artifact_set_sha256": canonical_sha256(artifact_rows),
        "findings": tuple(findings),
        "status": "valid" if not findings else "invalid",
    }
    return ProgramPreviewValidationReceiptV1(
        **body, receipt_sha256=canonical_sha256(body)
    )


def _validate_native_input(expectation, paths):
    suffixes = (
        (".com", ".gjf") if expectation.program == "gaussian" else (".inp",)
    )
    candidates = [path for path in paths if path.suffix.lower() in suffixes]
    if not candidates:
        return [_missing("native_input", suffixes, "workspace")]
    settings_cls = _settings_class(expectation.program)
    expected_settings = dict(expectation.settings)
    expected_settings.update(
        {
            "charge": expectation.charge,
            "multiplicity": expectation.multiplicity,
        }
    )
    parsed_candidates = [
        (path, settings_cls.from_filepath(str(path))) for path in candidates
    ]
    if expectation.program == "gaussian" and expectation.jobtype == "link":
        return _validate_gaussian_link_input(
            expectation,
            parsed_candidates,
            expected_settings=expected_settings,
        )
    if expectation.program == "gaussian" and expectation.jobtype == "td":
        return _validate_gaussian_td_input(
            expectation,
            parsed_candidates,
            expected_settings=expected_settings,
        )
    if expectation.program == "gaussian" and expectation.jobtype == "irc":
        return _validate_gaussian_irc_bundle(
            expectation,
            parsed_candidates,
            expected_settings=expected_settings,
        )
    matches = []
    for path, parsed in parsed_candidates:
        matches.append(
            _settings_match(
                parsed,
                expected_settings,
                native_input=path if expectation.program == "orca" else None,
            )
        )
    findings = [] if any(not item for item in matches) else matches[0]
    if findings:
        return findings
    if not _geometry_sets_equal(expectation.input_artifact.path, candidates):
        return [
            PreviewValidationFindingV1(
                "preview.geometry.mismatch",
                "geometry",
                expectation.input_artifact.sha256,
                "generated geometry differs",
                "generated:native_input",
            )
        ]
    return []


def _validate_gaussian_td_input(
    expectation,
    parsed_candidates,
    *,
    expected_settings,
):
    """Validate Gaussian TD semantics from the generated native route."""

    if len(parsed_candidates) != 1:
        return [
            _mismatch(
                "native_input_count",
                1,
                len(parsed_candidates),
                "generated:native_input_bundle",
            )
        ]
    import re

    path, parsed = parsed_candidates[0]
    route = str(getattr(parsed, "route_string", "") or "").casefold()
    match = re.search(r"\btd\s*\(([^)]*)\)", route)
    if match is None:
        return [_missing("td_route", "TD(...) route", path.name)]
    tokens = tuple(
        item.strip().casefold()
        for item in match.group(1).split(",")
        if item.strip()
    )
    values = {}
    flags = set()
    for token in tokens:
        if "=" in token:
            key, value = token.split("=", 1)
            values[key.strip()] = value.strip()
        else:
            flags.add(token)
    findings = []
    states = str(expected_settings.get("states") or "").strip().casefold()
    if states and states not in flags:
        findings.append(
            _mismatch("states", states, tuple(sorted(flags)), path.name)
        )
    for field in ("nstates", "root"):
        expected = expected_settings.get(field)
        if expected is not None and values.get(field) != str(expected):
            findings.append(
                _mismatch(field, int(expected), values.get(field), path.name)
            )
    eqsolv = str(expected_settings.get("eqsolv") or "").strip().casefold()
    if eqsolv and eqsolv not in flags:
        findings.append(
            _mismatch("eqsolv", eqsolv, tuple(sorted(flags)), path.name)
        )
    target_settings = dict(expected_settings)
    for field in ("eqsolv", "nstates", "root", "states"):
        target_settings.pop(field, None)
    # Gaussian's base parser calls a TD route an SP; the TD leaf itself is
    # established by the explicit route semantics above.
    target_settings.pop("jobtype", None)
    findings.extend(_settings_match(parsed, target_settings))
    if not _geometry_sets_equal(expectation.input_artifact.path, [path]):
        findings.append(
            PreviewValidationFindingV1(
                "preview.geometry.mismatch",
                "geometry",
                expectation.input_artifact.sha256,
                "generated geometry differs",
                "generated:native_input",
            )
        )
    return findings


def _validate_gaussian_link_input(
    expectation,
    parsed_candidates,
    *,
    expected_settings,
):
    """Validate both scientific route sections of one Gaussian link input."""

    if len(parsed_candidates) != 1:
        return [
            _mismatch(
                "native_input_count",
                1,
                len(parsed_candidates),
                "generated:native_input_bundle",
            )
        ]
    path, parsed_target = parsed_candidates[0]
    text = path.read_text(encoding="utf-8", errors="replace")
    if "--Link1--" not in text:
        return [
            _mismatch(
                "link_sections",
                "stability route followed by --Link1-- target route",
                "single route",
                path.name,
            )
        ]
    first_text, second_text = text.split("--Link1--", 1)

    def _route_block(value):
        lines = value.splitlines()
        for index, line in enumerate(lines):
            if not line.lstrip().startswith("#"):
                continue
            block = [line.strip()]
            for continuation in lines[index + 1 :]:
                if not continuation.strip():
                    break
                block.append(continuation.strip())
            return " ".join(" ".join(block).casefold().split())
        return ""

    first_route = _route_block(first_text)
    second_route = _route_block(second_text)
    findings = []
    stable = str(expected_settings.get("stable") or "").strip().casefold()
    if stable and f"stable={stable}" not in first_route:
        findings.append(
            _mismatch(
                "stable",
                stable,
                first_route or "missing route",
                path.name,
            )
        )
    guess = str(expected_settings.get("guess") or "").strip().casefold()
    normalized_guess = guess.strip("()")
    if guess and not any(
        token in first_route
        for token in (
            f"guess={guess}",
            f"guess=({normalized_guess})",
        )
    ):
        findings.append(
            _mismatch(
                "guess",
                guess,
                first_route or "missing route",
                path.name,
            )
        )
    link_route = str(expected_settings.get("link_route") or "").strip()
    if link_route:
        required_tokens = {
            token.casefold() for token in link_route.split() if token.strip()
        }
        observed_tokens = set(second_route.split())
        if not required_tokens.issubset(observed_tokens):
            findings.append(
                _mismatch(
                    "link_route",
                    tuple(sorted(required_tokens)),
                    tuple(sorted(observed_tokens)),
                    path.name,
                )
            )

    target_settings = dict(expected_settings)
    for field in ("guess", "link_route", "stable"):
        target_settings.pop(field, None)
    findings.extend(_settings_match(parsed_target, target_settings))
    if not _geometry_sets_equal(expectation.input_artifact.path, [path]):
        findings.append(
            PreviewValidationFindingV1(
                "preview.geometry.mismatch",
                "geometry",
                expectation.input_artifact.sha256,
                "generated geometry differs",
                "generated:native_input",
            )
        )
    return findings


def _validate_gaussian_irc_bundle(
    expectation,
    parsed_candidates,
    *,
    expected_settings,
):
    """Validate one canonical IRC node as forward and reverse native inputs."""

    directions = tuple(
        sorted(
            getattr(parsed, "jobtype", "")
            for _path, parsed in parsed_candidates
        )
    )
    if directions != ("ircf", "ircr"):
        return [
            _mismatch(
                "jobtype_bundle",
                ("ircf", "ircr"),
                directions,
                "generated:native_input_bundle",
            )
        ]

    shared_expected = dict(expected_settings)
    shared_expected.pop("jobtype", None)
    findings = []
    for path, parsed in parsed_candidates:
        findings.extend(_settings_match(parsed, shared_expected))
        if not _geometry_sets_equal(expectation.input_artifact.path, [path]):
            findings.append(
                PreviewValidationFindingV1(
                    "preview.geometry.mismatch",
                    "geometry",
                    expectation.input_artifact.sha256,
                    f"generated geometry differs in {path.name}",
                    path.name,
                )
            )
    return findings


def _validate_xtb_preview(expectation, paths):
    receipts = []
    for path in paths:
        if path.name != "xtb-preview-receipt-v1.json":
            continue
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            continue
        if payload.get("schema_version") == "chemsmart.xtb-preview.v1":
            receipts.append((path, payload))
    if len(receipts) != 1:
        return [
            _missing("preview_receipt", "one hash-bound receipt", "workspace")
        ]
    receipt_path, payload = receipts[0]
    from chemsmart.jobs.xtb.validation import canonical_sha256 as xtb_sha256

    receipt_body = dict(payload)
    observed_receipt_sha256 = str(receipt_body.pop("receipt_sha256", ""))
    findings = []
    if observed_receipt_sha256 != xtb_sha256(receipt_body):
        findings.append(
            _mismatch(
                "preview_receipt.digest",
                "valid",
                observed_receipt_sha256,
                receipt_path.name,
            )
        )
    for field, expected_value in (
        ("state", "previewed"),
        ("executed", False),
        ("execution_ready", False),
        ("program", "xtb"),
        ("jobtype", expectation.jobtype),
        ("required_version", "6.7.1"),
    ):
        if payload.get(field) != expected_value:
            findings.append(
                _mismatch(
                    field,
                    expected_value,
                    payload.get(field),
                    receipt_path.name,
                )
            )
    input_record = payload.get("input_artifact")
    generated_xyz = None
    if isinstance(input_record, dict):
        candidates = [
            path
            for path in paths
            if path.name == input_record.get("path")
            and path.suffix.lower() == ".xyz"
        ]
        generated_xyz = candidates[0] if len(candidates) == 1 else None
    if generated_xyz is None:
        findings.append(
            _missing("geometry", "receipt-bound XYZ", receipt_path.name)
        )
    else:
        if input_record.get(
            "size"
        ) != generated_xyz.stat().st_size or input_record.get(
            "sha256"
        ) != file_sha256(generated_xyz):
            findings.append(
                _mismatch(
                    "input_artifact",
                    "matching size and SHA-256",
                    input_record,
                    receipt_path.name,
                )
            )
        elif not _geometry_sets_equal(
            expectation.input_artifact.path, [generated_xyz]
        ):
            findings.append(
                _mismatch(
                    "geometry",
                    expectation.input_artifact.sha256,
                    input_record.get("sha256"),
                    receipt_path.name,
                )
            )
    expected = dict(expectation.settings)
    expected.update(
        {
            "charge": expectation.charge,
            "multiplicity": expectation.multiplicity,
            "jobtype": expectation.jobtype,
        }
    )
    from chemsmart.jobs.xtb.runner import XTBJobRunner
    from chemsmart.jobs.xtb.settings import XTBJobSettings

    xtb_settings = {
        key: value
        for key, value in expected.items()
        if key in XTBJobSettings.FIELDS
    }
    required = XTBJobRunner._settings_args(
        XTBJobSettings.from_dict(xtb_settings)
    )
    observed_argv = payload.get("canonical_argv")
    if not isinstance(observed_argv, list) or not observed_argv:
        findings.append(
            _missing("canonical_argv", "receipt argv", receipt_path.name)
        )
    elif (
        observed_argv[0] != "xtb"
        or generated_xyz is None
        or len(observed_argv) < 2
        or Path(observed_argv[1]).resolve() != generated_xyz.resolve()
        or not _contains_subsequence(observed_argv, required)
    ):
        findings.append(
            PreviewValidationFindingV1(
                "preview.xtb.argv_mismatch",
                "xtb_argv",
                tuple(required),
                tuple(observed_argv),
                receipt_path.name,
            )
        )
    expected_settings = {
        key: canonical_data(value)
        for key, value in sorted(xtb_settings.items())
    }
    observed_settings = payload.get("settings")
    for key, value in expected_settings.items():
        if (
            not isinstance(observed_settings, dict)
            or canonical_data(observed_settings.get(key)) != value
        ):
            findings.append(
                _mismatch(
                    "settings." + key,
                    value,
                    (
                        observed_settings.get(key)
                        if isinstance(observed_settings, dict)
                        else "missing"
                    ),
                    receipt_path.name,
                )
            )
    if any(path.suffix.lower() in {".out", ".err"} for path in paths):
        findings.append(
            _mismatch(
                "completion_artifacts",
                "absent in preview",
                "present",
                "workspace",
            )
        )
    return findings


def _validate_pyscf_preview(expectation, paths):
    receipts = []
    for path in paths:
        if path.suffix.lower() != ".json":
            continue
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            continue
        if isinstance(payload, dict):
            receipts.append((path, payload))
    run = next(
        (
            item
            for item in receipts
            if item[1].get("schema_version") == "chemsmart.pyscf-run.v1"
        ),
        None,
    )
    if run is None:
        return [_missing("run_receipt", "chemsmart.pyscf-run.v1", "workspace")]
    path, payload = run
    observed_digest = payload.pop("receipt_sha256", "")
    findings = []
    if observed_digest != canonical_sha256(payload):
        findings.append(
            _mismatch("run_receipt.digest", "valid", "invalid", path.name)
        )
    if payload.get("state") != "previewed" or payload.get("fake") is not True:
        findings.append(
            _mismatch(
                "run_receipt.state",
                "previewed fake",
                payload.get("state"),
                path.name,
            )
        )
    if payload.get("project_yaml_sha256") != expectation.project_sha256:
        findings.append(
            _mismatch(
                "project_yaml_sha256",
                expectation.project_sha256,
                payload.get("project_yaml_sha256"),
                path.name,
            )
        )
    if payload.get("findings"):
        findings.append(
            _mismatch(
                "run_receipt.findings", [], payload.get("findings"), path.name
            )
        )
    result = next(
        (item for item in paths if item.suffix.lower() == ".h5"), None
    )
    if result is None or payload.get("result_sha256") != file_sha256(result):
        findings.append(
            _mismatch(
                "result_sha256",
                payload.get("result_sha256"),
                "missing_or_mismatch",
                path.name,
            )
        )
    return findings


def _settings_class(program):
    if program == "gaussian":
        from chemsmart.jobs.gaussian.settings import GaussianJobSettings

        return GaussianJobSettings
    from chemsmart.jobs.orca.settings import ORCAJobSettings

    return ORCAJobSettings


def _settings_match(parsed, expected, *, native_input=None):
    findings = []
    route_tokens: set[str] = set()
    from chemsmart.jobs.gaussian.settings import GaussianJobSettings
    from chemsmart.jobs.orca.settings import (
        ORCAJobSettings,
        _normalize_orca_functional,
    )

    is_gaussian = isinstance(parsed, GaussianJobSettings)
    is_orca = isinstance(parsed, ORCAJobSettings)
    if native_input is not None:
        from chemsmart.io.orca.input import ORCAInput

        route_tokens = set(
            ORCAInput(str(native_input)).route_object.route_keywords
        )
    for field, value in sorted(expected.items()):
        if value is None:
            continue
        if native_input is not None and field == "basis":
            if str(value).strip().lower() in route_tokens:
                continue
        if native_input is not None and field == "additional_route_parameters":
            required = {
                item.lower() for item in str(value).split() if item.strip()
            }
            if required.issubset(route_tokens):
                continue
        if not hasattr(parsed, field):
            findings.append(
                _mismatch(
                    field,
                    value,
                    "missing_from_parsed_native_input",
                    "generated:native_input",
                )
            )
            continue
        observed = getattr(parsed, field)
        if is_orca and field == "functional":
            expected_functional = _normalize_orca_functional(value)
            observed_functional = _normalize_orca_functional(observed)
            if (
                str(expected_functional).strip().casefold()
                == str(observed_functional).strip().casefold()
            ):
                continue
        if is_orca and field == "scf_tol":

            def _orca_scf_preset(item):
                preset = str(item or "").strip().casefold()
                return preset[:-3] if preset.endswith("scf") else preset

            if _orca_scf_preset(value) == _orca_scf_preset(observed):
                continue
        if is_gaussian and field == "dispersion":
            from chemsmart.io.gaussian.route import (
                normalize_gaussian_dispersion,
            )

            if normalize_gaussian_dispersion(
                value
            ) == normalize_gaussian_dispersion(observed):
                continue
        if is_gaussian and field in {
            "basis",
            "high_level_basis",
            "medium_level_basis",
            "low_level_basis",
        }:
            from chemsmart.jobs.gaussian.settings import (
                gaussian_native_basis_token,
            )

            expected_basis = gaussian_native_basis_token(value)
            observed_basis = gaussian_native_basis_token(observed)
            if (
                str(expected_basis).strip().casefold()
                == str(observed_basis).strip().casefold()
            ):
                continue
        if is_gaussian and field == "additional_route_parameters":
            from chemsmart.io.gaussian.route import (
                normalize_gaussian_dispersion,
                split_gaussian_dispersion_tokens,
            )

            expected_route, expected_dispersion = (
                split_gaussian_dispersion_tokens(value)
            )
            observed_route, observed_dispersion = (
                split_gaussian_dispersion_tokens(observed)
            )
            parsed_dispersion = normalize_gaussian_dispersion(
                getattr(parsed, "dispersion", None)
            )
            semantic_dispersion = observed_dispersion or parsed_dispersion
            if (
                expected_dispersion is not None
                and expected_dispersion != semantic_dispersion
            ):
                findings.append(
                    _mismatch(
                        "dispersion",
                        expected_dispersion,
                        semantic_dispersion,
                        "generated:native_input",
                    )
                )
                continue
            value = expected_route or None
            observed = observed_route or None
            if value is None and observed is None:
                continue
        if is_gaussian and field == "functional" and isinstance(value, str):
            from chemsmart.io.gaussian.route import (
                gaussian_functional_without_dispersion_shorthand,
                normalize_gaussian_dispersion,
                split_gaussian_dispersion_tokens,
                split_gaussian_functional_dispersion_shorthand,
            )

            expected_functional, embedded_dispersion = (
                split_gaussian_dispersion_tokens(value)
            )
            expected_functional, shorthand_dispersion = (
                split_gaussian_functional_dispersion_shorthand(
                    expected_functional
                )
            )
            parsed_dispersion = normalize_gaussian_dispersion(
                getattr(parsed, "dispersion", None)
            )
            expected_dispersions = {
                item
                for item in (embedded_dispersion, shorthand_dispersion)
                if item is not None
            }
            if len(expected_dispersions) > 1:
                findings.append(
                    _mismatch(
                        "dispersion",
                        tuple(sorted(expected_dispersions)),
                        parsed_dispersion,
                        "generated:native_input",
                    )
                )
                continue
            expected_dispersion = next(iter(expected_dispersions), None)
            if (
                expected_dispersion is not None
                and expected_dispersion != parsed_dispersion
            ):
                findings.append(
                    _mismatch(
                        "dispersion",
                        expected_dispersion,
                        parsed_dispersion,
                        "generated:native_input",
                    )
                )
                continue
            expected_base = gaussian_functional_without_dispersion_shorthand(
                expected_functional,
                parsed_dispersion,
            )
            observed_base = gaussian_functional_without_dispersion_shorthand(
                observed,
                parsed_dispersion,
            )
            if (
                str(expected_base).strip().casefold()
                == str(observed_base).strip().casefold()
            ):
                continue
        # A method or basis name means the same chemistry in any case, in
        # every program ChemSmart drives -- and the observed side is read
        # back by ChemSmart's own parser, which lowercases. Comparing case
        # therefore tested the parser's normalisation, not the chemistry, and
        # no spelling could pass: a live session wrote the standard
        # ``6-31G(d)``, ChemSmart's writer emitted it, its reader returned
        # ``6-31g(d)``, and the preview was refused as
        # ``preview.semantic.mismatch``. The session recompiled six times and
        # abandoned Gaussian and ORCA. This rule was already applied to ORCA
        # alone, which is why only Gaussian's failure looked mysterious.
        if isinstance(value, str) and isinstance(observed, str):
            if value.strip().casefold() == observed.strip().casefold():
                continue
        if native_input is not None:
            # An absent opt-in boolean is the native representation of
            # ``false`` in ORCA's simple input.
            if value is False and observed in {None, False}:
                continue
        if canonical_data(observed) != canonical_data(value):
            findings.append(
                _mismatch(field, value, observed, "generated:native_input")
            )
    return findings


def _geometry_sets_equal(source, generated):
    from chemsmart.io.molecules.structure import Molecule

    source_molecules = Molecule.from_filepath(
        source, index=":", return_list=True
    )
    source_rows = sorted(_molecule_identity(item) for item in source_molecules)
    generated_rows = []
    for path in generated:
        molecules = Molecule.from_filepath(
            str(path), index=":", return_list=True
        )
        generated_rows.extend(_molecule_identity(item) for item in molecules)
    return source_rows == sorted(generated_rows)


def _molecule_identity(molecule):
    positions = getattr(molecule, "positions")
    return canonical_sha256(
        {
            "symbols": list(getattr(molecule, "symbols")),
            "positions": [
                [round(float(value), 10) for value in row] for row in positions
            ],
        }
    )


def _contains_subsequence(values, required):
    if len(required) > len(values):
        return False
    return any(
        values[index : index + len(required)] == required
        for index in range(len(values) - len(required) + 1)
    )


def _missing(field, expected, evidence):
    return PreviewValidationFindingV1(
        "preview.artifact.missing", field, expected, "missing", evidence
    )


def _mismatch(field, expected, observed, evidence):
    return PreviewValidationFindingV1(
        "preview.semantic.mismatch",
        field,
        canonical_data(expected),
        canonical_data(observed),
        evidence,
    )


__all__ = [
    "PROGRAM_PREVIEW_VERIFIERS",
    "ProgramPreviewExpectationV1",
    "ProgramPreviewValidationReceiptV1",
    "ProgramSupportStatusV1",
    "ProgramVerifierAdapterV1",
    "build_preview_expectation",
    "program_support_status_matrix",
    "validate_preview_workspace",
]
