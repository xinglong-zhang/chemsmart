"""Project validation must expose loader overrides of scientific settings."""

from __future__ import annotations

import hashlib
import re
from types import SimpleNamespace
from unittest.mock import patch

from chemsmart.agent._contracts import TrustedArtifactRefV1
from chemsmart.agent.capabilities import CapabilityQueryStatus
from chemsmart.agent.projects import (
    project_document,
    project_section_application_observation,
)
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1
from chemsmart.io.gaussian.route import GaussianRoute
from chemsmart.settings.gaussian import GaussianProjectSettings
from chemsmart.settings.capabilities import PROGRAM_CAPABILITIES


def test_loader_case_normalization_is_not_reported_as_scientific_override():
    document = project_document(
        program="xtb",
        sections={
            "opt": {
                "gfn_version": " GFN2 ",
                "optimization_level": "VTIGHT",
            }
        },
    )

    observation = project_section_application_observation(
        document,
        jobtype="opt",
        applied_settings={
            "gfn_version": "gfn2",
            "optimization_level": "vtight",
        },
    )

    assert observation["status"] == "effective_settings_present"
    assert "overridden_settings" not in observation
    assert "next_action" not in observation


def _artifact(path) -> TrustedArtifactRefV1:
    payload = path.read_bytes()
    return TrustedArtifactRefV1(
        artifact_id="project.gaussian.freq",
        kind="project_yaml",
        sha256=hashlib.sha256(payload).hexdigest(),
        size_bytes=len(payload),
        path=str(path.resolve()),
        cli_value=str(path.resolve()),
    )


def _capability(jobtype: str):
    return SimpleNamespace(
        query=SimpleNamespace(program="gaussian", jobtype=jobtype),
        status=CapabilityQueryStatus.SUPPORTED,
        capability=PROGRAM_CAPABILITIES["gaussian"],
        receipt_sha256=("a" if jobtype == "sp" else "b") * 64,
    )


def _validate_for(artifact, jobtype: str):
    capability = _capability(jobtype)
    host = object.__new__(CommandCompiledToolHostV1)
    host.artifacts = {artifact.artifact_id: artifact}
    host.capabilities = {capability.receipt_sha256: capability}
    host.project_validations = {}
    host.functional_resolutions = {}
    host.project_promotions = {}
    host._emit = lambda *_args, **_kwargs: None
    host._project_workflow_binding_observation = lambda _artifact_id: {
        "status": "no_scientific_workflow_planned"
    }
    result = host._validate_project_yaml(
        "turn",
        {
            "project_artifact_id": artifact.artifact_id,
            "capability_receipt_sha256": capability.receipt_sha256,
        },
    )
    return result, result["section_application"]


def test_gas_frequency_is_applied_to_opt_but_not_silently_to_sp(tmp_path):
    path = tmp_path / "gaussian.yaml"
    path.write_text(
        "gas:\n  functional: b3lyp\n  basis: def2-tzvp\n  freq: true\n",
        encoding="utf-8",
    )
    artifact = _artifact(path)

    opt_receipt, opt_observation = _validate_for(artifact, "opt")
    sp_receipt, sp_observation = _validate_for(artifact, "sp")

    assert opt_receipt["status"] == sp_receipt["status"] == "valid"
    assert dict(opt_receipt["settings"])["freq"] is True
    assert opt_observation["status"] == "effective_settings_present"
    assert "overridden_settings" not in opt_observation

    assert dict(sp_receipt["settings"])["freq"] is False
    assert sp_observation["status"] == "declared_settings_overridden"
    assert sp_observation["overridden_settings"] == (
        {
            "setting_name": "freq",
            "declared_value": True,
            "applied_value": False,
        },
    )
    assert "explicit 'sp:' section" in sp_observation["next_action"]


def test_td_does_not_inherit_the_historical_frequency_default(tmp_path):
    path = tmp_path / "gaussian-td.yaml"
    path.write_text(
        "td:\n"
        "  functional: CAM-B3LYP\n"
        "  basis: aug-cc-pVDZ\n"
        "  states: singlets\n"
        "  nstates: 5\n",
        encoding="utf-8",
    )

    settings = GaussianProjectSettings.from_project(str(path)).td_settings()

    assert settings.freq is False
    assert _frequency_route_tokens(settings.route_string) == ()
    assert "td(singlets,nstates=5" in settings.route_string.lower()


def test_explicit_td_frequency_is_not_silently_removed(tmp_path):
    path = tmp_path / "gaussian-td-frequency.yaml"
    path.write_text(
        "td:\n  functional: CAM-B3LYP\n  basis: aug-cc-pVDZ\n  freq: true\n",
        encoding="utf-8",
    )

    settings = GaussianProjectSettings.from_project(str(path)).td_settings()

    assert settings.freq is True
    assert _frequency_route_tokens(settings.route_string) == ("freq",)


def test_explicit_sp_frequency_remains_applied_and_needs_no_repair(tmp_path):
    path = tmp_path / "gaussian-explicit-sp.yaml"
    path.write_text(
        "gas:\n"
        "  functional: b3lyp\n"
        "  basis: def2-tzvp\n"
        "  freq: true\n"
        "sp:\n"
        "  freq: true\n",
        encoding="utf-8",
    )
    artifact = _artifact(path)

    receipt, observation = _validate_for(artifact, "sp")

    assert receipt["status"] == "valid"
    assert dict(receipt["settings"])["freq"] is True
    assert observation["status"] == "effective_settings_present"
    assert "overridden_settings" not in observation


def _frequency_route_tokens(route: str) -> tuple[str, ...]:
    return tuple(
        re.findall(r"(?i)(?<![a-z0-9_])freq(?:=numer)?(?![a-z0-9_])", route)
    )


def test_explicit_sp_frequency_reaches_native_route_and_result_validation(
    tmp_path,
):
    explicit_path = tmp_path / "gaussian-explicit-sp.yaml"
    explicit_path.write_text(
        "gas:\n"
        "  functional: b3lyp\n"
        "  basis: def2-tzvp\n"
        "  freq: true\n"
        "sp:\n"
        "  freq: true\n",
        encoding="utf-8",
    )
    ordinary_path = tmp_path / "gaussian-ordinary-sp.yaml"
    ordinary_path.write_text(
        "gas:\n  functional: b3lyp\n  basis: def2-tzvp\n  freq: true\n",
        encoding="utf-8",
    )

    explicit_settings = GaussianProjectSettings.from_project(
        str(explicit_path)
    ).sp_settings()
    ordinary_settings = GaussianProjectSettings.from_project(
        str(ordinary_path)
    ).sp_settings()
    explicit_route = explicit_settings.route_string
    ordinary_route = ordinary_settings.route_string

    assert explicit_settings.freq is True
    assert _frequency_route_tokens(explicit_route) == ("freq",)
    assert GaussianRoute(explicit_route).jobtype == "sp"
    assert ordinary_settings.freq is False
    assert _frequency_route_tokens(ordinary_route) == ()
    assert GaussianRoute(ordinary_route).jobtype == "sp"

    output_path = tmp_path / "fixed-geometry-frequency.log"
    output_path.write_text("bound Gaussian result", encoding="utf-8")
    output_artifact = TrustedArtifactRefV1(
        artifact_id="result.gaussian.fixed-geometry-frequency",
        kind="gaussian_output",
        sha256=hashlib.sha256(output_path.read_bytes()).hexdigest(),
        size_bytes=output_path.stat().st_size,
        path=str(output_path.resolve()),
        cli_value=str(output_path.resolve()),
    )
    parsed_output = SimpleNamespace(
        energies=(-204.0,),
        vibrational_frequencies=(750.0, 1320.0, 1610.0),
        tddft_transitions=(),
        wavefunction_stability_history=(),
        jobtype=GaussianRoute(explicit_route).jobtype,
        contents=("Normal termination of Gaussian",),
        convergence_criterion_not_met=False,
        normal_termination=True,
        charge=0,
        multiplicity=2,
        method="b3lyp",
        basis="def2-tzvp",
    )
    with patch(
        "chemsmart.io.gaussian.output.Gaussian16Output",
        return_value=parsed_output,
    ):
        evaluation = CommandCompiledToolHostV1._evaluate_execution_outputs(
            program="gaussian",
            jobtype="sp",
            charge=0,
            multiplicity=2,
            expected_settings={
                "jobtype": "sp",
                "functional": "b3lyp",
                "basis": "def2-tzvp",
                "freq": True,
            },
            output_artifacts=(output_artifact,),
            exit_status=0,
        )

    assert evaluation.validated is True
    assert evaluation.findings == ()
    assert evaluation.observations["gaussian"]["outputs"][0]["jobtype"] == "sp"
