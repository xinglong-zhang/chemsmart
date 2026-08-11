from __future__ import annotations

import yaml

from chemsmart.agent._contracts import canonical_sha256
from chemsmart.agent.projects import (
    ProjectValidationReceiptV1,
    project_document,
    project_scientific_materializations,
    render_project_yaml,
)
from chemsmart.agent.tool_runtime import (
    _scientific_decision_binding_requirement,
)


def _validation_receipt(settings):
    body = {
        "schema_version": "chemsmart.project-validation-receipt.v1",
        "project_artifact_id": "project.test",
        "project_sha256": "a" * 64,
        "capability_receipt_sha256": "b" * 64,
        "program": "pyscf",
        "jobtype": "sp",
        "loader_id": "chemsmart.settings.pyscf.YamlPySCFProjectSettings",
        "settings_sha256": canonical_sha256(dict(settings)),
        "settings": tuple(sorted(settings.items())),
        "status": "valid",
        "error_class": "",
        "diagnostic": "",
        "rule_ids": ("project.loader.valid",),
    }
    return ProjectValidationReceiptV1(
        **body, receipt_sha256=canonical_sha256(body)
    )


def test_agent_render_materializes_host_owned_pyscf_stage_defaults():
    document = project_document(
        program="pyscf",
        sections={
            "opt": {"functional": "b3lyp", "basis": "def2-svp"},
        },
    )

    receipt = render_project_yaml(document)
    rendered = yaml.safe_load(receipt.rendered_yaml)

    assert rendered == {
        "opt": {
            "basis": "def2-svp",
            "density_fit": False,
            "engine": "cpu",
            "freq": False,
            "functional": "b3lyp",
            "opt_maxsteps": 100,
            "opt_solver": "geometric",
            "scf_maxiter": None,
            "scf_tol": None,
        }
    }
    assert "project.render.host_effective_settings" in receipt.rule_ids
    assert receipt.status == "candidate_rendered"


def test_project_validation_exposes_resolution_without_changing_settings():
    validation = _validation_receipt(
        {"basis": "def2-svp", "functional": "b3lyp"}
    )

    (resolution,) = project_scientific_materializations(validation)

    assert dict(validation.settings)["functional"] == "b3lyp"
    assert resolution.requested_literal == "b3lyp"
    assert resolution.applied_xc == "b3lypg"
    assert resolution.correlation_convention == "vwn3_gaussian"
    assert resolution.evidence_ref == (
        f"functional_resolution:{resolution.receipt_sha256}"
    )
    assert resolution.project_validation_receipt_sha256 == (
        validation.receipt_sha256
    )

    binding = _scientific_decision_binding_requirement((resolution,))
    assert binding["next_tool"] == "record_scientific_decision"
    assert binding["evidence_refs"] == (resolution.evidence_ref,)
    assert binding["status"] == (
        "required_if_rendering_implementation_semantics"
    )


def test_hf_resolution_is_not_a_fictitious_xc_literal():
    validation = _validation_receipt(
        {"ab_initio": "hf", "basis": "def2-svp", "functional": None}
    )

    (resolution,) = project_scientific_materializations(validation)

    assert resolution.requested_method_kind == "hf"
    assert resolution.applied_xc is None
    assert resolution.correlation_convention == "not_applicable"
