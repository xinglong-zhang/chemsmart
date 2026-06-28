from __future__ import annotations

from chemsmart.agent.harness import spec_invariants as SI


def _workflow(kind: str, settings: dict | None = None) -> dict:
    job = {
        "id": 1,
        "kind": kind,
        "file": "mol.xyz",
        "charge": 0,
        "mult": 1,
    }
    if settings is not None:
        job["settings"] = settings
    return {"intent": "workflow", "jobs": [job]}


def test_noncanonical_opt_freq_kind_rejects():
    issues = SI.check_spec(_workflow("gaussian.opt+freq"), "opt+freq mol.xyz")

    assert any(issue.rule_id == "spec.kind.canonical" for issue in issues)


def test_ts_runtime_route_extra_rejects():
    issues = SI.check_spec(
        _workflow(
            "gaussian.ts",
            {"additional_opt_options_in_route": ["calcfc/noeigentest"]},
        ),
        "TS search in Gaussian from mol.xyz, calcfc/noeigentest",
    )

    assert any(issue.rule_id == "spec.ts.bad_extra" for issue in issues)


def test_ts_true_extra_passes_spec_gate():
    issues = SI.check_spec(
        _workflow("gaussian.ts", {"additional_opt_options_in_route": ["maxstep=8"]}),
        "TS search in Gaussian from mol.xyz with maxstep 8",
    )

    assert not [issue for issue in issues if issue.severity == "reject"]


def test_missing_modred_atoms_requires_decline():
    issues = SI.check_spec(
        _workflow("gaussian.modred", {"modred": [[2, 3]]}),
        "Run Gaussian modredundant job for mol.xyz",
    )

    assert any(issue.rule_id == "spec.decline_contract" for issue in issues)


def test_fragment_indices_with_underscore_counts_as_structural_slot():
    issues = SI.check_spec(
        _workflow("gaussian.dias", {"fragment_indices": [[1], [2]]}),
        "Need a Gaussian DIAS job for mol.xyz; fragment_indices [[1], [2]].",
    )

    assert not [issue for issue in issues if issue.rule_id == "spec.decline_contract"]
