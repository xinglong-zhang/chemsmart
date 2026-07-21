from __future__ import annotations

from chemsmart.agent.harness import spec_invariants as SI


def _workflow(
    kind: str,
    settings: dict | None = None,
    **job_overrides,
) -> dict:
    job = {
        "id": 1,
        "kind": kind,
        "file": "mol.xyz",
        "charge": 0,
        "mult": 1,
    }
    job.update(job_overrides)
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
        _workflow(
            "gaussian.ts", {"additional_opt_options_in_route": ["maxstep=8"]}
        ),
        "TS search in Gaussian from mol.xyz with maxstep 8",
    )

    assert not [issue for issue in issues if issue.severity == "reject"]


def test_gaussian_scan_rejects_constraint_in_route_options():
    issues = SI.check_spec(
        _workflow(
            "gaussian.scan",
            {
                "scan_definition": "B 1 2 S 10 0.05",
                "additional_opt_options_in_route": "B 1 3 F",
            },
        ),
        "Scan bond 1-2 and freeze bond 1-3 in mol.xyz",
    )

    assert any(
        issue.rule_id == "spec.scan.coordinate_in_route" for issue in issues
    )


def test_missing_modred_atoms_requires_decline():
    issues = SI.check_spec(
        _workflow("gaussian.modred", {"modred": [[2, 3]]}),
        "Run Gaussian modredundant job for mol.xyz",
    )

    assert any(issue.rule_id == "spec.decline_contract" for issue in issues)


def test_job_project_is_runtime_owned():
    spec = _workflow("orca.modred", {"modred": [[1, 2]]})
    spec["jobs"][0]["project"] = "model-chosen-project"

    issues = SI.check_spec(spec, "Freeze bond 1-2 in mol.xyz")

    assert any(
        issue.rule_id == "spec.project.runtime_owned" for issue in issues
    )


def test_fragment_indices_with_underscore_counts_as_structural_slot():
    issues = SI.check_spec(
        _workflow("gaussian.dias", {"fragment_indices": [[1], [2]]}),
        "Need a Gaussian DIAS job for mol.xyz; fragment_indices [[1], [2]].",
    )

    assert not [
        issue for issue in issues if issue.rule_id == "spec.decline_contract"
    ]


def test_chinese_atom_enumeration_satisfies_decline_contract():
    issues = SI.check_spec(
        _workflow("orca.modred", {"modred": [[1, 2, 3]]}),
        "Pd、ipso-C 和 Br 分别是原子 1、2、3。保持当前角度。",
    )

    assert not [
        issue for issue in issues if issue.rule_id == "spec.decline_contract"
    ]


def test_db_record_selector_passes_for_db_file():
    issues = SI.check_spec(
        _workflow(
            "gaussian.sp",
            file="results.db",
            record_index=1,
            structure_index="2",
        ),
        "Run a Gaussian single point from structure 2 of record 1 in results.db",
    )

    assert not [
        issue for issue in issues if issue.rule_id.startswith("spec.db.")
    ]


def test_db_selector_cardinality_rejects_missing_selector_for_db_file():
    issues = SI.check_spec(
        _workflow("orca.sp", file="results.db"),
        "Run an ORCA single point from results.db",
    )

    assert any(
        issue.rule_id == "spec.db.selector_cardinality" for issue in issues
    )


def test_db_molecule_id_rejects_for_job_submission():
    issues = SI.check_spec(
        _workflow("orca.sp", file="results.db", molecule_id="mol-abc"),
        "Run an ORCA single point from molecule mol-abc in results.db",
    )

    assert any(issue.rule_id == "spec.db.molecule_id_job" for issue in issues)


def test_db_structure_index_without_record_rejects():
    issues = SI.check_spec(
        _workflow("gaussian.sp", file="results.db", structure_index="2"),
        "Run Gaussian on structure 2 from results.db",
    )

    assert any(
        issue.rule_id == "spec.db.structure_index_requires_record"
        for issue in issues
    )


def test_db_selector_without_db_file_rejects():
    issues = SI.check_spec(
        _workflow("gaussian.sp", record_index=1),
        "Run Gaussian single point for mol.xyz",
    )

    assert any(
        issue.rule_id == "spec.db.selector_without_db" for issue in issues
    )
