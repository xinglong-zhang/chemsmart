"""Every archived PySCF result reaches its recorded verdict, through the host.

Round 2 shipped a ``td`` consumer the host evaluator refused with
``geometry_invalid`` while the result validator, called directly, passed
it: the evaluator handed the expected input geometry for two of its
three fixed-geometry job types and the third was added to the set later.
Every test that touched the validator was green, because none of them
went through the organ the goal driver actually calls.

So this drives that organ, over every archived fixture rather than the
ones a test author chose, and compares its verdict with the verdict the
run itself recorded. The parametrisation is the directory listing: a
fixture generated for a new job type joins by existing.

The input artifact is resolved through the host's own canonical
geometry identity, never by the path in ``.input.json`` -- those are the
campaign host's paths and do not exist here -- and a fixture launched
from another result's archive is handed that archive rather than its
own.

The second test is the general form of the loss: for every job type held
to the geometry it was handed, evaluating the same result *without* that
geometry must refuse. Round 2's evaluator passed the expected geometry
for two of its three fixed-geometry job types and every td consumer
under the goal driver died on the third.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from chemsmart.agent._contracts import TrustedArtifactRefV1, file_sha256
from chemsmart.agent.tool_runtime import (
    CommandCompiledToolHostV1,
    _output_artifact_kind,
    _pyscf_input_geometry_sha256,
)
from chemsmart.jobs.pyscf.validation import FIXED_GEOMETRY_JOBTYPES

FIXTURES = (
    Path(__file__).resolve().parents[1] / "data" / "PySCFTests" / "outputs"
)
#: Directories that hold inputs, project YAML or another program's
#: differential rather than one archived PySCF run.
_NOT_A_RUN = frozenset({"inputs", "projects", "orca_differential"})

#: Fixtures whose archived receipt records the runner's verdict and the
#: host then adds its own. The receipt is the program's word about its
#: own artifact; the evaluation applies the program-neutral rules above
#: it, and where those two differ the difference is named here with the
#: physics behind it rather than smoothed away.
_HOST_RULE_ABOVE_THE_RECEIPT = {
    # An exactly planar ammonia relaxes onto its D3h inversion saddle
    # and its Hessian carries one imaginary mode at -830 cm-1. PySCF's
    # own validator has nothing to say about which stationary point a
    # Hessian describes; the host's order rule does, and it types the
    # result failed_wrong_stationary_point with the anomaly recorded.
    "nh3_planar_hess": {"result.stationary_point_order"},
}

pytestmark = pytest.mark.capability("program_jobtype:pyscf:cpu:*")


def _cases() -> tuple[str, ...]:
    return tuple(
        sorted(
            path.name
            for path in FIXTURES.iterdir()
            if path.is_dir() and path.name not in _NOT_A_RUN
        )
    )


def _artifact(path: Path, kind: str, prefix: str) -> TrustedArtifactRefV1:
    return TrustedArtifactRefV1(
        artifact_id=f"{prefix}.{path.name}",
        kind=kind,
        sha256=file_sha256(path),
        size_bytes=path.stat().st_size,
        path=str(path),
        cli_value=str(path),
    )


def _outputs(case: str, label: str) -> tuple[TrustedArtifactRefV1, ...]:
    found = []
    for path in sorted((FIXTURES / case).glob(f"{label}*")):
        if path.suffix == ".json" and "reference" in path.name:
            continue
        kind = (
            "pyscf_hdf5"
            if path.suffix == ".h5"
            else _output_artifact_kind("pyscf", path)
        )
        found.append(_artifact(path, kind, "result"))
    return tuple(found)


def _input_for(spec: dict) -> TrustedArtifactRefV1 | None:
    """The geometry this run was handed, found by the identity it recorded.

    The writer stamps ``input_geometry_sha256`` over the symbols, the
    positions, the unit and the bound state, which is the same identity
    the host computes for an approved artifact, so the archive can be
    joined to the file it came from without trusting a stale path.
    """

    candidates = [
        (path, "geometry_xyz")
        for path in sorted((FIXTURES / "inputs").glob("*.xyz"))
    ]
    candidates += [
        (path, "pyscf_hdf5")
        for path in sorted(FIXTURES.glob("*/*.h5"))
        if path.parent.name not in _NOT_A_RUN
    ]
    # A run handed another result's archive recorded that file's digest,
    # and only those exact bytes satisfy the run receipt: an equivalent
    # geometry in another container is the same molecule and the wrong
    # artifact.
    recorded = str(spec.get("input_artifact_sha256") or "")
    if recorded:
        for path, kind in candidates:
            if file_sha256(path) == recorded:
                return _artifact(
                    path,
                    str(spec.get("input_artifact_kind") or kind),
                    "geometry",
                )
        return None
    wanted = str(spec.get("input_geometry_sha256") or "")
    if not wanted:
        return None
    charge = int(spec["charge"])
    multiplicity = int(spec["multiplicity"])
    for path, kind in candidates:
        artifact = _artifact(path, kind, "geometry")
        identity = _pyscf_input_geometry_sha256(
            artifact, charge=charge, multiplicity=multiplicity
        )
        if identity and identity == wanted:
            return artifact
    return None


#: The settings an approved node carries into the evaluator. The
#: archive records what was applied, so the evaluation runs against the
#: same scientific request the run was launched with rather than against
#: a bare job type, which derives the wrong stage list for every
#: correlated or excited-root result.
_EVALUATED_SETTINGS = (
    "ab_initio",
    "excited_state_root",
    "frozen_core",
    "response_method",
    "state_manifold",
    "nstates",
)


def _settings(spec: dict) -> dict:
    return {
        "jobtype": spec["jobtype"],
        **{
            field: spec[field]
            for field in _EVALUATED_SETTINGS
            if spec.get(field) is not None
        },
    }


def _spec(case: str, label: str) -> dict:
    from chemsmart.io.pyscf.output import read_pyscf_h5

    spec, *_rest = read_pyscf_h5(str(FIXTURES / case / f"{label}.h5"))
    return dict(spec)


@pytest.mark.parametrize("case", _cases())
def test_the_host_evaluator_agrees_with_the_recorded_run(case):
    directory = FIXTURES / case
    archives = sorted(directory.glob("*.h5"))
    assert len(archives) == 1, f"{case} holds {len(archives)} archives"
    label = archives[0].stem
    receipt = json.loads((directory / f"{label}.receipt.json").read_text())
    spec = _spec(case, label)
    expected_input = _input_for(spec)
    assert expected_input is not None, (
        f"{case}: no archived geometry matches the identity the run "
        "recorded, so this fixture cannot be evaluated as it ran"
    )

    evaluation = CommandCompiledToolHostV1._evaluate_execution_outputs(
        program="pyscf",
        jobtype=str(spec["jobtype"]),
        charge=int(spec["charge"]),
        multiplicity=int(spec["multiplicity"]),
        expected_settings=_settings(spec),
        expected_input_artifact=expected_input,
        output_artifacts=_outputs(case, label),
        exit_status=int(receipt.get("child_returncode", 0)),
    )

    recorded = str(receipt["state"])
    expected_above = _HOST_RULE_ABOVE_THE_RECEIPT.get(case, set())
    if recorded == "validated" and not evaluation.validated:
        assert set(evaluation.findings) == expected_above, (
            f"{case}: the run recorded validated and the host evaluator "
            f"refuses with {sorted(evaluation.findings)}; a host rule "
            "above the receipt is named in this module, anything else "
            "is a disagreement between two organs"
        )
    else:
        assert expected_above == set(), (
            f"{case}: a host rule above the receipt is recorded here and "
            "no longer fires"
        )
        assert evaluation.validated is (recorded == "validated"), (
            f"{case}: the run recorded {recorded!r} and the host "
            f"evaluator says validated={evaluation.validated} with "
            f"findings {sorted(evaluation.findings)}"
        )
    recorded_findings = {
        item if isinstance(item, str) else str(item.get("rule_id", item))
        for item in (receipt.get("findings") or ())
    }
    missing = recorded_findings - set(evaluation.findings)
    assert (
        not missing
    ), f"{case}: findings the evaluator lost: {sorted(missing)}"


@pytest.mark.parametrize("case", _cases())
def test_a_fixed_geometry_result_is_refused_without_the_geometry(case):
    """The general form of round 2's F1 loss.

    A job type held to the geometry it was handed is only held if the
    caller hands it over. Evaluating the same archived result with no
    expected geometry must refuse for exactly those job types, so a job
    type added to the validator's set and not to the caller's cannot
    pass unnoticed.
    """

    directory = FIXTURES / case
    label = sorted(directory.glob("*.h5"))[0].stem
    receipt = json.loads((directory / f"{label}.receipt.json").read_text())
    spec = _spec(case, label)
    if str(receipt["state"]) != "validated":
        pytest.skip("only a validated result can lose its verdict")

    evaluation = CommandCompiledToolHostV1._evaluate_execution_outputs(
        program="pyscf",
        jobtype=str(spec["jobtype"]),
        charge=int(spec["charge"]),
        multiplicity=int(spec["multiplicity"]),
        expected_settings=_settings(spec),
        expected_input_artifact=None,
        output_artifacts=_outputs(case, label),
        exit_status=0,
    )
    assert evaluation.validated is False, (
        f"{case}: the evaluator validated a result without the geometry "
        "the run was handed"
    )
    # A job type held to that geometry says so by name, rather than
    # failing for some other reason that would survive the set changing.
    if str(spec["jobtype"]) in FIXED_GEOMETRY_JOBTYPES:
        assert any(
            "geometry" in finding for finding in evaluation.findings
        ), f"{case}: {sorted(evaluation.findings)}"


def test_the_parametrisation_is_the_directory_listing():
    """A fixture joins this test by existing, not by being remembered."""

    assert len(_cases()) >= 25
    assert "inputs" not in _cases()
