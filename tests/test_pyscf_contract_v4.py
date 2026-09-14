"""Result contract v4: what the PySCF driver now says about its own result.

Each fact here was missing on a real artifact and cost something a
scientist needed: the gradient at a Hessian geometry (a projected spectrum
is all-real at a non-stationary point, so 0 imaginary modes proved
nothing), per-atom spin populations (one call away, discarded), the mass
table behind the frequencies, the tolerance behind the point group, the
optimiser's own convergence standard, and the fact that the final SCF of
an optimisation continues from the optimiser's density.  Two Hessians
PySCF 2.14 cannot compute are refused at preflight instead of dying in
the engine after the call was spent.
"""

import re
from types import SimpleNamespace

import h5py
import numpy as np
import pytest

from chemsmart.io.pyscf.output import PySCFOutput
from chemsmart.jobs.pyscf.settings import PySCFJobSettings
from chemsmart.jobs.pyscf.validation import (
    RULE_HESSIAN_NLC_OPEN_SHELL,
    RULE_HESSIAN_REFERENCE,
    _result_contract_validation,
    preflight,
)
from chemsmart.jobs.pyscf.writer import (
    APPLIED_SPEC_FIELDS,
    APPLIED_SPEC_FIELDS_V4,
    APPLIED_SPEC_FIELDS_V5,
    LEGACY_APPLIED_SPEC_FIELDS,
    PREVIOUS_RESULT_CONTRACT_VERSIONS,
    RESULT_CONTRACT_VERSION,
    RESULT_UNITS,
    SUPPORTED_RESULT_CONTRACT_VERSIONS,
    PySCFScriptWriter,
    applied_pyscf_spec_fields,
    write_pyscf_h5,
)


def _molecule(symbols, positions, charge, multiplicity):
    return SimpleNamespace(
        chemical_symbols=list(symbols),
        positions=[list(row) for row in positions],
        charge=charge,
        multiplicity=multiplicity,
        info={},
    )


def _water(multiplicity=1):
    return _molecule(
        ["O", "H", "H"],
        [[0.0, 0.0, 0.0], [0.7586, 0.0, 0.5043], [-0.7586, 0.0, 0.5043]],
        0,
        multiplicity,
    )


def _hess_settings(**overrides):
    values = dict(
        jobtype="hess",
        freq=True,
        basis="def2-svp",
        charge=0,
        multiplicity=1,
        engine="cpu",
    )
    values.update(overrides)
    return PySCFJobSettings(**values)


# ----------------------------------------------------------------------
# the contract version is a supported set, not one string
# ----------------------------------------------------------------------


@pytest.mark.capability("program_jobtype:pyscf:cpu:hess")
def test_a_previous_supported_contract_is_evidence_not_a_downgrade():
    # Every version this ChemSmart reads, with the current one last and
    # named nowhere but the module that declares it: what is pinned is
    # the relation between the versions, not which round we are in.
    assert SUPPORTED_RESULT_CONTRACT_VERSIONS == (
        PREVIOUS_RESULT_CONTRACT_VERSIONS + (RESULT_CONTRACT_VERSION,)
    )
    assert RESULT_CONTRACT_VERSION not in PREVIOUS_RESULT_CONTRACT_VERSIONS

    # The applied-settings digest of an archived artifact is
    # reconstructed from the vocabulary of *its* contract, so each
    # version's tuple is frozen once it has artifacts on disk and the
    # current one only ever extends. Extending a tuple in place would
    # mark every archived artifact of that version tampered.
    frozen = {
        "chemsmart.pyscf-result-contract.v3": APPLIED_SPEC_FIELDS_V4,
        "chemsmart.pyscf-result-contract.v4": APPLIED_SPEC_FIELDS_V4,
        "chemsmart.pyscf-result-contract.v5": APPLIED_SPEC_FIELDS_V5,
    }
    for version in PREVIOUS_RESULT_CONTRACT_VERSIONS:
        vocabulary = applied_pyscf_spec_fields(
            {"result_contract_version": version}
        )
        assert vocabulary == frozen[version], version
        # A previous vocabulary is a prefix of the current one: the
        # current contract adds fields and moves none.
        assert APPLIED_SPEC_FIELDS[: len(vocabulary)] == vocabulary, version
    assert (
        applied_pyscf_spec_fields(
            {"result_contract_version": RESULT_CONTRACT_VERSION}
        )
        == APPLIED_SPEC_FIELDS
    )
    # What v6 adds: the electronic surface a result's geometry and total
    # energy belong to, and how a Hessian's second derivative was taken.
    assert APPLIED_SPEC_FIELDS[len(APPLIED_SPEC_FIELDS_V5) :] == (
        "surface",
        "hessian_derivative",
        "fd_step_angstrom",
    )

    complete_spec = {
        "reference_family": "rks",
        "charge": 0,
        "multiplicity": 1,
        "spin": 0,
        "num_electrons": 10,
        "nelec": [5, 5],
    }
    complete_status = {
        "stages": {"scf": {"converged": True}},
        "engine_complete": True,
        "normal_termination": True,
        "failure": None,
        "properties": {},
    }
    for version in SUPPORTED_RESULT_CONTRACT_VERSIONS:
        observation, current, findings, _advisories = (
            _result_contract_validation(
                {"result_contract_version": version, **complete_spec},
                complete_status,
            )
        )
        assert findings == []
        assert observation["new_execution_admissible"] is True
        assert observation["state"] == (
            "current" if version == RESULT_CONTRACT_VERSION else "supported"
        )
        assert current is True, "every supported contract is a strict record"

    assert (
        applied_pyscf_spec_fields({}) == LEGACY_APPLIED_SPEC_FIELDS
    ), "a contract-less artifact keeps the legacy digest vocabulary"
    observation, _current, findings, _advisories = _result_contract_validation(
        {"result_contract_version": "chemsmart.pyscf-result-contract.v9"},
        {},
    )
    assert observation["state"] == "unsupported"
    assert observation["new_execution_admissible"] is False
    assert [item.field for item in findings] == [
        "spec.result_contract_version"
    ], "an unknown contract is refused by name, once"


# ----------------------------------------------------------------------
# the generated driver
# ----------------------------------------------------------------------


def _driver_source():
    return PySCFScriptWriter.render(
        {"schema_version": "2.0", "label": "contract-v4"}
    )


def _function_text(source, name):
    """The text of one driver function, or '' when it has none."""

    marker = "def %s(" % name
    if marker not in source:
        return ""
    start = source.index(marker)
    remainder = source[start + len(marker) :]
    offsets = [
        remainder.index(nxt)
        for nxt in ("\ndef ", "\nclass ")
        if nxt in remainder
    ]
    end = start + len(marker) + (min(offsets) if offsets else len(remainder))
    return source[start:end]


def _stage_function(source, stage):
    """``_run_<stage>`` and the driver helpers it calls.

    A stage's work is where the stage put it, and a stage with two ways
    of doing its job puts some of it in a helper: the Hessian stage
    takes the gradient of whichever surface it differentiated. What
    these tests pin is the work, so the reader follows the calls rather
    than the layout.
    """

    text = _function_text(source, "_run_%s" % stage)
    if not text:
        return ""
    called = {
        name
        for name in re.findall(r"\b(_[A-Za-z0-9_]+)\s*\(", text)
        if name != "_run_%s" % stage
    }
    for name in sorted(called):
        text += _function_text(source, name)
    return text


def _stage_branch(source, stage):
    """The text of one stage's work, wherever the driver keeps it.

    A stage that grew a second way of doing its job is extracted into
    ``_run_<stage>`` beside the others, so the branch is one call and the
    work is in the function. What these tests pin is the work -- which
    stage launches a gradient, which one states its mass convention --
    and following the delegation keeps them pinning that rather than the
    layout it happened to have.
    """

    order = ("scf", "opt", "td", "corr", "hess")
    start = source.index(
        ('if stage == "%s":' if stage == "scf" else 'elif stage == "%s":')
        % stage
    )
    following = order[order.index(stage) + 1 :]
    # The last branch ends where the stage loop refuses an unknown stage.
    end = min(
        [
            source.index('elif stage == "%s":' % name)
            for name in following
            if ('elif stage == "%s":' % name) in source
        ]
        + [source.index('raise ValueError("Unknown stage')]
    )
    return source[start:end] + _stage_function(source, stage)


@pytest.mark.capability("program_jobtype:pyscf:cpu:hess")
def test_the_gradient_is_computed_inside_the_hess_stage_only():
    """A single point still pays for no undeclared gradient.

    The former test pinned the text ``nuc_grad_method`` out of the whole
    driver; the invariant it protected is that the ``scf`` stage launches
    no gradient of its own, which is what is pinned now.  Contract v5's
    ``opt`` branch builds the gradient *scanner* an excited-root or a
    correlated optimisation walks on -- that is the stage's declared work,
    and every such call is a scanner construction, never a free gradient
    -- and the ``td`` and ``corr`` stages compute none.
    """

    source = _driver_source()
    compile(source, "<pyscf-driver>", "exec")
    assert "nuc_grad_method" in _stage_branch(source, "hess")
    for stage in ("scf", "td", "corr"):
        assert "nuc_grad_method" not in _stage_branch(source, stage)
    opt = _stage_branch(source, "opt")
    assert opt.count("nuc_grad_method") == opt.count(
        "nuc_grad_method().as_scanner("
    )
    hess = _stage_branch(source, "hess")
    assert 'results["forces"] = -gradient' in hess
    assert '"max_abs_gradient_eh_per_bohr"' in hess
    assert '"gradient_computed"' in hess
    assert RESULT_UNITS["forces"] == "Eh/Bohr"


@pytest.mark.capability("program_jobtype:pyscf:cpu:hess")
def test_the_hess_stage_states_its_mass_convention():
    hess = _stage_branch(_driver_source(), "hess")
    assert '"mass_convention": "isotope_averaged"' in hess
    assert "atom_mass_list(isotope_avg=True)" in hess


@pytest.mark.capability("program_jobtype:pyscf:cpu:opt")
def test_the_opt_stage_states_its_standard_and_its_continuation():
    source = _driver_source()
    opt = _stage_branch(source, "opt")
    assert '"final_scf_from_optimizer_density"' in opt
    assert "_optimizer_criteria(" in opt
    assert "geometric.params import OptParams" in source


@pytest.mark.capability("program_jobtype:pyscf:cpu:sp")
def test_spin_populations_are_written_for_open_shells_only():
    source = _driver_source()
    assert "mulliken_spin_pop" in source
    assert 'int(config.get("spin") or 0) == 0' in source
    assert '"not_applicable"' in source
    assert RESULT_UNITS["mulliken_spin_populations"] == "electron"


@pytest.mark.capability("program_jobtype:pyscf:cpu:sp")
def test_the_point_group_carries_the_tolerance_that_decided_it():
    source = _driver_source()
    assert '"tolerance_bohr": float(symm_geom.TOLERANCE)' in source
    assert '"source": "pyscf.symm.detect_symm"' in source


# ----------------------------------------------------------------------
# the host-side writer and reader agree on the new datasets
# ----------------------------------------------------------------------


def _write(tmp_path, *, spin_populations):
    path = tmp_path / "v4.h5"
    spec = {
        "program": "pyscf",
        "jobtype": "sp",
        "symbols": ["O", "H"],
        "positions": [[0.0, 0.0, 0.0], [0.96, 0.0, 0.0]],
        "unit": "Angstrom",
        "charge": 0,
        "multiplicity": 2,
        "spin": 1,
        "result_contract_version": RESULT_CONTRACT_VERSION,
    }
    status = {
        "normal_termination": True,
        "engine_complete": True,
        "failure": None,
        "stages": {"scf": {"converged": True}},
        "properties": {
            "mulliken_spin_populations": (
                {"status": "ok"}
                if spin_populations is not None
                else {
                    "status": "not_applicable",
                    "reason": "closed-shell reference carries no spin",
                }
            )
        },
    }
    results = {
        "energies": np.array([-75.0]),
        "positions": np.array(spec["positions"], dtype=float),
        "mo_energy": np.array([-1.0, 0.1]),
        "mo_occ": np.array([1.0, 0.0]),
    }
    if spin_populations is not None:
        results["mulliken_spin_populations"] = np.asarray(spin_populations)
    write_pyscf_h5(
        path,
        spec=spec,
        provenance={"engine": "cpu"},
        status=status,
        results=results,
    )
    return path


@pytest.mark.capability("program_jobtype:pyscf:cpu:sp")
def test_the_supplied_geometry_carries_its_unit_on_the_dataset(tmp_path):
    path = _write(tmp_path, spin_populations=[1.04, -0.04])
    with h5py.File(path, "r") as handle:
        assert handle["spec/positions"].attrs["unit"] == "Angstrom"
        assert (
            handle["results/mulliken_spin_populations"].attrs["unit"]
            == "electron"
        )


@pytest.mark.capability("program_jobtype:pyscf:cpu:sp")
def test_spin_populations_read_back_and_a_closed_shell_is_not_a_failure(
    tmp_path,
):
    open_shell = PySCFOutput(_write(tmp_path, spin_populations=[1.04, -0.04]))
    assert open_shell.mulliken_atomic_spin_populations == [1.04, -0.04]
    assert open_shell.property_failures == {}

    closed = PySCFOutput(
        _write(tmp_path / "closed", spin_populations=None)
        if (tmp_path / "closed").mkdir() is None
        else None
    )
    assert closed.mulliken_atomic_spin_populations is None
    assert (
        closed.property_failures == {}
    ), "not_applicable is neither an attempted failure nor an omission"


# ----------------------------------------------------------------------
# two Hessians PySCF cannot compute are refused before the engine
# ----------------------------------------------------------------------


@pytest.mark.capability("program_jobtype:pyscf:cpu:hess")
def test_a_one_electron_hessian_is_refused_at_preflight():
    hydrogen = _molecule(["H"], [[0.0, 0.0, 0.0]], 0, 2)
    findings = preflight(
        _hess_settings(ab_initio="hf", multiplicity=2), hydrogen, {}
    )
    rules = {item.rule_id for item in findings}
    assert RULE_HESSIAN_REFERENCE in rules
    refusal = next(
        item for item in findings if item.rule_id == RULE_HESSIAN_REFERENCE
    )
    assert refusal.observed["reference_family"] == "rohf"
    assert "UKS" in refusal.observed["reason"], "the refusal names the route"

    # The same atom under a DFT functional runs UKS and is admitted.
    findings = preflight(
        _hess_settings(functional="pbe", multiplicity=2), hydrogen, {}
    )
    assert RULE_HESSIAN_REFERENCE not in {item.rule_id for item in findings}


@pytest.mark.capability("program_jobtype:pyscf:cpu:hess")
def test_an_open_shell_nlc_hessian_is_refused_when_the_environment_says_so():
    hydroxyl = _molecule(["O", "H"], [[0.0, 0.0, 0.0], [0.96, 0.0, 0.0]], 0, 2)
    evidence = {"functional_metadata": {"wb97m-v": {"nlc": True}}}
    findings = preflight(
        _hess_settings(functional="wb97m-v", multiplicity=2),
        hydroxyl,
        evidence,
    )
    assert RULE_HESSIAN_NLC_OPEN_SHELL in {item.rule_id for item in findings}

    # A closed shell under the same functional has a Hessian.
    findings = preflight(
        _hess_settings(functional="wb97m-v"), _water(), evidence
    )
    assert RULE_HESSIAN_NLC_OPEN_SHELL not in {
        item.rule_id for item in findings
    }
    # Without the environment's word (a preview) nothing is asserted.
    findings = preflight(
        _hess_settings(functional="wb97m-v", multiplicity=2), hydroxyl, {}
    )
    assert RULE_HESSIAN_NLC_OPEN_SHELL not in {
        item.rule_id for item in findings
    }
    # A closed-shell water Hessian raises neither rule.
    findings = preflight(_hess_settings(functional="b3lyp"), _water(), {})
    assert not (
        {RULE_HESSIAN_REFERENCE, RULE_HESSIAN_NLC_OPEN_SHELL}
        & {item.rule_id for item in findings}
    )
