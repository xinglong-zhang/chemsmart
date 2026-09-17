"""A declared structural state must be the state the accessor reads.

D2 lets a selector say which molecular state its value belongs to, so a
consumer can ask for the role it needs -- ``build_reached_geometry``
wants ``as_reached`` -- and be refused rather than silently served
another. The objection to that design is its own recursion: **nothing
stops someone declaring ``as_reached`` on an accessor that reads the
thermochemistry block.** A declaration that is merely asserted buys
nothing, which is the defect class this whole round exists to end.

So the declaration needs an oracle, and one exists that requires no
chemistry and no knowledge of which accessor is *right*: a metamorphic
relation over a result that carries more than one structure.

For a completed optimisation with a multi-frame trajectory:

- a selector declared ``as_reached`` returns the last structure;
- a selector declared ``as_supplied`` or ``thermochemistry_reference``
  does not return that structure when the run moved;
- two selectors declaring **different** states do not return identical
  coordinates, because then one of the two declarations is false.

The fixture is the live unconverged saddle search from ``po3-r18``:
100 printed structures, the supplied and reached geometries
**1.231588 A** apart, and their energies **182.2 kcal/mol** apart under
one artifact and one receipt. That is the exact state in which the old
code answered "the structure orca reached" with the seed.
"""

import glob
import pathlib

import numpy as np
import pytest

from chemsmart.analysis.result_readers import (
    RESULT_READERS,
    STRUCTURAL_STATES,
    reader_for,
)

#: A real unconverged ORCA OptTS: Hessian at step 0, 100 printed
#: structures, normal termination, not converged.
FIXTURE = (
    "/home/chemsmart/agent-campaigns/ax41-refine-100/novel-round-7"
    "/workspaces/po3-r18/nodes/ts-ester-c4"
)

GEOMETRY_SELECTORS = ("positions", "reached_positions", "supplied_positions")

#: Archived PySCF optimisations that moved (tests/data): one converged
#: from a distorted water (0.083 A), one stopped after a single step
#: (0.066 A). The artifact carries the supplied and the final structure.
PYSCF_FIXTURES = (
    "tests/data/PySCFTests/outputs/water_opt/water_opt_gas_phase.h5",
    "tests/data/PySCFTests/outputs/water_opt_maxsteps1/"
    "water_opt_maxsteps1_gas_phase.h5",
)


def _fixture_output():
    hits = sorted(glob.glob(f"{FIXTURE}/*_optts_optts.out"))
    if not hits or not pathlib.Path(hits[0]).is_file():
        pytest.skip("the archived multi-structure ORCA result is absent")
    return reader_for("orca"), hits[0]


def _moved_outputs():
    """Every archived result whose run moved, with its reader."""

    found = []
    hits = sorted(glob.glob(f"{FIXTURE}/*_optts_optts.out"))
    if hits and pathlib.Path(hits[0]).is_file():
        found.append(("orca", reader_for("orca"), hits[0]))
    root = pathlib.Path(__file__).resolve().parents[2]
    for relative in PYSCF_FIXTURES:
        path = root / relative
        if path.is_file():
            found.append(("pyscf", reader_for("pyscf"), str(path)))
    if not found:
        pytest.skip("no archived multi-structure result is present")
    return found


def _coordinates(reader, handle, selector):
    value = reader.accessors[selector](handle)
    return np.asarray(value, dtype=float)


@pytest.mark.capability("selector:orca:ts:positions")
def test_every_selector_declares_a_structural_state():
    """``stateless`` is a choice; silence is not.

    Exhaustive rather than heuristic: ORCA exposes 72 selectors and a
    keyword scan over-includes -- ``symbols`` is structure-invariant,
    ``solvation_model`` is a string -- so guessing which are
    state-bearing is the wrong instrument. Every selector answers, and
    adding one forces the decision at declaration time.
    """

    undeclared = {}
    for program, reader in sorted(RESULT_READERS.items()):
        if not reader.selector_structural_states:
            continue
        missing = sorted(
            name
            for name in reader.accessors
            if reader.structural_state(name) == "stateless"
            and name in {"positions", "energy", "gibbs_free_energy"}
        )
        if missing:
            undeclared[program] = missing
    assert not undeclared, (
        "these load-bearing selectors fall back to 'stateless', which "
        f"for a structure-dependent quantity is silence: {undeclared}"
    )


@pytest.mark.capability("selector:orca:ts:positions")
def test_declared_states_are_from_the_vocabulary():
    for program, reader in sorted(RESULT_READERS.items()):
        for selector, state in reader.selector_structural_states:
            assert state in STRUCTURAL_STATES, (
                f"{program}:{selector} declares {state!r}, which is not "
                f"one of {STRUCTURAL_STATES}"
            )


@pytest.mark.capability("selector:pyscf:td:excitation_energies")
def test_electronic_provenance_is_declared_from_one_vocabulary_and_resolved():
    """A structural state identifies a geometry, not a density.

    The second axis says whose density or method a value belongs to, from
    one vocabulary, declared per reader and resolved per artifact: the
    symbolic ``computed_surface`` a reader declares for ``energy`` never
    reaches a consumer, because the resolver turns it into the reference,
    the followed root or the correlated method this result's own record
    supports.  Every reader that declares the axis must declare it for the
    excitation set and the mean-field set alike, or a session reads "the
    S1 dipole" off a ground-state number under a true host word.
    """

    from chemsmart.analysis.result_readers import ELECTRONIC_PROVENANCES

    declaring = {
        program: reader
        for program, reader in RESULT_READERS.items()
        if reader.selector_electronic_provenance
    }
    assert {"orca", "pyscf"} <= set(declaring)
    for program, reader in sorted(declaring.items()):
        for selector, word in reader.selector_electronic_provenance:
            assert word in ELECTRONIC_PROVENANCES, (program, selector, word)
        for selector in ("excitation_energies", "oscillator_strengths"):
            assert reader.electronic_provenance(selector) == "excited_root"
        for selector in ("dipole_moment", "scf_energy"):
            assert reader.electronic_provenance(selector) == "reference"
        assert reader.electronic_provenance("energy") == "computed_surface"
        assert reader.resolve_electronic_provenance is not None
    root = pathlib.Path(__file__).resolve().parents[2]
    reader = reader_for("pyscf")
    for relative, expected in (
        ("water_sp/water_sp_gas_phase.h5", "reference"),
        ("water_mp2_sp/water_mp2_sp_gas_phase.h5", "correlated"),
        (
            "formaldehyde_s1_opt/formaldehyde_s1_opt_gas_phase.h5",
            "excited_root",
        ),
        ("water_td_singlet/water_td_singlet_gas_phase.h5", "reference"),
    ):
        path = root / "tests" / "data" / "PySCFTests" / "outputs" / relative
        if not path.is_file():
            pytest.skip("an archived PySCF result is absent")
        handle = reader.open_output(str(path))
        resolved = reader.electronic_provenance_for_output(handle, "energy")
        assert resolved == expected, (relative, resolved)
        assert resolved != "computed_surface"


@pytest.mark.capability("selector:orca:ts:reached_positions")
@pytest.mark.capability("selector:pyscf:opt:reached_positions")
def test_a_reached_selector_returns_the_last_structure():
    """The relation that needs no oracle: reached means the last one."""

    for _program, reader, path in _moved_outputs():
        handle = reader.open_output(path)
        reached = _coordinates(reader, handle, "reached_positions")

        molecule = handle.molecule
        if isinstance(molecule, (list, tuple)):
            molecule = molecule[-1]
        last = np.asarray(molecule.positions, dtype=float)
        assert reached.shape == last.shape
        assert np.allclose(reached, last, atol=1e-8), (
            "a selector declared 'as_reached' does not return the last "
            f"structure the result printed ({path})"
        )


@pytest.mark.capability("selector:orca:ts:positions")
@pytest.mark.capability("selector:pyscf:opt:supplied_positions")
def test_selectors_declaring_different_states_do_not_agree():
    """If two states return the same bytes, a declaration is false.

    On the ORCA fixture the run moved 1.23 A, so ``positions``
    (``thermochemistry_reference``) and ``reached_positions``
    (``as_reached``) must differ; on the PySCF fixtures the run moved
    0.08 and 0.07 A, so ``supplied_positions`` (``as_supplied``) and
    ``reached_positions`` must differ. Before D2 both questions were
    answered by one accessor and this relation could not be stated, let
    alone checked.
    """

    for _program, reader, path in _moved_outputs():
        handle = reader.open_output(path)
        seen: dict[str, np.ndarray] = {}
        for selector in GEOMETRY_SELECTORS:
            if selector not in reader.accessors:
                continue
            if selector not in (
                reader.selectors_for_jobtype(str(handle.jobtype)) or ()
            ):
                continue
            seen[reader.structural_state(selector)] = _coordinates(
                reader, handle, selector
            )
        assert len(seen) >= 2, (
            "the fixture exercises fewer than two structural states, so the "
            f"relation is untested: {sorted(seen)} ({path})"
        )
        states = sorted(seen)
        for index, first in enumerate(states):
            for second in states[index + 1 :]:
                left, right = seen[first], seen[second]
                if left.shape != right.shape:
                    continue
                rmsd = float(np.sqrt(((left - right) ** 2).sum(axis=1).mean()))
                assert rmsd > 1e-6, (
                    f"{first!r} and {second!r} are declared as different "
                    f"structural states and return identical coordinates "
                    f"(RMSD {rmsd:.3e} A), so one declaration is false"
                )


@pytest.mark.capability("selector:pyscf:hess:supplied_positions")
def test_a_fixed_geometry_result_reaches_what_it_was_handed():
    """The fourth relation, which PySCF makes checkable: on a stage that
    moves no atom the supplied and the final structure are one, and the
    runner's own fixed-geometry invariant is read back through the
    selector plane."""

    root = pathlib.Path(__file__).resolve().parents[2]
    reader = reader_for("pyscf")
    for relative in (
        "tests/data/PySCFTests/outputs/water_hess/water_hess_gas_phase.h5",
        "tests/data/PySCFTests/outputs/water_sp/water_sp_gas_phase.h5",
    ):
        path = root / relative
        if not path.is_file():
            pytest.skip("archived PySCF fixed-geometry result is absent")
        handle = reader.open_output(str(path))
        supplied = _coordinates(reader, handle, "supplied_positions")
        final = _coordinates(reader, handle, "positions")
        assert np.allclose(supplied, final, atol=1e-8)
        assert reader.structural_state("supplied_positions") == "as_supplied"
        assert reader.structural_state("positions") == "as_reached"


@pytest.mark.capability("tool:bind_reached_geometry")
def test_the_recovery_route_asks_for_the_role_not_the_accessor():
    """``build_reached_geometry`` must select by declared state.

    Pinned at the source because the defect was invisible to every
    behavioural test: the tool returned coordinates, they were
    well-typed and digest-bound, and they were the input. The
    behavioural half now exists -- the witness bank drives the public
    tool over this same archived result and measures the bytes -- and
    this pin stays for the direction a witness cannot see: that the
    accessor is never reached by name again.
    """

    source = pathlib.Path("chemsmart/agent/execution.py").read_text(
        encoding="utf-8"
    )
    marker = 'selectors_in_state_for_output(output, "as_reached")'
    assert marker in source, (
        "build_reached_geometry no longer selects the geometry by its "
        "declared structural state"
    )
    head = source.index("def build_reached_geometry")
    body = source[head : head + 6000]
    assert 'accessors["positions"](' not in body, (
        "build_reached_geometry reads the 'positions' accessor directly "
        "again; for ORCA that is the thermochemistry geometry, which is "
        "the structure the run was handed"
    )


@pytest.mark.capability("selector:orca:ts:reached_positions")
def test_the_role_is_resolved_against_this_results_own_jobtype():
    """A role is served only by a selector *this jobtype* declares.

    The program-level state map and the per-jobtype declaration are two
    questions, and a host organ asking only the first gets a selector
    nobody audited for the job that ran.
    """

    reader, path = _fixture_output()
    handle = reader.open_output(path)
    assert handle.jobtype == "ts"
    served = reader.selectors_in_state_for_output(handle, "as_reached")
    assert "reached_positions" in served, (
        "the jobtype that motivated the repair does not serve the role: "
        f"{served}"
    )


@pytest.mark.capability("selector:orca:irc:irc_direction")
def test_a_jobtype_with_one_printed_structure_declares_no_state():
    """The charter's IRC restriction, made computable.

    ORCA writes the reaction path to an XYZ sidecar; the log's only
    printed structure is where the path started, so its printed energy
    differs from the true endpoint by the entire barrier. The charter
    says every state-dependent selector is therefore undeclared for the
    jobtype -- and until this test, nothing checked it. A selector
    declared here in any state but ``stateless`` is the IRC defect
    returning, and it would now arrive through the role resolver that
    the reached-geometry route asks.
    """

    reader = RESULT_READERS["orca"]
    declared = reader.selectors_for_jobtype("irc")
    assert declared, "orca no longer declares an irc jobtype"
    state_bearing = {
        name: reader.structural_state(name)
        for name in declared
        if reader.structural_state(name) != "stateless"
    }
    assert not state_bearing, (
        "orca irc declares selectors that belong to a molecular state, "
        "but its log prints only the structure the path started from: "
        f"{state_bearing}"
    )


# ----------------------------------------------------------------------
# the inspection reply carries every declared axis beside the level
# ----------------------------------------------------------------------


@pytest.mark.capability("tool:inspect_run")
def test_the_inspection_reply_names_the_level_beside_each_selectors_axes(
    tmp_path,
):
    """``inspect_run`` on one artifact answers three declared questions per
    selector -- whether the job type declares it, which structure it
    belongs to, whose density it is -- and, from the artifact's own
    record, the level the result computed at: the response and the root
    for an excited-surface optimisation, the frozen-core count for a
    correlated method.  A session names the level beside the number it
    delivers instead of inferring it from a project it may not hold."""

    from chemsmart.agent._contracts import TrustedArtifactRefV1, file_sha256
    from chemsmart.agent.runtime.event_store import RuntimeEventStore
    from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1

    root = pathlib.Path(__file__).resolve().parents[2]
    outputs = root / "tests" / "data" / "PySCFTests" / "outputs"
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    host = CommandCompiledToolHostV1(
        event_store=RuntimeEventStore(
            tmp_path / "events.jsonl", session_id="s1"
        ),
        task_spec_sha256s=("a" * 64,),
        approved_workspace=workspace,
    )
    for artifact_id, relative in (
        ("s1-opt", "formaldehyde_s1_opt/formaldehyde_s1_opt_gas_phase.h5"),
        ("ccsdt-sp", "water_ccsdt_sp/water_ccsdt_sp_gas_phase.h5"),
    ):
        path = outputs / relative
        host.artifacts[artifact_id] = TrustedArtifactRefV1(
            artifact_id=artifact_id,
            kind="pyscf_hdf5",
            sha256=file_sha256(path),
            size_bytes=path.stat().st_size,
            path=str(path),
            cli_value=str(path),
        )

    excited = host._inspect_run(
        "t1", {"program": "pyscf", "artifact_id": "s1-opt"}
    )
    assert excited["jobtype"] == "opt"
    assert excited["level"] == {
        "functional": "b3lyp",
        "basis": "def2-svp",
        "response_method": "tda",
        "state_manifold": "singlet",
        "nstates": 1,
        "excited_state_root": 1,
    }
    assert "excited_state_followed_root" in excited["requestable_selectors"]
    assert excited["structural_states"]["reached_positions"] == "as_reached"
    assert excited["electronic_provenance"]["energy"] == "excited_root"
    assert excited["electronic_provenance"]["dipole_moment"] == "reference"

    correlated = host._inspect_run(
        "t2", {"program": "pyscf", "artifact_id": "ccsdt-sp"}
    )
    assert correlated["level"] == {
        "ab_initio": "ccsd(t)",
        "basis": "def2-svp",
        "frozen_core": 1,
    }
    assert correlated["electronic_provenance"]["energy"] == "correlated"
    assert correlated["electronic_provenance"]["scf_energy"] == "reference"
    assert "triples_correction" in correlated["requestable_selectors"]


@pytest.mark.capability("tool:inspect_run")
@pytest.mark.capability("selector:xtb:opt:xtb_scc_atomic_charges")
def test_an_xtb_population_reaches_the_agent_with_its_own_scheme(tmp_path):
    """The inspection reply must not call an SCC population generic charge."""

    from chemsmart.agent._contracts import TrustedArtifactRefV1, file_sha256
    from chemsmart.agent.runtime.event_store import RuntimeEventStore
    from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1

    path = (
        pathlib.Path(__file__).resolve().parents[2]
        / "tests/data/XTBTests/outputs/co2_ohess/co2_ohess.out"
    )
    host = CommandCompiledToolHostV1(
        event_store=RuntimeEventStore(
            tmp_path / "events.jsonl", session_id="xtb-populations"
        ),
        task_spec_sha256s=("a" * 64,),
        approved_workspace=tmp_path / "workspace",
    )
    host.artifacts["co2-opt"] = TrustedArtifactRefV1(
        artifact_id="co2-opt",
        kind="xtb_output",
        sha256=file_sha256(path),
        size_bytes=path.stat().st_size,
        path=str(path),
        cli_value=str(path),
    )

    inspected = host._inspect_run(
        "t1", {"program": "xtb", "artifact_id": "co2-opt"}
    )
    assert "xtb_scc_atomic_charges" in inspected["requestable_selectors"]
    assert "mulliken_atomic_charges" not in inspected["requestable_selectors"]
    assert inspected["atom_resolved_metadata"] == {
        "xtb_scc_atomic_charges": {
            "semantic_quantity": "atomic_partial_charge",
            "population_scheme": "xTB self-consistent-charge population",
            "atom_order": "zero-based molecular atom order",
        }
    }
    assert inspected["level"] == {"method": "GFN2-xTB"}
