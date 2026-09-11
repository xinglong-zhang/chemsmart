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

GEOMETRY_SELECTORS = ("positions", "reached_positions")


def _fixture_output():
    hits = sorted(glob.glob(f"{FIXTURE}/*_optts_optts.out"))
    if not hits or not pathlib.Path(hits[0]).is_file():
        pytest.skip("the archived multi-structure ORCA result is absent")
    return reader_for("orca"), hits[0]


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


@pytest.mark.capability("selector:orca:ts:reached_positions")
def test_a_reached_selector_returns_the_last_structure():
    """The relation that needs no oracle: reached means the last one."""

    reader, path = _fixture_output()
    handle = reader.open_output(path)
    reached = _coordinates(reader, handle, "reached_positions")

    molecule = handle.molecule
    if isinstance(molecule, (list, tuple)):
        molecule = molecule[-1]
    last = np.asarray(molecule.positions, dtype=float)
    assert reached.shape == last.shape
    assert np.allclose(reached, last, atol=1e-8), (
        "a selector declared 'as_reached' does not return the last "
        "structure the result printed"
    )


@pytest.mark.capability("selector:orca:ts:positions")
def test_selectors_declaring_different_states_do_not_agree():
    """If two states return the same bytes, a declaration is false.

    On this fixture the run moved 1.23 A, so ``positions``
    (``thermochemistry_reference``) and ``reached_positions``
    (``as_reached``) must differ. Before D2 both questions were answered
    by one accessor and this relation could not be stated, let alone
    checked.
    """

    reader, path = _fixture_output()
    handle = reader.open_output(path)
    seen: dict[str, np.ndarray] = {}
    for selector in GEOMETRY_SELECTORS:
        if selector not in reader.accessors:
            continue
        seen[reader.structural_state(selector)] = _coordinates(
            reader, handle, selector
        )
    assert len(seen) >= 2, (
        "the fixture exercises fewer than two structural states, so the "
        f"relation is untested: {sorted(seen)}"
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
