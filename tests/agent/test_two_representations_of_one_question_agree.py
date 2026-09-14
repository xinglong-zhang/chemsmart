"""Two host representations of one question must not disagree in silence.

Every oracle the previous round built asks *does the code do what the code
says?* -- declaration against implementation, name against structural
state, parameter against round trip, resolver against resolver. Each is a
**self-consistency** oracle, and a fully declared, correctly provenanced
quantity can still be false about the world. `connectivity` carries a
name, a unit, a dimension, a per-program structural state, fourteen
jobtype declarations, a ladder rung and tests; every one of those
declarations is true; and on a converged formaldehyde it reported that
neither hydrogen is bonded to the carbon (``sm1-formaldehyde``,
2026-09-11).

This is the missing class: a **referential** oracle. It needs no chemistry
table and no list of molecules. The host answers "which atoms are bonded"
in more than one place, and the RED condition is that two answers the host
presents as interchangeable disagree on some input.

The generator is derived, not enumerated:

1. take the element pairs the radii table serves;
2. take each registered policy's decision boundary for that pair;
3. probe just below, at, and just above each boundary, and between two
   policies' boundaries;
4. drive the real consumers and compare.

That third step is what finds this defect class, because a threshold
defect lives at a boundary and nowhere else. All three of the cases the
owner named -- H2 at 0.7410, H2CO's C-H at 1.1215, a hydrogen-bonded O-H
at 1.0300 -- lie *between* two of this host's own boundaries, so they are
generated rather than remembered.

A policy that is deliberately **not** claimed interchangeable (the
conformer grouper's own buffer, the rdkit wrapper's) is declared legacy
and excluded here: disagreement between two differently-named conventions
is a scientific observation, not a defect. Only silence is the defect.
"""

import itertools
import json

import numpy as np
import pytest

from chemsmart.io.molecules import get_covalent_radius
from chemsmart.io.molecules.structure import Molecule

#: Element pairs to generate over. Drawn from what the host's own readers
#: and test corpus actually contain, plus hydrogen against each -- not a
#: chemistry selection, an inventory one.
ELEMENTS = ("H", "C", "N", "O", "F", "P", "S", "Cl", "Si", "Ti")

#: How far off a boundary to probe. Small enough that nothing but the
#: boundary decision can change, large enough to survive float noise.
EPSILON = 5.0e-4


def _radii_sum(first: str, second: str) -> float:
    return get_covalent_radius(first) + get_covalent_radius(second)


def _interchangeable_perceivers():
    """The perceivers the host presents as answering one question.

    Both are agent-reachable and neither is documented as a distinct
    convention, so a session reading one and a sensor reading the other
    are entitled to the same answer.
    """

    from chemsmart.agent.execution import _molecule_graph
    from chemsmart.analysis.result_readers import _connectivity_matrix

    def by_selector(molecule):
        matrix = _connectivity_matrix(molecule)
        return {
            (i, j)
            for i in range(len(matrix))
            for j in range(i + 1, len(matrix))
            if matrix[i][j]
        }

    def by_execution_plane(molecule):
        graph = _molecule_graph(molecule)
        # ``_molecule_graph`` numbers atoms from one.
        return {(min(a, b) - 1, max(a, b) - 1) for a, b in graph.edges}

    return {
        "analysis.connectivity_selector": by_selector,
        "execution.molecule_graph": by_execution_plane,
    }


def _boundaries(first: str, second: str) -> tuple[float, ...]:
    """Every decision boundary this host has for one element pair.

    Read from the perceivers themselves rather than from a list, so a
    policy change moves the probes with it.
    """

    from chemsmart.io.molecules import get_bond_cutoff

    total = _radii_sum(first, second)
    seen = {
        # the additive conventions the tree has used
        get_bond_cutoff(first, second, 0.3),
        get_bond_cutoff(first, second, 0.05),
        get_bond_cutoff(first, second, 0.12),
        # a multiplicative one, and the two-regime form
        1.30 * total,
        min(1.30 * total, total + 0.45),
    }
    return tuple(sorted(value for value in seen if value > 0.2))


def _probe_distances(first: str, second: str) -> tuple[float, ...]:
    marks = _boundaries(first, second)
    points: set[float] = set()
    for mark in marks:
        points.update({mark - EPSILON, mark, mark + EPSILON})
    # and midway between consecutive boundaries, which is where a case
    # that one policy calls bonded and another does not actually lives
    for low, high in zip(marks, marks[1:]):
        points.add((low + high) / 2.0)
    return tuple(sorted(round(value, 6) for value in points if value > 0.2))


def _diatomic(first: str, second: str, distance: float) -> Molecule:
    return Molecule(
        symbols=[first, second],
        positions=np.array([[0.0, 0.0, 0.0], [0.0, 0.0, distance]]),
    )


def _disagreements():
    perceivers = _interchangeable_perceivers()
    names = sorted(perceivers)
    found = []
    for first, second in itertools.combinations_with_replacement(ELEMENTS, 2):
        for distance in _probe_distances(first, second):
            molecule = _diatomic(first, second, distance)
            answers = {}
            for name in names:
                try:
                    answers[name] = perceivers[name](molecule)
                except Exception as error:  # noqa: BLE001
                    answers[name] = f"raised {type(error).__name__}: {error}"
            distinct = {repr(value) for value in answers.values()}
            if len(distinct) > 1:
                found.append(
                    (
                        f"{first}-{second}",
                        distance,
                        {
                            name: (
                                sorted(value)
                                if isinstance(value, set)
                                else value
                            )
                            for name, value in answers.items()
                        },
                    )
                )
    return found


@pytest.mark.capability("policy:bond_perception")
def test_the_generator_probes_the_boundaries_it_should():
    """A generator that probes nothing proves nothing."""

    assert _probe_distances("C", "H"), "no probe points for C-H"
    # the three cases the owner named must be generated, not listed
    for first, second, named in (
        ("H", "H", 0.7410),
        ("C", "H", 1.1215),
        ("O", "H", 1.0300),
    ):
        marks = _boundaries(first, second)
        assert any(
            low <= named < high for low, high in zip(marks, marks[1:])
        ), (
            f"{first}-{second} at {named} does not lie between two of this "
            f"host's own boundaries {marks}, so the generator cannot reach "
            "it"
        )


@pytest.mark.capability("policy:bond_perception")
def test_no_two_interchangeable_perceivers_disagree():
    """The referential oracle.

    Red on 57926c6b: the analysis plane shrinks every X-H tolerance to
    0.05 A while the execution plane does not, so the two disagree across
    a 0.25 A band for every X-H pair -- which is where H2, formaldehyde's
    C-H and a hydrogen-bonded O-H all sit.
    """

    found = _disagreements()
    assert not found, (
        f"{len(found)} input(s) on which two host perceivers of one "
        "question return different answers; the host presents them as "
        "interchangeable, so a session and a sensor reading the same "
        "structure are told different things. First five: "
        + "; ".join(
            f"{pair} at {distance}: {answers}"
            for pair, distance, answers in found[:5]
        )
    )


@pytest.mark.capability("policy:bond_perception")
def test_no_agent_reachable_module_reads_a_distance_derived_bond_order():
    """A bond order is not a distance, and the host no longer claims one.

    `determine_bond_order` derives order from the same cutoff the
    adjacency uses, so it moves whenever the cutoff moves. Measured at an
    additive 0.3 A buffer: ethane's single C-C reads **2.0**, benzene
    reads **3.0**, ethene **3.0**, and formaldehyde's C=O **3.0**. One
    earlier repair widened that buffer to fix a false connectivity and
    traded it for a false bond order, and nothing noticed because no
    oracle read the result.

    Where the electrons are is the scientist's question, or a program's.
    The agent plane therefore reads no distance-derived order at all; the
    human rdkit export keeps its own perception, quarantined and declared
    legacy.
    """

    import pathlib

    offenders = []
    for area in ("chemsmart/agent", "chemsmart/analysis"):
        for path in sorted(pathlib.Path(area).rglob("*.py")):
            text = path.read_text(encoding="utf-8")
            for number, line in enumerate(text.splitlines(), start=1):
                stripped = line.strip()
                if stripped.startswith("#") or "bond_order" not in stripped:
                    continue
                offenders.append(f"{path}:{number}: {stripped[:70]}")
    assert not offenders, (
        "these agent-reachable modules read a distance-derived bond "
        f"order, which the host does not claim: {offenders}"
    )


@pytest.mark.capability("policy:bond_perception")
def test_the_perception_owner_emits_no_bond_order():
    """The owner's own surface carries adjacency and margins, nothing more."""

    from chemsmart.io.molecules import perception

    surface = set(perception.__all__)
    forbidden = {
        name
        for name in surface
        if any(
            token in name.lower()
            for token in ("order", "aromatic", "valence", "isomer", "stereo")
        )
    }
    assert not forbidden, (
        "the perception owner exports names that assert chemistry it "
        f"cannot see from distances: {sorted(forbidden)}"
    )
    policy = perception.declared_policy()
    assert "bond order" in policy["excludes"]


# The family marker is deliberate and its limit is stated: this test
# exercises the *registration* of every host policy -- that each is
# declared, owned by a module, and family-tagged -- and not the science
# of any one threshold. Family coverage does not lift the rung, so the
# six policies with no oracle of their own stay visible at `advertised`
# on the ladder rather than being quietly promoted.
@pytest.mark.capability("policy:*")
@pytest.mark.capability("policy:bond_perception")
def test_every_host_policy_is_declared_and_owned():
    """A number the host decides science through is a declared policy.

    Seven such numbers existed as bare module constants with no registry,
    no rung and no oracle. The ladder could not report that one of them
    was wrong because it had no kind that could hold the thing: its only
    numeric kind, `constant`, is the *literature* registry -- values read
    from a source text with provenance -- and a threshold is a decision
    the host makes.
    """

    from pathlib import Path

    from chemsmart.agent.capability_registry import build_capability_registry
    from chemsmart.agent.rules import HOST_POLICIES

    registry = build_capability_registry(tests_root=Path("tests"))
    rows = (
        list(registry.values())
        if isinstance(registry, dict)
        else list(registry)
    )
    policies = {row.id: row for row in rows if row.kind == "policy"}
    assert set(policies) == {item[0] for item in HOST_POLICIES}
    unwired = sorted(pid for pid, row in policies.items() if not row.wired_by)
    assert not unwired, f"declared policies that own no module: {unwired}"
    # The shared policies are the ones presented as interchangeable; a
    # legacy policy is advertised as legacy on purpose, so that its
    # disagreement is an observation rather than a silent contradiction.
    families = {row.family for row in policies.values()}
    assert families == {"shared", "legacy"}


@pytest.mark.capability("tool:plan_scientific_workflow")
def test_every_declared_hessian_role_reaches_the_sentence_the_model_reads():
    """A role the host admits is a role the model can be told about.

    The two Hessian roles were written out seven times, and one of the
    seven is the tool description the planning session actually reads.
    A role added to the table and not to that sentence is admitted by
    the host and unknown to the only party that can ask for it, which is
    the affordance-that-never-joins class in its quietest form.
    """

    from chemsmart.agent.execution import HESSIAN_CONSUMER_ROLES
    from chemsmart.agent.tool_specs import (
        build_command_compiled_tool_surface,
    )

    rendered = json.dumps(
        build_command_compiled_tool_surface().tool_definitions
    )
    for role in HESSIAN_CONSUMER_ROLES.values():
        assert role.consumer_input_id in rendered, role.consumer_input_id
        assert role.artifact_class in rendered, role.artifact_class
        for stage in role.consumer_stages:
            assert stage in rendered, stage
