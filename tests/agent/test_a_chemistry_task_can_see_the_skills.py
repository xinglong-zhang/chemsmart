"""Knowledge the model is told to use must be knowledge it was shown.

`consult_domain_skill` describes itself as reading "a skill listed in the
system prompt", so the list is not a convenience -- it is the whole of what the
model knows exists. That list came from `activated_skill_documents`, which
gated every skill behind a plain lowercase substring match against each pack's
`activation_terms`. Those terms are program brand names: `orca`, `gaussian`,
`xtb`, `pyscf`, `gfn2`, `cpcm`, `dlpno`.

The consequence was invisible from inside the code. A task that describes its
chemistry correctly without naming a product -- "a DFT setup and a cheaper
semi-empirical one" -- matched nothing, so the prompt listed no skills, and a
tool that says to choose from a list was pointed at an empty one. Three
recorded observations of a frozen task were read as evidence about how the
model behaves; they were measuring this gate.

None of the three installed skills is program-specific, so the gate never
discriminated between programs; it only decided whether general chemistry
knowledge was mentioned at all.

The task text quoted below is the opening of that frozen task, kept verbatim
because a paraphrase naming a program would pass while the real one failed.
"""

from __future__ import annotations

from chemsmart.agent.knowledge_packs import (
    BUILTIN_PROGRAM_PACKS,
    activate_program_knowledge,
)
from chemsmart.agent.live_session import (
    _system_prompt,
    activated_skill_documents,
    cross_program_skill_ids,
)

#: Verbatim from the frozen cycle-027 control task.
_NAMES_NO_PROGRAM = (
    "We are checking how well our usual setup handles alkane branching. "
    "Four gas-phase calculations have finished and their outputs are in the "
    "workspace: straight-chain hexane and the branched isomer "
    "2,2-dimethylbutane, each done twice, once with a DFT setup and once "
    "with a cheaper semi-empirical one."
)


def test_the_control_task_names_no_activation_term():
    """If this ever fails the fixture stopped reproducing the real case."""

    lowered = _NAMES_NO_PROGRAM.lower()
    for pack in BUILTIN_PROGRAM_PACKS:
        for term in pack.activation_terms:
            assert term.lower() not in lowered, (
                f"the task mentions {term!r}, so it would have activated "
                f"{pack.pack_id} even before this fix and proves nothing"
            )


def test_a_chemistry_task_that_names_no_program_still_sees_the_skills():
    _, documents = activated_skill_documents(_NAMES_NO_PROGRAM)

    assert {doc.skill_id for doc in documents} == set(cross_program_skill_ids())
    assert documents, "the model would be told to choose from an empty list"


def test_the_listing_reaches_the_prompt_the_tool_points_at():
    """The production shape: the prompt, not the resolver's return value."""

    _, documents = activated_skill_documents(_NAMES_NO_PROGRAM)
    prompt = _system_prompt(
        None,
        bounded_review_requested=False,
        skill_index=tuple(doc.index_entry() for doc in documents),
    )

    assert "consult_domain_skill" in prompt
    for skill_id in cross_program_skill_ids():
        assert skill_id in prompt


def test_availability_is_not_a_claim_that_a_pack_fired():
    """Making a skill readable must not fabricate activation provenance."""

    pack_sha256s, documents = activated_skill_documents(_NAMES_NO_PROGRAM)

    assert documents
    assert pack_sha256s == ()


def test_program_specific_activation_still_discriminates():
    """The keyword gate keeps the job it is genuinely for."""

    silent = activate_program_knowledge(
        request=_NAMES_NO_PROGRAM, program="orca", engine="cpu"
    )
    named = activate_program_knowledge(
        request="Run this in ORCA with CPCM water", program="orca", engine="cpu"
    )

    assert silent.activated_pack_sha256s == ()
    assert named.activated_pack_sha256s != ()


def test_no_installed_skill_is_program_specific():
    """Explains why gating them by program name discriminated nothing.

    A future program-specific skill is free to exist; this records that today
    every skill applies to any chemistry, which is what makes unconditional
    availability correct rather than merely convenient.
    """

    shared = set(cross_program_skill_ids())
    for pack in BUILTIN_PROGRAM_PACKS:
        assert set(pack.skill_ids) == shared, (
            f"{pack.pack_id} carries a program-specific skill; availability "
            f"must then be reconsidered rather than left unconditional"
        )
