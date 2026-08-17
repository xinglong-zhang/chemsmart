"""Domain-knowledge skill resolution, authority limits, and skills-off parity."""

import re

import pytest

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.knowledge_packs import (
    BUILTIN_PROGRAM_PACKS,
    activate_program_knowledge,
    skills_for_activation,
)
from chemsmart.agent.skills import (
    SkillDocumentV1,
    available_skill_ids,
    resolve_skill,
    skills_enabled,
)
from chemsmart.agent.skills.conventions import (
    BUILTIN_CONVENTION_RULES,
    conventions_for_scope,
)
from chemsmart.agent.skills.rules import (
    DeterministicConventionRuleV1,
    build_convention_rule,
)

_CH2_TASK = (
    "Determine the adiabatic singlet-triplet splitting T0 of CH2 using PySCF "
    "through ChemSmart."
)


def _write_overlay(root, skill_id, description, body):
    target = root / skill_id
    target.mkdir(parents=True)
    (target / "SKILL.md").write_text(
        f"---\nname: {skill_id}\nversion: 9.9.9\ndescription: {description}\n"
        f"---\n\n{body}\n",
        encoding="utf-8",
    )
    return target


def test_builtin_conventions_skill_resolves():
    document = resolve_skill("scientific-conventions")
    assert document is not None
    assert document.origin == "builtin"
    assert document.skill_id == "scientific-conventions"
    assert document.body.strip()
    assert "scientific-conventions" in available_skill_ids()


def test_overlay_takes_precedence_over_builtin(tmp_path):
    _write_overlay(
        tmp_path,
        "scientific-conventions",
        "Overlaid description",
        "Overlaid body.",
    )
    document = resolve_skill("scientific-conventions", overlay_root=tmp_path)
    assert document.origin == "overlay"
    assert document.skill_version == "9.9.9"
    assert document.body == "Overlaid body."


def test_unknown_skill_resolves_to_none(tmp_path):
    assert resolve_skill("no-such-skill", overlay_root=tmp_path) is None


def test_skill_document_digest_is_body_bound():
    document = resolve_skill("scientific-conventions")
    with pytest.raises(ContractError):
        SkillDocumentV1(
            schema_version=document.schema_version,
            skill_id=document.skill_id,
            skill_version=document.skill_version,
            description=document.description,
            body=document.body + " tampered",
            origin=document.origin,
            body_sha256=document.body_sha256,
            document_sha256=document.document_sha256,
        )


def test_convention_rule_refuses_readiness_and_accuracy_authority():
    rule = BUILTIN_CONVENTION_RULES[0]
    body = rule.as_dict()
    for field in ("readiness_authority", "accuracy_authority"):
        with pytest.raises(ContractError):
            DeterministicConventionRuleV1(**{**body, field: True})


def test_convention_rule_digest_is_content_bound():
    rule = BUILTIN_CONVENTION_RULES[0]
    with pytest.raises(ContractError):
        DeterministicConventionRuleV1(
            **{**rule.as_dict(), "convention": "something else"}
        )


def test_convention_rule_rejects_unknown_scope():
    with pytest.raises(ContractError):
        build_convention_rule(
            rule_id="bad-scope",
            skill_id="scientific-conventions",
            scope="readiness",
            subject="x",
            convention="y",
        )


def test_conventions_are_partitioned_by_scope():
    thermo = conventions_for_scope("thermochemistry")
    results = conventions_for_scope("result_validation")
    assert thermo and results
    assert not set(thermo) & set(results)
    assert all(item.scope == "thermochemistry" for item in thermo)


def test_packs_surface_the_cross_program_skills():
    receipt = activate_program_knowledge(
        request=_CH2_TASK, program="pyscf", engine="cpu"
    )
    assert receipt.advisory_only is True
    assert skills_for_activation(receipt) == (
        "method-adequacy",
        "scientific-conventions",
        "typed-analysis-contract",
    )


def test_pack_activation_still_respects_exclusions():
    receipt = activate_program_knowledge(
        request="A periodic PySCF hessian calculation",
        program="pyscf",
        engine="cpu",
    )
    assert receipt.activated_pack_sha256s == ()
    assert skills_for_activation(receipt) == ()


def test_packs_carry_convention_rule_digests():
    digests = {
        digest
        for pack in BUILTIN_PROGRAM_PACKS
        for digest in pack.convention_rule_sha256s
    }
    assert digests == {item.rule_sha256 for item in BUILTIN_CONVENTION_RULES}


@pytest.mark.parametrize(
    "value,expected", [("0", False), ("false", False), ("1", True)]
)
def test_skills_enabled_toggle(monkeypatch, value, expected):
    monkeypatch.setenv("CHEMSMART_AGENT_SKILLS", value)
    assert skills_enabled() is expected


def test_skills_off_restores_historical_prompt_and_tool_surface(monkeypatch):
    from chemsmart.agent.live_session import (
        _system_prompt,
        activated_skill_documents,
    )
    from chemsmart.agent.tool_specs import build_command_compiled_tool_surface

    monkeypatch.delenv("CHEMSMART_AGENT_SKILLS", raising=False)
    _, documents = activated_skill_documents(_CH2_TASK)
    assert documents
    enabled_surface = build_command_compiled_tool_surface()
    assert "consult_domain_skill" in {
        item["function"]["name"] for item in enabled_surface.tool_definitions
    }

    monkeypatch.setenv("CHEMSMART_AGENT_SKILLS", "0")
    _, disabled_documents = activated_skill_documents(_CH2_TASK)
    assert disabled_documents == ()
    disabled_surface = build_command_compiled_tool_surface()
    assert "consult_domain_skill" not in {
        item["function"]["name"] for item in disabled_surface.tool_definitions
    }
    # A skills-off session is byte-identical to the pre-skill surface.
    assert _system_prompt({}, skill_index=()) == _system_prompt({})
    assert (
        disabled_surface.tool_schema_sha256
        != enabled_surface.tool_schema_sha256
    )


#: Systems used as benchmark probes while authoring the skill.  The skill must
#: state general principles; naming a probe system inside it is the signature of
#: overfitting the knowledge to its own test case.
_BENCHMARK_PROBE_TERMS = (
    "ch2",
    "ch₂",
    "methylene",
    "h2o",
    "h₂o",
    "water",
    "o2",
    "o₂",
    "dioxygen",
    "carbene",
)


#: Every skill is held to the authoring norm, not just the first one written.
_GENERAL_SKILLS = (
    "method-adequacy",
    "scientific-conventions",
    "typed-analysis-contract",
)


@pytest.mark.parametrize("skill_id", _GENERAL_SKILLS)
def test_skill_states_general_principles_not_benchmark_systems(skill_id):
    body = resolve_skill(skill_id).body.lower()
    named = [term for term in _BENCHMARK_PROBE_TERMS if term in body]
    assert not named, (
        f"{skill_id} must generalise: it names the benchmark probe "
        f"system(s) {named}. State the principle that generates the fact "
        "instead of the fact for one system."
    )


@pytest.mark.parametrize("skill_id", _GENERAL_SKILLS)
def test_skill_encodes_no_expected_value_or_tolerance(skill_id):
    """A skill states principles; a number in one is an answer in disguise.

    Cutoffs, tolerances and expected magnitudes belong to the validators and
    the project settings, which own them and can be checked against a result.
    A skill asserting one would be a second, unfalsifiable control plane.
    """

    body = resolve_skill(skill_id).body
    for pattern in (
        r"\d+\s*(?:kcal|kj)\b",
        r"\bwithin\s+\d",
        r"\btolerance of\b",
        r"\bcutoff of\b",
        r"\bmust be (?:below|above|less than|greater than)\b",
    ):
        assert not re.search(pattern, body, re.IGNORECASE), (
            f"{skill_id} states a numeric expectation matching {pattern!r}; "
            "skills carry principles, not thresholds or answers"
        )


def test_method_adequacy_covers_the_gaps_it_exists_for():
    """The two failures that motivated it, plus the neighbouring ones."""

    body = resolve_skill("method-adequacy").body.lower()
    for concept in (
        "diffuse",
        "dispersion",
        "superposition",
        "cancel",
        "conformer",
        "uncertainty",
        "multireference",
        "spin",
    ):
        assert concept in body, f"method-adequacy omits: {concept}"


def test_the_analysis_contract_skill_does_not_restate_a_stale_demand():
    """It must describe the contract as implemented, not as it once was.

    The semantic role is derived by the host now.  A skill that revived the
    old instruction to author one per input would reintroduce the failure the
    prompt and the argument description were repaired to stop.
    """

    body = resolve_skill("typed-analysis-contract").body.lower()
    assert "optional" in body
    assert "unique within its request" in body or "unique within" in body


def test_skill_covers_more_than_one_quantity_class():
    """A general skill serves the whole agent, not one benchmark case."""

    body = resolve_skill("scientific-conventions").body.lower()
    for concept in (
        "reaction",
        "barrier",
        "interaction",
        "vertical",
        "adiabatic",
        "symmetry number",
        "standard state",
        "unit",
    ):
        assert concept in body, f"skill omits the general concept: {concept}"


def test_thermochemistry_discloses_its_conventions():
    """The analysis surfaces the conventions in force without changing values."""

    from chemsmart.analysis.thermochemistry import Thermochemistry

    disclosed = Thermochemistry.applied_conventions.fget(None)
    subjects = {item["subject"] for item in disclosed}
    assert "standard state pressure" in subjects
    assert "rotational symmetry number" in subjects
    assert all(item["readiness_authority"] is False for item in disclosed)
    assert all(item["accuracy_authority"] is False for item in disclosed)


def test_result_validation_conventions_are_disclosure_only():
    from chemsmart.jobs.pyscf.validation import _result_validation_conventions

    rules = _result_validation_conventions()
    assert rules
    assert all(item.scope == "result_validation" for item in rules)
    assert all(item.readiness_authority is False for item in rules)


def test_skill_index_entry_is_one_line():
    document = resolve_skill("scientific-conventions")
    entry = document.index_entry()
    assert "\n" not in entry
    assert entry.startswith("scientific-conventions:")
