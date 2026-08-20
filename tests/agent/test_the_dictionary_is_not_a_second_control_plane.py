"""The vocabulary tells the truth early, and never outranks the program.

A session authored MN15 into an ORCA project. The YAML validated, the CLI
compiled, the native input was generated, and only the safe preview's
membership test refused it -- after the whole stage was built, and on the
authority of a parser list nothing had ever verified. A campaign probe then
asked the installed ORCA 6.1.1 itself: MN15 is genuinely not there, fourteen
functionals the binary DOES accept were missing from the list, and the
repository shipped an ORCA template using mn15 that its own preview would
refuse.

The dictionary that came out of that has three states, and only one refuses:

- a value in this program's declared domain passes silently;
- a value another program's domain declares is refused at authoring time,
  naming the implementing programs and the nearest names here -- the truth
  MN15 cost a whole stage to learn now arrives before anything is built;
- a value no declared domain knows passes with an advisory, because the
  dictionary is not the authority on chemistry -- the program validator and
  the safe preview are, and a name the dictionary merely does not know must
  never be blocked on the dictionary's word alone.

These checks pin the boundary between those states, in both directions, and
the single-source rule that keeps the parser, the verifier, and the receipt
reading one artifact.
"""

from __future__ import annotations

import pytest

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.projects import project_document, render_project_yaml


def _render(program, section):
    return render_project_yaml(
        project_document(program=program, sections={"gas": section})
    )


def test_a_gaussian_functional_in_an_orca_project_is_refused_with_the_truth():
    with pytest.raises(ContractError) as refusal:
        _render("orca", {"functional": "mn15", "basis": "def2-tzvp"})

    message = str(refusal.value)
    assert "not in orca's functional vocabulary" in message
    assert "implemented by: gaussian" in message
    assert "probe-verified" in message


def test_the_redirect_works_in_the_other_direction_too():
    with pytest.raises(ContractError) as refusal:
        _render("gaussian", {"functional": "wb97m-v", "basis": "def2tzvp"})

    message = str(refusal.value)
    assert "implemented by: orca" in message
    assert "Gaussian 16" in message


def test_a_probe_added_functional_authors_cleanly():
    """wb97m-v is real ORCA 6 and was refused at preview before the probe."""

    receipt = _render("orca", {"functional": "wb97m-v", "basis": "def2-tzvp"})

    assert receipt.status == "candidate_rendered"
    assert receipt.vocabulary_advisories == ()


def test_an_unknown_name_passes_with_an_advisory_not_a_refusal():
    receipt = _render(
        "orca", {"functional": "superduperfunc2029", "basis": "def2-tzvp"}
    )

    assert receipt.status == "candidate_rendered"
    assert len(receipt.vocabulary_advisories) == 1
    advisory = receipt.vocabulary_advisories[0]
    assert "outside every declared" in advisory
    assert "safe preview remain the authority" in advisory
    assert "project.render.vocabulary_advisory" in receipt.rule_ids


def test_a_clean_render_is_byte_identical_to_before_the_dictionary():
    """Old receipts keep their identity: no advisories, no new body field."""

    receipt = _render("orca", {"functional": "b3lyp", "basis": "def2-svp"})

    assert receipt.vocabulary_advisories == ()
    assert receipt.rule_ids == ("project.render.candidate_only",)


def test_the_domains_are_projections_of_the_io_tables():
    """Single source: parser, verifier and receipt cannot fork."""

    from chemsmart.io.orca import (
        ORCA_ALL_BASIS_SETS,
        ORCA_ALL_FUNCTIONALS,
        ORCA_ALL_SOLVENTS,
    )
    from chemsmart.settings.capabilities import PROGRAM_CAPABILITIES

    domains = dict(PROGRAM_CAPABILITIES["orca"].project_parameter_domains)

    def normalized(values):
        return tuple(sorted({str(v).strip().lower() for v in values}))

    assert domains["functional"] == normalized(ORCA_ALL_FUNCTIONALS)
    assert domains["basis"] == normalized(ORCA_ALL_BASIS_SETS)
    assert domains["solvent_id"] == normalized(ORCA_ALL_SOLVENTS)


def test_the_probe_facts_are_in_the_domains():
    from chemsmart.settings.capabilities import PROGRAM_CAPABILITIES

    orca = dict(PROGRAM_CAPABILITIES["orca"].project_parameter_domains)
    gaussian = dict(PROGRAM_CAPABILITIES["gaussian"].project_parameter_domains)

    assert "mn15" not in orca["functional"]
    assert "mn15" in gaussian["functional"]
    assert "wb97m-v" in orca["functional"]
    assert "xalpha" not in orca["functional"]


def test_the_shipped_template_no_longer_contradicts_the_binary():
    """The repo shipped an ORCA template its own preview refused."""

    from pathlib import Path

    template = Path(
        "chemsmart/settings/templates/.chemsmart/orca/test_qmmm.yaml"
    ).read_text(encoding="utf-8")

    assert "functional: mn15" not in template
    assert 'high_level_functional: "mn15"' not in template


def test_the_vocabulary_reaches_the_capability_receipt():
    """Delivery: the model's first tool call in every recorded session."""

    import tempfile
    from pathlib import Path

    from chemsmart.agent.runtime.event_store import RuntimeEventStore
    from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1

    with tempfile.TemporaryDirectory() as tmp:
        host = CommandCompiledToolHostV1(
            event_store=RuntimeEventStore(
                Path(tmp) / "events.jsonl", session_id="s1"
            ),
            task_spec_sha256s=("a" * 64,),
            approved_workspace=Path(tmp) / "workspace",
        )
        payload = host.dispatch(
            turn_id="t1",
            tool_name="inspect_program_capability",
            arguments={"program": "orca", "jobtype": "opt", "engine": "cpu"},
        )["result"]

    domains = dict(
        tuple(item)
        for item in payload["capability"]["project_parameter_domains"]
    )
    assert "wb97m-v" in domains["functional"]
    assert "mn15" not in domains["functional"]
    assert len(domains["basis"]) == 356
