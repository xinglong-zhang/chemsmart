"""What PySCF may take over is read from the registry, not from a list.

The Gaussian-to-PySCF transfer matrix carried two hand-written sets of
job families, and one of them forbade ``td`` for a release in which
PySCF ``td`` executes under sealed live goals.  Nothing was red, because
this gate had no test at all: the words and the capability they claimed
were never compared.

The mapping from a Gaussian job family to the PySCF job types that would
carry it stays a judgement about chemistry.  Whether those job types can
run is the registry's word, per engine, and these are the assertions
that hold the two together.
"""

from __future__ import annotations

import pytest

from chemsmart.agent.capabilities import (
    ProgramCapabilityQueryV1,
    query_capability,
)
from chemsmart.agent.knowledge import (
    _PYSCF_SUBSTITUTION_JOB_TYPES,
    assess_typed_program_substitution,
    build_program_substitution_request,
    pyscf_executable_jobtypes,
)

pytestmark = pytest.mark.capability("tool:assess_program_candidate")

_FAMILY_RULES = frozenset(
    {
        "program.substitution.job_family_unsupported",
        "program.substitution.job_family_forbidden",
    }
)


def _family_rules(*, families, engine):
    capability = query_capability(
        ProgramCapabilityQueryV1("pyscf", "sp", engine)
    )
    request = build_program_substitution_request(
        request_id="substitution-probe",
        requested_program="gaussian",
        selected_program="pyscf",
        requested_engine="cpu",
        selected_engine=engine,
        job_families=families,
        method_family="dft",
        method_name="b3lyp",
    )
    receipt = assess_typed_program_substitution(request, capability)
    return _FAMILY_RULES & set(receipt.rule_ids)


def test_a_family_whose_job_types_execute_is_not_refused_as_a_family():
    executable = pyscf_executable_jobtypes("cpu")
    transferable = sorted(
        family
        for family, jobtypes in _PYSCF_SUBSTITUTION_JOB_TYPES.items()
        if jobtypes and set(jobtypes).issubset(executable)
    )
    # The join: every family the registry can actually carry passes the
    # family rules, and the two can no longer disagree by one word.
    assert transferable, "the registry executes no PySCF job type at all"
    assert _family_rules(families=transferable, engine="cpu") == set()
    for family in transferable:
        assert _family_rules(families=(family,), engine="cpu") == set()


def test_the_release_that_runs_td_no_longer_forbids_substituting_it():
    assert "td" in pyscf_executable_jobtypes("cpu")
    assert _family_rules(families=("td",), engine="cpu") == set()


def test_a_job_type_the_engine_cannot_execute_is_refused_on_that_engine():
    # Every GPU4PySCF cell is preview-only, so the same family that
    # transfers on the CPU engine is refused on the GPU one.
    assert pyscf_executable_jobtypes("gpu") == frozenset()
    assert _family_rules(families=("sp",), engine="gpu") == {
        "program.substitution.job_family_forbidden"
    }


def test_a_family_with_no_pyscf_job_type_behind_it_is_unsupported():
    for family in ("neb", "qmmm", "link", "not_a_job_family"):
        assert _family_rules(families=(family,), engine="cpu") == {
            "program.substitution.job_family_unsupported"
        }
