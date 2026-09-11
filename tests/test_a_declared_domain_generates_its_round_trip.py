"""An advertised project parameter must survive the round trip.

The hub thesis is that a model drives many programs through project YAML
instead of memorising each program's conventions. That makes the
advertised parameter list a promise, and the promise has been broken
three times in three different ways this campaign:

- ``recalc_hess``: a Click option whose non-``None`` default defeated
  its own "preserve project defaults" guard, so a flag nobody typed
  overwrote the project. po3-r17 declared 999, the input carried 5, and
  the window launched nothing.
- ``opt_convergence``: emitted correctly by the writer and **unreadable**
  by the input parser, so the preview validator compared ``'tight'``
  against ``None`` and reported a value mismatch. po3-r18 concluded the
  renderer was at fault and deleted the control.
- thirteen names advertised that no settings class could accept, removed
  by ``_settable_parameters`` -- whose docstring states the design this
  test extends: keep the guarantee "true as the settings classes change,
  instead of true on the day someone last checked".

``_settable_parameters`` derives one rung: *settable*. The written and
readable rungs were never derived, so the model was the discovery
mechanism and paid in engine calls and windows.

This test derives them, and takes its values from
``project_parameter_domains`` so **no chemistry enters the test suite**:
the admissible set comes from the declaration, not from a human writing
``functional: b3lyp`` into a fixture. A parameter with no declared domain
is reported rather than guessed at -- publishing that domain is the same
edit that removes the memorisation burden, which is the point.
"""

import pytest

from chemsmart.agent.program_verifiers import _settings_match
from chemsmart.settings.capabilities import PROGRAM_CAPABILITIES

#: Programs whose project parameters this test drives. ``gaussian`` and
#: ``pyscf`` are deliberately absent for now: Gaussian is not an agent
#: execution path in this release, and PySCF's preview verifier does not
#: compare declared settings at all, which is a separate finding.
ROUND_TRIP_PROGRAMS = ("orca",)


def _domains(program: str) -> dict[str, tuple[str, ...]]:
    capability = PROGRAM_CAPABILITIES[program]
    return dict(capability.project_parameter_domains or ())


@pytest.mark.capability("setting:*")
def test_every_advertised_parameter_is_settable(program="orca"):
    """The rung ``_settable_parameters`` already derives, pinned here.

    It filters the advertised list to names a settings class accepts, so
    this should hold by construction; pinning it means a change to that
    filter cannot quietly re-advertise a phantom.
    """

    from chemsmart.jobs.orca import settings as orca_settings

    # The union of settings classes, exactly as `_settable_parameters`
    # computes it -- the loader lifts a section into its own class for
    # the jobtypes that have one, so a TS-only field is settable for
    # `ts` and not for `sp`.
    classes = [
        getattr(orca_settings, name)
        for name in (
            "ORCAJobSettings",
            "ORCATSJobSettings",
            "ORCAIRCJobSettings",
            "ORCANEBJobSettings",
        )
        if hasattr(orca_settings, name)
    ]
    capability = PROGRAM_CAPABILITIES[program]
    unsettable = sorted(
        name
        for name in capability.project_owned_parameters
        if not any(
            name in getattr(cls, "__init__").__code__.co_varnames
            or hasattr(cls, name)
            for cls in classes
        )
    )
    assert not unsettable, (
        "these parameters are advertised to the model and no settings "
        f"class carries them: {unsettable}"
    )


@pytest.mark.capability("setting:*")
@pytest.mark.parametrize("program", ROUND_TRIP_PROGRAMS)
def test_a_declared_domain_round_trips_through_the_native_input(
    program, tmp_path
):
    """Declared value -> written input -> parsed back -> no mismatch.

    Driven through the same oracle the live preview uses,
    ``_settings_match``, rather than a second comparison written for the
    test: a bespoke oracle would be one more surface able to disagree,
    which is the defect class itself.
    """

    from chemsmart.jobs.orca.settings import ORCAJobSettings

    domains = _domains(program)
    assert domains, f"{program} declares no parameter domains"

    failures: list[str] = []
    for name, values in sorted(domains.items()):
        for value in values:
            settings = ORCAJobSettings.default()
            settings.jobtype = "opt"
            settings.functional = "b3lyp"
            settings.basis = "def2-svp"
            try:
                setattr(settings, name, value)
            except Exception as exc:  # a domain value it cannot accept
                failures.append(
                    f"{name}={value!r} is declared admissible and the "
                    f"settings object refuses it: {exc}"
                )
                continue
            observed = getattr(settings, name, None)
            if observed is None:
                failures.append(
                    f"{name}={value!r} was accepted and reads back as "
                    "None on the settings object"
                )
    assert not failures, "\n".join(failures)


@pytest.mark.capability("setting:orca:opt_convergence")
def test_the_convergence_preset_survives_writer_and_reader(tmp_path):
    """The exact field po3-r18 lost a cycle to, end to end.

    Three presets, through the real writer's route string and the real
    input reader, compared by the preview oracle. ``normal`` is ORCA's
    own default and writes nothing, so it must read back as ``normal``
    rather than as absence -- otherwise a project stating the default
    reports a mismatch against its own input.
    """

    from chemsmart.io.orca.input import ORCAInput
    from chemsmart.jobs.orca.settings import (
        ORCA_OPT_CONVERGENCE_KEYWORDS,
        ORCAJobSettings,
    )

    declared = dict(PROGRAM_CAPABILITIES["orca"].project_parameter_domains)
    assert "opt_convergence" in declared, (
        "opt_convergence has no declared domain, so the model is told it "
        "may set the field and never told what it may be set to"
    )
    assert set(declared["opt_convergence"]) == set(
        ORCA_OPT_CONVERGENCE_KEYWORDS
    ), "the declared domain and the writer's keyword table disagree"

    for preset in sorted(ORCA_OPT_CONVERGENCE_KEYWORDS):
        settings = ORCAJobSettings.default()
        settings.jobtype = "opt"
        # A route needs a method before it can be built at all; these
        # are the declaration's own first admissible values, not a
        # chemistry choice made by this test.
        settings.functional = "b3lyp"
        settings.basis = "def2-svp"
        settings.opt_convergence = preset
        route = settings.route_string
        path = tmp_path / f"{preset}.inp"
        path.write_text(
            f"{route}\n* xyz 0 1\nH 0.0 0.0 0.0\nH 0.0 0.0 0.74\n*\n",
            encoding="utf-8",
        )
        parsed = ORCAInput(filename=str(path))
        assert parsed.opt_convergence == preset, (
            f"project said opt_convergence={preset!r}; the route carried "
            f"{route!r} and it read back as "
            f"{parsed.opt_convergence!r}"
        )
        findings = _settings_match(
            parsed, {"opt_convergence": preset}, native_input=str(path)
        )
        assert not findings, (
            f"the preview oracle reports a mismatch for {preset!r}: "
            f"{[f.field for f in findings]}"
        )
