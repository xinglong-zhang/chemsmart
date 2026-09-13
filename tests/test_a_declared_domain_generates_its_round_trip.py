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

#: Programs whose project parameters the settings-object round trip drives.
#: ``gaussian`` is deliberately absent: it is not an agent execution path in
#: this release. ``pyscf`` is driven by its own, stronger round trip below
#: -- through the public ``run --fake`` command and the live preview
#: verifier -- because a PySCF value is written into a driver CONFIG and
#: read back from an artifact, not parsed from a route string.
ROUND_TRIP_PROGRAMS = ("orca",)
PYSCF_ROUND_TRIP_PROGRAM = "pyscf"


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


# ----------------------------------------------------------------------
# PySCF: declared value -> run --fake -> fake artifact spec -> no mismatch
# ----------------------------------------------------------------------

_WATER_XYZ = (
    "3\nwater\nO 0.0 0.0 0.1173\nH 0.0 0.7572 -0.4692\nH 0.0 -0.7572 -0.4692\n"
)
_HYDROXYL_XYZ = "2\nhydroxyl radical\nO 0.0 0.0 0.0\nH 0.0 0.0 0.97\n"


def _pyscf_case(parameter: str, value: str):
    """The stage, section and molecule that make one declared value legal.

    Only the declaration's own words appear: a method or basis literal is
    the same minimal spelling the PySCF hardening tests use, never a
    chemistry choice made by this test.  The ``unrestricted`` manifold is
    the one an open-shell reference has, so it rides a doublet input.
    """

    base = {"basis": "def2-svp"}
    dft = {**base, "functional": "b3lyp"}
    response = {**dft, "response_method": "tda", "nstates": 2}
    if parameter == "ab_initio":
        return "sp", {**base, "ab_initio": value}, _WATER_XYZ, (0, 1)
    if parameter == "frozen_core":
        return (
            "sp",
            {**base, "ab_initio": "mp2", "frozen_core": value},
            _WATER_XYZ,
            (0, 1),
        )
    if parameter == "defgrid":
        return "sp", {**dft, "defgrid": value}, _WATER_XYZ, (0, 1)
    if parameter == "opt_solver":
        return "opt", {**dft, "opt_solver": value}, _WATER_XYZ, (0, 1)
    if parameter == "solvent_model":
        return (
            "sp",
            {**dft, "solvent_model": value, "solvent_id": "water"},
            _WATER_XYZ,
            (0, 1),
        )
    if parameter == "response_method":
        return (
            "td",
            {
                **response,
                "response_method": value,
                "state_manifold": "singlet",
            },
            _WATER_XYZ,
            (0, 1),
        )
    if parameter == "state_manifold":
        if value == "unrestricted":
            return (
                "td",
                {**response, "state_manifold": value},
                _HYDROXYL_XYZ,
                (0, 2),
            )
        return "td", {**response, "state_manifold": value}, _WATER_XYZ, (0, 1)
    raise AssertionError(
        f"a PySCF domain was declared for {parameter!r} and this round trip "
        "does not know which stage makes it legal; extend _pyscf_case"
    )


def _pyscf_domain_cases():
    return [
        (parameter, value)
        for parameter, values in sorted(
            _domains(PYSCF_ROUND_TRIP_PROGRAM).items()
        )
        for value in values
    ]


def _pyscf_fake_preview(tmp_path, parameter, value):
    """Run the public ``run --fake`` command and validate its workspace.

    Returns the verifier's receipt for the expectation built from the
    loader's own settings rows -- the chain the planning session walks:
    ``validate_project_yaml`` -> ``build_preview_expectation`` ->
    ``validate_preview_workspace`` -- and the expectation itself, so a
    caller can also show the instrument red.
    """

    import yaml
    from click.testing import CliRunner

    from chemsmart.agent._contracts import TrustedArtifactRefV1, file_sha256
    from chemsmart.agent.capabilities import (
        ProgramCapabilityQueryV1,
        build_command_compiled_preview_overlay,
        build_program_component_conformance_receipt,
        load_program_capabilities,
        query_capability,
    )
    from chemsmart.agent.cli_schema import build_live_click_schema
    from chemsmart.agent.live_session import _preview_server_profile
    from chemsmart.agent.program_verifiers import (
        build_preview_expectation,
        validate_preview_workspace,
    )
    from chemsmart.agent.projects import validate_project_yaml
    from chemsmart.cli.main import entry_point

    stage, section, xyz_text, (charge, multiplicity) = _pyscf_case(
        parameter, value
    )
    xyz = tmp_path / "input.xyz"
    xyz.write_text(xyz_text, encoding="utf-8")
    project = tmp_path / "project.yaml"
    project.write_text(yaml.safe_dump({stage: section}), encoding="utf-8")
    server = tmp_path / "preview-server.yaml"
    server.write_text(_preview_server_profile(), encoding="utf-8")
    workspace = tmp_path / "workspace"
    workspace.mkdir()

    argv = [
        "run",
        "--server",
        str(server),
        "--fake",
        "--no-scratch",
        "pyscf",
        "--project",
        str(project),
        "--filename",
        str(xyz),
        "--charge",
        str(charge),
        "--multiplicity",
        str(multiplicity),
        "--no-gpu",
        stage,
    ]
    runner = CliRunner()
    with runner.isolated_filesystem(temp_dir=workspace) as cwd:
        result = runner.invoke(entry_point, argv)
        preview_dir = cwd
    assert result.exit_code == 0, (
        f"{parameter}={value!r} is declared admissible and the public "
        f"command refused it: {result.output[-600:]} {result.exception!r}"
    )

    registry = load_program_capabilities()
    live_schema = build_live_click_schema()
    conformance = build_program_component_conformance_receipt(
        program="pyscf",
        registry_sha256=registry.registry_sha256,
        live_cli_schema_sha256=live_schema.schema_sha256,
        fixture_bundle_sha256="1" * 64,
        covered_jobtypes=("hess", "opt", "sp", "td"),
        covered_engines=("cpu",),
        covered_engine_job_pairs=(
            ("cpu", "hess"),
            ("cpu", "opt"),
            ("cpu", "sp"),
            ("cpu", "td"),
        ),
        compiler_receipt_sha256="2" * 64,
        preview_receipt_sha256="3" * 64,
        preflight_receipt_sha256="4" * 64,
        verifier_receipt_sha256="5" * 64,
        compiler_status="passed",
        preview_status="passed",
        preflight_status="passed",
        verifier_status="passed",
    )
    overlay = build_command_compiled_preview_overlay(
        registry, conformance_receipts=(conformance,), live_schema=live_schema
    )
    capability = query_capability(
        ProgramCapabilityQueryV1("pyscf", stage, "cpu"),
        registry=registry,
        live_schema=live_schema,
        overlay=overlay,
    )

    def _artifact(path, kind):
        return TrustedArtifactRefV1(
            artifact_id=path.stem,
            kind=kind,
            sha256=file_sha256(path),
            size_bytes=path.stat().st_size,
            path=str(path),
            cli_value=str(path),
        )

    validation = validate_project_yaml(
        _artifact(project, "project_yaml"), capability=capability
    )
    assert validation.status == "valid", validation.diagnostic
    expectation = build_preview_expectation(
        program="pyscf",
        jobtype=stage,
        input_artifact=_artifact(xyz, "geometry_xyz"),
        project=validation,
        charge=charge,
        multiplicity=multiplicity,
    )
    return validate_preview_workspace(expectation, preview_dir), validation


@pytest.mark.capability("setting:pyscf:*")
@pytest.mark.parametrize(("parameter", "value"), _pyscf_domain_cases())
def test_a_declared_pyscf_domain_round_trips_through_the_fake_preview(
    tmp_path, parameter, value
):
    """Declared value -> public run --fake -> artifact spec -> no mismatch.

    The verifier compares the loader's own settings rows with the applied
    spec the driver wrote into the fake artifact through the live
    provenance oracle (``verify_provenance``), so a value the writer drops
    or rewrites is a ``preview.semantic.mismatch`` here, before any engine.
    """

    receipt, _validation = _pyscf_fake_preview(tmp_path, parameter, value)
    assert receipt.status == "valid", [
        (item.field, item.expected, item.observed) for item in receipt.findings
    ]


@pytest.mark.capability("setting:pyscf:*")
def test_the_pyscf_preview_verifier_can_report_red(tmp_path):
    """An instrument that cannot report red witnesses nothing.

    The same workspace, validated against an expectation whose declared
    values differ from what the driver wrote, names the field.
    """

    from chemsmart.agent._contracts import canonical_sha256
    from chemsmart.agent.program_verifiers import (
        ProgramPreviewExpectationV1,
        validate_preview_workspace,
    )

    receipt, validation = _pyscf_fake_preview(tmp_path, "ab_initio", "mp2")
    assert receipt.status == "valid"
    preview_dir = next(tmp_path.joinpath("workspace").iterdir())
    rows = dict(validation.settings)
    rows["ab_initio"] = "ccsd"
    body = {
        "schema_version": "chemsmart.program-preview-expectation.v1",
        "program": "pyscf",
        "jobtype": "sp",
        "input_artifact": None,
        "project_receipt_sha256": validation.receipt_sha256,
        "project_sha256": validation.project_sha256,
        "settings": tuple(sorted(rows.items())),
        "charge": 0,
        "multiplicity": 1,
    }
    from chemsmart.agent._contracts import TrustedArtifactRefV1, file_sha256

    xyz = tmp_path / "input.xyz"
    body["input_artifact"] = TrustedArtifactRefV1(
        artifact_id="input",
        kind="geometry_xyz",
        sha256=file_sha256(xyz),
        size_bytes=xyz.stat().st_size,
        path=str(xyz),
        cli_value=str(xyz),
    )
    wrong = ProgramPreviewExpectationV1(
        **body, expectation_sha256=canonical_sha256(body)
    )
    red = validate_preview_workspace(wrong, preview_dir)
    assert red.status == "invalid"
    assert {item.field for item in red.findings} >= {
        "settings.ab_initio",
        "settings.method",
    }, [item.field for item in red.findings]
