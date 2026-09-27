"""A setting a project states reaches the input the program reads.

CHEMSMART owns the translation from a project's typed settings to each
program's native input. A translation that drops what the scientist said
runs a second, hidden program, and nothing downstream can see it: the
engine is honest about the input it was given.

Census R10 Q31 (the public ``run --fake`` over every program, job type
and settings field, each written at a non-default value and compared with
the input written without it; 783 cases on b96e63ee) found the class open
in five shapes. This file pins the first: a command option that shares a
project setting's name and has a default -- a value, not None -- is
applied whether or not anyone typed it, so it replaces what the project
stated. ORCA's ``--tssearch-type`` defaulted to ``optts`` and ran OptTS
for a project asking for ScanTS; five ORCA IRC booleans defaulted to
False and were applied unconditionally, so no project value of theirs
ever reached ``%irc``; ``orca opt`` wrote ``invert_constraints`` from its
default; the ORCA group merged ``forces=False`` into every job. R10 Q28
paid for the same shape on Gaussian's IRC (every IRC ran maxpoints=512).
"""

from __future__ import annotations

import click
import pytest
import yaml

from chemsmart.settings.capabilities import AGENT_PROGRAM_JOBTYPES

#: The smallest level each program's loader accepts in every section.
_LEVEL = {
    "gaussian": {"functional": "b3lyp", "basis": "def2svp"},
    "orca": {"functional": "b3lyp", "basis": "def2-svp"},
    "pyscf": {"functional": "b3lyp", "basis": "def2-svp"},
    "xtb": {"gfn_version": "gfn2"},
}


def _project_fields(program, tmp_path):
    """Every field a project section of *program* can state.

    Asked of the live loader, not listed here: one project carrying a
    section for each job type the Agent can run is loaded, and each job
    type's settings object -- in the class the loader lifts that section
    into -- names the fields a project can set for it.
    """

    import importlib

    from chemsmart.agent.projects import _yaml_project_loader

    jobtypes = AGENT_PROGRAM_JOBTYPES[program]
    sections = {jobtype: dict(_LEVEL[program]) for jobtype in jobtypes}
    if program in {"gaussian", "orca"}:
        sections["gas"] = dict(_LEVEL[program])
    if program == "gaussian" and "link" in sections:
        sections["link"]["jobtype"] = "opt"
    if program in {"orca", "pyscf"} and "td" in sections:
        sections["td"].update(
            response_method="tddft", nstates=3, state_manifold="singlet"
        )
    if program == "pyscf" and "irc" in sections:
        sections["irc"]["irc_direction"] = "forward"
    if program == "orca" and "neb" in sections:
        sections["neb"] = {"semiempirical": "XTB2", "joboption": "NEB-CI"}
        sections["neb"]["nimages"] = 4
    path = tmp_path / f"{program}-all-sections.yaml"
    path.write_text(yaml.safe_dump(sections), encoding="utf-8")
    loader = _yaml_project_loader(
        importlib.import_module(f"chemsmart.settings.{program}")
    )
    project = loader.from_yaml(str(path))
    fields = set()
    for jobtype in jobtypes:
        settings = getattr(project, f"{jobtype}_settings")()
        fields.update(
            name for name in vars(settings) if not name.startswith("_")
        )
    assert fields, f"the {program} loader named no settings field"
    return fields


def _commands(command, ctx, path):
    yield path, command
    if isinstance(command, click.Group):
        for name in command.list_commands(ctx):
            child = command.get_command(ctx, name)
            if child is not None:
                yield from _commands(
                    child, click.Context(child, parent=ctx), (*path, name)
                )


@pytest.mark.capability(
    "setting:gaussian:*", "setting:orca:*", "setting:pyscf:*", "setting:xtb:*"
)
def test_no_command_option_named_for_a_project_setting_has_a_default(
    tmp_path,
):
    """An option a person did not type leaves the project's value alone.

    The contract every program command states -- "use the project value
    unless the command line gives one" -- holds only when "not given" is
    None. Any other default is indistinguishable from a typed value, so
    it wins over the project on every run that did not type it.
    """

    from chemsmart.cli.run import run

    root = click.Context(run)
    offenders = []
    for program in sorted(AGENT_PROGRAM_JOBTYPES):
        fields = _project_fields(program, tmp_path)
        group = run.get_command(root, program)
        for path, command in _commands(
            group, click.Context(group, parent=root), (program,)
        ):
            for option in command.params:
                if not isinstance(option, click.Option):
                    continue
                if option.name in fields and option.default is not None:
                    offenders.append(
                        f"run {' '.join(path)} --{option.name.replace('_', '-')}"
                        f" defaults to {option.default!r}"
                    )
    assert not offenders, (
        "command options that replace a project's setting when not typed:\n"
        + "\n".join(offenders)
    )


# ----------------------------------------------------------------------
# The census: every setting, every executable stage, through run --fake
# ----------------------------------------------------------------------

_WATER = (
    "3\nwater\nO 0.0 0.0 0.1173\nH 0.0 0.7572 -0.4692\nH 0.0 -0.7572 -0.4692\n"
)

#: One stated value per setting, different from what the stage writes
#: without it, and the section values that make it admissible beside it
#: (None removes a base value). This is the program's own grammar -- the
#: oracle the writer is held to -- so it is written here rather than read
#: from the writer's tables. Every setting the capability advertises must
#: appear: a setting added to the capability without a row fails
#: ``test_every_advertised_setting_has_a_stated_value``.
_STATED = {
    "gaussian": {
        "ab_initio": ("mp2", {"functional": None}),
        "additional_opt_options_in_route": ("maxstep=5", {}),
        "additional_route_parameters": ("nosymm", {}),
        "additional_solvent_options": (
            "iterative",
            {"solvent_model": "smd", "solvent_id": "water"},
        ),
        "basis": ("def2tzvp", {}),
        "broken_symmetry": (True, {}),
        "custom_solvent": (
            "eps=4.0\nepsinf=2.0",
            {"solvent_model": "pcm", "solvent_id": "generic,read"},
        ),
        "defgrid": ("superfinegrid", {}),
        "dieze_tag": ("#p", {}),
        "direction": ("reverse", {}),
        "dispersion": ("gd3bj", {}),
        "eqsolv": ("eqsolv", {"solvent_model": "smd", "solvent_id": "water"}),
        "flat_irc": (True, {}),
        "forces": (True, {"freq": False}),
        "freq": ("TOGGLE", {}),
        "functional": ("pbe0", {}),
        "geom_maxiter": (7, {}),
        "guess": ("mix,always", {}),
        "heavy_elements_basis": ("def2tzvp", {}),
        "jobtype": ("ts", {}),
        "link_route": ("opt=tight freq", {}),
        "maxcycles": (31, {}),
        "maxpoints": (23, {}),
        "nstates": (7, {}),
        "numfreq": (True, {"freq": False}),
        "predictor": ("lqa", {"recorrect": "never"}),
        "recalc_step": (4, {}),
        "recorrect": ("never", {"predictor": "hpc"}),
        "response_method": ("tda", {}),
        "root": (2, {}),
        "scf_convergence": ("tight", {}),
        "semiempirical": ("pm6", {"functional": None, "basis": None}),
        "solvent_id": ("toluene", {"solvent_model": "smd"}),
        "solvent_model": ("smd", {"solvent_id": "water"}),
        "stable": ("opt,qrhf", {}),
        "state_manifold": ("triplet", {}),
        "states": ("triplets", {}),
        "stepsize": (7, {}),
    },
    "orca": {
        "ab_initio": ("mp2", {"functional": None}),
        "additional_route_parameters": ("NoUseSym", {}),
        "additional_solvent_options": (
            "epsilon 80.0",
            {"solvent_model": "cpcm", "solvent_id": "water"},
        ),
        "aux_basis": ("def2/j", {}),
        "basis": ("def2-tzvp", {}),
        "broken_symmetry": (True, {}),
        "custom_solvent": (
            "epsilon 80.0\nrefrac 1.33",
            {"solvent_model": "cpcm"},
        ),
        "defgrid": ("defgrid3", {}),
        "dipole": (True, {}),
        "direction": ("forward", {}),
        "dispersion": ("d3bj", {}),
        "extrapolation_basis": ("Extrapolate(2/3,def2)", {}),
        "forces": (True, {}),
        "freq": (True, {}),
        "frozen_core": ("fc_none", {"functional": None, "ab_initio": "mp2"}),
        "frozen_core_electrons": (
            2,
            {
                "functional": None,
                "ab_initio": "mp2",
                "frozen_core": "fc_electrons",
            },
        ),
        "full_scan": (True, {}),
        "functional": ("pbe0", {}),
        "gbw": (True, {}),
        "geom_maxiter": (7, {}),
        "heavy_elements": (["O"], {"heavy_elements_basis": "def2-tzvp"}),
        "heavy_elements_basis": ("def2-tzvp", {"heavy_elements": ["O"]}),
        "hessmode": (1, {}),
        "inithess": ("calc_numfreq", {}),
        "joboption": ("NEB-TS", {}),
        "light_elements_basis": ("def2-svp", {}),
        "mdci_cutoff": (
            "tight",
            {
                "functional": None,
                "ab_initio": "dlpno-ccsd(t)",
                "aux_basis": "def2-svp/c",
            },
        ),
        "mdci_density": (
            "unrelaxed",
            {
                "functional": None,
                "ab_initio": "dlpno-ccsd(t)",
                "aux_basis": "def2-svp/c",
            },
        ),
        "nimages": (6, {}),
        "nstates": (7, {}),
        "numfreq": (True, {}),
        "numhess": (True, {}),
        "opt_convergence": ("tight", {}),
        "preopt_ends": (True, {}),
        "quadrupole": (True, {}),
        "recalc_hess": (3, {}),
        "reference": ("uhf", {}),
        "relativistic": ("zora", {"basis": "zora-def2-svp"}),
        "response_method": ("tda", {}),
        "ri_approximation": ("rijcosx", {"aux_basis": "def2/j"}),
        "scf_algorithm": ("kdiis", {}),
        "scf_convergence": ("verytight", {}),
        "scf_maxiter": (77, {}),
        "semiempirical": ("pm3", {"functional": None, "basis": None}),
        # On the census's singlet water a flip has no high-spin state to
        # start from, so it ends refused at compile with its sentence.
        "site_spin_flip": ({"atoms": [1], "final_ms": 0}, {}),
        "solvent_id": ("toluene", {"solvent_model": "cpcm"}),
        "solvent_model": ("smd", {"solvent_id": "water"}),
        "solventfilename": ("census.cosmorsxyz", {"solvent_model": "cosmors"}),
        "state_manifold": ("singlet_triplet", {}),
        "trust_radius": (0.2, {}),
        "tssearch_type": (
            "scants",
            {
                "scants_modred": {
                    "coords": [[1, 2]],
                    "dist_start": 0.9,
                    "dist_end": 1.1,
                    "num_steps": 3,
                }
            },
        ),
        "vpt2": (True, {}),
        "vpt2_anharmonic_displacement": (0.4, {"vpt2": True}),
        "vpt2_hessian_cutoff": (50.0, {"vpt2": True}),
        # Not advertised; each was dropped by the command's default (R10
        # Q31), so they stay in the census a person's project goes through.
        "adapt_scale_displ": (False, {}),
        "do_sd_corr": (False, {}),
        "interpolate_only": (False, {}),
        "monitor_internals": (True, {"internal_modred": [[1, 2]]}),
        "sd_corr_parabolicfit": (False, {}),
        "sd_parabolicfit": (False, {}),
    },
    "pyscf": {
        "ab_initio": ("mp2", {"functional": None}),
        "aux_basis": ("def2-universal-jkfit", {"density_fit": True}),
        "basis": ("def2-tzvp", {}),
        "broken_symmetry": (True, {}),
        "cc_max_cycle": (77, {"functional": None, "ab_initio": "ccsd"}),
        "defgrid": ("defgrid3", {}),
        "density_fit": (True, {}),
        "dispersion": ("d3bj", {}),
        "excited_state_root": (
            2,
            {
                "response_method": "tddft",
                "nstates": 3,
                "state_manifold": "singlet",
            },
        ),
        "fd_step_angstrom": (
            0.01,
            {"hessian_derivative": "finite_difference"},
        ),
        "freq": (True, {}),
        "frozen_core": ("auto", {"functional": None, "ab_initio": "mp2"}),
        "functional": ("pbe0", {}),
        "hessian_derivative": ("finite_difference", {}),
        "irc_direction": ("backward", {}),
        "nstates": (5, {}),
        "opt_maxsteps": (33, {}),
        "opt_solver": ("berny", {}),
        "response_method": ("tda", {}),
        "scf_maxiter": (77, {}),
        "scf_stability": (True, {}),
        "scf_tol": (1e-10, {}),
        "solvent_id": ("toluene", {"solvent_model": "pcm"}),
        "solvent_model": ("pcm", {"solvent_id": "water"}),
        "state_manifold": ("triplet", {}),
        "td_max_cycle": (55, {}),
    },
    "xtb": {
        "charge": (1, {"multiplicity": 2}),
        "gfn_version": ("gfn1", {}),
        "grad": (True, {}),
        "jobtype": ("opt", {}),
        "multiplicity": (3, {}),
        "optimization_level": ("loose", {}),
        "solvent_id": ("toluene", {"solvent_model": "alpb"}),
        "solvent_model": ("alpb", {"solvent_id": "water"}),
    },
}

#: The stages whose project is not a level in the phase (or own) section.
_BASE = {
    ("gaussian", "td"): {"td": dict(_LEVEL["gaussian"])},
    ("orca", "td"): {
        "td": {
            **_LEVEL["orca"],
            "response_method": "tddft",
            "nstates": 3,
            "state_manifold": "singlet",
        }
    },
    ("pyscf", "irc"): {"irc": {**_LEVEL["pyscf"], "irc_direction": "forward"}},
    ("pyscf", "td"): {
        "td": {
            **_LEVEL["pyscf"],
            "response_method": "tddft",
            "nstates": 3,
            "state_manifold": "singlet",
        }
    },
}

#: The driven or held coordinate a stage cannot be compiled without.
_JOB_ARGUMENTS = {
    ("gaussian", "scan"): (
        "--coordinates",
        "[[1,2]]",
        "--step-size",
        "0.05",
        "--num-steps",
        "3",
    ),
    ("gaussian", "modred"): ("--coordinates", "[[1,2]]"),
    ("orca", "scan"): (
        "--coordinates",
        "[[1,2]]",
        "--dist-start",
        "0.95",
        "--dist-end",
        "1.05",
        "--num-steps",
        "3",
    ),
    ("orca", "modred"): ("--coordinates", "[[1,2]]"),
}


def _executable_stages():
    from chemsmart.settings.capabilities import PROGRAM_CAPABILITIES

    return sorted(
        (program, item.jobtype)
        for program in AGENT_PROGRAM_JOBTYPES
        for item in PROGRAM_CAPABILITIES[
            program
        ].resolved_engine_job_capabilities
        if item.execution_supported and item.engine == "cpu"
    )


def _base(program, jobtype):
    if (program, jobtype) in _BASE:
        return _BASE[(program, jobtype)]
    if program in {"gaussian", "orca"}:
        return {"gas": dict(_LEVEL[program])}
    return {jobtype: dict(_LEVEL[program])}


def _stated_sections(program, jobtype, field, value, extra):
    sections = yaml.safe_load(yaml.safe_dump(_base(program, jobtype)))
    section = sections.setdefault(jobtype, {})
    for key, item in {**extra, field: value}.items():
        section[key] = item
    return sections


def _stage_fields(program, jobtype, tmp_path):
    """The fields the loader's settings object for *jobtype* carries."""

    import importlib

    from chemsmart.agent.projects import _yaml_project_loader

    path = tmp_path / f"{program}-{jobtype}-stage.yaml"
    path.write_text(yaml.safe_dump(_base(program, jobtype)), encoding="utf-8")
    loader = _yaml_project_loader(
        importlib.import_module(f"chemsmart.settings.{program}")
    )
    project = loader.from_yaml(str(path))
    return set(vars(getattr(project, f"{jobtype}_settings")()))


def _written_input(program, workspace):
    """What the program reads: the native input, or its canonical stand-in.

    Gaussian and ORCA read their input file; PySCF reads the CONFIG its
    driver carries (run identifiers and the project digest removed); xTB
    reads its command line (the geometry path removed).
    """

    import json
    import re

    written = {}
    for path in sorted(workspace.rglob("*")):
        if not path.is_file():
            continue
        if path.suffix in {".com", ".gjf", ".inp"}:
            written[path.name] = path.read_text(errors="replace")
        elif program == "pyscf" and path.suffix == ".py":
            match = re.search(
                r'CONFIG = json.loads\("""(.*?)"""\)', path.read_text(), re.S
            )
            if match:
                config = json.loads(match.group(1))
                for volatile in ("run_id", "run_nonce", "project_yaml_digest"):
                    config.pop(volatile, None)
                written[path.name] = json.dumps(config, sort_keys=True)
        elif path.name == "xtb-preview-receipt-v1.json":
            argv = json.loads(path.read_text()).get("canonical_argv") or []
            written["xtb-argv"] = " ".join(argv[2:])
    return written


def _stage_run(program, jobtype, tmp_path, sections):
    """Validate as the session does, then write through ``run --fake``.

    Returns ``(validation, result, written, receipt)``: ``receipt`` is the
    preview verifier's over the workspace the command wrote, and
    ``written`` and ``receipt`` are None when it wrote nothing.
    """

    from pathlib import Path

    from click.testing import CliRunner

    from chemsmart.agent.live_session import _preview_server_profile
    from chemsmart.agent.program_verifiers import (
        build_preview_expectation,
        validate_preview_workspace,
    )
    from chemsmart.cli.main import entry_point

    from .gaussian_fake_preview import artifact, validate

    tmp_path.mkdir(parents=True, exist_ok=True)
    xyz = tmp_path / "input.xyz"
    xyz.write_text(_WATER, encoding="utf-8")
    project, validation = validate(tmp_path, program, sections, jobtype)
    if validation.status != "valid":
        return validation, None, None, None
    applied = dict(validation.settings)
    charge = int(applied.get("charge") or 0) if program == "xtb" else 0
    multiplicity = (
        int(applied.get("multiplicity") or 1) if program == "xtb" else 1
    )
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
        program,
        "--project",
        str(project),
        "--filename",
        str(xyz),
        "--charge",
        str(charge),
        "--multiplicity",
        str(multiplicity),
        jobtype,
        *_JOB_ARGUMENTS.get((program, jobtype), ()),
    ]
    runner = CliRunner()
    with runner.isolated_filesystem(temp_dir=workspace) as cwd:
        result = runner.invoke(entry_point, argv)
        if result.exit_code != 0:
            return validation, result, None, None
        expectation = build_preview_expectation(
            program=program,
            jobtype=jobtype,
            input_artifact=artifact(xyz, "geometry_xyz"),
            project=validation,
            charge=charge,
            multiplicity=multiplicity,
        )
        receipt = validate_preview_workspace(expectation, Path(cwd))
        written = _written_input(program, Path(cwd))
    return validation, result, written, receipt


@pytest.mark.capability(
    "setting:gaussian:*", "setting:orca:*", "setting:pyscf:*", "setting:xtb:*"
)
@pytest.mark.parametrize("program", sorted(_STATED))
def test_every_advertised_setting_has_a_stated_value(program):
    """A setting the capability offers is in the census before it ships."""

    from chemsmart.settings.capabilities import PROGRAM_CAPABILITIES

    advertised = set(PROGRAM_CAPABILITIES[program].project_owned_parameters)
    missing = sorted(advertised - set(_STATED[program]))
    assert not missing, (
        f"{program} advertises settings this census never states: "
        f"{missing}; add a row to _STATED with a value the program reads "
        "differently from its default"
    )


@pytest.mark.capability(
    "setting:gaussian:*", "setting:orca:*", "setting:pyscf:*", "setting:xtb:*"
)
@pytest.mark.parametrize(("program", "jobtype"), _executable_stages())
def test_every_stated_setting_reaches_the_input_or_says_why_not(
    tmp_path, program, jobtype
):
    """Stated -> validated -> written -> read back, for one stage.

    Each census setting the stage's settings class carries is stated alone
    at a value that is not its default, through the session's own
    validation and the public ``run --fake``. It must end one of four ways:
    - refused, at validation or at compile, with a sentence (never a crash);
    - reported overridden: the loader applied another value, which the
      session is told (``declared_settings_overridden``);
    - written for no stage of this job type by the writer's own table
      (``settings_not_written_for``), and so not demanded of its input;
    - or reaching the input: the preview verifier, reading the written
      input back, is green -- and an input identical to the one written
      without the setting passes only if the verifier compared it.
    A setting validated, never written and compared by nothing is the class
    R10 Q31's census found open in five shapes.
    """

    import importlib

    from chemsmart.agent._contracts import canonical_data

    module = importlib.import_module(f"chemsmart.jobs.{program}.settings")
    ask = getattr(module, "settings_not_written_for", lambda _jobtype: ())
    not_written = set(ask(jobtype))
    stage_fields = _stage_fields(program, jobtype, tmp_path)
    base_validation, base_result, baseline, base_receipt = _stage_run(
        program, jobtype, tmp_path / "baseline", _base(program, jobtype)
    )
    assert baseline is not None, (base_validation.diagnostic, base_result)
    assert base_receipt.status == "valid", base_receipt.findings

    failures = []
    for field, (value, extra) in sorted(_STATED[program].items()):
        if field not in stage_fields:
            continue
        if value == "TOGGLE":
            value = not bool(dict(base_validation.settings).get(field))
        sections = _stated_sections(program, jobtype, field, value, extra)
        validation, result, written, receipt = _stage_run(
            program, jobtype, tmp_path / field, sections
        )
        if validation.status != "valid":
            if not validation.diagnostic:
                failures.append(f"{field}: refused without a sentence")
            continue
        if written is None:
            exception = result.exception
            if not isinstance(exception, (ValueError, SystemExit)):
                failures.append(
                    f"{field}={value!r}: validated, then the command "
                    f"crashed: {type(exception).__name__}: {exception}"
                )
            continue
        applied = dict(validation.settings)
        overridden = field in applied and canonical_data(
            applied[field]
        ) != canonical_data(value)
        if overridden:
            continue
        if receipt.status != "valid":
            findings = [
                (item.field, item.expected, item.observed)
                for item in receipt.findings
            ]
            failures.append(
                f"{field}={value!r}: the preview is red on the input it "
                f"wrote: {findings}"
            )
            continue
        if (
            written == baseline
            and field not in applied
            and field not in not_written
        ):
            failures.append(
                f"{field}={value!r}: validated, never written, and "
                "compared by nothing"
            )
    assert not failures, f"{program} {jobtype}:\n" + "\n".join(failures)
