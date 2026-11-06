"""Plain test helpers: stub/fake classes and builder functions that need
no pytest fixtures. Import them directly, e.g.
``from tests.helpers import ConcreteJobRunner``.

Fixtures live in ``tests/fixtures/`` (loaded by ``tests/conftest.py``).
"""

import importlib
import json
import logging
import os
import shutil
from pathlib import Path
from textwrap import dedent
from unittest.mock import MagicMock, patch

import click
import numpy as np
import pytest
import yaml
from click.testing import CliRunner
from rdkit import Chem
from rdkit.Chem import AllChem
from requests.exceptions import HTTPError as RequestsHTTPError

from chemsmart.cli.gaussian.gaussian import gaussian
from chemsmart.cli.job import click_folder_options
from chemsmart.cli.nciplot.nciplot import nciplot
from chemsmart.cli.run import run
from chemsmart.cli.thermochemistry.thermochemistry import thermochemistry
from chemsmart.io.molecules.structure import Molecule
from chemsmart.io.xtb.input import XTBInput
from chemsmart.jobs.crest.settings import CRESTJobSettings
from chemsmart.jobs.gaussian.settings import (
    GaussianJobSettings,
    GaussianpKaJobSettings,
)
from chemsmart.jobs.grouper.base import MatrixGrouper, MoleculeGrouper
from chemsmart.jobs.job import Job
from chemsmart.jobs.runner import JobRunner
from chemsmart.jobs.xtb.settings import XTBJobSettings
from chemsmart.settings.gaussian import GaussianProjectSettings
from chemsmart.settings.xtb import XTBProjectSettings
from chemsmart.utils.mixins import RegistryMixin

_UNSET = object()
mol_cli_module = importlib.import_module("chemsmart.cli.mol.mol")
# ``chemsmart.cli.nciplot`` re-exports the ``nciplot`` Click command under
# the submodule's own name, so look the submodule up explicitly to patch
# ``NCIPLOTJob`` on it.
nciplot_cli_module = importlib.import_module("chemsmart.cli.nciplot.nciplot")
thermochemistry_cli_module = importlib.import_module(
    "chemsmart.cli.thermochemistry.thermochemistry"
)


class FakePKaSettings(GaussianpKaJobSettings):
    """Settings double with deterministic sub-job wiring.

    Overrides the pair/settings-generation methods to return plain
    sentinel strings instead of doing real molecule-copying or
    basis/functional resolution, since ``GaussianOptJob`` and
    ``GaussianSinglePointJob`` are mocked out in these tests and never
    actually inspect these values.
    """

    def conjugate_pair_molecules(self, molecule):
        return "PROT_MOL", "CONJ_MOL"

    def _create_gas_phase_job_settings(self, molecule):
        return "PROT_OPT_SETTINGS", "CONJ_OPT_SETTINGS"

    def reference_pair_molecules(self):
        return "REF_ACID_MOL", "REF_CONJ_MOL"

    def reference_pair_job_settings(self):
        return "REF_ACID_OPT_SETTINGS", "REF_CONJ_OPT_SETTINGS"

    def _create_solution_phase_sp_settings(self, molecule):
        return "PROT_SP_SETTINGS", "CONJ_SP_SETTINGS"

    def reference_pair_sp_job_settings(self):
        return "REF_ACID_SP_SETTINGS", "REF_CONJ_SP_SETTINGS"


def make_pka_sub_job_factory():
    """Return a side_effect factory for GaussianOptJob/GaussianSinglePointJob.

    Records every created mock job in ``factory.created`` (in creation
    order) so tests can reach in and configure ``is_complete``/
    ``_output``/``outputfile`` per job.
    """

    created = []

    def factory(*, molecule, settings, label, jobrunner, skip_completed):
        job = MagicMock()
        job.molecule = molecule
        job.settings = settings
        job.label = label
        job.jobrunner = jobrunner
        job.skip_completed = skip_completed
        job.is_complete.return_value = False
        job._output.return_value = None
        job.outputfile = f"{label}.log"

        # By default, simulate a successful run: once `.run()` is
        # called the job reports itself complete. Tests that need to
        # simulate a job staying incomplete after running should clear
        # this with `job.run.side_effect = None`.
        def _mark_complete_on_run():
            job.is_complete.return_value = True

        job.run.side_effect = _mark_complete_on_run
        created.append(job)
        return job

    factory.created = created
    return factory


class DummyJob:
    """Minimal stand-in for a Job, exposing only the attributes the job
    runners read. File paths are derived from the *current* label, so
    runners that rename the job (e.g. appending ``_fake``) see the
    updated paths."""

    execution_name = Job.execution_name

    def __init__(self, program, input_ext, folder, label, settings, molecule):
        self.PROGRAM = program
        self._input_ext = input_ext
        self.folder = str(folder)
        self.label = label
        self.settings = settings
        self.molecule = molecule if molecule is not None else MagicMock()
        # InputWriter.__init__ reads job.jobrunner; runners never do.
        self.jobrunner = None

    def _path(self, ext):
        return str(Path(self.folder) / f"{self.label}{ext}")

    @property
    def inputfile(self):
        return self._path(self._input_ext)

    @property
    def outputfile(self):
        return self._path(".out")

    @property
    def errfile(self):
        return self._path(".err")

    @property
    def xyzfile(self):
        return self._path(".xyz")

    @property
    def chkfile(self):
        return self._path(".chk")

    @property
    def gbwfile(self):
        return self._path(".gbw")

    def is_complete(self):
        return False


# program key -> (Job.PROGRAM, input file extension, default settings)
DUMMY_JOB_SPECS = {
    "gaussian": ("Gaussian", ".com", lambda: None),
    "orca": ("ORCA", ".inp", lambda: None),
    "crest": (
        "crest",
        ".inp",
        lambda: CRESTJobSettings(charge=0, multiplicity=1),
    ),
    "xtb": (
        "xtb",
        ".inp",
        lambda: XTBJobSettings(charge=0, multiplicity=1),
    ),
}


class ConcreteJobRunner(JobRunner):
    """Minimal concrete :class:`JobRunner` for exercising base-class
    behavior. Defining it registers it as the runner for the
    ``concrete_runner_type`` jobtype (used by ``JobRunner.from_job`` tests).
    """

    JOBTYPES = ["concrete_runner_type"]
    PROGRAM = "concrete"
    FAKE = False
    SCRATCH = False

    @property
    def executable(self):
        return None

    def _get_command(self, job):
        return "echo hi"

    def _create_process(self, job, command, env):
        return MagicMock()


class ConcreteMatrixGrouper(MatrixGrouper):
    """Matrix-capable concrete grouper so ``center``/``top3`` are allowed.

    ``group()`` just returns singleton groups; tests use it to exercise
    the shared ``MatrixGrouper``/``MoleculeGrouper`` helpers directly.
    """

    def group(self):
        return self._build_singleton_groups()


class ConcreteMoleculeGrouper(MoleculeGrouper):
    """Minimal concrete subclass so the abstract base can be instantiated.

    ``group()`` returns every molecule in a single group.
    """

    def group(self):
        indices = list(range(len(self.molecules)))
        return [self.molecules], [indices]


def cleanup_pymol_grouper(grouper):
    """Neutralise a ``PymolRMSDGrouper`` so ``__del__`` can't call a real
    ``cmd.quit()``: PyMOL's embedded session is process-wide, and quitting
    it hangs every later PyMOL test."""
    if getattr(grouper, "_temp_dir", None):
        shutil.rmtree(grouper._temp_dir, ignore_errors=True)
        grouper._temp_dir = None
    grouper.cmd = None


class BaseRegistry(RegistryMixin):
    pass


class SubRegistry1(BaseRegistry):
    pass


class SubRegistry2(BaseRegistry):
    pass


class FakeIterateQueue:
    """Stand-in for a ``multiprocessing.Manager().Queue()``.

    ``get_nowait_plan`` is a list consumed one entry per call to
    ``get_nowait``: ``None`` means "return the next queued item (or
    raise ``queue.Empty`` if none)", and any exception instance/class
    means "raise this instead". Once the plan is exhausted, falls back
    to the "return next item or raise Empty" behavior.
    """

    def __init__(self, get_nowait_plan=None):
        self._items = []
        self._plan = list(get_nowait_plan or [])

    def put(self, item):
        self._items.append(item)

    def get_nowait(self):
        import queue

        if self._plan:
            effect = self._plan.pop(0)
            if effect is not None:
                raise effect
        if self._items:
            return self._items.pop(0)
        raise queue.Empty()


class FakeIterateProcess:
    """Stand-in for a ``multiprocessing.Process``.

    ``alive`` controls every ``is_alive()`` call for this process
    (kept simple: a process is either always "still running" or
    already "finished", which is all the watchdog logic needs).
    """

    def __init__(self, target=None, args=(), daemon=False, alive=False):
        self.target = target
        self.args = args
        self.daemon = daemon
        self.pid = 4242
        self._alive = alive
        self.terminated = False
        self.killed = False
        self.joined_timeouts = []

    def start(self):
        pass

    def is_alive(self):
        return self._alive

    def join(self, timeout=None):
        self.joined_timeouts.append(timeout)

    def terminate(self):
        self.terminated = True

    def kill(self):
        self.killed = True


class ConcreteJob(Job):
    """Minimal concrete subclass for exercising base ``Job`` behavior.

    Overrides every hook, with ``complete``/``output`` controlling
    ``_job_is_complete()``/``_output()`` and ``ran`` recording ``_run()``.
    """

    TYPE = "concrete"

    def __init__(self, complete=False, output=None, **kwargs):
        super().__init__(**kwargs)
        self._complete = complete
        self._output_obj = output
        self.ran = False

    def _run(self, **kwargs):
        self.ran = True

    def _job_is_complete(self):
        return self._complete

    def _output(self):
        return self._output_obj

    def _backup_files(self, **kwargs):
        pass


class DefaultBehaviorJob(Job):
    """Job subclass that keeps the base class's default ``_run``,
    ``_backup_files`` and ``_job_is_complete``, overriding only ``_output``
    (returns ``output``)."""

    TYPE = "default_behavior"

    def __init__(self, output=None, **kwargs):
        super().__init__(**kwargs)
        self._output_obj = output

    def _run(self, **kwargs):
        return super()._run(**kwargs)

    def _output(self):
        return self._output_obj


# ── Factories and CLI-invocation helpers (formerly fixtures) ──
def make_cli_ctx_obj(jobrunner, **extra):
    """Factory for the minimal Click context object.

    Extra keyword arguments (e.g. ``molecules``, ``label``) are merged
    into the returned dict.
    """
    return {"jobrunner": jobrunner, **extra}


def make_cli_ctx(molecules, molecule_indices, label, jobrunner=None):
    """Factory for a lightweight stand-in for a Click context.

    Returns an object exposing only ``.obj``, which is sufficient for
    helpers such as ``chemsmart.cli.utils.build_jobs`` that read state
    from ``ctx.obj``. A ``MagicMock`` jobrunner is used when none is given.
    """
    from types import SimpleNamespace

    if jobrunner is None:
        jobrunner = MagicMock(name="jobrunner")
    return SimpleNamespace(
        obj=make_cli_ctx_obj(
            jobrunner,
            molecules=molecules,
            molecule_indices=molecule_indices,
            label=label,
        )
    )


def invoke_config_server(args=None):
    """Invoke ``chemsmart config server`` via Click's CliRunner::

    result = invoke_config_server()
    assert result.exit_code == 0
    """
    from chemsmart.cli.config import config

    runner = CliRunner()
    return runner.invoke(config, ["server"] + (args or []))


def write_server_yaml(config_root, name, content):
    """Write a server YAML file at ``<config_root>/server/<name>``.

    The content is dedented and left-stripped before being written, so
    callers can pass indented triple-quoted strings::

        path = write_server_yaml(
            tmp_path, "A.yaml", "SERVER:\n    USER_VALUE: preserved\n"
        )
    """
    server_dir = Path(config_root) / "server"
    server_dir.mkdir(parents=True, exist_ok=True)
    path = server_dir / name
    path.write_text(dedent(content).lstrip(), encoding="utf-8")
    return path


def invoke_update_configs(config_root, args=None):
    """Invoke ``chemsmart update configs`` against ``config_root`` via
    Click's CliRunner::

        result = invoke_update_configs(tmp_path, ["-s", "A"])
        assert result.exit_code == 0
    """
    from chemsmart.cli.update import update

    return CliRunner().invoke(
        update,
        ["configs"] + (args or []),
        env={"CHEMSMART_CONFIG_DIR": str(config_root)},
    )


def read_yaml_file(path):
    """Read and parse the YAML file at ``path``."""
    return yaml.safe_load(Path(path).read_text(encoding="utf-8"))


def write_text_file(path, content):
    """Write ``content`` to ``path`` (creating parent directories as needed)
    and return the path as a :class:`~pathlib.Path`::

        path = write_text_file(tmp_path / "a" / "b.yaml", "content\n")
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    return path


def invoke_update_projects(config_root):
    """Invoke ``chemsmart update projects`` against ``config_root`` via
    Click's CliRunner::

        result = invoke_update_projects(tmp_path)
        assert result.exit_code == 0
    """
    from chemsmart.cli.update import update

    return CliRunner().invoke(
        update,
        ["projects"],
        env={"CHEMSMART_CONFIG_DIR": str(config_root)},
    )


@click.command()
@click_folder_options
@click.pass_context
def _invoke_folder_command_folder_cmd(ctx, directory, filetype, program):
    """Minimal command used only to inspect registered options."""
    ctx.ensure_object(dict)
    ctx.obj["directory"] = directory
    ctx.obj["filetype"] = filetype
    ctx.obj["program"] = program


def invoke_folder_command(args=None):
    """Invoke a minimal Click command decorated with ``click_folder_options``
    with ``args``; the parsed options are left in the result's context
    ``obj`` (``directory``, ``filetype``, ``program``)."""
    runner = CliRunner()
    result = runner.invoke(_invoke_folder_command_folder_cmd, args or [])
    return result


def run_thermochemistry_and_capture_settings(
    extra_args=None,
    ctx_obj=None,
    filename="dummy.log",
    detected_program="gaussian",
):
    """Run the thermochemistry CLI with mocked job construction."""
    runner = CliRunner()
    captured_settings = None
    mock_job = MagicMock()

    base_args = ["-f", filename, "-T", "298.15"]
    cli_args = base_args + (extra_args or [])

    with (
        patch.object(
            thermochemistry_cli_module,
            "get_program_type_from_file",
            return_value=detected_program,
        ),
        patch.object(
            thermochemistry_cli_module.ThermochemistryJob,
            "from_filename",
            return_value=mock_job,
        ) as mock_from_filename,
    ):
        result = runner.invoke(
            thermochemistry,
            cli_args,
            obj=ctx_obj or {},
            catch_exceptions=False,
        )
        if mock_from_filename.call_args is not None:
            captured_settings = mock_from_filename.call_args[1].get("settings")

    return result, captured_settings


def run_thermochemistry_with_directory(extra_args, mock_files=None):
    """Invoke the thermochemistry CLI with directory options and a mocked
    folder.

    Patches ``BaseFolder`` so that
    ``get_all_output_files_in_current_folder_by_program``,
    ``get_all_files_in_current_folder_by_suffix``, and
    ``get_all_files_in_current_folder_by_program_and_suffix`` all return
    ``mock_files``.  Also patches ``ThermochemistryJob.from_filename`` to
    avoid real job execution::

        result, mock_from_filename = run_thermochemistry_with_directory(
            ["-d", str(tmp_path), "-p", "gaussian", "-T", "298.15"],
            mock_files=["/fake/a.log", "/fake/b.log"],
        )
        assert result.exit_code == 0
        assert mock_from_filename.call_count == 2
    """
    if mock_files is None:
        mock_files = ["/fake/dir/mol1.log"]

    runner = CliRunner()
    mock_job = MagicMock()
    mock_job.label = "mol1"

    with (
        patch.object(
            thermochemistry_cli_module, "BaseFolder"
        ) as mock_folder_cls,
        patch.object(
            thermochemistry_cli_module.ThermochemistryJob, "from_filename"
        ) as mock_from_filename,
    ):
        mock_folder = MagicMock()
        # Configure every discovery method to return the caller-supplied list
        mock_folder.get_all_output_files_in_current_folder_by_program.return_value = (
            mock_files
        )
        mock_folder.get_all_output_files_in_current_folder_and_subfolders_by_program.return_value = (
            mock_files
        )
        mock_folder.get_all_files_in_current_folder_by_suffix.return_value = (
            mock_files
        )
        mock_folder.get_all_files_in_current_folder_by_program_and_suffix.return_value = (
            mock_files
        )
        mock_folder_cls.return_value = mock_folder
        mock_from_filename.return_value = mock_job

        result = runner.invoke(thermochemistry, extra_args)
        return result, mock_from_filename


def run_gaussian_and_capture_settings(
    job_class_path,
    cli_args,
    ctx_obj,
    *,
    job_cls=None,
    catch_exceptions=False,
):
    """Run the gaussian CLI with a patched job class and capture settings.

    ``job_class_path`` is either a dotted path for :func:`patch` or a
    ``(target, attribute)`` tuple for :func:`patch.object` (needed when
    the dotted path is shadowed, e.g. ``chemsmart.cli.gaussian.traj``).

    Pass ``job_cls`` to supply the mock yourself and inspect its calls
    (labels, molecules, call counts) afterwards. Pass
    ``catch_exceptions=True`` to assert on ``result.exception`` instead
    of letting errors propagate.
    """
    runner = CliRunner()
    captured_settings = None
    if job_cls is None:
        job_cls = MagicMock(name="job_cls")
        job_cls.return_value = MagicMock()

    if isinstance(job_class_path, tuple):
        patcher = patch.object(*job_class_path, new=job_cls)
    else:
        patcher = patch(job_class_path, new=job_cls)

    with patcher as mock_job_cls:
        result = runner.invoke(
            gaussian,
            cli_args,
            obj=ctx_obj,
            catch_exceptions=catch_exceptions,
        )
        if mock_job_cls.call_args is not None:
            captured_settings = mock_job_cls.call_args[1].get("settings")

    return result, captured_settings


def run_gaussian_qmmm_callback(**ctx_obj_overrides):
    """Invoke the ``opt qmmm`` callback directly and return the job mock.

    Bypasses the real "opt"/"ts"/... parent commands: every real parent
    unconditionally sets ``ctx.obj["parent_settings"]`` to a real settings
    object (see e.g. ``chemsmart/cli/gaussian/opt.py``), so
    ``parent_settings`` being None, ``jobtype`` being None (no
    ``ctx.parent``), and absent ``parent_jobtype``/``high_level_atoms``
    are branches no genuine CLI invocation reaches. Building ``ctx.obj``
    by hand reaches them directly. Keyword arguments override entries in
    the default ``ctx.obj``.
    """
    opt_mod = importlib.import_module("chemsmart.cli.gaussian.opt")
    qmmm_cmd = opt_mod.opt.commands["qmmm"]

    project_settings = MagicMock()
    project_settings.qmmm_settings.return_value = None
    molecule = MagicMock()
    molecule.frozen_atoms = None

    ctx_obj = {
        "jobrunner": MagicMock(),
        "project_settings": project_settings,
        "label": "test",
        "molecules": [molecule],
    }
    ctx_obj.update(ctx_obj_overrides)

    ctx = click.Context(qmmm_cmd)
    ctx.obj = ctx_obj

    kwargs = dict(
        high_level_functional="b3lyp",
        high_level_basis="6-31g",
        high_level_force_field=None,
        medium_level_functional=None,
        medium_level_basis=None,
        medium_level_force_field=None,
        low_level_functional=None,
        low_level_basis=None,
        low_level_force_field=None,
        charge_total=None,
        mult_total=None,
        charge_intermediate=None,
        mult_intermediate=None,
        charge_high=None,
        mult_high=None,
        high_level_atoms=None,
        medium_level_atoms=None,
        low_level_atoms=None,
        bonded_atoms=None,
        scale_factors=None,
        mm_atom_info=None,
        mm_parameters_file=None,
    )

    with patch("chemsmart.jobs.gaussian.qmmm.GaussianQMMMJob") as mock_job:
        mock_job.return_value = MagicMock()
        with ctx:
            qmmm_cmd.callback(**kwargs)
        return mock_job


def make_gaussian_writer(settings, molecule=None, label="testjob"):
    """Factory for a ``GaussianInputWriter`` around a lightweight fake job.

    Bypasses full Job construction so private writer methods can be
    tested directly. The fake job exposes only ``settings``,
    ``jobrunner`` (``num_cores=4``, ``mem_gb=8``), ``label`` and
    ``molecule``.
    """
    from types import SimpleNamespace

    from chemsmart.jobs.gaussian.writer import GaussianInputWriter

    job = SimpleNamespace(
        settings=settings,
        jobrunner=SimpleNamespace(num_cores=4, mem_gb=8),
        label=label,
        molecule=molecule,
    )
    return GaussianInputWriter(job=job)


def make_orca_writer(settings, molecule=None):
    """Factory for an ``ORCAInputWriter`` around a lightweight fake job.

    Bypasses full Job construction so individual writer methods can be
    tested directly. The fake job exposes ``settings``, ``molecule``
    (default: a ``MagicMock`` with ``frozen_atoms=None``), ``jobrunner``
    (``num_cores=4``, ``mem_gb=8``), ``label`` and ``folder``.
    """
    from types import SimpleNamespace

    from chemsmart.jobs.orca.writer import ORCAInputWriter

    job = SimpleNamespace(
        settings=settings,
        molecule=(
            molecule if molecule is not None else MagicMock(frozen_atoms=None)
        ),
        jobrunner=MagicMock(num_cores=4, mem_gb=8),
        label="test_job",
        folder="/tmp",
    )
    return ORCAInputWriter(job=job)


def make_uvvis_job(settings=None):
    """Factory for a ``GaussianUVVISJob`` built without ``__init__``.

    ``GaussianUVVISJob.__init__`` always raises (BUGS_FOUND.md #11), so
    the instance is created with ``__new__`` and its attributes are set
    by hand.
    """
    from chemsmart.jobs.gaussian.uvvis import GaussianUVVISJob

    job = GaussianUVVISJob.__new__(GaussianUVVISJob)
    job.folder = "/tmp/uvvis_job"
    job.molecule = MagicMock(name="molecule")
    job.settings = settings or GaussianJobSettings.default()
    return job


def make_matrix_grouper(molecules, **kwargs):
    """Build a minimal concrete :class:`MatrixGrouper`::

    grouper = make_matrix_grouper(
        methanol_molecules, representative_strategy="center"
    )
    """
    return ConcreteMatrixGrouper(molecules, **kwargs)


def make_molecule_grouper(molecules, **kwargs):
    """Factory for a minimal concrete :class:`MoleculeGrouper`, e.g.
    ``make_molecule_grouper(molecules=[water_molecule], num_procs=2)``.
    """
    return ConcreteMoleculeGrouper(molecules=molecules, **kwargs)


def make_formaldehyde_molecule(x_offset=0.0, swap_hydrogens=False):
    """Factory for planar formaldehyde (H₂CO), optionally translated along
    x, or with its two symmetric H atoms swapped (``swap_hydrogens``) to
    force an atom-assignment search.

    One C, one O, two H: exercises the Hungarian assignment's
    single-atom-per-element shortcut (C and O) alongside the real
    ``linear_sum_assignment`` path (the two H atoms).
    """
    positions = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.2, 0.0, 0.0],
            [-0.6, 0.9, 0.0],
            [-0.6, -0.9, 0.0],
        ]
    ) + np.array([x_offset, 0.0, 0.0])
    if swap_hydrogens:
        positions = positions[[0, 1, 3, 2]]
    return Molecule(symbols=["C", "O", "H", "H"], positions=positions)


def make_mock_molecules(n=3):
    """Factory for ``n`` distinct ``MagicMock`` molecules, for job/runner
    tests that never inspect the molecules themselves.
    """
    return [MagicMock(name=f"molecule{i}") for i in range(n)]


def make_thermo_with_output(**output_attrs):
    """Factory for a mock thermochemistry object whose ``file_object``
    carries the given attributes, e.g.
    ``make_thermo_with_output(energies=[1.0, 2.0])``.
    """
    thermo = MagicMock()
    output = MagicMock()
    for key, value in output_attrs.items():
        setattr(output, key, value)
    thermo.file_object = output
    return thermo


def make_scripted_minimize(*script, default=None):
    """Factory for a ``minimize`` ``side_effect`` that follows a per-call
    script, e.g. ``make_scripted_minimize(None, RuntimeError("boom"))``.

    Each positional entry handles one call, in order: ``None`` delegates to
    the real ``scipy.optimize.minimize``, an exception is raised, and any
    other object is returned as the result. Once the script is used up,
    every later call is handled by ``default`` (same rules).
    """
    from scipy.optimize import minimize as real_minimize

    remaining = list(script)

    def _side_effect(*args, **kwargs):
        step = remaining.pop(0) if remaining else default
        if step is None:
            return real_minimize(*args, **kwargs)
        if isinstance(step, BaseException):
            raise step
        return step

    return _side_effect


def make_iterate_combination(**overrides):
    """Factory for an ``IterateCombination`` pairing skeleton ``s`` and
    substituent ``sub`` at link index 1 (label ``s_1_sub_1``); any field
    can be overridden by keyword.
    """
    from chemsmart.jobs.iterate.runner import IterateCombination

    fields = dict(
        skeleton_idx=0,
        skeleton_label="s",
        skeleton_link_index=1,
        skeleton_indices=None,
        substituent_idx=0,
        substituent_label="sub",
        substituent_link_index=1,
    )
    fields.update(overrides)
    return IterateCombination(**fields)


def dummy_job(program, folder, label, settings=None, molecule=None):
    """Factory for lightweight jobs used to unit-test job runners, e.g.
    ``dummy_job("xtb", folder=tmp_path, label="myjob")``.
    """
    program_name, input_ext, default_settings = DUMMY_JOB_SPECS[program]
    return DummyJob(
        program=program_name,
        input_ext=input_ext,
        folder=folder,
        label=label,
        settings=settings if settings is not None else default_settings(),
        molecule=molecule,
    )


def make_concrete_job(label="x", molecule=None, jobrunner=None, **kwargs):
    """Factory for :class:`ConcreteJob` with no molecule or jobrunner, e.g.
    ``make_concrete_job(label="j1", complete=True)``.
    """
    return ConcreteJob(
        molecule=molecule, label=label, jobrunner=jobrunner, **kwargs
    )


def make_default_behavior_job(
    label="myjob", molecule=None, jobrunner=None, **kwargs
):
    """Factory for :class:`DefaultBehaviorJob` with no molecule or
    jobrunner, e.g. ``make_default_behavior_job(output=fake_output)``.
    """
    return DefaultBehaviorJob(
        molecule=molecule, label=label, jobrunner=jobrunner, **kwargs
    )


def make_log_record(msg, args=None):
    """Build an INFO-level :class:`logging.LogRecord` for ``msg``/``args``."""
    return logging.LogRecord(
        name="test",
        level=logging.INFO,
        pathname="",
        lineno=0,
        msg=msg,
        args=args,
        exc_info=None,
    )


def make_mixin_instance(
    mixin_cls, *, contents=_UNSET, route_string=_UNSET, **attrs
):
    """Factory for a stand-alone instance of a mixin class, for testing
    the mixin without a real parser behind it.

    ``contents`` and ``route_string`` (if given, ``None`` included) are
    overridden as fixed read-only properties on a throwaway subclass; all
    other keyword arguments become instance attributes. The mixin's own
    ``__init__`` (if any) is not run::

        stub = make_mixin_instance(
            ORCAFileMixin, contents=["%cpcm", "%end"], filename="cpcm.inp"
        )
    """
    namespace = {}
    if contents is not _UNSET:
        namespace["contents"] = property(lambda self: contents)
    if route_string is not _UNSET:
        namespace["route_string"] = property(lambda self: route_string)
    stub_cls = type(f"_Stub{mixin_cls.__name__}", (mixin_cls,), namespace)
    instance = object.__new__(stub_cls)
    for name, value in attrs.items():
        setattr(instance, name, value)
    return instance


def make_dummy_file(filename, **attrs):
    """Factory for a :class:`FileMixin` reading ``filename`` from disk, with
    two energies, no forces and no input coordinates block (overridable
    by keyword).
    """
    from chemsmart.utils.mixins import FileMixin

    defaults = dict(
        forces=None, energies=[1.0, 2.0], input_coordinates_block=None
    )
    defaults.update(attrs)
    return make_mixin_instance(FileMixin, filename=filename, **defaults)


def make_dummy_gaussian_file(
    route_string="modred", filename="test.gjf", **attrs
):
    """Factory for a :class:`GaussianFileMixin` with ``%chk``/``%mem=32GB``/
    ``%nproc=8`` link-0 lines, a modredundant group mixing a frozen and a
    scan coordinate, and the given ``route_string`` (default ``"modred"``).
    """
    from chemsmart.utils.mixins import GaussianFileMixin

    return make_mixin_instance(
        GaussianFileMixin,
        contents=[
            "%chk=test.chk",
            "%mem=32GB",
            "%nproc=8",
            "#p opt freq",
            "modred",
        ],
        route_string=route_string,
        filename=filename,
        modredundant_group=["B 2 12 F", "B 9 2 S 10 0.05"],
        **attrs,
    )


def make_dummy_folder(folder):
    """Factory for a :class:`FolderMixin` over ``folder``."""
    from chemsmart.utils.mixins import FolderMixin

    return make_mixin_instance(FolderMixin, folder=folder)


def make_dummy_folder_output(main_out=None, secondary_out=None):
    """Factory for a :class:`FolderOutputMixin` delegating, in order, to
    the ``main_out`` and ``secondary_out`` parsers (either may be
    ``None``).
    """
    from chemsmart.utils.mixins import FolderOutputMixin

    return make_mixin_instance(
        FolderOutputMixin,
        FILE_PARSERS=("main_out", "secondary_out"),
        main_out=main_out,
        secondary_out=secondary_out,
    )


def invoke_orca_cli(cli_args, ctx_obj=None, *, catch_exceptions=False):
    """Invoke the ``orca`` CLI group with ``cli_args`` (``ctx_obj``
    defaults to ``{}``) and return the Click result. Pass
    ``catch_exceptions=True`` to assert on ``result.exception`` /
    ``result.output`` instead of letting errors propagate."""
    from chemsmart.cli.orca.orca import orca as orca_cli

    return CliRunner().invoke(
        orca_cli,
        cli_args,
        obj={} if ctx_obj is None else ctx_obj,
        catch_exceptions=catch_exceptions,
    )


def run_orca_and_capture_settings(
    job_class_path,
    cli_args,
    ctx_obj=None,
    *,
    job_cls=None,
    catch_exceptions=False,
    autospec=False,
):
    """Run the orca CLI with a patched job class and capture settings.

    Returns ``(result, settings)``; ``settings`` is ``None`` when the job
    class was never instantiated (e.g. the CLI rejected the arguments).

    Pass ``autospec=True`` to patch with a signature-checked mock, so the
    test fails if the CLI calls the job with arguments its constructor
    doesn't accept. Pass ``job_cls`` to supply the mock yourself and
    inspect its calls (labels, call counts) afterwards. Pass
    ``catch_exceptions=True`` to assert on error exits instead of letting
    exceptions propagate.
    """
    captured_settings = None

    if job_cls is not None:
        patcher = patch(job_class_path, new=job_cls)
    else:
        patcher = patch(job_class_path, autospec=autospec)

    with patcher as mock_job_cls:
        if job_cls is None:
            mock_job_cls.return_value = MagicMock()
        result = invoke_orca_cli(
            cli_args, ctx_obj, catch_exceptions=catch_exceptions
        )
        if mock_job_cls.call_args is not None:
            captured_settings = mock_job_cls.call_args[1].get("settings")

    return result, captured_settings


def run_grouper_and_capture_job_kwargs(
    cli_args, *, job_cls=None, catch_exceptions=False
):
    """Run the grouper CLI with a patched ``GrouperJob`` and capture the
    keyword arguments it was constructed with.

    Grouper jobs have no settings object, so the grouper analogue of
    ``run_gaussian_and_capture_settings`` captures the ``GrouperJob``
    constructor kwargs instead. Returns ``(result, job_kwargs)``;
    ``job_kwargs`` is ``None`` when ``GrouperJob`` was never instantiated
    (e.g. the CLI rejected the arguments).

    Pass ``job_cls`` to supply the mock yourself and inspect its calls
    (e.g. call counts) afterwards. Pass ``catch_exceptions=True`` to
    assert on ``result.exception`` instead of letting errors propagate.
    """
    from chemsmart.cli.grouper import grouper as grouper_cli

    captured_kwargs = None
    if job_cls is None:
        job_cls = MagicMock(name="grouper_job_cls")

    with patch("chemsmart.jobs.grouper.GrouperJob", new=job_cls):
        result = CliRunner().invoke(
            grouper_cli, cli_args, catch_exceptions=catch_exceptions
        )
    if job_cls.call_args is not None:
        captured_kwargs = job_cls.call_args.kwargs

    return result, captured_kwargs


def run_crest_and_capture_settings(job_class_path, cli_args, ctx_obj=None):
    """Run the crest CLI with a patched job class and capture settings.

    Unlike the other ``run_*_and_capture_settings`` helpers, this also
    returns the mock job class itself, since several crest CLI tests need
    to inspect other constructor kwargs (e.g. ``label``) in addition to
    ``settings``::

        result, settings, mock_job_cls = run_crest_and_capture_settings(
            "chemsmart.jobs.crest.conformers.CRESTConformerSearchJob",
            ["-p", "test", "-f", filename, "-c", "0", "-m", "1",
             "conformers"],
        )
        assert result.exit_code == 0
    """
    from chemsmart.cli.crest.crest import crest as crest_cli

    if ctx_obj is None:
        ctx_obj = {"jobrunner": MagicMock()}
    runner = CliRunner()
    captured_settings = None

    with patch(job_class_path) as mock_job_cls:
        mock_job_cls.return_value = MagicMock()
        result = runner.invoke(
            crest_cli,
            cli_args,
            obj=ctx_obj,
            catch_exceptions=False,
        )
        if mock_job_cls.call_args is not None:
            captured_settings = mock_job_cls.call_args[1].get("settings")

    return result, captured_settings, mock_job_cls


def invoke_mol_with_visualize(cli_args, ctx_obj=None):
    """Invoke ``mol … visualize`` with all PyMOL job execution mocked out."""
    from chemsmart.cli.mol.mol import mol as mol_group

    runner = CliRunner()
    if ctx_obj is None:
        ctx_obj = {}

    @mol_group.command("_test_noop")
    @click.pass_context
    def _noop(ctx):
        pass

    try:
        with (
            patch.object(
                mol_cli_module,
                "Molecule",
                MagicMock(),
            ) as mock_molecule_cls,
            patch.object(
                mol_cli_module,
                "BaseFolder",
                MagicMock(),
            ) as mock_folder_cls,
        ):
            # Setup folder mock to return mock files
            mock_folder_instance = MagicMock()
            mock_folder_instance.get_all_output_files_in_current_folder_by_program.return_value = [
                "/fake/dir/mol.log"
            ]
            mock_folder_cls.return_value = mock_folder_instance

            # Setup molecule mock
            mock_molecule_cls.from_filepath.return_value = MagicMock()

            result = runner.invoke(
                mol_group,
                cli_args + ["_test_noop"],
                obj=ctx_obj,
                catch_exceptions=False,
            )
    finally:
        mol_group.commands.pop("_test_noop", None)

    return result


def invoke_mol_cli(cli_args, ctx_obj=None):
    """Invoke ``mol`` CLI with provided args."""
    from chemsmart.cli.mol.mol import mol as mol_group

    runner = CliRunner()
    if ctx_obj is None:
        ctx_obj = {}
    return runner.invoke(
        mol_group, cli_args, obj=ctx_obj, catch_exceptions=False
    )


def run_mol_and_capture_kwargs(job_class_path, cli_args):
    """Run the ``mol`` CLI with ``job_class_path`` patched and capture the
    keyword arguments the job class was constructed with (``None`` if it
    never was)::

        result, kwargs = run_mol_and_capture_kwargs(
            "chemsmart.jobs.mol.nci.PyMOLNCIJob",
            ["-f", xyz_file, "nci"],
        )
    """
    with patch(job_class_path) as mock_job_cls:
        mock_job_cls.return_value = MagicMock()
        result = invoke_mol_cli(cli_args)
        call_kwargs = None
        if mock_job_cls.call_args is not None:
            call_kwargs = mock_job_cls.call_args.kwargs
    return result, call_kwargs


def run_nciplot_and_capture_kwargs(cli_args):
    """Run the ``nciplot`` CLI with ``NCIPLOTJob`` patched and capture the
    keyword arguments it was constructed with (``None`` if it never was)::

        result, kwargs = run_nciplot_and_capture_kwargs(["-f", xyz_file])
    """
    with patch.object(nciplot_cli_module, "NCIPLOTJob") as mock_job_cls:
        mock_job_cls.return_value = MagicMock()
        result = CliRunner().invoke(nciplot, cli_args, catch_exceptions=False)
        call_kwargs = None
        if mock_job_cls.call_args is not None:
            call_kwargs = mock_job_cls.call_args.kwargs
    return result, call_kwargs


def make_thermochemistry_mock(
    vibrational_frequencies,
    jobtype="opt",
    check_imaginary_frequencies=True,
    s_freq_cutoff_cm=None,
    h_freq_cutoff_cm=None,
    rotational_mode="gaussian",
):
    """Build a ``MagicMock`` mimicking a :class:`Thermochemistry` instance,
    configured with the attributes ``Thermochemistry.cleaned_frequencies``
    reads, so tests stay free of direct ``MagicMock`` construction::

        mock = make_thermochemistry_mock(
            vibrational_frequencies=[-50.0, 100.0],
            jobtype="opt",
            check_imaginary_frequencies=False,
        )
        result = Thermochemistry.cleaned_frequencies.fget(mock)
        assert result == [100.0, 100.0]
    """
    from chemsmart.analysis.thermochemistry import Thermochemistry

    mock = MagicMock(spec=Thermochemistry)
    mock.vibrational_frequencies = vibrational_frequencies
    mock.jobtype = jobtype
    mock.check_imaginary_frequencies = check_imaginary_frequencies
    mock.s_freq_cutoff_cm = s_freq_cutoff_cm
    mock.h_freq_cutoff_cm = h_freq_cutoff_cm
    mock.rotational_mode = rotational_mode
    mock.filename = "dummy.log"
    mock.target = "dummy.log"
    return mock


def bare_xtb_input(content_groups):
    """Factory for an XTBInput that skips file parsing and uses the
    given content_groups directly (seeding the cached_property).
    """
    inp = XTBInput.__new__(XTBInput)
    inp.__dict__["content_groups"] = content_groups
    return inp


# ── pKa CLI helpers ──
def write_program_signature_file(path: Path, program: str):
    signatures = {
        "gaussian": "Gaussian, Inc.\n",
        "orca": "* O   R   C   A *\n",
        "unknown": "Some random text\n",
    }
    path.write_text(signatures[program])


def build_pka_outputs(tmp_path: Path, program: str):
    names = [
        "ha.log",
        "a.log",
        "hb.log",
        "b.log",
        "has.log",
        "as.log",
        "hbs.log",
        "bs.log",
    ]
    files = {}
    for n in names:
        p = tmp_path / n
        write_program_signature_file(p, program)
        files[n] = str(p)
    return files


def invoke_pka_direct(runner, files, delta_g_proton=-265.9):
    return runner.invoke(
        run,
        [
            "pka",
            "-s",
            "direct",
            "-dG",
            str(delta_g_proton),
            "analyze",
            "-ha",
            files["ha.log"],
            "-a",
            files["a.log"],
            "-has",
            files["has.log"],
            "-as",
            files["as.log"],
        ],
    )


def invoke_pka(runner, files):
    return runner.invoke(
        run,
        [
            "pka",
            "analyze",
            "-ha",
            files["ha.log"],
            "-a",
            files["a.log"],
            "-hr",
            files["hb.log"],
            "-r",
            files["b.log"],
            "-has",
            files["has.log"],
            "-as",
            files["as.log"],
            "--href-solv",
            files["hbs.log"],
            "--ref-solv",
            files["bs.log"],
            "-rp",
            "6.75",
        ],
    )


def require_backend_pka_subcommand(command_group, backend):
    runner = CliRunner()
    result = runner.invoke(command_group, [backend, "--help"])
    assert result.exit_code == 0, result.output
    if "\n  pka" not in result.output:
        pytest.skip(
            f"{backend} backend pka subcommand is not registered in this build."
        )


class FakeThermochemistry:
    def __init__(self, filename, **kwargs):
        self.filename = filename
        self.electronic_energy = -627.0
        self.qrrho_gibbs_free_energy = -628.0


def write_test_backend_project(tmp_path, backend):
    config_root = tmp_path / "chemsmart_cfg"
    backend_cfg_dir = config_root / backend
    backend_cfg_dir.mkdir(parents=True)
    (backend_cfg_dir / "test.yaml").write_text(
        "gas:\n"
        "  functional: B3LYP\n"
        "  basis: def2-SVP\n"
        "solv:\n"
        "  functional: B3LYP\n"
        "  basis: def2-SVP\n"
        "  freq: false\n"
        "  solvent_model: smd\n"
        "  solvent_id: water\n"
    )
    return config_root


def build_pka_batch_table(tmp_path):
    acid1 = tmp_path / "acid1.xyz"
    acid1.write_text("2\nacid1\nC 0.0 0.0 0.0\nH 0.0 0.0 1.0\n")
    acid2 = tmp_path / "acid2.xyz"
    acid2.write_text("2\nacid2\nN 0.0 0.0 0.0\nH 0.0 0.0 1.0\n")

    table = tmp_path / "pka_scale.csv"
    table.write_text(
        "filepath,proton_index,charge,multiplicity\n"
        f"{acid1},2,0,1\n"
        f"{acid2},2,1,2\n"
    )
    return table


def build_pka_shared(**overrides):
    """Minimal 'pka_shared' dict matching the keys the `pka` group builds."""
    shared = dict(
        scheme="direct",
        reference=None,
        reference_proton_index=None,
        reference_color_code=None,
        reference_charge=None,
        reference_multiplicity=None,
        reference_conjugate_base_charge=None,
        reference_conjugate_base_multiplicity=None,
        delta_g_proton=None,
        conjugate_base_charge=None,
        conjugate_base_multiplicity=None,
        solvent_model=None,
        solvent_id=None,
        temperature=None,
        concentration=None,
        pressure=None,
        cutoff_entropy_grimme=None,
        cutoff_enthalpy=None,
        entropy_method=None,
        skip_completed=False,
    )
    shared.update(overrides)
    return shared


def invoke_sub_process_pipeline_directly(
    cli_tokens,
    batch_entry=None,
    subcommand_list=None,
    _omit_subcommand_key=False,
    **result_kwargs,
):
    """Directly invoke ``sub``'s result callback with a stubbed
    ``CtxObjArguments.reconstruct_command_line`` output.

    This bypasses the need to build a fully realistic
    ``ctx.obj["subcommand"]`` chain (which only ``MyGroup``/``MyCommand``
    invocation normally populates) so that ``_replace_batch_table_tokens``
    can be driven with exact, hand-picked token lists to reach its less
    common rewrite branches.
    """
    from types import SimpleNamespace
    from unittest.mock import MagicMock, patch

    import click

    import chemsmart.cli.sub as sub_module
    from chemsmart.settings.server import Server

    class _StubArgs:
        def __init__(self, commands, entry_point):
            pass

        def reconstruct_command_line(self):
            return ["placeholder"] + list(cli_tokens)

    fake_server = Server(name="dummy")
    captured = {"cli_args": None}
    fake_server.submit = (
        lambda job, test=False, cli_args=None, **kw: captured.update(
            cli_args=cli_args
        )
    )

    job = SimpleNamespace()
    if batch_entry is not None:
        job._batch_entry = batch_entry

    if subcommand_list is None:
        subcommand_list = [{"name": "sub", "kwargs": {}}]

    kwargs = {"test": False, "print_command": False, "server": "dummy"}
    kwargs.update(result_kwargs)

    with (
        patch.object(sub_module, "CtxObjArguments", _StubArgs),
        patch(
            "chemsmart.settings.server.Server.from_servername",
            lambda _name: fake_server,
        ),
    ):
        ctx = click.Context(sub_module.sub)
        mock_jobrunner = MagicMock()
        mock_jobrunner.server = fake_server
        ctx.obj = {"jobrunner": mock_jobrunner}
        if not _omit_subcommand_key:
            ctx.obj["subcommand"] = subcommand_list
        with ctx:
            sub_module.sub._result_callback(job, **kwargs)

    return captured["cli_args"]


# ── pubchem HTTP fakes ──
class MockHTTPResponse:
    def __init__(self, text=None, status_code=200):
        self.text = text
        self.status_code = status_code

    def raise_for_status(self):
        pass


def make_http_response(text=None, status_code=200, raise_exc=None):
    response = MockHTTPResponse(text=text, status_code=status_code)
    if raise_exc is not None:
        response.raise_for_status = lambda: (_ for _ in ()).throw(raise_exc)
    return response


def make_http_error(status_code):
    error = RequestsHTTPError(f"HTTP {status_code}")
    error.response = MockHTTPResponse(status_code=status_code)
    return error


# ── Gaussian project-settings stub ──
class ConcreteGaussianProjectSettings(GaussianProjectSettings):
    """A hardcoded (non-YAML) project settings subclass, exercising
    the base class's own default method implementations."""

    PROJECT_NAME = "concrete_test_project"
    functional = "b3lyp"
    small_basis = "sto-3g"
    large_basis = "def2tzvp"


# ── xTB project-settings stub ──
class MyXTBProjectSettings(XTBProjectSettings):
    PROJECT_NAME = "my_project"
    gfn_version = "gfnff"


# ── CtxObjArguments option kwargs ──
def make_ctx_option_kwargs(
    value,
    nargs=1,
    is_multiple=False,
    type_name="string",
    is_flag=False,
    secondary_opts=None,
):
    return {
        "value": value,
        "nargs": nargs,
        "is_multiple": is_multiple,
        "type": MagicMock(name=type_name),
        "is_flag": is_flag,
        "secondary_opts": secondary_opts or [],
    }


# ── script-analysis mocks ──
def make_fukui_output(shift):
    output = MagicMock()
    output.energies = [-100.0 + shift]
    output.mulliken_atomic_charges = {"1C": 0.1 + shift}
    output.natural_charges = {"1C": 0.2 + shift}
    output.hirshfeld_charges = {"1C": 0.3 + shift}
    output.hirshfeld_cm5_charges = {"1C": 0.4 + shift}
    return output


def make_thermochemistry_results_mock():
    thermo = MagicMock()
    thermo.electronic_energy = -100.0
    thermo.zero_point_energy = 0.05
    thermo.enthalpy = -99.9
    thermo.qrrho_enthalpy = -99.9
    thermo.entropy_times_temperature = -0.02
    thermo.qrrho_entropy_times_temperature = -0.02
    thermo.gibbs_free_energy = -99.92
    thermo.qrrho_gibbs_free_energy = -99.92
    thermo.qrrho_gibbs_free_energy_qh = -99.92
    thermo.qrrho_gibbs_free_energy_qs = -99.92
    thermo.imaginary_frequencies = []
    thermo.jobtype = "opt"
    thermo.vibrational_frequencies = [100.0]
    return thermo


# ── mock job for submitter tests ──
def make_mock_job(
    label="job1", program="Gaussian", complete=False, folder="."
):
    job = MagicMock()
    job.label = label
    job.PROGRAM = program
    job.execution_name = Job.execution_name.fget(job)
    job.folder = folder
    job.is_complete.return_value = complete
    return job


# ── structure (InChI / CXSMILES) test helpers ──
def load_last_molecule(filepath):
    """Read ``filepath`` with ``Molecule.from_filepath``; if it yields a list
    (multi-structure file), return its last molecule."""
    mol = Molecule.from_filepath(filepath)
    if isinstance(mol, list):
        mol = mol[-1]
    return mol


def cxsmiles_core(cxsmiles):
    """Extract the SMILES part before the CX coordinate extension."""
    return cxsmiles.split(" |")[0]


def read_stripped_text(filepath):
    """Return the contents of ``filepath`` with surrounding whitespace
    stripped (e.g. an expected-output reference file)."""
    with open(filepath, "r") as f:
        return f.read().strip()


# ── thermochemistry / Boltzmann CLI helpers ──
BOLTZMANN_JOB = "chemsmart.jobs.thermochemistry.boltzmann.BoltzmannAverageThermochemistryJob"


def run_boltzmann_and_capture_kwargs(cli_args, obj=None):
    """Run the ``thermochemistry ... boltzmann`` CLI with a patched job class
    and capture the constructor call kwargs, plus the mock job instance (to
    check whether ``compute_boltzmann_averages``/``show_results`` were
    invoked). ``obj`` (default ``{}``) is the Click context object.
    Returns ``(result, call_kwargs, mock_instance)``."""
    obj = {} if obj is None else obj
    with patch(BOLTZMANN_JOB) as mock_job_cls:
        mock_instance = MagicMock()
        mock_job_cls.return_value = mock_instance
        result = CliRunner().invoke(
            thermochemistry,
            cli_args,
            obj=obj,
            catch_exceptions=True,
        )
        call_kwargs = None
        if mock_job_cls.call_args is not None:
            call_kwargs = mock_job_cls.call_args.kwargs
    return result, call_kwargs, mock_instance


def run_thermochemistry_directory_mode(cli_args, obj=None):
    """Invoke ``thermochemistry <cli_args> boltzmann`` with a patched
    Boltzmann job class, so the group's own directory/filename-resolution
    logic runs (populating ``obj``, the context object) without reaching
    the real per-job computation. Returns ``(result, obj)``."""
    obj = {} if obj is None else obj
    result, _, _ = run_boltzmann_and_capture_kwargs(
        [*cli_args, "boltzmann"], obj
    )
    return result, obj


def make_spec_thermochemistry_mock(**attrs):
    """Build a ``MagicMock(spec=Thermochemistry)`` with the given attributes,
    for calling ``Thermochemistry`` property/method implementations on a
    controlled stand-in."""
    from chemsmart.analysis.thermochemistry import Thermochemistry

    mock = MagicMock(spec=Thermochemistry)
    for key, value in attrs.items():
        setattr(mock, key, value)
    return mock


def thermochemistry_cached_value(name, mock):
    """Evaluate ``Thermochemistry``'s cached property ``name`` on ``mock``
    (calls the underlying function, bypassing the cache)."""
    from chemsmart.analysis.thermochemistry import Thermochemistry

    return Thermochemistry.__dict__[name].func(mock)


# ── update / configuration-updater helpers ──
class FakeYamlParser:
    """Stand-in for the real ruamel.yaml.YAML() parser instance returned
    by ConfigurationUpdater._yaml(), for exercising load()'s error and
    not-a-Mapping branches without needing malformed real files."""

    def __init__(self, load_result=None, load_error=None):
        self._result = load_result
        self._error = load_error

    def load(self, text):
        if self._error is not None:
            raise self._error
        return self._result


# ── utils test helpers ──
def make_executable_file(path):
    """Write a minimal ``#!/bin/sh`` script at ``path`` and mark it
    executable (``0o755``), e.g. a stand-in for an external program on
    ``PATH``."""
    write_text_file(path, "#!/bin/sh\n")
    os.chmod(str(path), 0o755)


# ── xTB CLI helpers ──
def invoke_xtb_cli(cli_args, ctx_obj=None, *, catch_exceptions=False):
    """Invoke the ``xtb`` CLI group with ``cli_args`` and return the Click
    result. ``ctx_obj`` defaults to ``{"jobrunner": MagicMock()}``; pass
    ``catch_exceptions=True`` to assert on ``result.exception`` /
    ``result.output`` instead of letting errors propagate."""
    from chemsmart.cli.xtb.xtb import xtb as xtb_cli

    return CliRunner().invoke(
        xtb_cli,
        cli_args,
        obj={"jobrunner": MagicMock()} if ctx_obj is None else ctx_obj,
        catch_exceptions=catch_exceptions,
    )


def run_xtb_and_capture_settings(job_path, args, *, catch_exceptions=False):
    """Run the ``xtb`` CLI with the job class at ``job_path`` patched and
    capture the ``settings`` it was constructed with (``None`` if it never
    was). Returns ``(result, settings, mock_job_cls)``; inspect
    ``mock_job_cls`` for labels, call counts, etc."""
    mock_job = MagicMock()
    with patch(job_path, return_value=mock_job) as mock_job_cls:
        result = invoke_xtb_cli(args, catch_exceptions=catch_exceptions)
        settings = (
            mock_job_cls.call_args[1]["settings"]
            if mock_job_cls.call_args is not None
            else None
        )
    return result, settings, mock_job_cls


def assert_runner_was_configured(job, options):
    """Assert a job factory passed its runner options to the runner only.

    ``options`` is the ``server``, ``scratch`` and ``fake`` mapping given to
    the factory, as in the ``runner_options`` fixture.
    """
    assert job.jobrunner is not None
    assert job.jobrunner.server is options["server"]
    assert job.jobrunner.scratch is False
    assert job.jobrunner.fake is True
    assert not set(options) & set(job.kwargs)


def invoke_sterimol(*args):
    """Run ``chemsmart run sterimol`` with ``args``.

    Names a server because ``run`` otherwise needs ``local.yaml`` in the user's
    ``~/.chemsmart``, which CI does not have.
    """
    server_yaml = Path(__file__).parent / "data" / "YAMLTests" / "server.yaml"
    return CliRunner().invoke(
        run, ["-s", str(server_yaml), "sterimol", *map(str, args)]
    )


def run_sterimol_to_json(tmp_path, *args):
    """Run ``chemsmart run sterimol`` successfully; return its JSON report."""
    output = tmp_path / "report.json"
    result = invoke_sterimol(*args, "--output", output)
    assert result.exit_code == 0, result.output
    return json.loads(output.read_text())


def assert_sterimol_matches(report, reference):
    """Compare a report with reference (L, B1, B5) values.

    L and B5 are analytical. B1 is sampled more coarsely than the morfeus
    references, so it has a wider tolerance.
    """
    length, b1, b5 = reference
    assert report["L"] == pytest.approx(length, abs=1e-4)
    assert report["B1"] == pytest.approx(b1, abs=1e-3)
    assert report["B5"] == pytest.approx(b5, abs=1e-4)


def write_xyz_file(path, symbols, frames):
    """Write the frames (each an (N, 3) array) to ``path``; return the path."""
    lines = []
    for positions in frames:
        lines += [str(len(symbols)), "generated for a test"]
        lines += [
            f"{symbol} {x:.6f} {y:.6f} {z:.6f}"
            for symbol, (x, y, z) in zip(symbols, positions)
        ]
    path.write_text("\n".join(lines) + "\n")
    return str(path)


def rdkit_embed_and_optimize(smiles, seed=1):
    """Return an RDKit Mol with one UFF-optimized 3D conformer (Hs added)."""
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    # fixed seed for deterministic embedding
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    AllChem.EmbedMolecule(mol, params)
    AllChem.UFFOptimizeMolecule(mol)
    return mol
