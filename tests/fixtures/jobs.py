"""Job, job-settings, mock and file/server fixtures."""

from shutil import copyfile
from unittest.mock import MagicMock

import click
import pytest
from pytest_mock import MockerFixture

from chemsmart.cli.run import run
from chemsmart.cli.update import DependencyUpdater, VersionUpdater
from chemsmart.database.database import Database
from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.crest.settings import CRESTJobSettings
from chemsmart.jobs.gaussian.settings import (
    GaussianIRCJobSettings,
    GaussianJobSettings,
    GaussianLinkJobSettings,
)
from chemsmart.jobs.nciplot.job import NCIPLOTJob
from chemsmart.jobs.nciplot.settings import NCIPLOTJobSettings
from chemsmart.jobs.orca.job import ORCAJob
from chemsmart.jobs.orca.settings import ORCAJobSettings
from chemsmart.jobs.runner import JobRunner
from chemsmart.jobs.thermochemistry.job import ThermochemistryJob
from chemsmart.jobs.xtb.settings import XTBJobSettings
from chemsmart.settings.crest import CRESTProjectSettings
from chemsmart.settings.server import Server
from chemsmart.settings.submitters import PBSSubmitter
from chemsmart.settings.user import CHEMSMARTUserSettings
from tests.helpers import make_concrete_job, make_matrix_grouper, make_mock_job


############ Job Settings Fixtures ##################
@pytest.fixture()
def crest_settings():
    return CRESTJobSettings.default()


@pytest.fixture()
def my_crest_project_settings():
    class MyCRESTProjectSettings(CRESTProjectSettings):
        PROJECT_NAME = "my_project"
        gfn_version = "gfnff"

    return MyCRESTProjectSettings()


@pytest.fixture()
def gaussian_irc_settings():
    return GaussianIRCJobSettings.default()


@pytest.fixture()
def gaussian_settings():
    return GaussianJobSettings.default()


@pytest.fixture()
def orca_settings():
    return ORCAJobSettings.default()


@pytest.fixture()
def xtb_settings():
    return XTBJobSettings.default()


@pytest.fixture()
def nciplot_settings():
    return NCIPLOTJobSettings()


@pytest.fixture()
def orca_qrc_settings():
    settings = ORCAJobSettings.default()
    settings.jobtype = "qrc"
    return settings


@pytest.fixture
def gaussian_qrc_settings():
    # Use a real settings object to pass isinstance checks
    settings = GaussianJobSettings()
    settings.jobtype = "qrc"
    return settings


@pytest.fixture()
def gaussian_link_irc_settings():
    return GaussianLinkJobSettings(jobtype="irc", direction=None)


@pytest.fixture()
def gaussian_link_opt_settings():
    return GaussianLinkJobSettings(jobtype="opt")


############ Mock Fixtures ##################
@pytest.fixture()
def mock_molecule():
    mol = MagicMock(spec=Molecule)
    mol.has_vibrations = True
    mol.copy.return_value = mol
    mol.vibrationally_displaced.return_value = mol
    mol.get_chemical_formula.return_value = "C1H4"
    return mol


@pytest.fixture()
def mock_jobrunner():
    return MagicMock(spec=JobRunner)


@pytest.fixture()
def mock_popen(mocker):
    """Fixture to mock subprocess.Popen."""
    return mocker.patch("subprocess.Popen")


@pytest.fixture(scope="session")
def session_mocker(pytestconfig):
    """Session-scoped mocker fixture for patching during the test session."""
    from unittest.mock import MagicMock

    mocker = MockerFixture(pytestconfig)
    mock = MagicMock()
    mocker.patch = mock.patch
    mocker.patch.object = mock.patch.object
    yield mocker
    mocker.resetall()


############ File, Database and Server Fixtures ##################
@pytest.fixture
def raw_rows_db(tmp_path):
    """Database with two records and two molecules sharing ID prefixes,
    and no structures linked to either record."""
    db = Database(str(tmp_path / "raw_rows.db"))
    db.create()
    conn = db.get_connection()
    for record_id, program in (("rec111", "Gaussian"), ("rec222", "ORCA")):
        conn.execute(
            "INSERT INTO records (record_id, program) VALUES (?, ?)",
            (record_id, program),
        )
    for molecule_id in ("mol111", "mol222"):
        conn.execute(
            "INSERT INTO molecules (molecule_id) VALUES (?)",
            (molecule_id,),
        )
    conn.commit()
    conn.close()
    return db


@pytest.fixture()
def water_sp_inp_copy(tmp_path, water_sp_input_path):
    target = tmp_path / "water_sp.inp"
    copyfile(water_sp_input_path, target)
    return str(target)


@pytest.fixture()
def valid_iterate_config_toml(tmp_path):
    path = tmp_path / "config.toml"
    path.write_text("[[skeletons]]\nfile = 'a.xyz'\n")
    return str(path)


@pytest.fixture
def make_dummy_files(tmp_path):
    def _make(*names):
        for name in names:
            (tmp_path / name).write_text("dummy content\n")
        return tmp_path

    return _make


@pytest.fixture()
def pbs_test_server():
    return Server(
        "custom-pbs",
        SCHEDULER="PBS",
        NUM_CORES=8,
        MEM_GB=24,
        NUM_GPUS=0,
        NUM_HOURS=24,
        QUEUE_NAME="normal",
    )


@pytest.fixture()
def make_pbs_submitter_with_executable(pbs_test_server, monkeypatch):
    """Factory for a ``PBSSubmitter`` (mock job, ``pbs_test_server``) whose
    ``executable`` is a mock with ``conda_env``/``modules``/``scripts``/
    ``envars`` set to ``None`` unless overridden by keyword. The patch is
    undone automatically at teardown::

        submitter = make_pbs_submitter_with_executable(modules=["module load g16\n"])
    """

    def _make(**exec_kwargs):
        submitter = PBSSubmitter(job=make_mock_job(), server=pbs_test_server)
        mock_exec = MagicMock(
            conda_env=None, modules=None, scripts=None, envars=None
        )
        for key, value in exec_kwargs.items():
            setattr(mock_exec, key, value)
        monkeypatch.setattr(type(submitter), "executable", mock_exec)
        return submitter

    return _make


@pytest.fixture()
def slurm_test_server():
    return Server(
        "custom-slurm",
        SCHEDULER="SLURM",
        NUM_CORES=8,
        MEM_GB=24,
        NUM_GPUS=0,
        NUM_HOURS=24,
        QUEUE_NAME="normal",
    )


@pytest.fixture()
def concrete_job(tmp_path, monkeypatch):
    """A ``skip_completed`` :class:`ConcreteJob` labelled ``myjob``, created
    from inside ``tmp_path`` (so its folder is ``tmp_path``)."""
    monkeypatch.chdir(tmp_path)
    return make_concrete_job(label="myjob", skip_completed=True)


@pytest.fixture()
def nciplot_job_from_files(
    single_molecule_xyz_file, nciplot_settings, tmp_path, monkeypatch
):
    monkeypatch.chdir(tmp_path)
    return NCIPLOTJob(
        filenames=[single_molecule_xyz_file],
        settings=nciplot_settings,
        label="ncitest",
    )


@pytest.fixture()
def orca_job(single_molecule_xyz_file, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    molecule = Molecule.from_filepath(single_molecule_xyz_file)
    settings = ORCAJobSettings.default()
    return ORCAJob(molecule=molecule, settings=settings, label="orcatest")


@pytest.fixture()
def thermo_job(gaussian_co2_opt_outfile, tmp_path):
    job = ThermochemistryJob(
        filename=gaussian_co2_opt_outfile, label="co2_thermo"
    )
    job.folder = str(tmp_path)
    return job


@pytest.fixture()
def run_ctx():
    server = Server(name="dummy")
    ctx = click.Context(run)
    ctx.ensure_object(dict)
    ctx.obj["jobrunner"] = JobRunner(server=server, fake=True)
    return ctx


@pytest.fixture
def user_settings(tmp_path, monkeypatch):
    monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(tmp_path))
    return CHEMSMARTUserSettings()


@pytest.fixture()
def dependency_updater(tmp_path):
    u = DependencyUpdater()
    u._package_path = tmp_path
    u._pyproject_path = tmp_path / "pyproject.toml"
    return u


@pytest.fixture()
def version_updater(tmp_path):
    u = VersionUpdater("0.0.0")
    u._pyproject_path = tmp_path / "pyproject.toml"
    u._version_file_path = tmp_path / "VERSION"
    u._docs_conf_file_path = tmp_path / "conf.py"
    return u


@pytest.fixture
def center_matrix_grouper(methanol_molecules):
    return make_matrix_grouper(
        methanol_molecules, representative_strategy="center"
    )
