"""Shared pytest configuration.

Fixtures are defined in ``tests/fixtures/<area>.py`` and registered here
via ``pytest_plugins``, so every test can request them by name without
importing. Plain stub classes and helper functions live in
``tests/helpers.py`` and are imported directly.
"""

import numpy as np
import pytest

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.gaussian.mecp import GaussianMECPJob
from chemsmart.jobs.gaussian.settings import GaussianMECPJobSettings

pytest_plugins = [
    "tests.fixtures.io",
    "tests.fixtures.thermochemistry",
    "tests.fixtures.cli",
    "tests.fixtures.gaussian",
    "tests.fixtures.orca",
    "tests.fixtures.common",
    "tests.fixtures.xtb",
    "tests.fixtures.crest",
    "tests.fixtures.structures",
    "tests.fixtures.jobrunners",
    "tests.fixtures.molecules",
    "tests.fixtures.grouper",
    "tests.fixtures.mixins",
    "tests.fixtures.logs",
    "tests.fixtures.iterate",
    "tests.fixtures.database",
    "tests.fixtures.jobs",
    "tests.fixtures.pymol",
    "tests.fixtures.pka",
]


MECP_DRIVER_SETTINGS = dict(
    step_size=0.1,
    trust_radius=0.3,
    adaptive_step_size=False,
    step_size_method="bb",
    restart=False,
    mecp_numfreq=False,
    energy_diff_tol=1e-4,
    force_max_tol=1e-3,
    force_rms_tol=1e-3,
    disp_max_tol=1e-3,
    disp_rms_tol=1e-3,
)


@pytest.fixture()
def make_mecp_job(
    tmp_path,
    gaussian_full_gen_outfile,
    gaussian_jobrunner_no_scratch,
):
    """Build a real GaussianMECPJob rooted in ``tmp_path``.

    The molecule defaults to the final geometry of the
    bromochloromethane Gaussian output fixture; keyword arguments are
    forwarded to :class:`GaussianMECPJobSettings`.
    """

    def _make(label="crossing", molecule=None, **settings_kwargs):
        if molecule is None:
            molecule = Molecule.from_filepath(gaussian_full_gen_outfile)
        job = GaussianMECPJob(
            molecule=molecule,
            settings=GaussianMECPJobSettings(**settings_kwargs),
            label=label,
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        job.set_folder(str(tmp_path))
        return job

    return _make


@pytest.fixture()
def mecp_job(make_mecp_job):
    return make_mecp_job()


@pytest.fixture()
def make_driver_job(make_mecp_job):
    """MECP job on a single H atom with fast optimizer settings."""

    def _make(max_steps=2, **settings_kwargs):
        settings = {**MECP_DRIVER_SETTINGS, **settings_kwargs}
        job = make_mecp_job(
            label="driver",
            molecule=Molecule(symbols=["H"], positions=np.zeros((1, 3))),
            max_steps=max_steps,
            **settings,
        )
        job._export_final_gaussian_logs = lambda sources: None
        return job

    return _make


@pytest.fixture()
def adaptive_step_job(make_mecp_job):
    return make_mecp_job(
        step_size=0.1,
        step_size_grow=1.2,
        step_size_shrink=0.5,
        step_size_min=1.0e-4,
        step_size_max=1.0,
    )


@pytest.fixture
def make_cli_ctx_obj():
    from tests.helpers import make_cli_ctx_obj as helper

    return helper


@pytest.fixture
def run_gaussian_and_capture_settings():
    from tests.helpers import run_gaussian_and_capture_settings as helper

    return helper


@pytest.fixture
def run_thermochemistry_and_capture_settings():
    from tests.helpers import (
        run_thermochemistry_and_capture_settings as helper,
    )

    return helper


@pytest.fixture
def run_orca_and_capture_settings():
    from tests.helpers import run_orca_and_capture_settings as helper

    return helper
