"""Server and job-runner fixtures."""

import pytest

from chemsmart.jobs.crest.runner import FakeCRESTJobRunner
from chemsmart.jobs.gaussian.runner import FakeGaussianJobRunner
from chemsmart.jobs.iterate.runner import IterateJobRunner
from chemsmart.jobs.mol.runner import (
    PyMOLAlignJobRunner,
    PyMOLHybridVisualizationJobRunner,
    PyMOLIRCMovieJobRunner,
    PyMOLMOJobRunner,
    PyMOLMovieJobRunner,
    PyMOLScientificStyleVisualizationJobRunner,
    PyMOLVisualizationJobRunner,
)
from chemsmart.jobs.nciplot.runner import (
    FakeNCIPLOTJobRunner,
    NCIPLOTJobRunner,
)
from chemsmart.jobs.orca.runner import FakeORCAJobRunner
from chemsmart.jobs.thermochemistry.runner import ThermochemistryJobRunner
from chemsmart.jobs.xtb.runner import FakeXTBJobRunner
from chemsmart.settings.server import Server
from tests.helpers import (
    ConcreteJobRunner,
)


@pytest.fixture()
def clear_scheduler_caches():
    """Clear ``Server``'s cached scheduler detection and scheduler-type
    lookup before and after the test, so results cached by one test (or
    from the real host) can't leak into another. Apply with
    ``@pytest.mark.usefixtures("clear_scheduler_caches")``."""
    Server.detect_server_scheduler.cache_clear()
    Server.from_scheduler_type.cache_clear()
    yield
    Server.detect_server_scheduler.cache_clear()
    Server.from_scheduler_type.cache_clear()


@pytest.fixture()
def pbs_server(server_yaml_file):
    return Server.from_yaml(server_yaml_file)


@pytest.fixture()
def runner_options(pbs_server):
    """Options a job factory uses to select a fake runner on ``pbs_server``."""
    return {"server": pbs_server, "scratch": False, "fake": True}


@pytest.fixture()
def gaussian_jobrunner_no_scratch(pbs_server):
    return FakeGaussianJobRunner(server=pbs_server, scratch=False, fake=True)


@pytest.fixture()
def gaussian_jobrunner_scratch(tmpdir, pbs_server):
    return FakeGaussianJobRunner(
        scratch_dir=tmpdir, server=pbs_server, scratch=True, fake=True
    )


@pytest.fixture()
def orca_jobrunner_no_scratch(pbs_server):
    return FakeORCAJobRunner(server=pbs_server, scratch=False, fake=True)


@pytest.fixture()
def orca_jobrunner_scratch(tmpdir, pbs_server):
    return FakeORCAJobRunner(
        scratch_dir=tmpdir, server=pbs_server, scratch=True, fake=True
    )


@pytest.fixture()
def xtb_jobrunner_no_scratch(pbs_server):
    return FakeXTBJobRunner(server=pbs_server, scratch=False, fake=True)


@pytest.fixture()
def xtb_jobrunner_scratch(tmpdir, pbs_server):
    return FakeXTBJobRunner(
        scratch_dir=tmpdir, server=pbs_server, scratch=True, fake=True
    )


@pytest.fixture()
def crest_jobrunner_no_scratch(pbs_server):
    return FakeCRESTJobRunner(server=pbs_server, scratch=False, fake=True)


@pytest.fixture()
def crest_jobrunner_scratch(tmpdir, pbs_server):
    return FakeCRESTJobRunner(
        scratch_dir=tmpdir, server=pbs_server, scratch=True, fake=True
    )


@pytest.fixture()
def make_concrete_jobrunner(pbs_server):
    """Factory for :class:`ConcreteJobRunner`, on ``pbs_server`` without
    scratch by default, e.g. ``make_concrete_jobrunner(num_cores=99)`` or
    ``make_concrete_jobrunner(server="myserver")``."""

    def _make(server=pbs_server, scratch=False, **kwargs):
        return ConcreteJobRunner(server=server, scratch=scratch, **kwargs)

    return _make


@pytest.fixture()
def concrete_jobrunner(make_concrete_jobrunner):
    """A :class:`ConcreteJobRunner` on ``pbs_server`` without scratch."""
    return make_concrete_jobrunner()


@pytest.fixture()
def make_gaussian_com_file(tmp_path):
    """Factory writing a minimal CH Gaussian input file into ``tmp_path``
    and returning its path, e.g.
    ``make_gaussian_com_file(route="# opt freq pm6", mult=2)``."""

    def _make(filename="test.com", route="# opt freq b3lyp/6-31g(d)", mult=1):
        lines = [
            "%chk=test.chk",
            "%mem=4GB",
            route,
            "",
            "title",
            "",
            f"0 {mult}",
            "C 0.0 0.0 0.0",
            "H 0.0 0.0 1.0",
            "",
        ]
        path = tmp_path / filename
        path.write_text("\n".join(lines))
        return str(path)

    return _make


@pytest.fixture()
def pymol_visualization_jobrunner(pbs_server):
    return PyMOLVisualizationJobRunner(server=pbs_server, scratch=False)


@pytest.fixture()
def pymol_hybrid_visualization_jobrunner(pbs_server):
    return PyMOLHybridVisualizationJobRunner(server=pbs_server, scratch=False)


@pytest.fixture()
def pymol_scientific_style_visualization_jobrunner(pbs_server):
    return PyMOLScientificStyleVisualizationJobRunner(
        server=pbs_server, scratch=False
    )


@pytest.fixture()
def pymol_movie_jobrunner(pbs_server):
    return PyMOLMovieJobRunner(server=pbs_server, scratch=False)


@pytest.fixture()
def nciplot_jobrunner(pbs_server):
    """Real (non-fake) :class:`NCIPLOTJobRunner` on ``pbs_server``, no
    scratch."""
    return NCIPLOTJobRunner(server=pbs_server, scratch=False)


@pytest.fixture()
def nciplot_jobrunner_no_scratch(pbs_server):
    return FakeNCIPLOTJobRunner(server=pbs_server, scratch=False, fake=True)


@pytest.fixture()
def nciplot_jobrunner_scratch(tmpdir, pbs_server):
    return FakeNCIPLOTJobRunner(
        scratch_dir=tmpdir, server=pbs_server, scratch=True, fake=True
    )


@pytest.fixture()
def pymol_align_jobrunner(pbs_server):
    return PyMOLAlignJobRunner(server=pbs_server, scratch=False)


@pytest.fixture()
def pymol_ircmovie_jobrunner(pbs_server):
    return PyMOLIRCMovieJobRunner(server=pbs_server, scratch=False)


@pytest.fixture()
def pymol_mo_jobrunner(pbs_server):
    return PyMOLMOJobRunner(server=pbs_server, scratch=False)


@pytest.fixture()
def thermochemistry_jobrunner(pbs_server):
    """:class:`ThermochemistryJobRunner` on ``pbs_server``, no scratch."""
    return ThermochemistryJobRunner(server=pbs_server, scratch=False)


@pytest.fixture()
def iterate_jobrunner(pbs_server):
    return IterateJobRunner(server=pbs_server, scratch=False)


@pytest.fixture()
def fake_iterate_jobrunner(pbs_server):
    return IterateJobRunner(server=pbs_server, scratch=False, fake=True)
