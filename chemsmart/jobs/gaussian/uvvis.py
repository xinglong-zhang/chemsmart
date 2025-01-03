from copy import deepcopy

from pyatoms.jobs.gaussian.geomopt import GaussianGeomOptJob
from pyatoms.jobs.gaussian.job import GaussianJob
from pyatoms.jobs.gaussian.singlepoint import GaussianSinglePointJob


class GaussianUVVISJob(GaussianJob):
    # TODO: incomplete job
    TYPE = 'g16uvvis'

    def __init__(self, folder, atoms, settings):
        super().__init__(folder=folder, atoms=atoms, settings=settings)
        self.settings = settings

    def _gs_geom_opt(self):
        # Ground state geometry optimization and frequencies
        # add equilibrium solvation if necessary
        gs_geom_opt_settings = deepcopy(self.settings)
        return GaussianGeomOptJob(folder=self.folder, atoms=self.atoms, settings=gs_geom_opt_settings)

    def _run_gs_geom_opt(self, jobrunner, queue_manager=None):
        self._gs_geom_opt().run(jobrunner=jobrunner, queue_manager=queue_manager)

    def _gs_geom_opt_complete(self):
        return self._gs_geom_opt().is_complete()

    def _vertical_excitation(self):
        # Vertical excitation with linear response solvation
        # single-point TD-DFT of vertical excitation (with default solvation), at  ground state equilibrium geometry
        # i.e., linear response, non-equilibrium.
        if self._gs_geom_opt_complete():
            sp_settings = deepcopy(self.settings)
            sp_settings.job_type = 'sp'
            sp_settings.freq = False
            return GaussianSinglePointJob(folder=self.folder, atoms=self.atoms, settings=sp_settings)
        return None

    def _run_vertical_excitation(self, jobrunner, queue_manager=None):
        self._vertical_excitation().run(jobrunner=jobrunner, queue_manager=queue_manager)
        return self._vertical_excitation().is_complete()

    def _run(self, jobrunner, queue_manager=None, **kwargs):
        self._run_gs_geom_opt(jobrunner=jobrunner, queue_manager=queue_manager)
        self._run_vertical_excitation(jobrunner=jobrunner, queue_manager=queue_manager)
