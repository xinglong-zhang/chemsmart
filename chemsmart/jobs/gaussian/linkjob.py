from pyatoms.jobs.gaussian.job import GaussianJob
from pyatoms.jobs.gaussian.settings import GaussianLinkJobSettings


class GaussianLinkJob(GaussianJob):
    TYPE = 'g16link'
    _SETTINGS_CLS = GaussianLinkJobSettings

    def __init__(self, folder, atoms, settings, **kwargs):
        super().__init__(folder=folder, atoms=atoms, settings=settings, **kwargs)
