from pyatoms.jobs.gaussian.job import GaussianJob
from pyatoms.jobs.gaussian.settings import GaussianTDDFTJobSettings


class GaussianTDDFTJob(GaussianJob):
    TYPE = 'g16td'
    _SETTINGS_CLS = GaussianTDDFTJobSettings

    def __init__(self, folder, atoms, settings, **kwargs):
        super().__init__(folder=folder, atoms=atoms, settings=settings, **kwargs)
