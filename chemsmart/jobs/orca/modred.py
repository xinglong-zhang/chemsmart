from pyatoms.jobs.orca.job import ORCAJob


class ORCAModredJob(ORCAJob):
    TYPE = "orcamodred"

    def __init__(self, folder, atoms, settings, **kwargs):
        super().__init__(
            folder=folder, atoms=atoms, settings=settings, **kwargs
        )
