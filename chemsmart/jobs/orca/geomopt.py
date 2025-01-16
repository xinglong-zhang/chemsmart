from pyatoms.jobs.orca.job import ORCAJob


class ORCAGeomOptJob(ORCAJob):
    TYPE = "orcaopt"

    def __init__(self, folder, atoms, settings, **kwargs):
        super().__init__(
            folder=folder, atoms=atoms, settings=settings, **kwargs
        )
