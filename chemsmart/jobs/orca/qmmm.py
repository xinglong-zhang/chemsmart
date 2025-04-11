from chemsmart.jobs.orca.job import ORCAJob


class ORCAQMMMJob(ORCAJob):
    TYPE = "orcaqmmm"

    def __init__(self, molecule, settings, label, **kwargs):
        super().__init__(
            molecule=molecule, settings=settings, label=label, **kwargs
        )
