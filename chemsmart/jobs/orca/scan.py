from chemsmart.jobs.orca.job import ORCAJob


class ORCAScanJob(ORCAJob):
    TYPE = "orcascan"

    def __init__(self, molecule, settings, label, **kwargs):
        super().__init__(
            molecule=molecule, settings=settings, label=label, **kwargs
        )

