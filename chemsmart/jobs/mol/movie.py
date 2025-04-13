from chemsmart.jobs.mol.job import PyMOLJob


class PyMOLMovieJob(PyMOLJob):
    TYPE = "pymol_movie"

    def __init__(
        self,
        molecule,
        label,
        overwrite=False,
        **kwargs,
    ):
        super().__init__(
            molecule=molecule,
            label=label,
            **kwargs,
        )
        self.overwrite = overwrite
