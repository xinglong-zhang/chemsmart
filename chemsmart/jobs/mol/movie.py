from chemsmart.jobs.mol.job import PyMOLJob


class PyMOLMovieJob(PyMOLJob):
    TYPE = "pymol_movie"

    def __init__(
        self,
        molecule,
        label,
        pymol_script=None,
        render=None,
        vdw=None,
        quiet_mode=True,
        command_line_only=True,
        trace=True,
        overwrite=False,
        **kwargs,
    ):
        super().__init__(
            molecule=molecule,
            label=label,
            pymol_script=pymol_script,
            render=render,
            trace=trace,
            vdw=vdw,
            quiet_mode=quiet_mode,
            command_line_only=command_line_only,
            **kwargs,
        )
        self.overwrite = overwrite
