from chemsmart.jobs.mol.job import PyMOLJob


class PyMOLVisualizationJob(PyMOLJob):
    TYPE = "pymol_visualization"

    def __init__(
        self,
        molecule,
        label,
        pymol_script=None,
        render_style=None,
        vdw=None,
        quiet_mode=True,
        command_line_only=True,
        **kwargs,
    ):
        super().__init__(
            molecule=molecule,
            label=label,
            pymol_script=pymol_script,
            render_style=render_style,
            vdw=vdw,
            quiet_mode=quiet_mode,
            command_line_only=command_line_only,
            **kwargs,
        )
