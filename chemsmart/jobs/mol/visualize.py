from chemsmart.jobs.mol.job import PyMOLJob


class PyMOLVisualizationJob(PyMOLJob):
    TYPE = "pymol_visualization"

    def __init__(
        self,
        molecule,
        label,
        **kwargs,
    ):
        super().__init__(
            molecule=molecule,
            label=label,
            **kwargs,
        )
