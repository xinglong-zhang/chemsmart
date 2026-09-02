"""PyMOL electrostatic potential (ESP) visualization jobs."""

from chemsmart.jobs.mol.job import PyMOLJob


class PyMOLESPJob(PyMOLJob):
    """PyMOL job for electrostatic potential surface visualization."""

    TYPE = "pymol_esp"

    def __init__(
        self,
        molecule,
        label,
        esp_basename=None,
        color_range=None,
        npts="-2",
        isosurface_value=None,
        transparency_value=None,
        **kwargs,
    ):
        super().__init__(
            molecule=molecule,
            label=label,
            **kwargs,
        )
        # set defaults
        if color_range is None:
            color_range = 0.04
        if npts is None:
            npts = "-2"
        if isosurface_value is None:
            isosurface_value = 0.001
        if transparency_value is None:
            transparency_value = 0.5
        self.color_range = color_range
        self.npts = npts
        self.isosurface_value = isosurface_value
        self.transparency_value = transparency_value

        if esp_basename is None:
            esp_basename = f"{self.label}_ESP"

        self.esp_basename = esp_basename

    def _get_job_basename(self):
        return self.esp_basename
