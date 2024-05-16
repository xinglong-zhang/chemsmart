from chemsmart.utils.mixins import FileMixin
from chemsmart.utils.periodictable import PeriodicTable as p

class GaussianCubeFile(FileMixin):
    def __init__(self, filename):
        self.filename = filename

    @property
    def cube_job_title(self):
        """ First line of the cube file gives the title of the job/cube."""
        return self.contents[0]

    @property
    def cube_job_description(self):
        """ Second line of the cube file gives the description of the job/cube."""
        return self.contents[1]

    @property
    def num_atoms(self):
        """ Get number of atoms for the system."""
        return int(self.contents[2].split()[0])

    @property
    def coordinate_origin(self):
        """ Get the  x-,y-,z-coordinates of origin."""
        line_elements = self.contents[2].split()
        x = float(line_elements[1])
        y = float(line_elements[2])
        z = float(line_elements[3])
        return (x, y, z)

    @property
    def structure(self):
        pass




