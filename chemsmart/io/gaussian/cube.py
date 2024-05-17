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
    def grid_points(self):
        """ Get grid points information."""
        # gridpoints information will be specified before molecular structure
        grid_points_data = []
        for line in self.contents[3:6]:
            grid_points_data.append(int(line.split()[0]))
        return tuple(grid_points_data)

    @property
    def grid_increment_vector(self):
        """ Get increment vector of the grids."""
        grid_vectors = []
        for line in self.contents[3:6]:
            grid_vector = (float(line.split()[1]), float(line.split()[2]), float(line.split()[3]))
            grid_vectors.append(grid_vector)
        return tuple(grid_vectors)

    @property
    def structure(self):
        pass


class CubeFileOperator:
    """ Class to operate on two cube files.
    Args: cubefile1: file path for cube file 1.
          cubefile2: file path for cube file 2.
    """
    def __init__(self, cubefile1, cubefile2):
        self.cubefile1 = cubefile1
        self.cubefile2 = cubefile2
        self.cube1 = GaussianCubeFile(filename=self.cubefile1)
        self.cube2 = GaussianCubeFile(filename=self.cubefile2)

    def _check_natoms_matched(self):
        return self.cube1.num_atoms == self.cube2.num_atoms

    def _check_coordinate_origin_matched(self):
        return self.cube1.coordinate_origin == self.cube2.coordinate_origin







