import os
import numpy as np
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

    @property
    def values_by_lines(self):
        """ Return a numpy array of list of lists where each member list is a list of floats
        corresponding to the line of values"""
        lines_of_values = []
        for line in self.contents[6+self.num_atoms:]:
            # skip first 2 header lines + next 1 num atoms line + 3 grid vectors lines
            line_of_floats_as_list = [float(x) for x in line.split()]
            lines_of_values.append(line_of_floats_as_list)
        return np.array(lines_of_values)


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

    def _check_grid_points_matched(self):
        grid_points_matched = self.cube1.grid_points == self.cube2.grid_points
        grid_increment_vectors_matched = self.cube1.grid_increment_vector == self.cube2.grid_increment_vector
        return grid_points_matched and grid_increment_vectors_matched

    def _all_checked(self):
        return self._check_natoms_matched() and self._check_coordinate_origin_matched() and self._check_grid_points_matched()

    def add_values(self):
        cube1_values_by_lines = self.cube1.values_by_lines
        cube2_values_by_lines = self.cube2.values_by_lines
        resultant_values_by_lines = np.add(cube1_values_by_lines, cube2_values_by_lines)  # Element-wise addition
        return resultant_values_by_lines

    def subtract_values(self):
        cube1_values_by_lines = self.cube1.values_by_lines
        cube2_values_by_lines = self.cube2.values_by_lines
        resultant_values_by_lines = np.subtract(cube1_values_by_lines, cube2_values_by_lines)  # Element-wise subtraction
        return resultant_values_by_lines

    def write_results(self, operation='subtract', output_cubefile=None):
        if output_cubefile is None:
            output_cubefile = os.path.join(self.cube1.filepath_directory,
                                  f'{self.cube1.base_filename_with_extension}_{operation}.cube')
        if operation.lower() == 'subtract':
            resultant_values = self.subtract_values()
        elif operation.lower() == 'add':
            resultant_values = self.add_values()

        with open(output_cubefile, 'w') as f:
            self._write_header()  # 2 lines
            self._write_num_atoms()  # next 1 line
            self._write_grids()  # next 3 lines
            self._write_geometry()
            self._write_values_by_line(resultant_values)
        f.close()
