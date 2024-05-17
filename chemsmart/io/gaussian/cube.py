import os
import numpy as np
from chemsmart.utils.mixins import FileMixin

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
        """ Third line: Get number of atoms for the system."""
        return int(self.contents[2].split()[0])

    @property
    def coordinate_origin(self):
        """ Third line: Get the  x-,y-,z-coordinates of origin."""
        line_elements = self.contents[2].split()
        x = float(line_elements[1])
        y = float(line_elements[2])
        z = float(line_elements[3])
        return (x, y, z)

    @property
    def grid_points(self):
        """ Fourth to sixth line: Get grid points information."""
        # gridpoints information will be specified before molecular structure
        grid_points_data = []
        for line in self.contents[3:6]:
            grid_points_data.append(int(line.split()[0]))
        return tuple(grid_points_data)

    @property
    def grid_increment_vector(self):
        """ Fourth to sixth line: Get increment vector of the grids."""
        grid_vectors = []
        for line in self.contents[3:6]:
            grid_vector = (float(line.split()[1]), float(line.split()[2]), float(line.split()[3]))
            grid_vectors.append(grid_vector)
        return tuple(grid_vectors)

    @property
    def coordinate_block_as_list(self):
        """ 7th to (7+num_atoms)th line: coordinate block."""
        return self.contents[6:6+self.num_atoms]

    @property
    def coordinate_block_object(self):
        from chemsmart.utils.utils import CoordinateBlock
        return CoordinateBlock(coordinate_block=self.coordinate_block_as_list)

    @property
    def chemical_symbols(self):
        return self.coordinate_block_object.chemical_symbols

    @property
    def positions(self):
        return self.coordinate_block_object.positions

    @property
    def structure(self):
        return self.coordinate_block_object.molecule

    @property
    def values_by_lines(self):
        """ Return a list of lists where each member list is a list of floats
        corresponding to the line of values"""
        lines_of_values = []
        for line in self.contents[6+self.num_atoms:]:
            # skip first 2 header lines + next 1 num atoms line + 3 grid vectors lines
            line_of_floats_as_list = [float(x) for x in line.split()]
            lines_of_values.append(line_of_floats_as_list)
        return lines_of_values


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

    def _check_geometries_matched(self):
        """ Method to check that the geometries given in the two cube files are the same."""
        cube1_contents = self.cube1.contents
        cube2_contents = self.cube2.contents
        for i, line in enumerate(cube1_contents[6:]):
            #TODO
            pass

    def _all_checked(self):
        return (self._check_natoms_matched() and self._check_coordinate_origin_matched()
                and self._check_grid_points_matched() and self._check_geometries_matched())

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

    def _write_header(self):
        # write the first two lines

        pass

    def _write_geometry(self):
        """ Write the geometry that is the same. Instead of using structure.write(), we will simply write
        the geometry from the input cube file. This is due to the non-standard geometry specification in
        cube file, where the atomic_number is repeated as a float in the second element."""
        pass
