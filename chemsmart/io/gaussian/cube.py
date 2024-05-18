import os
import numpy as np
from chemsmart.utils.mixins import FileMixin


class GaussianCubeFile(FileMixin):
    def __init__(self, filename):
        self.filename = filename

    @property
    def cube_job_title(self):
        """First line of the cube file gives the title of the job/cube."""
        return self.contents[0]

    @property
    def cube_job_description(self):
        """Second line of the cube file gives the description of the job/cube."""
        return self.contents[1]

    @property
    def num_atoms(self):
        """Third line: Get number of atoms for the system."""
        return int(self.contents[2].split()[0])

    @property
    def coordinate_origin(self):
        """Third line: Get the  x-,y-,z-coordinates of origin."""
        line_elements = self.contents[2].split()
        x = float(line_elements[1])
        y = float(line_elements[2])
        z = float(line_elements[3])
        return (x, y, z)

    @property
    def grid_points(self):
        """Fourth to sixth line: Get grid points information."""
        # gridpoints information will be specified before molecular structure
        grid_points_data = []
        for line in self.contents[3:6]:
            grid_points_data.append(int(line.split()[0]))
        return tuple(grid_points_data)

    @property
    def grid_increment_vector(self):
        """Fourth to sixth line: Get increment vector of the grids."""
        grid_vectors = []
        for line in self.contents[3:6]:
            grid_vector = (
                float(line.split()[1]),
                float(line.split()[2]),
                float(line.split()[3]),
            )
            grid_vectors.append(grid_vector)
        return tuple(grid_vectors)

    @property
    def coordinate_block_as_list(self):
        """7th to (7+num_atoms)th line: coordinate block."""
        return self.contents[6 : 6 + self.num_atoms]

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
        """Return a list of lists where each member list is a list of floats
        corresponding to the line of values"""
        lines_of_values = []
        for line in self.contents[6 + self.num_atoms :]:
            # skip first 2 header lines + next 1 num atoms line + 3 grid vectors lines
            line_of_floats_as_list = [float(x) for x in line.split()]
            lines_of_values.append(line_of_floats_as_list)
        return lines_of_values


class CubeFileOperator:
    """Class to operate on two cube files.
    Args: cubefile1: file path for cube file 1.
          cubefile2: file path for cube file 2.
    """

    def __init__(
        self, cubefile1, cubefile2, operation="subtract", output_cubefile=None
    ):
        self.cubefile1 = cubefile1
        self.cubefile2 = cubefile2
        self.cube1 = GaussianCubeFile(filename=self.cubefile1)
        self.cube2 = GaussianCubeFile(filename=self.cubefile2)
        self.operation = operation
        if output_cubefile is None:
            output_cubefile = os.path.join(
                self.cube1.filepath_directory,
                f"{self.cube1.basename}_{operation}.cube",
            )
        self.output_cubefile = output_cubefile

    def _check_natoms_matched(self):
        return self.cube1.num_atoms == self.cube2.num_atoms

    def _check_coordinate_origin_matched(self):
        return all(
            np.isclose(
                self.cube1.coordinate_origin,
                self.cube2.coordinate_origin,
                atol=1e-4,
            )
        )

    def _check_grid_points_matched(self):
        grid_points_matched = self.cube1.grid_points == self.cube2.grid_points
        grid_increment_vectors_matched = (
            self.cube1.grid_increment_vector
            == self.cube2.grid_increment_vector
        )
        return grid_points_matched and grid_increment_vectors_matched

    def _check_geometries_matched(self):
        """Method to check that the geometries given in the two cube files are the same."""
        symbols_matched = (
            self.cube1.chemical_symbols == self.cube2.chemical_symbols
        )
        positions_are_close = np.isclose(
            self.cube1.positions, self.cube2.positions, atol=1e-3
        )
        positions_all_close = np.all(positions_are_close)

        if symbols_matched and positions_all_close:
            return True
        else:
            return False

    def _all_checked(self):
        return (
            self._check_natoms_matched()
            and self._check_coordinate_origin_matched()
            and self._check_grid_points_matched()
            and self._check_geometries_matched()
        )

    def _check_value_lines(self):
        if len(self.cube1.values_by_lines) != len(self.cube2.values_by_lines):
            raise ValueError(
                f"The two cube files do not have the same number of lines.\n"
                f"Cube 1 has {len(self.cube1.values_by_lines)} lines and"
                f"Cube 2 has {len(self.cube2.values_by_lines)} lines."
            )

        for line1, line2 in zip(
            self.cube1.values_by_lines, self.cube2.values_by_lines
        ):
            if len(line1) != len(line2):
                raise ValueError(
                    "The two cube files have lines with different lengths."
                )

    def add_values(self):
        cube1_values_by_lines = self.cube1.values_by_lines
        cube2_values_by_lines = self.cube2.values_by_lines
        print(cube1_values_by_lines)

        # Ensure both lists have the same structure
        self._check_value_lines()
        resultant_values_by_lines = [
            [val1 + val2 for val1, val2 in zip(line1, line2)]
            for line1, line2 in zip(
                cube1_values_by_lines, cube2_values_by_lines
            )
        ]

        return resultant_values_by_lines

    def subtract_values(self):
        cube1_values_by_lines = self.cube1.values_by_lines
        cube2_values_by_lines = self.cube2.values_by_lines

        # Ensure both lists have the same structure
        self._check_value_lines()

        resultant_values_by_lines = [
            [val1 - val2 for val1, val2 in zip(line1, line2)]
            for line1, line2 in zip(
                cube1_values_by_lines, cube2_values_by_lines
            )
        ]

        return resultant_values_by_lines

    def write_results(self):
        if not self._all_checked():
            raise ValueError(
                "The geometries of the two cube files are different, please check carefully!"
            )
        if self.operation.lower() == "subtract":
            resultant_values = self.subtract_values()
        elif self.operation.lower() == "add":
            resultant_values = self.add_values()
        else:
            raise ValueError(f"Unknown operation given: {self.operation}!")

        with open(self.output_cubefile, "w") as f:
            self._write_header(f)  # 2 lines
            self._write_num_atoms(f)  # next 1 line
            self._write_grids(f)  # next 3 lines
            self._write_geometry(f)
            self._write_values_by_line(f, resultant_values)
        f.close()

    def _write_header(self, f):
        # write the first two lines
        f.write(self.cube1.cube_job_title + "\n")
        f.write(self.cube1.cube_job_description + "\n")

    def _write_num_atoms(self, f):
        # write the third line
        assert (
            self._check_natoms_matched()
            and self._check_coordinate_origin_matched()
        ), "Number of atoms and coordinate of the origin should be the same!"
        # f.write(f'{self.cube1.num_atoms:5} '
        #         f'{self.cube1.coordinate_origin[0]:12.6f} '
        #         f'{self.cube1.coordinate_origin[1]:12.6f} '
        #         f'{self.cube1.coordinate_origin[2]:12.6f} \n' )
        # instead of reconstructing as above, just write the 3rd line
        f.write(self.cube1.contents[2] + "\n")

    def _write_grids(self, f):
        assert (
            self._check_grid_points_matched()
        ), "Grid points should be the same for both cubes but are different!"
        for line in self.cube1.contents[3:6]:
            f.write(line + "\n")

    def _write_geometry(self, f):
        """Write the geometry that is the same. Instead of using structure.write(), we will simply write
        the geometry from the input cube file. This is due to the non-standard geometry specification in
        cube file, where the atomic_number is repeated as a float in the second element.
        """
        assert (
            self._check_geometries_matched()
        ), "Geometries should be the same for both cubes but are different!"
        for line in self.cube1.contents[6 : 6 + self.cube1.num_atoms]:
            f.write(line + "\n")

    def _write_values_by_line(self, f, resultant_values):
        for line in resultant_values:
            string = ""
            for x in line:
                string += f"{x:12.6}".upper()
                string += "  "
            f.write(f"{string}\n")
