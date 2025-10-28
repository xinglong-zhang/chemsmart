import os

import numpy as np

from chemsmart.utils.mixins import FileMixin


class GaussianCubeFile(FileMixin):
    """
    Parser for Gaussian cube files.
    """

    def __init__(self, filename):
        """
        Initialize cube file parser.
        """
        self.filename = filename

    @property
    def cube_job_title(self):
        """Get the job title from the first line of the cube file."""
        return self.contents[0]

    @property
    def cube_job_description(self):
        """
        Get the job description from the second line of the cube file.
        """
        return self.contents[1]

    @property
    def num_atoms(self):
        """
        Get number of atoms from the third line.
        """
        return int(self.contents[2].split()[0])

    @property
    def coordinate_origin(self):
        """
        Get the x-, y-, z-coordinates of the grid origin.
        """
        line_elements = self.contents[2].split()
        x = float(line_elements[1])
        y = float(line_elements[2])
        z = float(line_elements[3])
        return (x, y, z)

    @property
    def grid_points(self):
        """
        Get grid point counts for each dimension.
        """
        # Grid points information is specified in lines 4-6
        grid_points_data = []
        for line in self.contents[3:6]:
            grid_points_data.append(int(line.split()[0]))
        return tuple(grid_points_data)

    @property
    def grid_increment_vector(self):
        """
        Get increment vectors for the grid.
        """
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
        """
        Get coordinate block as list of strings.
        """
        return self.contents[6 : 6 + self.num_atoms]

    @property
    def coordinate_block_object(self):
        """
        Get coordinate block as CoordinateBlock object.
        """
        from chemsmart.io.molecules.structure import CoordinateBlock

        return CoordinateBlock(coordinate_block=self.coordinate_block_as_list)

    @property
    def chemical_symbols(self):
        """
        Get chemical symbols of all atoms.
        """
        return self.coordinate_block_object.chemical_symbols

    @property
    def positions(self):
        """
        Get atomic positions.
        """
        return self.coordinate_block_object.positions

    @property
    def structure(self):
        """
        Get molecular structure object.
        """
        return self.coordinate_block_object.molecule

    @property
    def values_by_lines(self):
        """
        Get volumetric data values organized by lines.
        """
        lines_of_values = []
        for line in self.contents[6 + self.num_atoms :]:
            # skip first 2 header lines + next 1 num atoms line + 3 grid vectors lines
            line_of_floats_as_list = [float(x) for x in line.split()]
            lines_of_values.append(line_of_floats_as_list)
        return lines_of_values


class CubeFileOperator:
    """Operator for mathematical operations between two cube files.

    Args:
        cubefile1 (str): Path to the first cube file
        cubefile2 (str): Path to the second cube file
        operation (str): Type of operation ('add' or 'subtract')
        output_cubefile (str, optional): Path for output file. If None,
            generates automatically based on first file and operation.
    """

    def __init__(
        self, cubefile1, cubefile2, operation="subtract", output_cubefile=None
    ):
        self.cubefile1 = cubefile1
        self.cubefile2 = cubefile2
        self.cube1 = GaussianCubeFile(filename=self.cubefile1)
        self.cube2 = GaussianCubeFile(filename=self.cubefile2)
        self.operation = operation

        # Generate output filename if not provided
        if output_cubefile is None:
            output_cubefile = os.path.join(
                self.cube1.filepath_directory,
                f"{self.cube1.basename}_{operation}.cube",
            )
        self.output_cubefile = output_cubefile

    def _check_natoms_matched(self):
        """
        Check if both cube files have the same number of atoms.
        """
        return self.cube1.num_atoms == self.cube2.num_atoms

    def _check_coordinate_origin_matched(self):
        """
        Check if both cube files have the same coordinate origin.
        """
        return all(
            np.isclose(
                self.cube1.coordinate_origin,
                self.cube2.coordinate_origin,
                atol=1e-4,
            )
        )

    def _check_grid_points_matched(self):
        """
        Check if both cube files have compatible grids.
        """
        grid_points_matched = self.cube1.grid_points == self.cube2.grid_points
        grid_increment_vectors_matched = (
            self.cube1.grid_increment_vector
            == self.cube2.grid_increment_vector
        )
        return grid_points_matched and grid_increment_vectors_matched

    def _check_geometries_matched(self):
        """
        Check if molecular geometries in both cube files are identical.
        """
        symbols_matched = (
            self.cube1.chemical_symbols == self.cube2.chemical_symbols
        )
        positions_are_close = np.isclose(
            self.cube1.positions, self.cube2.positions, atol=1e-3
        )
        positions_all_close = np.all(positions_are_close)

        return symbols_matched and positions_all_close

    def _all_checked(self):
        """
        Verify all compatibility checks pass.
        """
        return (
            self._check_natoms_matched()
            and self._check_coordinate_origin_matched()
            and self._check_grid_points_matched()
            and self._check_geometries_matched()
        )

    def _check_value_lines(self):
        """
        Validate that both cube files have compatible value data.
        """
        if len(self.cube1.values_by_lines) != len(self.cube2.values_by_lines):
            raise ValueError(
                f"The two cube files do not have the same number of lines.\n"
                f"Cube 1 has {len(self.cube1.values_by_lines)} lines and "
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
        """
        Perform element-wise addition of cube file values.
        """
        cube1_values_by_lines = self.cube1.values_by_lines
        cube2_values_by_lines = self.cube2.values_by_lines

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
        """
        Perform element-wise subtraction of cube file values.
        """
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
        """
        Write the result of the operation to a new cube file.
        """
        if not self._all_checked():
            raise ValueError(
                "The geometries of the two cube files are different, "
                "please check carefully!"
            )
        if self.operation.lower() == "subtract":
            resultant_values = self.subtract_values()
        elif self.operation.lower() == "add":
            resultant_values = self.add_values()
        else:
            raise ValueError(f"Unknown operation given: {self.operation}!")

        with open(self.output_cubefile, "w") as f:
            self._write_header(f)  # Write first 2 lines
            self._write_num_atoms(f)  # Write 3rd line
            self._write_grids(f)  # Write next 3 lines
            self._write_geometry(f)  # Write molecular geometry
            self._write_values_by_line(f, resultant_values)  # Write values
        f.close()

    def _write_header(self, f):
        """
        Write the header lines (title and description).
        """
        f.write(self.cube1.cube_job_title + "\n")
        f.write(self.cube1.cube_job_description + "\n")

    def _write_num_atoms(self, f):
        """
        Write the number of atoms and coordinate origin line.
        """
        assert (
            self._check_natoms_matched()
            and self._check_coordinate_origin_matched()
        ), (
            "Number of atoms and coordinate of the origin should be "
            "the same!"
        )

        # Write the original 3rd line instead of reconstructing
        f.write(self.cube1.contents[2] + "\n")

    def _write_grids(self, f):
        """
        Write the grid specification lines.
        """
        assert self._check_grid_points_matched(), (
            "Grid points should be the same for both cubes but are "
            "different!"
        )
        for line in self.cube1.contents[3:6]:
            f.write(line + "\n")

    def _write_geometry(self, f):
        """
        Write the molecular geometry section.

        This writes the geometry exactly as it appears in the input cube file
        to preserve the non-standard format where atomic number is repeated
        as a float in the second column.
        """
        assert self._check_geometries_matched(), (
            "Geometries should be the same for both cubes but are "
            "different!"
        )
        for line in self.cube1.contents[6 : 6 + self.cube1.num_atoms]:
            f.write(line + "\n")

    def _write_values_by_line(self, f, resultant_values):
        """
        Write the volumetric data values.
        """
        for line in resultant_values:
            string = ""
            for x in line:
                string += f"{x:12.6}".upper()
                string += "  "
            f.write(f"{string}\n")
