import logging
import os
import re

from chemsmart.io.molecules.structure import CoordinateBlock, Molecule
from chemsmart.io.orca.route import ORCARoute
from chemsmart.utils.mixins import ORCAFileMixin
from chemsmart.utils.repattern import (
    standard_coord_pattern,
    xyz_filename_pattern,
)

logger = logging.getLogger(__name__)


class ORCAInput(ORCAFileMixin):
    """
    Parser for ORCA quantum chemistry input files.

    This class provides comprehensive parsing capabilities for ORCA input files,
    extracting molecular coordinates, calculation parameters, job settings, and
    other computational directives. It handles both embedded coordinates and
    external coordinate file references.

    Args:
        filename (str): Path to the ORCA input file
    """

    def __init__(self, filename):
        """
        Initialize ORCA input file parser.

        Args:
            filename (str): Path to the ORCA input file to parse
        """
        self.filename = filename

        try:
            cb = CoordinateBlock(coordinate_block=self.coordinate_lines)
            self.cb = cb
        except ValueError as err:
            logger.error(f"Error parsing coordinate block: {err}")
            self.cb = None

    @property
    def route_string(self):
        """
        Extract and format the ORCA route string from input file.

        Combines all lines starting with '!' to form the complete route
        specification for the ORCA calculation.

        Returns:
            str: Complete route string in lowercase format
        """
        route_string_lines = []
        for line in self.contents:
            if line.startswith("!"):
                new_line = line.replace("!", "")
                route_string_lines.append(new_line)

        route_string = " ".join(route_string_lines)
        return f"! {route_string}".lower()

    @property
    def route_object(self):
        """
        Create ORCARoute object from the route string.

        Returns:
            ORCARoute: Parsed route object containing calculation methods
                      and keywords, or None if parsing fails
        """
        try:
            route_object = ORCARoute(route_string=self.route_string)
            return route_object
        except TypeError as err:
            print(err)

    @property
    def coordinate_type(self):
        """
        Extract coordinate specification type from input file.

        Returns:
            str: Coordinate type (e.g., 'xyz', 'xyzfile') or None if not found
        """
        for line in self.contents:
            if line.startswith("*") and len(line) > 1:
                line_elem = line.split()
                return line_elem[1]
        return None

    @property
    def charge(self):
        """
        Extract molecular charge from coordinate specification.

        Returns:
            int: Molecular charge or None if not specified
        """
        for line in self.contents:
            if line.startswith("*") and len(line) > 1:
                line_elem = line.split()
                return int(line_elem[2])
        return None

    @property
    def multiplicity(self):
        """
        Extract spin multiplicity from coordinate specification.

        Returns:
            int: Spin multiplicity or None if not specified
        """
        for line in self.contents:
            if line.startswith("*") and len(line) > 1:
                line_elem = line.split()
                return int(line_elem[3])
        return None

    @property
    def coordinate_lines(self):
        """
        Extract coordinate lines from input file.

        Uses regex pattern matching to identify lines containing atomic
        coordinates in standard format.

        Returns:
            list: Lines containing atomic coordinate specifications
        """
        pattern = re.compile(standard_coord_pattern)
        coordinate_lines = []
        for line in self.contents:
            if pattern.match(line):
                coordinate_lines.append(line)
        return coordinate_lines

    @property
    def molecule(self):
        """
        Create Molecule object from coordinate data.

        Attempts to create a Molecule object from embedded coordinates or
        external coordinate files. Sets charge and multiplicity from the
        input file specification.

        Returns:
            Molecule: Molecular structure object with coordinates, charge,
                     and spin multiplicity, or None if parsing fails

        Raises:
            FileNotFoundError: If external coordinate file is not found
        """
        molecule = None
        try:
            molecule = self.cb.molecule
        except ValueError as err:
            logger.debug(
                f"Error creating molecule from coordinate block: {err}"
            )
            for line in self.contents:
                if line.startswith("* xyzfile"):
                    xyz_file = line.strip().split()[-1]
                    xyz_filepath = os.path.join(
                        self.filepath_directory, xyz_file
                    )
                    if os.path.exists(xyz_filepath):
                        molecule = Molecule.from_filepath(
                            filepath=xyz_filepath
                        )
                    else:
                        raise FileNotFoundError(
                            f"Coordinate file {xyz_filepath} not found."
                        )
            # update charge and spin multiplicity
        if molecule:
            molecule.charge = self.charge
            molecule.spin_multiplicity = self.multiplicity
        return molecule

    @property
    def scf_maxiter(self):
        """
        Extract SCF maximum iteration count from input file.

        Searches for 'maxiter' keyword in the input file to determine
        the maximum number of SCF iterations allowed.

        Returns:
            int: Maximum SCF iterations or None if not specified
        """
        for i, raw_line in enumerate(self.contents):
            line = raw_line.lower()
            if "maxiter" in line:
                try:
                    # Find the index of "maxiter" in the same line
                    index = line.index("maxiter")
                    # Find the substring immediately following "maxiter"
                    num_maxiter = line[index + len("maxiter") :].split()[0]
                    return int(num_maxiter)
                except ValueError:
                    # If "maxiter" is not found, search the next line for it
                    next_lines = self.contents[i + 1 :]
                    for line in next_lines:
                        if "maxiter" in line:
                            # Find the index of "maxiter" in the same line
                            index = line.index("maxiter")
                            # Find the substring immediately following "maxiter"
                            num_maxiter = line[
                                index + len("maxiter") :
                            ].split()[0]
                            return int(num_maxiter)
        return None

    @property
    def scf_convergence(self):
        """
        Extract SCF convergence criteria from %scf block.

        Returns:
            str: SCF convergence criterion (e.g., 'tight', 'loose')
                 or None if not specified
        """
        for i, line in enumerate(self.contents):
            if "%scf" not in line:
                continue

            if "convergence" in line:
                # Find the index of "convergence" in the same line
                index = line.index("convergence")
                # Find the substring immediately following "convergence"
                return line[index + len("convergence") :].split()[0]

            # If "convergence" is not found, search the next line for it
            next_lines = self.contents[i + 1 :]
            for next_line in next_lines:
                if "convergence" in next_line:
                    # Find the index of "convergence" in the same line
                    index = next_line.index("convergence")
                    # Find the substring immediately following "convergence"
                    return next_line[index + len("convergence") :].split()[0]
        return None

    @property
    def dipole(self):
        """
        Extract dipole moment calculation specification from elprop block.

        Returns:
            str: Dipole calculation setting or None if not specified
        """
        # in elprop block
        for line in self.contents:
            if "dipole" in line:
                return line.split()[-1]
        return None

    @property
    def quadrupole(self):
        """
        Extract quadrupole moment calculation specification from elprop block.

        Returns:
            str: Quadrupole calculation setting or None if not specified
        """
        # in elprop block
        for line in self.contents:
            if "quadrupole" in line:
                return line.split()[-1]
        return None


class ORCANEBInput(ORCAInput):
    @property
    def starting_xyzfile(self):
        starting_xyzfile, _, _, _ = self._get_geometries()
        return starting_xyzfile

    @property
    def ending_xyzfile(self):
        _, ending_xyzfile, _, _ = self._get_geometries()
        return ending_xyzfile

    @property
    def ts_xyzfile(self):
        _, _, ts_xyzfile, _ = self._get_geometries()
        return ts_xyzfile

    @property
    def restarting_allxyzfile(self):
        _, _, _, restarting_allxyzfile = self._get_geometries()
        return restarting_allxyzfile

    @property
    def nimages(self):
        return self._get_number_of_images()

    @property
    def pre_optimization(self):
        return self._get_pre_optimization_status()

    def _get_geometries(self):
        neb_starting_xyz = neb_end_xyzile = neb_ts_xyzile = (
            restart_allxyzfile
        ) = None
        for line in self.contents:
            match = re.search(xyz_filename_pattern, line)
            if "neb_end_xyzfile" in line.lower():
                neb_end_xyzile = match.group(1)
            elif "neb_ts_xyzfile" in line.lower():
                neb_ts_xyzile = match.group(1)
            elif "restart_allxyzfile" in line.lower():
                restart_allxyzfile = match.group(1)
            elif line.startswith("* xyzfile"):
                neb_starting_xyz = match.group(1)
            elif self.coordinate_lines:
                neb_starting_xyz = self.coordinate_lines
        return (
            neb_starting_xyz,
            neb_end_xyzile,
            neb_ts_xyzile,
            restart_allxyzfile,
        )

    def _get_number_of_images(self):
        nimages = None
        for line in self.contents:
            if "nimages" in line.lower():
                nimages = int(line.split()[-1])
        return nimages

    def _get_pre_optimization_status(self):
        pre_opt = False
        for line in self.contents:
            line = line.lower()
            if "preopt_ends" in line and "true" in line:
                pre_opt = True
        return pre_opt
