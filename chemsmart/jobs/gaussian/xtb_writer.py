"""
Input file writer for Gaussian XTB jobs.

This module provides specialized writers for Gaussian jobs that use
xtb as an external calculator. Handles the unique input file format
requirements for QST2 and other multi-structure optimizations.
"""

import logging

from chemsmart.jobs.gaussian.writer import GaussianInputWriter

logger = logging.getLogger(__name__)


class GaussianXTBInputWriter(GaussianInputWriter):
    """
    Input file writer for Gaussian XTB external calculator jobs.

    Extends GaussianInputWriter to handle the specific requirements
    of xtb-based calculations, including proper external keyword
    formatting and suppression of basis set specifications.
    """

    def _write_all(self, f):
        """
        Write complete Gaussian XTB input file.

        Writes the input file without GenECP or basis set sections
        since these are handled by the external xtb calculator.

        Args:
            f (file): Open file object to write to.
        """
        logger.debug("Starting XTB input file generation")
        self._write_gaussian_header(f)
        self._write_route_section(f)
        self._write_gaussian_title(f)
        self._write_charge_and_multiplicity(f)
        self._write_cartesian_coordinates(f)

        # Skip GenECP and basis set sections for external calculator
        # The external calculator handles its own energy/gradient method

        self._append_other_additional_info(f)
        logger.debug("Completed XTB input file generation")


class GaussianXTBQST2Writer(GaussianInputWriter):
    """
    Input file writer for Gaussian QST2 calculations with xtb.

    Handles the special input format required for QST2 calculations
    which need both reactant and product structures specified in
    the same input file.
    """

    def _write_all(self, f):
        """
        Write complete QST2 input file with both structures.

        QST2 requires both reactant and product structures to be
        specified in sequence, each with its own title and
        charge/multiplicity specification.

        Args:
            f (file): Open file object to write to.
        """
        logger.debug("Starting QST2 input file generation")
        self._write_gaussian_header(f)
        self._write_route_section(f)

        # Write reactant section
        self._write_gaussian_title(f, title="Reactants")
        self._write_charge_and_multiplicity(f)
        self._write_cartesian_coordinates(f, molecule=self.job.reactant)

        # Write product section
        self._write_gaussian_title(f, title="Product")
        self._write_charge_and_multiplicity(f)
        self._write_cartesian_coordinates(f, molecule=self.job.product)

        # Optional: Write TS guess if provided (for QST3 jobs)
        ts_guess = getattr(self.job, 'ts_guess', None)
        if ts_guess is not None:
            self._write_gaussian_title(f, title="TS Guess")
            self._write_charge_and_multiplicity(f)
            self._write_cartesian_coordinates(f, molecule=ts_guess)

        self._append_other_additional_info(f)
        logger.debug("Completed QST2 input file generation")

    def _write_gaussian_title(self, f, title=None):
        """
        Write the job title section.

        Args:
            f (file): Open file object to write to.
            title (str, optional): Title to write. If None, uses settings.
        """
        if title is None:
            title = self.settings.title
        f.write(f"{title}\n")
        f.write("\n")

    def _write_cartesian_coordinates(self, f, molecule=None):
        """
        Write molecular Cartesian coordinates.

        Args:
            f (file): Open file object to write to.
            molecule (Molecule, optional): Molecule to write.
                If None, uses job's primary molecule.
        """
        if molecule is None:
            molecule = self.job.molecule

        logger.debug("Writing coordinates for molecule: %s", molecule)
        molecule.write_coordinates(f, program="gaussian")
        f.write("\n")
