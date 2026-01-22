"""
Gaussian input file writer for computational chemistry jobs.

This module provides the GaussianInputWriter class for generating
properly formatted Gaussian input files from job configurations.
It handles all aspects of input file creation including route sections,
molecular coordinates, basis set specifications, and advanced options
like modredundant coordinates and custom solvent parameters.

The writer supports both standard job types and link jobs with
multiple calculation steps.
"""

import logging
import os.path

from chemsmart.jobs.gaussian.settings import GaussianLinkJobSettings
from chemsmart.jobs.writer import InputWriter
from chemsmart.utils.utils import (
    get_prepend_string_list_from_modred_free_format,
)

logger = logging.getLogger(__name__)


class GaussianInputWriter(InputWriter):
    """
    Writer class for generating Gaussian input files from job settings.

    Creates properly formatted Gaussian input files including route
    sections, molecular coordinates, basis set specifications, and
    advanced calculation options. Supports both standard calculations
    and complex multi-step link jobs.

    Handles all aspects of input file generation including modredundant
    coordinates, custom solvent parameters, and GenECP specifications.
    """

    def write(self, **kwargs):
        """
        Write the Gaussian input file for the job.

        Main entry point for input file generation. Delegates to
        internal write method with optional keyword arguments.

        Args:
            **kwargs: Additional arguments passed to _write method.
        """
        self._write(**kwargs)

    def _write(self, target_directory=None):
        """
        Internal method to write the Gaussian input file.

        Creates the output directory if needed and writes the complete
        input file based on job settings and molecule configuration.

        Args:
            target_directory (str, optional): Directory to write file to.
                If None, uses job's default folder.
        """
        if target_directory is not None:
            if not os.path.exists(target_directory):
                logger.debug(f"Creating target directory: {target_directory}")
                os.makedirs(target_directory)
            folder = target_directory
        else:
            folder = self.job.folder
        job_inputfile = os.path.join(folder, f"{self.job.label}.com")
        logger.debug(f"Writing Gaussian input file: {job_inputfile}")

        f = open(job_inputfile, "w")
        if self.settings.input_string:
            # write the file itself for direct run
            self._write_self(f)
        else:
            # Write complete input file from settings and molecule
            logger.debug("Generating input file from settings and molecule")
            self._write_all(f)
        logger.info(f"Finished writing Gaussian input file: {job_inputfile}")
        f.close()

    def _write_all(self, f):
        """
        Write complete Gaussian input file with all sections.

        Orchestrates writing of all input file sections in the proper
        order including header, route, coordinates, basis sets, and
        additional job-specific information.

        Args:
            f (file): Open file object to write to.
        """
        logger.debug("Starting complete input file generation")
        self._write_gaussian_header(f)
        self._write_route_section(f)
        self._write_gaussian_title(f)
        self._write_charge_and_multiplicity(f)
        self._write_cartesian_coordinates(f)

        # Skip modredundant for link jobs (handled in link section)
        if not isinstance(self.settings, GaussianLinkJobSettings):
            self._append_modredundant(f)

        self._append_gen_genecp_basis(f)  # then write genecp info
        self._append_custom_solvent_parameters(
            f
        )  # followed by user defined solvent parameters
        self._append_job_specific_info(f)
        self._append_other_additional_info(f)

        # Write link section for multi-step calculations
        if isinstance(self.settings, GaussianLinkJobSettings):
            self._write_link_section(f)
        logger.debug("Completed input file generation")

    def _write_self(self, f):
        """
        Write the input file content directly from settings string.

        Used when the job has a pre-defined input string that should
        be written directly to the file without any processing.

        Args:
            f (file): Open file object to write to.
        """
        f.write(self.settings.input_string)

    def _write_gaussian_header(self, f):
        """
        Write the Gaussian header section with resource specifications.

        Writes checkpoint file specifications, processor count, and
        memory allocation directives to the input file header.

        Args:
            f (file): Open file object to write to.
        """
        logger.debug("Writing Gaussian header.")
        if self.settings.chk:
            logger.debug(f"Writing chk file: {self.job.label}.chk")
            f.write(f"%chk={self.job.label}.chk\n")

        # Set default values if jobrunner resources are not specified
        num_cores = self.jobrunner.num_cores if not None else 12
        mem_gb = self.jobrunner.mem_gb if not None else 16
        logger.debug(
            f"Writing resource specifications - cores: {num_cores}, "
            f"memory: {mem_gb}GB"
        )
        f.write(f"%nprocshared={num_cores}\n")
        f.write(f"%mem={mem_gb}GB\n")

    def _write_route_section(self, f):
        """
        Write the Gaussian route section with calculation keywords.

        Constructs and writes the route line (#-line) containing all
        calculation specifications. Handles basis set adjustments
        for mixed heavy/light element calculations and determines
        appropriate gen/genecp keyword based on elements present.

        Args:
            f (file): Open file object to write to.
        """
        logger.debug("Writing route section.")
        route_string = self.settings.route_string

        # Handle basis set replacement for structures without heavy elements
        # if project settings has heavy elements but molecule has no heavy
        # elements, then replace the basis set with light elements basis
        if self.settings.heavy_elements_basis is not None:
            heavy_elements_in_structure = self.settings.prune_heavy_elements(
                self.job.molecule
            )
            if (
                heavy_elements_in_structure is None
                # returns None if no heavy elements given in settings
                or len(heavy_elements_in_structure) == 0
                # returns empty list if no heavy elements found in structure
                # (heavy elements specified in settings)
            ):
                # first remove any '-' in light_elements_basis
                # this is because '-' is needed by bse package to get the right
                # basis set from https://www.basissetexchange.org/, eg, def2-SVP
                # but this is not needed in Gaussian input file
                # (will cause error when run), eg, def2svp
                light_elements_basis = (
                    self.settings.light_elements_basis.replace("-", "").lower()
                )

                route_string = route_string.replace(
                    self.settings.basis,
                    light_elements_basis,
                )
            else:
                # Determine the correct basis keyword (gen vs genecp) based on
                # the heavy elements actually present in the molecule
                determined_basis = self.settings.determine_basis_keyword(
                    self.job.molecule
                )
                if determined_basis != self.settings.basis:
                    logger.info(
                        f"Replacing basis keyword '{self.settings.basis}' with "
                        f"'{determined_basis}' based on heavy elements in molecule"
                    )
                    route_string = route_string.replace(
                        self.settings.basis,
                        determined_basis,
                    )
        f.write(route_string + "\n")
        f.write("\n")

    def _write_gaussian_title(self, f):
        """
        Write the job title section to the input file.

        Writes the descriptive title for the calculation which appears
        in the output file and helps identify the job purpose.

        Args:
            f (file): Open file object to write to.
        """
        logger.debug("Writing Gaussian title.")
        title = self.settings.title
        f.write(f"{title}\n")
        f.write("\n")

    def _write_charge_and_multiplicity(self, f):
        """
        Write molecular charge and spin multiplicity to input file.

        Writes the electronic configuration specification required
        for all Gaussian calculations. Validates that both values
        are properly defined.

        Args:
            f (file): Open file object to write to.

        Raises:
            AssertionError: If charge or multiplicity are not specified.
        """
        logger.debug("Writing charge and multiplicity.")
        charge = self.settings.charge
        multiplicity = self.settings.multiplicity
        assert (
            charge is not None and multiplicity is not None
        ), "Charge and multiplicity must be specified!"
        logger.debug(
            f"Molecular charge: {charge}, multiplicity: {multiplicity}"
        )
        f.write(f"{charge} {multiplicity}\n")

    def _write_cartesian_coordinates(self, f):
        """
        Write molecular Cartesian coordinates to the input file.

        Outputs the molecular geometry in Cartesian coordinate format
        suitable for Gaussian calculations. Validates that molecular
        structure is available.

        Args:
            f (file): Open file object to write to.

        Raises:
            AssertionError: If no molecular geometry is found.
        """
        logger.debug(
            f"Writing Cartesian coordinates of molecule: {self.job.molecule}."
        )
        assert self.job.molecule is not None, "No molecular geometry found!"

        # Log molecular information for debugging
        logger.debug(
            f"Molecule contains {self.job.molecule.num_atoms} atoms with formula: "
            f"{self.job.molecule.chemical_formula}"
        )

        self.job.molecule.write_coordinates(f, program="gaussian")
        f.write("\n")

    def _append_modredundant(self, f):
        """
        Write modredundant coordinate specifications to input file.

        Handles both constraint optimization (list format) and
        coordinate scanning (dictionary format) by writing appropriate
        modredundant sections. Supports various coordinate types.

        Args:
            f (file): Open file object to write to.

        Raises:
            ValueError: If modredundant format is invalid.
        """
        logger.debug("Writing modred section.")
        modredundant = self.settings.modred
        if modredundant is not None:
            if isinstance(modredundant, list):
                # for modredunant job
                # modred = [[1,2], [3,4]]
                logger.debug("Writing modredundant constraints")
                prepend_string_list = (
                    get_prepend_string_list_from_modred_free_format(
                        input_modred=modredundant
                    )
                )
                for prepend_string in prepend_string_list:
                    f.write(f"{prepend_string} F\n")
                f.write("\n")
            elif isinstance(modredundant, dict):
                # For coordinate scanning job
                # modred = {'num_steps': 10, 'step_size': 0.05,
                #           'coords': [[1,2], [3,4]]}
                logger.debug("Writing coordinate scan specifications")
                coords_list = modredundant["coords"]
                scan_prepend_string_list = (
                    get_prepend_string_list_from_modred_free_format(
                        input_modred=coords_list
                    )
                )
                for prepend_string, num_step, step_size in zip(
                    scan_prepend_string_list,
                    self.settings.modred["num_steps"],
                    self.settings.modred["step_size"],
                ):
                    f.write(f"{prepend_string} S {num_step} {step_size}\n")
                if "constrained_coordinates" in modredundant:
                    logger.debug(
                        "Writing modredundant constraints for scan job"
                    )
                    constrained_coordinates_list = modredundant[
                        "constrained_coordinates"
                    ]
                    prepend_string_list = (
                        get_prepend_string_list_from_modred_free_format(
                            input_modred=constrained_coordinates_list
                        )
                    )
                    for prepend_string in prepend_string_list:
                        f.write(f"{prepend_string} F\n")
                f.write("\n")
            else:
                logger.error(
                    f"Invalid modredundant format: {type(modredundant)}"
                )
                raise ValueError(
                    "modredundant must be a list or dictionary with "
                    "'num_steps', 'step_size', and 'coords'."
                )

    def _append_gen_genecp_basis(self, f):
        """
        Write GenECP basis set and pseudopotential information.

        Appends mixed basis set specifications and effective core
        potentials for calculations involving heavy elements. Handles
        automatic basis set selection when no heavy elements are present.

        Args:
            f (file): Open file object to write to.
        """
        logger.debug("Writing gen/genecp basis set information.")
        if self.settings.genecp:
            genecp_section = self.settings.get_genecp_section(
                molecule=self.job.molecule
            )

            # Write genecp section only if heavy elements are in the molecule
            heavy_elements_in_structure = self.settings.prune_heavy_elements(
                self.job.molecule
            )

            # For the given project setting, if heavy elements are not found
            # in the structure, and no user-specified gen_genecp_file is
            # supplied, then use light elements basis, e.g., organic reactant
            # molecules in TM catalysis
            if (
                heavy_elements_in_structure is None
                or len(heavy_elements_in_structure) == 0
            ) and self.settings.gen_genecp_file is None:
                # Replace gen or genecp keyword in route by light elements basis
                logger.debug(
                    "No heavy elements found - using light elements basis only"
                )
                pass

            else:
                logger.debug("Writing GenECP section for heavy elements")
                f.write(genecp_section.string)
                # Check that the last line of genecp_section.string is empty,
                # if not, add an empty line
                if genecp_section.string_list[-1] != "\n":
                    f.write("\n")

    def _append_custom_solvent_parameters(self, f):
        """
        Write custom solvent parameter specifications to input file.

        Appends user-defined solvent parameters for advanced solvation
        calculations that require custom dielectric constants or
        other solvent-specific properties.

        Args:
            f (file): Open file object to write to.
        """
        logger.debug("Writing custom solvent parameters.")
        custom_solvent = self.settings.custom_solvent
        if custom_solvent is not None:
            logger.debug("Adding custom solvent specification")
            f.write(custom_solvent)
            f.write("\n")

    def _append_job_specific_info(self, f):
        """
        Write job-type-specific information to the input file.

        Appends calculation-specific data required for specialized
        job types like NCI analysis, WBI calculations, and RESP
        charge fitting. Each job type has unique requirements.

        Args:
            f (file): Open file object to write to.
        """
        logger.debug("Writing job specific information.")
        job_type = self.settings.job_type
        job_label = self.job.label
        if job_type == "nci":
            # Appending wavefunction file specification for NCI job
            logger.debug("Adding NCI-specific wavefunction file")
            f.write(f"{job_label}.wfn\n")
            f.write("\n")
        elif job_type == "wbi":
            # Appending NBO directive for WBI job
            logger.debug("Adding WBI-specific NBO directive")
            f.write("$nbo bndidx $end\n")
            f.write("\n")
        elif job_type == "resp":
            # Appending electrostatic potential file for RESP job
            logger.debug("Adding RESP-specific potential file")
            f.write(f"{job_label}.gesp\n")
            f.write("\n")

    def _append_other_additional_info(self, f):
        """
        Write additional user-specified information to input file.

        Appends custom information that can be provided either as
        a file path or as direct string content. Supports flexible
        input file customization for advanced users.

        Args:
            f (file): Open file object to write to.
        """
        logger.debug("Writing other additional information.")
        append_additional_info = self.settings.append_additional_info
        if append_additional_info is not None:
            if isinstance(append_additional_info, str) and os.path.exists(
                os.path.expanduser(append_additional_info)
            ):
                # Path to the file for additional append info
                logger.debug(
                    f"Reading additional info from file: "
                    f"{append_additional_info}"
                )
                with open(append_additional_info) as g:
                    for line in g.readlines():
                        f.write(line)
            else:
                # Ensures that the additional append info can be supplied
                # in free string format
                logger.debug("Writing additional info from string")
                line_elem = append_additional_info.strip().split("\n")
                for line in line_elem:
                    f.write(f"{line}\n")

    def _write_link_section(self, f):
        """
        Write the complete link section for multi-step calculations.

        Creates the second part of a link job by writing the link
        separator and all necessary sections for the continuation
        calculation including headers, routes, and coordinates.

        Args:
            f (file): Open file object to write to.
        """
        logger.debug("Writing link section.")
        if self.settings.link:
            logger.debug("Adding link job separator and second calculation")
            f.write("--Link1--\n")
            self._write_gaussian_header(f)
            self._write_link_route(f)
            self._write_gaussian_title(f)
            self._write_charge_and_multiplicity(f)
            f.write("\n")
            self._append_modredundant(f)
            self._append_gen_genecp_basis(f)  # then write genecp info
            self._append_custom_solvent_parameters(
                f
            )  # followed by user defined solvent parameters
            self._append_job_specific_info(f)
            self._append_other_additional_info(f)

    def _write_link_route(self, f):
        """
        Write the route section for the link calculation step.

        Writes the route line for the second part of a link job,
        which typically includes geometry and orbital continuation
        from the previous calculation step.

        Args:
            f (file): Open file object to write to.
        """
        logger.debug("Writing link route section.")
        f.write(f"{self.settings.link_route_string}\n")
        f.write("\n")
