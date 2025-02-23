import logging
import os.path

from chemsmart.jobs.gaussian.settings import GaussianLinkJobSettings
from chemsmart.jobs.writer import InputWriter
from chemsmart.utils.utils import (
    get_prepend_string_list_from_modred_free_format,
)
from chemsmart.utils.utils import str_indices_range_to_list
logger = logging.getLogger(__name__)


class GaussianInputWriter(InputWriter):
    """Class that writes Gaussian input files for a job."""

    def write(self, **kwargs):
        self._write(**kwargs)

    def _write(self, target_directory=None):
        if target_directory is not None:
            if not os.path.exists(target_directory):
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
            self._write_all(f)
        logger.info(f"Finished writing Gaussian input file: {job_inputfile}")
        f.close()

    def _write_all(self, f):
        self._write_gaussian_header(f)
        self._write_route_section(f)
        self._write_gaussian_title(f)
        self._write_charge_and_multiplicity(f)
        self._write_cartesian_coordinates(f)
        if not isinstance(self.settings, GaussianLinkJobSettings):
            self._append_modredundant(f)
        self._append_gen_genecp_basis(f)  # then write genecp info
        self._append_custom_solvent_parameters(
            f
        )  # followed by user defined solvent parameters
        self._append_job_specific_info(f)
        self._append_other_additional_info(f)
        if isinstance(self.settings, GaussianLinkJobSettings):
            self._write_link_section(f)

    def _write_self(self, f):
        """Write the input file itself for direct run."""
        f.write(self.settings.input_string)

    def _write_gaussian_header(self, f):
        logger.debug("Writing Gaussian header.")
        if self.settings.chk:
            logger.debug(f"Writing chk file: {self.job.label}.chk")
            f.write(f"%chk={self.job.label}.chk\n")
        num_cores = self.jobrunner.num_cores if not None else 12
        mem_gb = self.jobrunner.mem_gb if not None else 16
        logger.debug(f"Writing nprocshared={num_cores} and mem={mem_gb}GB.")
        f.write(f"%nprocshared={num_cores}\n")
        f.write(f"%mem={mem_gb}GB\n")

    def _write_route_section(self, f):
        logger.debug("Writing route section.")
        route_string = self.settings.route_string
        # if project settings has heavy elements but molecule has no heavy elements,
        # then replace the basis set with light elements basis
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
        f.write(route_string + "\n")
        f.write("\n")

    def _write_gaussian_title(self, f):
        logger.debug("Writing Gaussian title.")
        title = self.settings.title
        f.write(f"{title}\n")
        f.write("\n")

    def _write_charge_and_multiplicity(self, f):
        logger.debug("Writing charge and multiplicity.")
        charge = self.settings.charge
        multiplicity = self.settings.multiplicity
        assert (
            charge is not None and multiplicity is not None
        ), "Charge and multiplicity must be specified!"
        f.write(f"{charge} {multiplicity}\n")

    def _write_cartesian_coordinates(self, f):
        logger.debug(
            f"Writing Cartesian coordinates of molecule: {self.job.molecule}."
        )
        assert self.job.molecule is not None, "No molecular geometry found!"
        self.job.molecule.write_coordinates(f, program="gaussian")
        f.write("\n")

    def _append_modredundant(self, f):
        """Write modred section if present in the job settings.
        Given a list or list of lists, write the modred section for the constraints.
        Given a dictionary with 'num_steps', 'step_size', and 'coords', write the scan section.
        """
        logger.debug("Writing modred section.")
        modredundant = self.settings.modred
        if modredundant is not None:
            if isinstance(modredundant, list):
                # for modredunant job
                # modred = [[1,2], [3,4]]
                prepend_string_list = (
                    get_prepend_string_list_from_modred_free_format(
                        input_modred=modredundant
                    )
                )
                for prepend_string in prepend_string_list:
                    f.write(f"{prepend_string} F\n")
                f.write("\n")
            elif isinstance(modredundant, dict):
                # for scanning job
                # modred = {'num_steps': 10, 'step_size': 0.05, 'coords': [[1,2], [3,4]]}
                coords_list = modredundant["coords"]
                prepend_string_list = (
                    get_prepend_string_list_from_modred_free_format(
                        input_modred=coords_list
                    )
                )
                for prepend_string in prepend_string_list:
                    f.write(
                        f"{prepend_string} S {self.settings.modred['num_steps']} {self.settings.modred['step_size']}\n"
                    )
                f.write("\n")
            else:
                raise ValueError(
                    "modredundant must be a list or dictionary with 'num_steps', 'step_size', and 'coords'."
                )

    def _append_gen_genecp_basis(self, f):
        """Write the genecp basis set information if present in the job settings."""
        logger.debug("Writing gen/genecp basis set information.")
        if self.settings.genecp:
            genecp_section = self.settings.get_genecp_section(
                molecule=self.job.molecule
            )

            # write genecp section only if heavy elements are in the molecule
            heavy_elements_in_structure = self.settings.prune_heavy_elements(
                self.job.molecule
            )

            # for the given project setting, if heavy elements are not found in the structure,
            # and no user-specified gen_genecp_file is supplied, then
            # use light elements basis, e.g., organic reactant molecules in TM catalysis
            if (
                heavy_elements_in_structure is None
                or len(heavy_elements_in_structure) == 0
            ) and self.settings.gen_genecp_file is None:
                # replace gen or genecp keyword in route by light elements basis
                pass

            else:
                f.write(genecp_section.string)
                # check that the last line of genecp_section.string is empty,
                # if not, add an empty line
                if genecp_section.string_list[-1] != "\n":
                    f.write("\n")

    def _append_custom_solvent_parameters(self, f):
        """Write the custom solvent parameters if present in the job settings."""
        logger.debug("Writing custom solvent parameters.")
        custom_solvent = self.settings.custom_solvent
        if custom_solvent is not None:
            f.write(custom_solvent)
            f.write("\n")

    def _append_job_specific_info(self, f):
        """Write any job specific information that needs to be appended to the input file."""
        logger.debug("Writing job specific information.")
        job_type = self.settings.job_type
        job_label = self.job.label
        if job_type == "nci":
            # appending for nci job
            f.write(f"{job_label}.wfn\n")
            f.write("\n")
        elif job_type == "wbi":
            # appending for wbi job
            f.write("$nbo bndidx $end\n")
            f.write("\n")
        elif job_type == "resp":
            # appending for resp job
            f.write(f"{job_label}.gesp\n")
            f.write("\n")

    def _append_other_additional_info(self, f):
        """Write any additional information that needs to be appended to the input file."""
        logger.debug("Writing other additional information.")
        append_additional_info = self.settings.append_additional_info
        if append_additional_info is not None:
            if isinstance(append_additional_info, str) and os.path.exists(
                os.path.expanduser(append_additional_info)
            ):
                # path to the file for additional append info
                with open(append_additional_info) as g:
                    for line in g.readlines():
                        f.write(line)
            else:
                # ensures that the additional append info can be supplied in free string format
                line_elem = append_additional_info.strip().split("\n")
                for line in line_elem:
                    f.write(f"{line}\n")

    def _write_link_section(self, f):
        """Write the link section for the input file."""
        logger.debug("Writing link section.")
        if self.settings.link:
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
        """Write the route section for the link job."""
        logger.debug("Writing link route section.")
        f.write(f"{self.settings.link_route_string}\n")
        f.write("\n")
