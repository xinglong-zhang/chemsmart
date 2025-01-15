import logging
import os.path

from chemsmart.jobs.gaussian.settings import GaussianLinkJobSettings
from chemsmart.jobs.writer import InputWriter
from chemsmart.io.gaussian.gengenecp import GenGenECPSection
from chemsmart.utils.utils import (
    get_prepend_string_list_from_modred_free_format,
    prune_list_of_elements,
)

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
        if self.job.settings.input_string:
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
        self._write_pbc(f)
        self._append_modredundant(f)
        self._append_gen_genecp_basis(f)  # then write genecp info
        self._append_custom_solvent_parameters(
            f
        )  # followed by user defined solvent parameters
        self._append_job_specific_info(f)
        self._append_other_additional_info(f)
        if isinstance(self.job.settings, GaussianLinkJobSettings):
            self._write_link_section(f)

    def _write_self(self, f):
        """Write the input file itself for direct run."""
        f.write(self.job.settings.input_string)

    def _write_gaussian_header(self, f):
        logger.debug("Writing Gaussian header.")
        if self.job.settings.chk:
            logger.debug(f"Writing chk file: {self.job.label}.chk")
            f.write(f"%chk={self.job.label}.chk\n")
        num_cores = self.jobrunner.num_cores if not None else 12
        mem_gb = self.jobrunner.mem_gb if not None else 16
        logger.debug(f"Writing nprocshared={num_cores} and mem={mem_gb}GB.")
        f.write(f"%nprocshared={num_cores}\n")
        f.write(f"%mem={mem_gb}GB\n")

    def _write_route_section(self, f):
        logger.debug("Writing route section.")
        route_string = self.job.settings.route_string
        f.write(route_string + "\n")
        f.write("\n")

    def _write_gaussian_title(self, f):
        logger.debug("Writing Gaussian title.")
        title = self.job.settings.title
        f.write(f"{title}\n")
        f.write("\n")

    def _write_charge_and_multiplicity(self, f):
        logger.debug("Writing charge and multiplicity.")
        charge = self.job.settings.charge
        multiplicity = self.job.settings.multiplicity
        f.write(f"{charge} {multiplicity}\n")

    def _write_cartesian_coordinates(self, f):
        logger.debug("Writing Cartesian coordinates.")
        assert self.job.molecule is not None, "No molecular geometry found!"
        self.job.molecule.write_coordinates(f)
        if not self.job.molecule.pbc_conditions:
            f.write("\n")

    def _write_pbc(self, f):
        """Write the periodic boundary conditions section if present in the job settings."""
        logger.debug("Writing periodic boundary conditions.")
        if self.job.molecule.pbc_conditions:
            for translation_vector in self.job.molecule.translation_vectors:
                f.write(
                    f"TV   "
                    f"{translation_vector[0]:15.10f} "
                    f"{translation_vector[1]:15.10f} "
                    f"{translation_vector[2]:15.10f}\n"
                )
            f.write("\n")

    def _append_modredundant(self, f):
        """Write modred section if present in the job settings.
        Given a list or list of lists, write the modred section for the constraints.
        Given a dictionary with 'num_steps', 'step_size', and 'coords', write the scan section.
        """
        logger.debug("Writing modred section.")
        modredundant = self.job.settings.modred
        if modredundant:
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
                        f"{prepend_string} S {self.job.settings.modred['num_steps']} {self.job.settings.modred['step_size']}\n"
                    )
                f.write("\n")

    def _append_gen_genecp_basis(self, f):
        """Write the genecp basis set information if present in the job settings."""
        logger.debug("Writing gen/genecp basis set information.")
        gen_genecp_file = self.job.settings.gen_genecp_file
        heavy_elements = self.job.settings.heavy_elements
        heavy_elements_basis = self.job.settings.heavy_elements_basis
        light_elements = self.job.settings.get_light_elements(
            self.job.molecule
        )
        light_elements_basis = self.job.settings.light_elements_basis
        if self.job.settings.genecp:
            if gen_genecp_file and heavy_elements and heavy_elements_basis:
                raise ValueError(
                    "Please provide either gen_genecp_file or heavy_elements and heavy_elements_basis."
                )
            if gen_genecp_file:
                assert os.path.exists(
                    gen_genecp_file
                ), f"File {gen_genecp_file} not found!"
                genecp_section = GenGenECPSection.from_genecp_path(
                    gen_genecp_file
                )
            elif heavy_elements and heavy_elements_basis:
                logger.info(
                    f"GENECP elements specified:\n"
                    f"Heavy elements: {heavy_elements}\n"
                    f"Heavy elements basis: {heavy_elements_basis}\n"
                    f"Light elements: {light_elements}\n"
                    f"Light elements basis: {light_elements_basis}\n"
                    "Using basis set exchange api to get gen/genecp basis set for heavy atoms.\n"
                )
                heavy_elements_in_structure = prune_list_of_elements(
                    heavy_elements, self.job.molecule
                )

                genecp_section = GenGenECPSection.from_bse_api(
                    light_elements=light_elements,
                    light_elements_basis=light_elements_basis,
                    heavy_elements=heavy_elements_in_structure,
                    heavy_elements_basis=heavy_elements_basis,
                )
            else:
                return
            f.write(genecp_section.string)
            f.write("\n")

    def _append_custom_solvent_parameters(self, f):
        """Write the custom solvent parameters if present in the job settings."""
        logger.debug("Writing custom solvent parameters.")
        custom_solvent = self.job.settings.custom_solvent
        if custom_solvent:
            f.write(custom_solvent)
            f.write("\n")

    def _append_job_specific_info(self, f):
        """Write any job specific information that needs to be appended to the input file."""
        logger.debug("Writing job specific information.")
        job_type = self.job.settings.job_type
        job_label = self.job.label
        if job_type == "nci":
            # appending for nci job
            f.write(f"{job_label}.wfn\n")
            f.write(f"\n")
        elif job_type == "wbi":
            # appending for wbi job
            f.write("$nbo bndidx $end\n")
            f.write(f"\n")
        elif job_type == "resp":
            # appending for resp job
            f.write(f"{job_label}.gesp\n")
            f.write(f"\n")

    def _append_other_additional_info(self, f):
        """Write any additional information that needs to be appended to the input file."""
        logger.debug("Writing other additional information.")
        append_additional_info = self.job.settings.append_additional_info
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
        if self.job.settings.link:
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
        f.write(f"{self.job.settings.link_route_string}\n")
        f.write("\n")
