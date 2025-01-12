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
    """Class that writes Gaussian input files for a job.

    Args:
        job (Job): The job for which the input file is being written.
    """

    def write(self):
        logger.info(f"Writing Gaussian input file: {self.job.inputfile}")
        f = open(self.job.inputfile, "w")
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
        f.close()

    def _write_gaussian_header(self, f):
        if self.job.settings.chk:
            f.write(f"%chk={self.job.label}.chk\n")
        f.write(f"%nprocshared={self.jobrunner.num_cores}\n")
        f.write(f"%mem={self.jobrunner.mem_gb}GB\n")

    def _write_route_section(self, f):
        route_string = self.job.settings.route_string
        f.write(route_string + "\n")
        f.write("\n")

    def _write_gaussian_title(self, f):
        title = self.job.settings.title
        f.write(f"{title}\n")
        f.write("\n")

    def _write_charge_and_multiplicity(self, f):
        charge = self.job.settings.charge
        multiplicity = self.job.settings.multiplicity
        f.write(f"{charge} {multiplicity}\n")
        f.write("\n")

    def _write_cartesian_coordinates(self, f):
        assert self.job.molecule is not None, "No molecular geometry found!"
        self.job.molecule.write_coordinates(f)

    def _write_pbc(self, f):
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
        """Write modredundant section if present in the job settings.
        Given a list or list of lists, write the modredundant section for the constraints.
        Given a dictionary with 'num_steps', 'step_size', and 'coords', write the scan section.
        """
        modredundant = self.job.settings.modredundant
        if modredundant:
            if isinstance(modredundant, list):
                # for modredunant job
                # modredundant = [[1,2], [3,4]]
                prepend_string_list = (
                    get_prepend_string_list_from_modred_free_format(
                        input_modred=modredundant
                    )
                )
                for prepend_string in prepend_string_list:
                    f.write(f"{prepend_string} F\n")

            elif isinstance(modredundant, dict):
                # for scanning job
                # modredundant = {'num_steps': 10, 'step_size': 0.05, 'coords': [[1,2], [3,4]]}
                coords_list = modredundant["coords"]
                prepend_string_list = (
                    get_prepend_string_list_from_modred_free_format(
                        input_modred=coords_list
                    )
                )
                for prepend_string in prepend_string_list:
                    f.write(
                        f"{prepend_string} S {self.modred['num_steps']} {self.modred['step_size']}\n"
                    )
                    f.write("\n")

    def _append_gen_genecp_basis(self, f):
        """Write the genecp basis set information if present in the job settings."""
        gen_genecp_file = self.job.settings.gen_genecp_file
        heavy_elements = self.job.settings.heavy_elements
        heavy_elements_basis = self.job.settings.heavy_elements_basis
        light_elements = self.job.settings.get_light_elements(
            self.job.molecule
        )
        light_elements_basis = self.job.settings.light_elements_basis
        if gen_genecp_file and heavy_elements and heavy_elements_basis:
            raise ValueError(
                "Please provide either gen_genecp_file or heavy_elements and heavy_elements_basis."
            )
        if gen_genecp_file:
            assert os.path.exists(
                gen_genecp_file
            ), f"File {gen_genecp_file} not found!"
            genecp_section = GenGenECPSection.from_genecp_path(gen_genecp_file)
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
        custom_solvent = self.job.settings.custom_solvent
        if custom_solvent:
            f.write(custom_solvent)
            f.write("\n")

    def _append_job_specific_info(self, f):
        """Write any job specific information that needs to be appended to the input file."""
        job_type = self.job.settings.job_type
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
        if self.job.settings.link:
            f.write("--Link1--\n")
            self._write_gaussian_header(f)
            f.write(f"# {self.job.settings.link_route_string}\n")
            f.write("\n")
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
