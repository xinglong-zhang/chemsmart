import logging
import os.path

from chemsmart.jobs.writer import InputWriter

logger = logging.getLogger(__name__)


class ORCAInputWriter(InputWriter):
    """Class that writes ORCA input files for a job."""

    def write(self, **kwargs):
        self._write(**kwargs)

    def _write(self, target_directory=None):
        if target_directory is not None:
            if not os.path.exists(target_directory):
                os.makedirs(target_directory)
            folder = target_directory
        else:
            folder = self.job.folder
        job_inputfile = os.path.join(folder, f"{self.job.label}.inp")
        logger.debug(f"Writing ORCA input file: {job_inputfile}")
        f = open(job_inputfile, "w")
        if self.job.settings.input_string:
            # write the file itself for direct run
            self._write_self(f)
        else:
            self._write_all(f)
        logger.info(f"Finished writing ORCA input file: {job_inputfile}")
        f.close()

    def _write_all(self, f):
        """Write the entire input file."""
        self._write_route_section(f)
        self._write_processors(f)
        self._write_memory(f)
        self._write_scf_iterations(f)
        self._write_mdci_block(f)
        self._write_elprop_block(f)
        self._write_charge_and_multiplicity(f)
        self._write_cartesian_coordinates(f)
        # other functionalities in ORCA input job may be implemented here
        # self._append_modredundant(f)
        # self._append_gen_genecp_basis(f)  # then write genecp info
        # self._append_custom_solvent_parameters(
        #     f
        # )  # followed by user defined solvent parameters
        # self._append_job_specific_info(f)
        # self._append_other_additional_info(f)

    def _write_self(self, f):
        """Write the input file itself for direct run."""
        f.write(self.job.settings.input_string)

    def _write_route_section(self, f):
        logger.debug("Writing ORCA route section.")
        route_string = self.settings.route_string
        f.write(route_string + "\n")
        f.write("\n")

    def _write_processors(self, f):
        logger.debug("Writing processors.")
        f.write("# Number of processors")
        f.write(f"%pal nprocs {self.jobrunner.num_cores} end\n")

    def _write_memory(self, f):
        logger.debug("Writing memory.")
        f.write("# Memory per core")
        mpc = self.jobrunner.mem_gb / self.jobrunner.num_cores
        f.write(f"%maxcore {mpc}\n")

    def _write_scf_iterations(self, f):
        logger.debug("Writing SCF iterations.")
        if self.settings.scf_maxiter is not None:
            f.write(f"%scf maxiter {self.settings.scf_maxiter} end\n")

    def _write_mdci_block(self, f):
        logger.debug("Writing MDCI block.")
        mdci_cutoff = self.settings.mdci_cutoff
        mdci_density = self.settings.mdci_density
        if mdci_cutoff is not None:
            # check that mdci_cutoff is one of the allowed values: ["loose", "normal", "tight"]
            assert mdci_cutoff in ["loose", "normal", "tight"], (
                "mdci_cutoff must be one of the allowed values: "
                "['loose', 'normal', 'tight']"
            )
            f.write("%mdci\n")
            if mdci_cutoff.lower() == "loose":
                f.write("  # loose cutoff\n")
                f.write("  TCutPairs 1e-3\n")
                f.write("  TCutPNO 1e-6\n")
                f.write("  TCutMKN 1e-3\n")
            elif mdci_cutoff.lower() == "normal":
                f.write("  # normal cutoff\n")
                f.write("  TCutPairs 1e-4\n")
                f.write("  TCutPNO 3.33e-7\n")
                f.write("  TCutMKN 1e-3\n")
            elif mdci_cutoff.lower() == "tight":
                f.write("  # tight cutoff\n")
                f.write("  TCutPairs 1e-5\n")
                f.write("  TCutPNO 1e-7\n")
                f.write("  TCutMKN 1e-4\n")
            if mdci_density is not None:
                # check that mdci_density is one of the allowed values: ["none", "unrelaxed", "relaxed"]
                assert mdci_density in ["none", "unrelaxed", "relaxed"], (
                    "mdci_density must be one of the allowed values: "
                    "['none', 'unrelaxed', 'relaxed']"
                )
                if mdci_density.lower() == "none":
                    f.write("  Density None  # no density\n")
                elif mdci_density.lower() == "unrelaxed":
                    f.write("  Density Unrelaxed  # unrelaxed density\n")
                elif mdci_density.lower() == "relaxed":
                    f.write("  Density Relaxed  # relaxed density\n")
            f.write("end\n")

    def _write_elprop_block(self, f):
        """Write the elprop block for the ORCA input file."""
        logger.debug("Writing elprop block.")
        dipole = self.settings.dipole
        quadrupole = self.settings.quadrupole
        if dipole is not None or quadrupole is not None:
            f.write("%elprop\n")
            if dipole is True:
                f.write("  Dipole True\n")
            else:
                f.write("  Dipole False\n")
            if quadrupole is True:
                f.write("  Quadrupole True\n")
            else:
                f.write("  Quadrupole False\n")
            f.write("end\n")

    def _write_charge_and_multiplicity(self, f):
        logger.debug("Writing charge and multiplicity.")
        charge = self.settings.charge
        multiplicity = self.settings.multiplicity
        assert (
            charge is not None and multiplicity is not None
        ), "Charge and multiplicity must be specified!"
        f.write(f"* xyz {charge} {multiplicity}\n")

    def _write_cartesian_coordinates(self, f):
        logger.debug(
            f"Writing Cartesian coordinates of molecule: {self.job.molecule}."
        )
        assert self.job.molecule is not None, "No molecular geometry found!"
        self.job.molecule.write_coordinates(f, program="orca")
        f.write("*\n")
