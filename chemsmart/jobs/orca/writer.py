"""
ORCA input file writer implementation.

This module contains the input file writer for ORCA quantum chemistry
calculations, handling the generation of properly formatted input files
based on job settings and molecular structures.
"""

import logging
import os.path

from chemsmart.jobs.orca.settings import (
    ORCAIRCJobSettings,
    ORCAQMMMJobSettings,
    ORCATSJobSettings,
)
from chemsmart.jobs.writer import InputWriter
from chemsmart.utils.io import remove_keyword
from chemsmart.utils.utils import (
    get_prepend_string_list_from_modred_free_format,
)

logger = logging.getLogger(__name__)


class ORCAInputWriter(InputWriter):
    """
    ORCA input file writer.

    This class handles the generation of ORCA input files with proper
    formatting for quantum chemistry calculations, including route sections,
    job parameters, and molecular geometry specifications.

    Attributes:
        job (ORCAJob): Target job for which the input is generated.
        settings (ORCAJobSettings): Settings used to generate the route and blocks.
        jobrunner (ORCAJobRunner): Runner providing cores/memory and paths.
    """

    def write(self, **kwargs):
        """
        Write the ORCA input file.

        Args:
            **kwargs: Additional arguments for writing
        """
        self._write(**kwargs)

    def _write(self, target_directory=None):
        """
        Internal method to write the input file.

        Args:
            target_directory: Directory to write the file to
        """
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
        """
        Write the complete input file with all sections.

        Args:
            f: File object to write to
        """
        logger.debug("Writing complete ORCA input file")
        self._write_route_section(f)
        self._write_processors(f)
        self._write_memory(f)
        self._write_scf_block(f)
        self._write_solvent_block(f)
        self._write_mdci_block(f)
        self._write_elprop_block(f)
        self._write_qmmm_block(f)
        self._write_modred_block(f)
        self._write_hessian_block(f)
        self._write_irc_block(f)
        self._write_constrained_atoms(f)
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
        """
        Write user-provided input string directly.

        Args:
            f: File object to write to
        """
        f.write(self.job.settings.input_string)

    def _write_route_section(self, f):
        """
        Write the ORCA route section (calculation keywords).

        Args:
            f: File object to write to
        """
        logger.debug("Writing ORCA route section")

        if self.job.molecule.is_monoatomic:
            # TODO: need to merge thermo branch into main branch
            logger.info(f"Molecule {self.job.molecule} is monoatomic.")
            logger.info(
                "Removing `opt` keyword from route string, since ORCA cannot run OPT on monoatomic molecule."
            )
            route_string = remove_keyword(self.settings.route_string, "opt")
        else:
            route_string = self.settings.route_string

        if isinstance(self.settings, ORCAQMMMJobSettings):
            logger.debug("Writing qmmm route string")
            route_string = self.settings.qmmm_route_string

        f.write(route_string + "\n")

    def _write_processors(self, f):
        """
        Write processor specification for parallel execution.

        Args:
            f: File object to write to
        """
        logger.debug("Writing processors.")
        f.write("# Number of processors\n")
        f.write(f"%pal nprocs {self.jobrunner.num_cores} end\n")

    def _write_memory(self, f):
        """
        Write memory specification per core.

        Args:
            f: File object to write to
        """
        logger.debug("Writing memory specification")
        f.write("# Memory per core\n")
        mpc = (self.jobrunner.mem_gb * 1000) / self.jobrunner.num_cores * 0.75
        # Safety Factor: Applies a 75% factor to the memory per core to
        # reduce the risk of out-of-memory errors, as recommended for ORCA.
        mpc = int(mpc)
        f.write(f"%maxcore {mpc}\n")

    def _write_scf_block(self, f):
        """
        Write SCF convergence parameters.

        Args:
            f: File object to write to
        """
        logger.debug("Writing SCF block")

        if self.settings.scf_convergence or self.settings.scf_maxiter:
            f.write("%scf\n")
            self._write_scf_maxiter(f)
            self._write_scf_convergence(f)
            f.write("end\n")

    def _write_scf_maxiter(self, f):
        """
        Write SCF maximum iterations parameter.

        Args:
            f: File object to write to
        """
        scf_maxiter = (
            self.settings.scf_maxiter
            if self.settings.scf_maxiter is not None
            else 200
        )
        f.write(f"  maxiter {scf_maxiter}\n")

    def _write_scf_convergence(self, f):
        """
        Write SCF convergence criteria.

        Args:
            f: File object to write to

        Raises:
            ValueError: If unsupported convergence option is specified
        """
        if self.settings.scf_convergence:
            from chemsmart.io.orca import ORCA_SCF_CONVERGENCE

            scf_conv = self.settings.scf_convergence.lower().strip()
            # check that the convergence is in the list given by ORCA
            if scf_conv.endswith("scf"):
                scf_conv = scf_conv[:-3]
                if scf_conv not in ORCA_SCF_CONVERGENCE:
                    raise ValueError(
                        f"Warning: SCF convergence {self.scf_convergence} is not supported by ORCA!\n"
                        f"Available SCF convergence options are: {ORCA_SCF_CONVERGENCE}"
                    )
            f.write(f"  convergence {scf_conv}\n")

    def _write_solvent_block(self, f):
        """
        Write solvent parameters block.

        Args:
            f: File object to write to

        Note:
            Currently placeholder for complex solvents specified via
            %cpcm, %cosmo, or %smd blocks that cannot be captured by route.
        """
        # to implement if there is more complex solvents to be specified via
        # %cpcm block, %cosmo block, or %smd block that cannot be capture by route
        pass

    def _write_mdci_block(self, f):
        """
        Write MDCI (domain-based local correlation) parameters.

        Args:
            f: File object to write to

        Raises:
            AssertionError: If invalid MDCI options are specified
        """
        mdci_cutoff = self.settings.mdci_cutoff
        mdci_density = self.settings.mdci_density

        if mdci_cutoff is not None:
            logger.debug("Writing MDCI block")
            # check that mdci_cutoff is one of the allowed values: ["loose", "normal", "tight"]
            assert mdci_cutoff.lower() in ["loose", "normal", "tight"], (
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
                assert mdci_density.lower() in [
                    "none",
                    "unrelaxed",
                    "relaxed",
                ], (
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
        """
        Write electronic properties calculation block.

        Args:
            f: File object to write to
        """
        dipole = self.settings.dipole
        quadrupole = self.settings.quadrupole

        if dipole or quadrupole:
            logger.debug("Writing elprop block")
            f.write("%elprop\n")
            if dipole:
                f.write("  Dipole True\n")
            else:
                f.write("  Dipole False\n")
            if quadrupole is True:
                f.write("  Quadrupole True\n")
            else:
                f.write("  Quadrupole False\n")
            f.write("end\n")

    def _write_qmmm_block(self, f):
        """
        Write QM/MM parameter block to ORCA input file.

        Writes the complete %qmmm block for multiscale calculations
        when using ORCAQMMMJobSettings. The block contains all necessary
        parameters for defining QM/MM partitioning and calculation setup.

        Args:
            f: File object to write to
        """
        if isinstance(self.settings, ORCAQMMMJobSettings):
            logger.debug("Writing qmmm block")
            f.write(f"{self.settings.qmmm_block}\n")

    def _write_modred_block(self, f):
        """
        Write modified redundant coordinates block.

        Args:
            f: File object to write to
        """
        if self.settings.modred:
            f.write("%geom\n")
            self._write_modred(f, modred=self.settings.modred)
            f.write("end\n")

    def _write_modred(self, f, modred):
        """
        Write modified redundant coordinates.

        Args:
            f: File object to write to
            modred: Modified redundant coordinates specification
        """
        if isinstance(modred, list):
            self._write_modred_if_list(f, modred)
        elif isinstance(modred, dict):
            self._write_modred_if_dict(f, modred)

    @staticmethod
    def _write_modred_if_list(f, modred):
        """
        Write modred as constraints (list format).

        Args:
            f: File object to write to
            modred: List of coordinate constraints
        """
        f.write("  Constraints\n")
        # append for modred jobs
        # 'self.modred' as list of lists, or a single list if only one fixed constraint
        prepend_string_list = get_prepend_string_list_from_modred_free_format(
            input_modred=modred, program="orca"
        )
        for prepend_string in prepend_string_list:
            f.write(f"  {{{prepend_string} C}}\n")
        # write 'end' for each modred specified
        f.write("  end\n")

    @staticmethod
    def _write_modred_if_dict(f, modred):
        """
        Write modred as scan coordinates (dictionary format).

        Args:
            f: File object to write to
            modred: Dictionary containing scan parameters
        """
        # append for additonal contraints for scan jobs
        if "constrained_coordinates" in modred:
            logger.debug("Writing modredundant constraints for scan job")
            f.write("  Constraints\n")
            constrained_coordinates_list = modred["constrained_coordinates"]
            prepend_string_list = (
                get_prepend_string_list_from_modred_free_format(
                    input_modred=constrained_coordinates_list,
                    program="orca",
                )
            )
            for prepend_string in prepend_string_list:
                f.write(f"  {{{prepend_string} C}}\n")
            # write 'end' for each modred specified
            f.write("  end\n")
            f.write("\n")

        f.write("  Scan\n")
        # append for scanning job
        coords_list = modred["coords"]
        scan_prepend_string_list = (
            get_prepend_string_list_from_modred_free_format(
                input_modred=coords_list, program="orca"
            )
        )
        for prepend_string, dist_start, dist_end, num_steps in zip(
            scan_prepend_string_list,
            modred["dist_start"],
            modred["dist_end"],
            modred["num_steps"],
        ):
            if prepend_string.lower().startswith("b"):
                scan_var = "bond distance"
                scan_unit = "Angstrom"
            elif prepend_string.lower().startswith("a"):
                scan_var = "angle"
                scan_unit = "degree"
            elif prepend_string.lower().startswith("d"):
                scan_var = "dihedral"
                scan_unit = "degree"
            else:
                scan_var = "variable"
                scan_unit = "unit"
            f.write(
                f"  {prepend_string} = {dist_start}, {dist_end}, "
                f"{num_steps}  # Scanning {scan_var} from {dist_start} {scan_unit} "
                f"to {dist_end} {scan_unit} in {num_steps} steps. \n"
            )
        # write 'end' for each modred specified
        f.write("  end\n")

    def _write_hessian_block(self, f):
        """
        Write Hessian calculation parameters for transition state jobs.

        Args:
            f: File object to write to
        """
        if isinstance(self.settings, ORCATSJobSettings):
            self._write_hessian_block_for_ts(f)

    def _write_hessian_block_for_ts(self, f):
        """
        Write ORCA Hessian block for transition state calculations.

        Args:
            f: File object to write to

        Raises:
            AssertionError: If required Hessian files or parameters are missing
        """
        # write orca block for hessian options
        f.write("%geom\n")

        # Read initial Hessian from file if desired
        if self.settings.inhess:
            f.write("  InHess Read  # Read Hessian from file\n")
            assert (
                self.settings.inhess_filename is not None
            ), "No Hessian file is given!"
            assert os.path.exists(
                self.settings.inhess_filename
            ), f"Hessian file {self.settings.inhess_filename} is not found!"
            f.write(
                f'  InHessName "{self.settings.inhess_filename}"  # Hessian file\n'
            )

        """Hybrid Hessian for speed up of TS search: TS mode is complicated and delocalized, 
        e.g. in a concerted proton transfer reaction, can use hybrid Hessian to calc 
        numerical second derivatives only for atoms involved in the TS mode"""
        if self.settings.hybrid_hess:
            assert (
                self.settings.hybrid_hess_atoms is not None
            ), "No atoms are specified for hybrid Hessian calculation!"
            hybrid_hess_atoms_string = " ".join(
                [str(i - 1) for i in self.settings.hybrid_hess_atoms]
                # using 1-indices from user, convert them to 0-indices for ORCA
            )
            f.write(
                f"  Hybrid_Hess {{{hybrid_hess_atoms_string}}} end  # Use hybrid Hessian\n"
            )

        # Hessian options
        f.write(
            "  Calc_Hess True  # calc initial Hessian\n"
        )  # for ts job, initial hessian is required
        f.write(
            f"  NumHess {self.settings.numhess}  # Request numerical Hessian (if analytical not available)\n"
        )
        f.write(
            f"  Recalc_Hess {self.settings.recalc_hess}   # Recalculate the Hessian every 5 step\n"
        )

        # trust radius update
        if self.settings.trust_radius is not None:
            f.write(f"  Trust {self.settings.trust_radius}\n")
            if self.settings.trust_radius < 0:
                f.write("  # use fixed trust radius (default: -0.3 au)\n")
            elif self.settings.trust_radius > 0:
                f.write("  # use trust radius update, i.e. 0.3 means:\n")
                f.write(
                    "  # start with trust radius 0.3 and use trust radius update\n"
                )

        # TS search type
        if self.settings.tssearch_type.lower() == "scants":
            # for scanTS, also include in the scanning coordinates
            assert (
                self.settings.scants_modred is not None
            ), "No modred is specified for scanTS!"
            self._write_modred_if_dict(f, self.settings.scants_modred)

        # full scan or not
        if self.settings.full_scan:
            f.write("  fullScan True\n")
        f.write("end\n")

    def _write_irc_block(self, f):
        """
        Write IRC calculation parameters.

        Args:
            f: File object to write to
        """
        if isinstance(self.settings, ORCAIRCJobSettings):
            self._write_irc_block_for_irc(f)

    def _write_irc_block_for_irc(self, f):
        """Writes the IRC block options.

        IRC block input example below:
        ! IRC
        %irc
            MaxIter    20
            PrintLevel 1
            Direction  both # both - default
                            # forward
                            # backward
                            # down
        # Initial displacement
            InitHess   read # by default ORCA uses the Hessian from AnFreq or NumFreq, or computes a new one
                            # read    - reads the Hessian that is defined via Hess_Filename
                            # calc_anfreq  - computes the analytic Hessian
                            # calc_numfreq - computes the numeric Hessian
            Hess_Filename "h2o.hess"  # Hessian for initial displacement, must be used together with InitHess = read
            hessMode   0  # Hessian mode that is used for the initial displacement. Default 0
            Init_Displ DE      # DE (default) - energy difference
                               # length       - step size
            Scale_Init_Displ 0.1 # step size for initial displacement from TS. Default 0.1 a.u.
            DE_Init_Displ    2.0 # energy difference that is expected for initial displacement
                                 #  based on provided Hessian (Default: 2 mEh)
        # Steps
            Follow_CoordType cartesian # default and only option
            Scale_Displ_SD    0.15  # Scaling factor for scaling the 1st SD step
            Adapt_Scale_Displ true  # modify Scale_Displ_SD when the step size becomes smaller or larger
            SD_ParabolicFit   true  # Do a parabolic fit for finding an optimal SD step length
            Interpolate_only  true  # Only allow interpolation for parabolic fit, not extrapolation
            Do_SD_Corr        true  # Apply a correction to the 1st SD step
            Scale_Displ_SD_Corr  0.333 # Scaling factor for scaling the correction step to the SD step.
                                       # It is multiplied by the length of the final 1st SD step
            SD_Corr_ParabolicFit true  # Do a parabolic fit for finding an optimal correction
                                       # step length
        # Convergence thresholds - similar to LooseOpt
            TolRMSG   5.e-4      # RMS gradient (a.u.)
            TolMaxG   2.e-3      # Max. element of gradient (a.u.)
        # Output options
            Monitor_Internals   # Up to three internal coordinates can be defined
                {B 0 1}         # for which the values are printed during the IRC run.
                {B 1 5}         # Possible are (B)onds, (A)ngles, (D)ihedrals and (I)mpropers
            end
        end.
        """
        irc_settings_keys = self.settings.__dict__.keys()
        from chemsmart.jobs.orca.settings import ORCAJobSettings

        parent_settings_keys = ORCAJobSettings().__dict__.keys()
        irc_specific_keys = set(irc_settings_keys) - set(parent_settings_keys)

        if not any(
            getattr(self.settings, key) is not None
            for key in irc_specific_keys
        ):
            return

        # write irc block if any option value is not None:
        f.write("%irc\n")
        for key in irc_specific_keys:
            value = getattr(self.settings, key)
            if value is None:
                continue  # ignore the rest of the code and go to next in the for loop
            # only write into IRC input if the value is not None
            if key == "internal_modred":
                pass  # internal_modred is not an option in ORCA IRC file
            elif key == "inithess":
                f.write(f"  {key} {value}\n")
                if value.lower() == "read":  # if initial hessian is to be read
                    assert (
                        self.settings.hess_filename is not None
                    ), "No Hessian file is given!"
                    assert os.path.exists(
                        self.settings.hess_filename
                    ), f"Hessian file {self.settings.hess_filename} is not found!"
                    f.write(
                        f'  Hess_Filename "{self.settings.hess_filename}"  # Hessian file\n'
                    )
            elif (
                key == "hess_filename"
            ):  # already used/written, if initial hessian is to be read
                pass
            elif key == "monitor_internals":
                if value is True:
                    f.write("  True\n")
                    assert (
                        self.settings.internal_modred is not None
                    ), 'No internal modred is specified for IRC job "monitor_intervals" option!'
                    prepend_string_list = (
                        get_prepend_string_list_from_modred_free_format(
                            self.settings.internal_modred, program="orca"
                        )
                    )
                    for prepend_string in prepend_string_list:
                        f.write(f"  {{ {prepend_string} }}\n")
                    f.write("  end\n")
            elif key == "adapt_scale_displ":
                if value is True:
                    f.write("  Adapt_Scale_Displ True\n")
            elif key == "sd_parabolicfit":
                if value is True:
                    f.write("  SD_ParabolicFit True\n")
            elif key == "interpolate_only":
                if value is True:
                    f.write("  Interpolate_only True\n")
            elif key == "do_sd_corr":
                if value is True:
                    f.write("  Do_SD_Corr True\n")
            elif key == "sd_corr_parabolicfit":
                if value is True:
                    f.write("  SD_Corr_ParabolicFit True\n")
            else:  # all other keys with given values
                f.write(f"  {key} {value}\n")
        f.write("end\n")

    def _write_constrained_atoms(self, f):
        """
        Write atomic constraints for frozen atoms.

        Args:
            f: File object to write to
        """
        molecule = self.job.molecule
        if molecule.frozen_atoms:
            f.write("%geom\n")
            for i, val in enumerate(molecule.frozen_atoms):
                if val == -1:
                    f.write(f"  {{ C {i} C }}\n")  # ORCA is 0-indexed
            if self.settings.invert_constraints:
                f.write("  InvertConstraints True\n")
            f.write("end\n")

    def _write_charge_and_multiplicity(self, f):
        """
        Write molecular charge and spin multiplicity.

        Args:
            f: File object to write to

        Raises:
            AssertionError: If charge or multiplicity is not specified
        """
        charge = getattr(self.settings, "charge", None)
        multiplicity = getattr(self.settings, "multiplicity", None)

        # If missing, attempt to populate from common QMMM-related fields.
        # Common names across settings: charge_qm, charge_medium, charge_total
        # and mult_qm, mult_medium, mult_total.
        if charge is None or multiplicity is None:
            # order of preference: medium (QM2) -> qm -> total
            candidate_charge_attrs = [
                "charge_medium",
                "charge_qm",
                "charge_total",
                "charge",
            ]
            candidate_mult_attrs = [
                "mult_medium",
                "mult_qm",
                "mult_total",
                "multiplicity",
            ]

            for attr in candidate_charge_attrs:
                if charge is None:
                    charge = getattr(self.settings, attr, None)
            for attr in candidate_mult_attrs:
                if multiplicity is None:
                    multiplicity = getattr(self.settings, attr, None)

            # if we found values, assign back to settings for downstream use
            if charge is not None:
                logger.debug(
                    f"Populating settings.charge from QMMM fields: {charge}"
                )
                self.settings.charge = charge
            if multiplicity is not None:
                logger.debug(
                    f"Populating settings.multiplicity from QMMM fields: {multiplicity}"
                )
                self.settings.multiplicity = multiplicity

        assert (
            charge is not None and multiplicity is not None
        ), "Charge and multiplicity must be specified!"
        f.write(f"* xyz {charge} {multiplicity}\n")

    def _write_cartesian_coordinates(self, f):
        """
        Write molecular Cartesian coordinates.

        Args:
            f: File object to write to

        Raises:
            AssertionError: If molecular geometry is not available
        """
        logger.debug(
            f"Writing Cartesian coordinates of molecule: {self.job.molecule}"
        )
        assert self.job.molecule is not None, "No molecular geometry found!"
        self.job.molecule.write_coordinates(f, program="orca")
        f.write("*\n")
