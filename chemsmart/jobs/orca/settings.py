"""
ORCA job settings implementation.

This module contains the settings classes for configuring ORCA quantum chemistry
calculations, including general settings and specialized settings for transition
state and IRC calculations.
"""

import copy
import logging
import os
import re

from chemsmart.io.orca import ORCA_ALL_SOLVENT_MODELS
from chemsmart.jobs.settings import MolecularJobSettings
from chemsmart.utils.utils import (
    get_list_from_string_range,
    get_prepend_string_list_from_modred_free_format,
)

logger = logging.getLogger(__name__)


class ORCAJobSettings(MolecularJobSettings):
    """
    Settings for ORCA quantum chemistry jobs.

    This class handles the configuration of ORCA calculations including
    method selection, basis sets, solvation, convergence criteria, and
    various calculation types.

    Attributes:
        ab_initio (str | None): Ab initio method (e.g., 'HF', 'MP2').
        functional (str | None): DFT functional (e.g., 'B3LYP', 'PBE0').
        dispersion (str | None): Dispersion correction (e.g., 'D3BJ').
        basis (str | None): Orbital basis set (e.g., 'def2-TZVP').
        aux_basis (str | None): Auxiliary basis set for RI/DF methods.
        extrapolation_basis (str | None): Basis for extrapolation schemes.
        defgrid (str | None): Grid quality keyword (e.g., 'defgrid2').
        scf_tol (str | float | None): SCF convergence tolerance.
        scf_algorithm (str | None): SCF algorithm.
        scf_maxiter (int | None): Maximum SCF iterations.
        scf_convergence (str | None): SCF convergence criteria.
        charge (int | None): Molecular charge.
        multiplicity (int | None): Spin multiplicity.
        gbw (bool): Write GBW file.
        freq (bool): Perform vibrational frequency calculation.
        numfreq (bool): Use numerical frequencies.
        dipole (bool): Compute dipole moment.
        quadrupole (bool): Compute quadrupole moment.
        mdci_cutoff (float | None): MDCI cutoff.
        mdci_density (float | None): MDCI density.
        jobtype (str | None): Calculation type (e.g., 'opt', 'sp').
        title (str | None): Job title string.
        solvent_model (str | None): Solvation model identifier.
        solvent_id (str | None): Solvent identifier.
        additional_route_parameters (str | None): Extra route parameters.
        route_to_be_written (str | None): Custom route string to write.
        modred (list | dict | None): Modredundant coordinates specification.
        gen_genecp_file (str | None): Path to gen/genecp file to include.
        heavy_elements (list | None): Heavy elements list for genecp.
        heavy_elements_basis (str | None): Basis for heavy elements.
        light_elements_basis (str | None): Basis for light elements.
        custom_solvent (str | None): Custom solvent parameters block.
        forces (bool): Calculate forces.
        input_string (str | None): Predefined input content to write directly.
        invert_constraints (bool): Invert modred constraints if True.
    """

    def __init__(
        self,
        ab_initio=None,
        functional=None,
        dispersion=None,
        basis=None,
        aux_basis=None,
        extrapolation_basis=None,
        defgrid=None,
        scf_tol=None,
        scf_algorithm=None,
        scf_maxiter=None,
        scf_convergence=None,
        charge=None,
        multiplicity=None,
        gbw=True,
        freq=False,
        numfreq=False,
        dipole=False,
        quadrupole=False,
        mdci_cutoff=None,
        mdci_density=None,
        jobtype=None,
        title=None,
        solvent_model=None,
        solvent_id=None,
        additional_route_parameters=None,
        route_to_be_written=None,
        modred=None,
        gen_genecp_file=None,
        heavy_elements=None,
        heavy_elements_basis=None,
        light_elements_basis=None,
        custom_solvent=None,
        forces=False,
        input_string=None,
        invert_constraints=False,
        **kwargs,
    ):
        """
        Initialize ORCA job settings.

        Args:
            ab_initio: Ab initio method (e.g., 'HF', 'MP2')
            functional: DFT functional (e.g., 'B3LYP', 'PBE0')
            dispersion: Dispersion correction (e.g., 'D3BJ')
            basis: Basis set (e.g., 'def2-TZVP')
            aux_basis: Auxiliary basis set
            extrapolation_basis: Extrapolation basis set
            defgrid: Grid quality (e.g., 'defgrid2')
            scf_tol: SCF convergence tolerance
            scf_algorithm: SCF algorithm
            scf_maxiter: Maximum SCF iterations
            scf_convergence: SCF convergence criteria
            charge: Molecular charge
            multiplicity: Spin multiplicity
            gbw: Whether to write GBW file
            freq: Whether to calculate frequencies
            numfreq: Whether to use numerical frequencies
            dipole: Whether to calculate dipole moment
            quadrupole: Whether to calculate quadrupole moment
            mdci_cutoff: MDCI cutoff
            mdci_density: MDCI density
            jobtype: Type of calculation
            title: Job title
            solvent_model: Solvation model
            solvent_id: Solvent identifier
            additional_route_parameters: Additional route parameters
            route_to_be_written: Custom route string
            modred: Modified redundant coordinates
            gen_genecp_file: General ECP file
            heavy_elements: Heavy elements list
            heavy_elements_basis: Basis for heavy elements
            light_elements_basis: Basis for light elements
            custom_solvent: Custom solvent parameters
            forces: Whether to calculate forces
            input_string: Custom input string
            invert_constraints: Whether to invert constraints
            **kwargs: Additional keyword arguments
        """
        super().__init__(
            ab_initio=ab_initio,
            functional=functional,
            dispersion=dispersion,
            basis=basis,
            defgrid=defgrid,
            charge=charge,
            multiplicity=multiplicity,
            freq=freq,
            numfreq=numfreq,
            jobtype=jobtype,
            title=title,
            solvent_model=solvent_model,
            solvent_id=solvent_id,
            additional_route_parameters=additional_route_parameters,
            route_to_be_written=route_to_be_written,
            modred=modred,
            gen_genecp_file=gen_genecp_file,
            heavy_elements=heavy_elements,
            heavy_elements_basis=heavy_elements_basis,
            light_elements_basis=light_elements_basis,
            custom_solvent=custom_solvent,
            forces=forces,
            input_string=input_string,
            **kwargs,
        )

        # ORCA-specific parameters
        self.aux_basis = aux_basis
        self.extrapolation_basis = extrapolation_basis
        self.scf_tol = scf_tol
        self.scf_algorithm = scf_algorithm
        self.scf_maxiter = scf_maxiter
        self.scf_convergence = scf_convergence
        self.gbw = gbw
        self.mdci_cutoff = mdci_cutoff
        self.mdci_density = mdci_density
        self.dipole = dipole
        self.quadrupole = quadrupole
        self.invert_constraints = invert_constraints

        # Validate frequency and force settings
        if forces is True and (freq is True or numfreq is True):
            raise ValueError(
                "Frequency and Force calculations cannot be performed by "
                "Orca at the same time!\n"
                'Such an input file will give "Illegal IType or MSType '
                'generated by parse." error.'
            )

    def merge(
        self, other, keywords=("charge", "multiplicity"), merge_all=False
    ):
        """
        Merge this settings object with another.

        Args:
            other: Settings object or dictionary to merge with
            keywords: Specific keywords to merge if merge_all is False
            merge_all: Whether to merge all attributes

        Returns:
            ORCAJobSettings: New merged settings object
        """

        other_dict = other if isinstance(other, dict) else other.__dict__

        if merge_all:
            # Update self with other for all
            merged_dict = self.__dict__.copy()
            merged_dict.update(other_dict)
            return type(self)(**merged_dict)

        if keywords is not None:
            other_dict = {
                k: other_dict[k] for k in keywords if k in other_dict
            }
        # Update self with other
        merged_dict = self.__dict__.copy()
        merged_dict.update(other_dict)
        return type(self)(**merged_dict)

    def copy(self):
        """
        Create a deep copy of the settings object.

        Returns:
            ORCAJobSettings: Deep copy of this settings object
        """
        return copy.deepcopy(self)

    def __getitem__(self, key):
        """
        Get settings attribute by key.

        Args:
            key: Attribute name to retrieve

        Returns:
            Value of the specified attribute
        """
        return self.__dict__[key]

    def __eq__(self, other):
        """
        Check equality with another settings object.

        Args:
            other: Another settings object to compare

        Returns:
            bool: True if settings are equal, False otherwise
        """
        if type(self) is not type(other):
            return NotImplemented

        # Exclude append_additional_info from the comparison
        self_dict = self.__dict__
        self_dict.pop("append_additional_info")

        other_dict = other.__dict__
        other_dict.pop("append_additional_info")

        return self_dict == other_dict

    @classmethod
    def from_comfile(cls, com_path):
        """
        Create ORCA job settings from a Gaussian .com file.

        Args:
            com_path: Path to the Gaussian .com file

        Returns:
            ORCAJobSettings: Settings object from Gaussian file
        """
        com_path = os.path.abspath(com_path)
        from chemsmart.io.gaussian.input import Gaussian16Input

        logger.info(f"Return Settings object from {com_path}")
        gaussian_settings_from_comfile = Gaussian16Input(
            filename=com_path
        ).read_settings()
        orca_default_settings = cls.default()
        return orca_default_settings.merge(
            gaussian_settings_from_comfile, merge_all=True
        )

    @classmethod
    def from_logfile(cls, log_path, **kwargs):
        """
        Create ORCA settings from Gaussian output .log file.

        Args:
            log_path: Path to the .log file
            **kwargs: Additional arguments for Gaussian16Output class

        Returns:
            ORCAJobSetting: Settings object from log file
        """
        log_path = os.path.abspath(log_path)
        from chemsmart.io.gaussian.output import Gaussian16Output

        logger.info(f"Return Settings object from {log_path}")
        gaussian_settings_from_logfile = Gaussian16Output(
            log_path
        ).read_settings()
        orca_default_settings = cls.default()
        orca_settings_from_logfile = orca_default_settings.merge(
            gaussian_settings_from_logfile, merge_all=True
        )
        # Convert def2 basis set naming
        if (
            "def2" in orca_settings_from_logfile.basis
            and "def2-" not in orca_settings_from_logfile.basis
        ):
            orca_settings_from_logfile.basis = (
                orca_settings_from_logfile.basis.replace("def2", "def2-")
            )
        return orca_settings_from_logfile

    @classmethod
    def from_inpfile(cls, inp_path):
        """
        Create ORCA job settings from ORCA .inp file.

        Args:
            inp_path: Path to the ORCA .inp file

        Returns:
            ORCAJobSettings: Settings object from input file
        """
        inp_path = os.path.abspath(inp_path)
        from chemsmart.io.orca.input import ORCAInput

        logger.info(f"Return settings object from {inp_path}")
        orca_settings_from_inpfile = ORCAInput(
            filename=inp_path
        ).read_settings()
        logger.info(f"with settings: {orca_settings_from_inpfile.__dict__}")
        return orca_settings_from_inpfile

    @classmethod
    def from_outfile(cls, out_path):
        """
        Create ORCA job settings from ORCA .out file.

        Args:
            out_path: Path to the ORCA .out file

        Returns:
            ORCAJobSettings: Settings object from output file
        """
        out_path = os.path.abspath(out_path)
        from chemsmart.io.orca.output import ORCAOutput

        logger.info(f"Return Settings object from {out_path}")
        return ORCAOutput(filename=out_path).read_settings()

    @classmethod
    def from_xyzfile(cls):
        """
        Create default ORCA job settings for .xyz file input.

        Returns:
            ORCAJobSettings: Default ORCA settings for xyz input
        """
        return ORCAJobSettings.default()

    @classmethod
    def from_filepath(cls, filepath, **kwargs):
        """
        Create settings from any supported file type.

        Args:
            filepath: Path to the input file
            **kwargs: Additional keyword arguments

        Returns:
            ORCAJobSettings or None: Settings object if file type supported
        """

        if ".com" in filepath or ".gjf" in filepath:
            return cls.from_comfile(filepath)

        if ".log" in filepath:
            return cls.from_logfile(filepath)

        if ".inp" in filepath:
            return cls.from_inpfile(filepath)

        if ".out" in filepath:
            return cls.from_outfile(filepath)

        if ".xyz" in filepath:
            return cls.from_xyzfile()

        return None

    @classmethod
    def default(cls):
        """
        Create default ORCA job settings.

        Returns:
            ORCAJobSettings: Default settings object
        """
        return cls(
            ab_initio=None,
            functional=None,
            dispersion=None,
            basis=None,
            aux_basis=None,
            extrapolation_basis=None,
            defgrid=None,
            scf_tol=None,
            scf_algorithm=None,
            scf_maxiter=None,
            scf_convergence=None,
            charge=None,
            multiplicity=None,
            gbw=True,
            freq=True,
            numfreq=False,
            mdci_cutoff=None,
            mdci_density=None,
            jobtype=None,
            title=None,
            solvent_model=None,
            solvent_id=None,
            additional_route_parameters=None,
            route_to_be_written=None,
            modred=None,
            gen_genecp_file=None,
            heavy_elements=None,
            heavy_elements_basis=None,
            light_elements_basis=None,
            custom_solvent=None,
            forces=False,
            input_string=None,
            invert_constraints=False,
        )

    @property
    def route_string(self):
        """
        Get the ORCA route string for the job.

        Returns:
            str: ORCA route string based on job settings
        """
        if self.route_to_be_written is not None:
            route_string = self._get_route_string_from_user_input()
        else:
            route_string = self._get_route_string_from_jobtype()
        logger.debug(f"Route for settings {self}: {route_string}")
        return route_string

    def _get_route_string_from_user_input(self):
        """
        Get route string from user-provided input.

        Returns:
            str: User-specified route string with '!' prefix if needed
        """
        route_string = self.route_to_be_written
        if not route_string.startswith("!"):
            route_string = f"! {route_string}"
        return route_string

    def _get_route_string_from_jobtype(self):
        """
        Generate ORCA route string based on job type and settings.

        Returns:
            str: Generated ORCA route string

        Raises:
            ValueError: If invalid combination of settings is specified
        """
        route_string = ""
        if not route_string.startswith("!"):
            route_string += "! "

        # route string depends on job type
        # determine if route string requires 'opt' keyword
        if self.jobtype in ("opt", "modred", "scan"):
            route_string += "Opt"
        elif self.jobtype == "ts":
            route_string += (
                "OptTS"  # Orca keyword for transition state optimization
            )
        elif self.jobtype == "irc":
            route_string += "IRC"
        elif self.jobtype == "sp":
            route_string += ""

        # add frequency calculation
        # not okay if both freq and numfreq are True
        if self.freq and self.numfreq:
            raise ValueError("Cannot specify both freq and numfreq!")

        if self.freq:
            route_string += " Freq"
        elif self.numfreq:
            route_string += " NumFreq"  # requires numerical frequency,
            # e.g., in SMD model where analytic Hessian is not available

        # write level of theory
        level_of_theory = self._get_level_of_theory()
        route_string += f" {level_of_theory}"

        # write grid information
        if self.defgrid is not None:
            route_string += (
                f" {self.defgrid}"  # default is 'defgrid2', if not specified
            )

        # write convergence criteria in simple input/route
        if self.scf_tol is not None:
            if not self.scf_tol.lower().endswith("scf"):
                self.scf_tol += "SCF"
            route_string += f" {self.scf_tol}"

        # write convergence algorithm if not default
        if self.scf_algorithm is not None:
            route_string += f" {self.scf_algorithm}"

        # write solvent if solvation is turned on
        if self.solvent_model is not None and self.solvent_id is not None:
            route_string += f" {self.solvent_model}({self.solvent_id})"
        elif self.solvent_model is not None and self.solvent_id is None:
            raise ValueError(
                "Warning: Solvent model is specified but solvent identity "
                "is missing!"
            )
        elif self.solvent_model is None and self.solvent_id is not None:
            logger.warning(
                "Warning: Solvent identity is specified but solvent model "
                "is missing!\nDefaulting to CPCM model."
            )
            route_string += f" CPCM({self.solvent_id})"
        else:
            pass

        return route_string

    def _get_level_of_theory(self):
        """
        Get the level of theory string.

        Returns:
            str: Level of theory specification

        Raises:
            ValueError: If invalid method combination or missing parameters
        """
        level_of_theory = ""

        # detect whether this is a QMMM-type job to relax some validations
        job_type_val = getattr(self, "job_type", None) or getattr(
            self, "jobtype", None
        )
        is_qmmm = False
        if job_type_val and (
            "qm" in str(job_type_val).lower()
            and "mm" in str(job_type_val).lower()
        ):
            try:
                is_qmmm = True
            except Exception:
                is_qmmm = False

        if self.ab_initio is not None and self.functional is not None:
            raise ValueError(
                "Warning: both ab initio and DFT are specified!\n"
                "Please specify only one method!"
            )

        if self.ab_initio is None and self.functional is None:
            # Allow missing ab_initio/functional for ORCA QMMM-type jobs.
            if not is_qmmm:
                raise ValueError(
                    "Warning: neither ab initio nor DFT is specified!\n"
                    "Please specify one method!"
                )
        if self.ab_initio is not None:
            level_of_theory += f"{self.ab_initio}"
        elif self.functional is not None:
            level_of_theory += f"{self.functional}"

        if self.basis is not None:
            level_of_theory += f" {self.basis}"
        elif self.basis is None:
            # allow missing basis for QMMM-type jobs where basis may be
            # provided per-layer (or omitted)
            if not is_qmmm:
                raise ValueError("Warning: basis is missing!")

        if self.aux_basis is not None:
            level_of_theory += f" {self.aux_basis}"

        if self.extrapolation_basis is not None:
            level_of_theory += f" {self.extrapolation_basis}"
        return level_of_theory

    def _write_geometry(self, f, atoms):
        """
        Write molecular geometry to input file.

        Args:
            f: File object to write to
            atoms: Atoms object containing molecular geometry

        Raises:
            AssertionError: If charge, multiplicity, or geometry is missing
        """
        # check that both charge and multiplicity are specified
        assert self.charge is not None, "No charge found!"
        assert self.multiplicity is not None, "No multiplicity found!"
        f.write(f"* xyz {self.charge} {self.multiplicity}\n")

        # check that a molecular geometry is given
        assert atoms is not None, "No molecular geometry found!"
        logger.info(f"Molecule given is: {atoms}")

        coordinates = ""
        for _i, (s, (x, y, z)) in enumerate(
            zip(atoms.symbols, atoms.positions, strict=False)
        ):
            string = f"{s:5} {float(x):15.10f} {float(y):15.10f} {float(z):15.10f}\n"
            coordinates += string
        f.write(coordinates)
        f.write("*\n")

    def _check_solvent(self, solvent_model):
        """
        Validate solvent model specification.

        Args:
            solvent_model: Solvent model name to validate

        Raises:
            ValueError: If solvent model is not supported
        """
        if solvent_model.lower() not in ORCA_ALL_SOLVENT_MODELS:
            raise ValueError(
                f"The specified solvent model {solvent_model} is not in \n"
                f"the available solvent models: {ORCA_ALL_SOLVENT_MODELS}"
            )


class ORCATSJobSettings(ORCAJobSettings):
    """
    Settings for ORCA transition state calculations.

    This class extends ORCAJobSettings with specialized options for
    transition state searches and optimizations.

    Attributes:
        inhess (bool): Read initial Hessian if True.
        inhess_filename (str | None): Path to initial Hessian file (used when inhess=True).
        hybrid_hess (bool): Use hybrid Hessian scheme.
        hybrid_hess_atoms (list[int] | None): 1‑based atom indices for hybrid Hessian region.
        numhess (bool): Use numerical Hessian.
        recalc_hess (int): Frequency (in cycles) to recalculate Hessian.
        trust_radius (float | None): Trust radius for optimization.
        tssearch_type (str): TS search method ('optts' or 'scants').
        scants_modred (list | dict | None): Modredundant coordinates for ScanTS.
        full_scan (bool): If True, do not abort ScanTS after highest point.
    """

    def __init__(
        self,
        inhess=False,
        inhess_filename=None,
        hybrid_hess=False,
        hybrid_hess_atoms=None,
        numhess=False,
        recalc_hess=5,
        trust_radius=None,
        tssearch_type="optts",
        scants_modred=None,
        full_scan=False,
        **kwargs,
    ):
        """
        Initialize ORCA transition state job settings.

        Args:
            inhess: Whether to read initial Hessian
            inhess_filename: Filename for initial Hessian
            hybrid_hess: Whether to use hybrid Hessian
            hybrid_hess_atoms: Atoms for hybrid Hessian (1-indexed)
            numhess: Whether to use numerical Hessian
            recalc_hess: Frequency for Hessian recalculation
            trust_radius: Trust radius for optimization
            tssearch_type: Type of TS search ('optts' or 'scants')
            scants_modred: Modified coordinates for ScanTS
            full_scan: Whether to do full scan or abort after highest point
            **kwargs: Additional keyword arguments
        """
        super().__init__(**kwargs)
        self.inhess = inhess
        self.inhess_filename = inhess_filename
        self.hybrid_hess = hybrid_hess
        self.hybrid_hess_atoms = (
            hybrid_hess_atoms  # supplied a list; 1-indexed as per user
        )
        self.numhess = numhess
        self.recalc_hess = recalc_hess
        self.trust_radius = trust_radius
        self.tssearch_type = (
            tssearch_type  # methods for TS search: OptTS, ScanTS
        )
        self.scants_modred = (
            scants_modred  # modred for scanTS (as in a scan job)
        )
        self.full_scan = full_scan  # full scan or not;  do or not abort scan after highest point is reached

    @property
    def route_string(self):
        """
        Get ORCA route string for transition state job.

        Overrides parent property to handle TS-specific keywords.

        Returns:
            str: ORCA route string for TS calculation
        """
        self.jobtype = "ts"
        route_string = self._get_route_string_from_jobtype()
        if self.tssearch_type.lower() == "scants":
            route_string = route_string.replace("OptTS", "ScanTS")
        return route_string


class ORCAIRCJobSettings(ORCAJobSettings):
    """
    Settings for ORCA intrinsic reaction coordinate calculations.

    This class extends ORCAJobSettings with specialized options for
    IRC path following calculations.

    Attributes:
        maxiter (int | None): Maximum number of IRC iterations.
        printlevel (int | None): Verbosity level for IRC output.
        direction (str | None): IRC direction ('both', 'forward', 'backward', 'down').
        inithess (str | None): Initial Hessian specification ('read', 'calc_anfreq', 'calc_numfreq').
        hess_filename (str | None): Hessian filename when inithess='read'.
        hessmode (int | None): Hessian mode index used for initial displacement.
        init_displ (str | None): Initial displacement type ('DE' or 'length').
        scale_init_displ (float | None): Step size for initial displacement.
        de_init_displ (float | None): Target energy difference for initial displacement (Eh or mEh as per ORCA).
        follow_coordtype (str | None): Coordinate type to follow (typically 'cartesian').
        scale_displ_sd (float | None): Scaling factor for the first SD step.
        adapt_scale_displ (bool | None): Adapt SD step scaling dynamically.
        sd_parabolicfit (bool | None): Use parabolic fit to optimize SD step length.
        interpolate_only (bool | None): Restrict parabolic fit to interpolation only.
        do_sd_corr (bool | None): Apply correction to the first SD step.
        scale_displ_sd_corr (float | None): Scaling factor for SD correction step.
        sd_corr_parabolicfit (bool | None): Use parabolic fit for the correction step.
        tolrmsg (float | None): RMS gradient tolerance (a.u.).
        tolmaxg (float | None): Max gradient element tolerance (a.u.).
        monitor_internals (bool | None): Print selected internal coordinates during IRC.
        internal_modred (list | dict | None): Internal coordinates (B/A/D/I) to monitor when monitor_internals=True.
    """

    def __init__(
        self,
        maxiter=None,
        printlevel=None,
        direction=None,
        inithess=None,
        hess_filename=None,
        hessmode=None,
        init_displ=None,
        scale_init_displ=None,
        de_init_displ=None,
        follow_coordtype=None,
        scale_displ_sd=None,
        adapt_scale_displ=None,
        sd_parabolicfit=None,
        interpolate_only=None,
        do_sd_corr=None,
        scale_displ_sd_corr=None,
        sd_corr_parabolicfit=None,
        tolrmsg=None,
        tolmaxg=None,
        monitor_internals=None,
        internal_modred=None,
        **kwargs,
    ):
        """
        Initialize ORCA IRC job settings.

        Args:
            maxiter: Maximum number of IRC iterations
            printlevel: Level of output detail
            direction: IRC direction ('both', 'forward', 'backward', 'down')
            inithess: Initial Hessian specification
            hess_filename: Filename for Hessian file
            hessmode: Hessian mode for initial displacement
            init_displ: Initial displacement type
            scale_init_displ: Scaling for initial displacement
            de_init_displ: Energy difference for initial displacement
            follow_coordtype: Coordinate type to follow
            scale_displ_sd: Scaling for steepest descent displacement
            adapt_scale_displ: Whether to adapt displacement scaling
            sd_parabolicfit: Whether to use parabolic fit for SD step
            interpolate_only: Whether to only allow interpolation
            do_sd_corr: Whether to apply SD correction
            scale_displ_sd_corr: Scaling for SD correction
            sd_corr_parabolicfit: Whether to use parabolic fit for correction
            tolrmsg: RMS gradient tolerance
            tolmaxg: Maximum gradient tolerance
            monitor_internals: Whether to monitor internal coordinates
            internal_modred: Internal coordinates to monitor
            **kwargs: Additional keyword arguments
        """
        super().__init__(**kwargs)
        self.maxiter = maxiter
        self.printlevel = printlevel
        self.direction = direction
        self.inithess = inithess
        self.hess_filename = hess_filename
        self.hessmode = hessmode
        self.init_displ = init_displ
        self.scale_init_displ = scale_init_displ
        self.de_init_displ = de_init_displ
        self.follow_coordtype = follow_coordtype
        self.scale_displ_sd = scale_displ_sd
        self.adapt_scale_displ = adapt_scale_displ
        self.sd_parabolicfit = sd_parabolicfit
        self.interpolate_only = interpolate_only
        self.do_sd_corr = do_sd_corr
        self.scale_displ_sd_corr = scale_displ_sd_corr
        self.sd_corr_parabolicfit = sd_corr_parabolicfit
        self.tolrmsg = tolrmsg
        self.tolmaxg = tolmaxg
        self.monitor_internals = monitor_internals
        self.internal_modred = internal_modred

    @property
    def route_string(self):
        """
        Get ORCA route string for IRC job.

        Overrides parent property, removes frequency calculation
        as it's not compatible with IRC.

        Returns:
            str: ORCA route string for IRC calculation
        """
        self.jobtype = "irc"
        route_string = self._get_route_string_from_jobtype()
        if "freq" in route_string.lower():
            route_string = re.sub(
                r"freq", "", route_string, flags=re.IGNORECASE
            )
        return route_string

    def _write_irc_block(self, f):
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
        irc_settings_keys = ORCAIRCJobSettings().__dict__.keys()
        parent_settings_keys = ORCAJobSettings().__dict__.keys()
        irc_specific_keys = set(irc_settings_keys) - set(parent_settings_keys)

        if not any(
            getattr(self, key) is not None for key in irc_specific_keys
        ):
            return

        # write irc block if any option value is not None:
        f.write("%irc\n")
        for key in irc_specific_keys:
            value = getattr(self, key)
            if value is None:
                continue  # ignore the rest of the code and go to next in the for loop
            # only write into IRC input if the value is not None
            if key == "internal_modred":
                pass  # internal_modred is not an option in ORCA IRC file
            elif key == "inithess":
                f.write(f"  {key} {value}\n")
                if value.lower() == "read":  # if initial hessian is to be read
                    assert (
                        self.hess_filename is not None
                    ), "No Hessian file is given!"
                    assert os.path.exists(
                        self.hess_filename
                    ), f"Hessian file {self.hess_filename} is not found!"
                    f.write(
                        f'  Hess_Filename "{self.hess_filename}"  # Hessian file\n'
                    )
            elif (
                key == "hess_filename"
            ):  # already used/written, if initial hessian is to be read
                pass
            elif key == "monitor_internals":
                if str(value).lower() == "true":
                    f.write(f"  {key}\n")
                    assert (
                        self.internal_modred is not None
                    ), 'No internal modred is specified for IRC job "monitor_intervals" option!'
                    prepend_string_list = (
                        get_prepend_string_list_from_modred_free_format(
                            self.internal_modred, program="orca"
                        )
                    )
                    for prepend_string in prepend_string_list:
                        f.write(f"  {{ {prepend_string} }}\n")
                    f.write("  end\n")
                else:  # monitor_internals has other value (false), then don't write it in input
                    pass
            else:  # all other keys with given values
                f.write(f"  {key} {value}\n")
        f.write("end\n")


class ORCAQMMMJobSettings(ORCAJobSettings):
    """
    Configuration for ORCA multiscale QM/MM calculations.

    Supports five multiscale calculation types:
    1. Additive QM/MM
    2. Subtractive QM/QM2 (2-layer ONIOM)
    3. Subtractive QM/QM2/MM (3-layer ONIOM)
    4. MOL-CRYSTAL-QMMM (molecular crystals)
    5. IONIC-CRYSTAL-QMMM (semiconductors/insulators)

    Attributes:
        jobtype (str): Multiscale calculation type
        high_level_functional (str): DFT functional for high-level (QM) region
        high_level_basis (str): Basis set for high-level (QM) region
        intermediate_level_functional (str): DFT functional for intermediate-level (QM2) region
        intermediate_level_basis (str): Basis set for intermediate-level (QM2) region
        intermediate_level_method (str): Built-in method for intermediate-level (XTB, HF-3C, etc.)
        low_level_force_field (str): Force field for low-level (MM) region
        high_level_atoms (list): Atom indices for high-level (QM) region
        intermediate_level_atoms (list): Atom indices for intermediate-level (QM2) region
        charge_total (int): Total system charge
        mult_total (int): Total system multiplicity
        charge_intermediate (int): Intermediate layer charge
        mult_intermediate (int): Intermediate layer multiplicity
        charge_high (int): High-level region charge
        mult_high (int): High-level region multiplicity
        intermediate_level_solvation (str): Solvation model for intermediate-level region
        active_atoms (list): Active atoms for optimization
        use_active_info_from_pbc (bool): Use PDB active atom info
        optregion_fixed_atoms (list): Fixed atoms in optimization
        high_level_h_bond_length (dict): Custom high-level-H bond distances
        delete_la_double_counting (bool): Remove bend/torsion double counting
        delete_la_bond_double_counting_atoms (bool): Remove bond double counting
        embedding_type (str): Electronic or mechanical embedding
        conv_charges (bool): Use converged charges for crystal QM/MM
        conv_charges_max_n_cycles (int): Max charge convergence cycles
        conv_charges_conv_thresh (float): Charge convergence threshold
        scale_formal_charge_mm_atom (float): MM charge scaling factor
        n_unit_cell_atoms (int): Atoms per unit cell (MOL-CRYSTAL-QMMM)
        ecp_layer_ecp (str): ECP type for boundary region
        ecp_layer (int): Number of ECP layers
        scale_formal_charge_ecp_atom (float): ECP charge scaling factor
    """

    # Class attribute: Built-in methods supported for intermediate level (QM2)
    INTERMEDIATE_LEVEL_BUILT_IN_METHODS = [
        "XTB",
        "XTB0",
        "XTB1",
        "HF-3C",
        "PBEH-3C",
        "R2SCAN-3C",
        "PM3",
        "AM1",
    ]

    def __init__(
        self,
        jobtype=None,  # corresponding to the 5 types of jobs mentioned above
        high_level_functional=None,
        high_level_basis=None,
        intermediate_level_functional=None,
        intermediate_level_basis=None,
        intermediate_level_method=None,
        low_level_force_field=None,  # level-of-theory for MM
        high_level_atoms=None,
        intermediate_level_atoms=None,
        charge_total=None,
        mult_total=None,
        charge_intermediate=None,
        mult_intermediate=None,
        charge_high=None,
        mult_high=None,
        intermediate_level_solvation=None,
        intermediate_solv_scheme=None,
        active_atoms=None,
        use_active_info_from_pbc=False,
        optregion_fixed_atoms=None,
        high_level_h_bond_length=None,  # similar to scale factors in Gaussian ONIOM jobs
        delete_la_double_counting=False,
        delete_la_bond_double_counting_atoms=False,
        embedding_type=None,  # optional
        # the followings are crystal QM/MM parameters
        conv_charges=True,
        conv_charges_max_n_cycles=None,
        conv_charges_conv_thresh=None,
        scale_formal_charge_mm_atom=None,
        n_unit_cell_atoms=None,  # for MOL-CRYSTAL-QMMM jobs
        # the followings are for INONIC-CRYSTAL-QMMM jobs
        ecp_layer_ecp=None,
        ecp_layer=None,
        scale_formal_charge_ecp_atom=None,
        parent_jobtype=None,
        **kwargs,
    ):
        """
        Initialize ORCA QM/MM job settings.

        Args:
            jobtype: Type of multiscale calculation (QMMM, QM/QM2, QM/QM2/MM, etc.)
            high_level_functional: DFT functional for high-level (QM) region
            high_level_basis: Basis set for high-level (QM) region
            intermediate_level_functional: DFT functional for intermediate-level (QM2) region
            intermediate_level_basis: Basis set for intermediate-level (QM2) region
            intermediate_level_method: Built-in method for intermediate-level (XTB, HF-3C, PBEH-3C, etc.)
            low_level_force_field: Force field for low-level (MM) region (MMFF, AMBER, CHARMM, etc.)
            high_level_atoms: Atom indices for high-level (QM) region
            intermediate_level_atoms: Atom indices for intermediate-level (QM2) region
            charge_total: Total system charge
            mult_total: Total system multiplicity
            charge_intermediate: Intermediate layer charge (QM2)
            mult_intermediate: Intermediate layer multiplicity (QM2)
            charge_high: High-level region charge
            mult_high: High-level region multiplicity
            intermediate_level_solvation: Solvation model for intermediate-level (CPCM, SMD, etc.)
            active_atoms: Active atom indices (default: whole system)
            optregion_fixed_atoms: Fixed atom indices in optimization
            high_level_h_bond_length: Custom bond lengths {(atom1, atom2): length}
            delete_la_double_counting: Remove bend/torsion double counting
            delete_la_bond_double_counting_atoms: Remove bond double counting
            embedding_type: Electronic (default) or mechanical embedding
            conv_charges: Use converged charges
            conv_charges_max_n_cycles: Max cycles for charge convergence
            conv_charges_conv_thresh: Convergence threshold for charges
            scale_formal_charge_mm_atom: MM atomic charge scaling factor
            n_unit_cell_atoms: Atoms per unit cell (required for MOL-CRYSTAL-QMMM)
            ecp_layer_ecp: ECP type for boundary region
            ecp_layer: Number of ECP layers around QM region
            scale_formal_charge_ecp_atom: ECP atomic charge scaling factor
        """
        super().__init__(**kwargs)
        self.intermediate_solv_scheme = intermediate_solv_scheme
        self.jobtype = jobtype
        self.parent_jobtype = parent_jobtype
        self.high_level_functional = high_level_functional
        self.high_level_basis = high_level_basis
        self.intermediate_level_functional = intermediate_level_functional
        self.intermediate_level_basis = intermediate_level_basis
        self.intermediate_level_method = intermediate_level_method
        self.low_level_force_field = low_level_force_field
        self.high_level_atoms = high_level_atoms
        self.intermediate_level_atoms = intermediate_level_atoms
        self.charge_total = charge_total
        self.mult_total = mult_total
        self.charge_intermediate = charge_intermediate
        self.mult_intermediate = mult_intermediate
        self.charge_high = charge_high
        self.mult_high = mult_high
        self.intermediate_level_solvation = intermediate_level_solvation
        self.active_atoms = active_atoms
        self.use_active_info_from_pbc = use_active_info_from_pbc
        self.optregion_fixed_atoms = optregion_fixed_atoms
        self.high_level_h_bond_length = high_level_h_bond_length
        self.delete_la_double_counting = delete_la_double_counting
        self.delete_la_bond_double_counting_atoms = (
            delete_la_bond_double_counting_atoms
        )
        self.embedding_type = embedding_type
        # Crystal QM/MM parameters
        self.conv_charges = conv_charges
        self.conv_charges_max_n_cycles = conv_charges_max_n_cycles
        self.conv_charges_conv_thresh = conv_charges_conv_thresh
        self.scale_formal_charge_mm_atom = scale_formal_charge_mm_atom
        self.n_unit_cell_atoms = n_unit_cell_atoms  # For MOL-CRYSTAL-QMMM
        # Ionic crystal QM/MM parameters
        self.ecp_layer_ecp = ecp_layer_ecp
        self.ecp_layer = ecp_layer
        self.scale_formal_charge_ecp_atom = scale_formal_charge_ecp_atom

        # Set parent class attributes from high-level (QM) region
        self.functional = self.high_level_functional
        self.basis = self.high_level_basis

        # Set charge/multiplicity with fallback priority: high -> intermediate -> total
        if self.charge_high is not None and self.mult_high is not None:
            self.charge = self.charge_high
            self.multiplicity = self.mult_high
        elif (
            self.charge_intermediate is not None
            and self.mult_intermediate is not None
        ):
            # the charge/multiplicity of the intermediate system corresponds to the
            # sum of the charge/multiplicity of the high level and low level regions
            self.charge = self.charge_intermediate
            self.multiplicity = self.mult_intermediate
        elif self.charge_total is not None and self.mult_total is not None:
            self.charge = self.charge_total
            self.multiplicity = self.mult_total
        else:
            self.charge = None
            self.multiplicity = None

        # Validate intermediate-level parameters match job type
        self._validate_intermediate_parameters()

    def re_init_and_validate(self):
        """Recompute derived fields and rerun validation after overrides."""
        # Update inherited level-of-theory aliases
        self.functional = self.high_level_functional
        self.basis = self.high_level_basis

        if self.charge_high is not None and self.mult_high is not None:
            self.charge = self.charge_high
            self.multiplicity = self.mult_high
        elif (
            self.charge_intermediate is not None
            and self.mult_intermediate is not None
        ):
            self.charge = self.charge_intermediate
            self.multiplicity = self.mult_intermediate
        elif self.charge_total is not None and self.mult_total is not None:
            self.charge = self.charge_total
            self.multiplicity = self.mult_total
        else:
            self.charge = None
            self.multiplicity = None

        self._validate_intermediate_parameters()

    def _validate_intermediate_parameters(self):
        """
        Validate that intermediate-level parameters are provided when needed and not provided when not needed.

        Raises:
            ValueError: If QM2 layer is required but intermediate parameters are missing,
                       or if intermediate parameters are provided but QM2 layer is not required.
        """
        # Check if intermediate parameters are provided
        has_intermediate_params = (
            self.intermediate_level_functional is not None
            and self.intermediate_level_basis is not None
        ) or self.intermediate_level_method is not None
        if has_intermediate_params:
            missing = [
                name
                for name, value in {
                    "intermediate_level_atoms": self.intermediate_level_atoms,
                    "charge_intermediate": self.charge_intermediate,
                    "mult_intermediate": self.mult_intermediate,
                }.items()
                if value is None
            ]

            if missing:
                raise ValueError(
                    "When intermediate-level theory is specified, the following "
                    f"parameters must also be provided: {', '.join(missing)}"
                )

        # Job types that require QM2 (intermediate) layer
        qm2_required_jobtypes = ["QM/QM2", "QM/QM2/MM"]

        # Check if this job type requires QM2 layer
        requires_qm2 = self.jobtype and self.jobtype.upper() in [
            jt.upper() for jt in qm2_required_jobtypes
        ]

        if requires_qm2 and not has_intermediate_params:
            raise ValueError(
                f"Job type '{self.jobtype}' requires QM2 (intermediate) layer, but no intermediate-level "
                f"parameters were provided. Please specify at least one of:\n"
                f"  - intermediate_level_functional and intermediate_level_basis\n"
                f"  - intermediate_level_method (e.g., 'XTB', 'HF-3C')\n"
                f"  - intermediate_level_atoms\n"
                f"  - charge_intermediate and mult_intermediate"
            )

        if not requires_qm2 and has_intermediate_params:
            provided_params = []
            if self.intermediate_level_functional is not None:
                provided_params.append("intermediate_level_functional")
            if self.intermediate_level_basis is not None:
                provided_params.append("intermediate_level_basis")
            if self.intermediate_level_method is not None:
                provided_params.append("intermediate_level_method")
            if self.intermediate_level_atoms is not None:
                provided_params.append("intermediate_level_atoms")
            if self.charge_intermediate is not None:
                provided_params.append("charge_intermediate")
            if self.mult_intermediate is not None:
                provided_params.append("mult_intermediate")
            if self.intermediate_level_solvation is not None:
                provided_params.append("intermediate_level_solvation")

            raise ValueError(
                f"Job type '{self.jobtype}' does not require QM2 (intermediate) layer, but "
                f"intermediate-level parameters were provided: {', '.join(provided_params)}.\n"
                f"QM2 layer is only used in job types: {', '.join(qm2_required_jobtypes)}.\n"
                f"Either:\n"
                f"  - Change job type to 'QM/QM2' or 'QM/QM2/MM', OR\n"
                f"  - Remove the intermediate-level parameters"
            )

    @property
    def qmmm_route_string(self):
        return self._get_level_of_theory_string()

    @property
    def qmmm_block(self):
        return self._write_qmmm_block()

    def validate_and_assign_level(
        self, functional, basis, built_in_method, level_name
    ):
        """
        Validate and assign level of theory for high/intermediate/low-level layers.

        Args:
            functional: DFT functional
            basis: Basis set
            built_in_method: Built-in ORCA method (XTB, HF-3C, etc.)
            level_name: Layer name (high_level, intermediate_level, low_level)

        Returns:
            str: Validated level of theory string

        Raises:
            ValueError: If incompatible options are specified
        """
        level_of_theory = ""
        if functional and basis and built_in_method:
            raise ValueError(
                f"For {level_name} level: specify either functional/basis OR built-in method, not both!"
            )
        if built_in_method:
            assert functional is None and basis is None, (
                f"Built-in method specified for {level_name} level - "
                f"functional and basis should not be provided!"
            )
            if (
                level_name == "intermediate_level"
                and built_in_method.upper()
                in self.INTERMEDIATE_LEVEL_BUILT_IN_METHODS
            ):
                level_of_theory = built_in_method
        elif functional and basis:
            if level_name == "intermediate_level":
                level_of_theory = "QM2"
            else:
                level_of_theory = f"{functional} {basis}"
        if level_name == "low_level":
            level_of_theory = "MM"
        logger.debug(f"Level of theory for {level_name}: {level_of_theory}")
        return level_of_theory

    def check_crystal_qmmm(self):
        """
        Validate crystal QM/MM job settings.

        Ensures required parameters are set for MOL-CRYSTAL-QMMM
        and IONIC-CRYSTAL-QMMM calculations.

        Raises:
            AssertionError: If required parameters are missing or invalid
        """
        jobtype = self.job_type.upper()
        if jobtype in ["IONIC-CRYSTAL-QMMM", "MOL-CRYSTAL-QMMM"]:
            assert (
                self.mult_high is None
                and self.mult_intermediate is None
                and self.mult_total is None
            ), f"Multiplicity should not be specified for {jobtype} job!"
            self.multiplicity = 0  # avoid conflicts from parent class
            if self.conv_charges is False:
                assert (
                    self.low_level_force_field is not None
                ), "Force field file containing convergence charges is not provided!"
            if jobtype == "MOL-CRYSTAL-QMMM":
                assert (
                    self.n_unit_cell_atoms
                ), f"The number of atoms per molecular subunit for {jobtype} job is not provided!"
            else:
                assert (
                    self.ecp_layer_ecp
                ), f"cECPs used for the boundary region for {jobtype} job must be specified! "
                assert (
                    self.n_unit_cell_atoms is None
                ), "The number of atoms per molecular subunit is only applicable to MOL-CRYSTAL-QMMM!"

    def _get_level_of_theory_string(self):
        """
        Generate ORCA route string for QM/MM calculations.

        Returns:
            str: Route string (e.g., '!QM/XTB', '!QM/HF-3C/MM', '!QMMM')
        """
        if (
            self.jobtype.upper() == "IONIC-CRYSTAL-QMMM"
            or self.jobtype.upper() == "MOL-CRYSTAL-QMMM"
        ):
            level_of_theory = f"! {self.jobtype.upper()}"
        else:
            level_of_theory = "!"
            parent_jobtype = self.parent_jobtype.lower()
            if parent_jobtype in ("opt", "modred", "scan"):
                level_of_theory += " Opt"
            elif parent_jobtype == "ts":
                level_of_theory += " OptTS"
            elif parent_jobtype == "irc":
                level_of_theory += " IRC"
            elif parent_jobtype == "sp":
                level_of_theory += ""
            level_of_theory += " QM"
            self.high_level_level_of_theory = self.validate_and_assign_level(
                self.high_level_functional,
                self.high_level_basis,
                None,
                level_name="high_level",
            )
            self.intermediate_level_level_of_theory = (
                self.validate_and_assign_level(
                    self.intermediate_level_functional,
                    self.intermediate_level_basis,
                    self.intermediate_level_method,
                    level_name="intermediate_level",
                )
            )
            self.low_level_level_of_theory = self.validate_and_assign_level(
                None, None, self.low_level_force_field, level_name="low_level"
            )
            # only "!QMMM" will be used for additive QMMM
            level_of_theory += f"/{self.intermediate_level_level_of_theory}"
            if self.low_level_level_of_theory is not None:
                if self.jobtype.upper() == "QMMM":
                    level_of_theory = "!QMMM"  # Additive QM/MM
                else:
                    level_of_theory += f"/{self.low_level_level_of_theory}"
            level_of_theory += f" {self.high_level_level_of_theory}"
            if self.solvent_model is not None:
                level_of_theory += f" {self.solvent_model}"
            self.intermediate_level_method = (
                self.intermediate_level_method or ""
            ).lower()
            if self.intermediate_level_solvation:
                if self.intermediate_level_solvation in [
                    "alpb(water)",
                    "ddcosmo(water)",
                    "cpcmx(water)",
                ]:
                    assert self.intermediate_level_method.lower() == "xtb", (
                        "The intermediate-level solvation models ALPB, DDCOSMO, CPCMX are only "
                        "compatible with XTB method!"
                    )
                    level_of_theory += f" {self.intermediate_level_solvation}"
                elif (
                    self.intermediate_level_solvation.lower() == "cpcm(water)"
                ):
                    level_of_theory += f" {self.intermediate_level_solvation}"
        return level_of_theory

    def _get_h_bond_length(self):
        """
        Generate custom high-level-H bond length specifications.

        Returns:
            str: Bond length specifications for %qmmm block

        Example output:
            Dist_C_HLA 1.09
            Dist_O_HLA 0.98
            or
            H_Dist_FileName "QM_H_dist.txt"
        """
        if isinstance(self.high_level_h_bond_length, dict):
            h_bond_length = ""
            for (
                atom_pair,
                bond_length,
            ) in self.high_level_h_bond_length.items():
                h_bond_length += (
                    f"Dist_{atom_pair[0]}_{atom_pair[1]} {bond_length}\n"
                )
            return h_bond_length
        elif isinstance(self.high_level_h_bond_length, str):
            # if the user provided a file with the d0_X-H values
            assert os.path.exists(
                self.high_level_h_bond_length
            ), f"File {self.high_level_h_bond_length} does not exist!"
            return f'H_Dist_FileName "{self.high_level_h_bond_length}"'

    def _get_embedding_type(self):
        """
        Generate embedding type specification for QM/MM calculations.

        Returns:
            str: Embedding directive (e.g., "Embedding Electronic")

        Raises:
            ValueError: If invalid embedding type specified
        """
        embedding_type = "Embedding "
        if self.embedding_type.lower() == "electronic":
            embedding_type += "Electronic"
        elif self.embedding_type.lower() == "mechanical":
            embedding_type += "Mechanical"
        else:
            raise ValueError(
                f"Invalid embedding type: {self.embedding_type}. "
                "Valid options are 'Electronic' or 'Mechanical'."
            )
        return embedding_type

    def _get_charge_and_multiplicity(self):
        """
        Generate charge and multiplicity lines for %qmmm block.

        Returns charge/multiplicity for total system (QM/MM) or intermediate system (QM/QM2/MM).
        ORCA specifies total/intermediate charge/multiplicity in %qmmm block while
        high-level charge/multiplicity are specified in the coordinate section.

        Returns:
            tuple: (charge_str, mult_str) for %qmmm block
        """
        if self.intermediate_level_atoms is not None:
            charge = f"Charge_Intermediate {self.charge_intermediate}"
            mult = f"Mult_Intermediate {self.mult_intermediate}"
        else:
            charge = f"Charge_Total {self.charge_total}"
            mult = f"Mult_Total {self.mult_total}"
        return charge, mult

    def _get_partition_string(self):
        """
        Generate atom partition specifications for high-level and intermediate-level regions.

        Returns:
            str: Partition block with QMAtoms and QM2Atoms directives

        Example:
            QMAtoms {1:10 15:20} end
            QM2Atoms {21:30} end
        """

        partition_string = ""
        high_fmt = self._get_formatted_partition_strings(self.high_level_atoms)
        if high_fmt is not None:
            partition_string += f"QMAtoms {{{high_fmt}}} end\n"
        intermediate_fmt = self._get_formatted_partition_strings(
            self.intermediate_level_atoms
        )
        if intermediate_fmt is not None:
            partition_string += f"QM2Atoms {{{intermediate_fmt}}} end\n"
        return partition_string

    def _get_formatted_partition_strings(self, atom_id):
        """
        Convert atom specifications to ORCA-formatted range strings.

        Normalizes various atom specification formats to sorted unique integers
        and compresses contiguous sequences into compact range notation.

        Args:
            atom_id: Atom specification (list/tuple of ints, string ranges, etc.)

        Returns:
            str: Formatted string with ranges and individual atoms

        Examples:
            Input: [1,2,3,5,7,8,9] → Output: "1:3 5 7:9"
            Input: "1-15,37,39" → Output: "1:15 37 39"
        """
        if atom_id is None:
            return None

        # If it's already a sequence of integers (list/tuple), coerce to ints
        if isinstance(atom_id, (list, tuple)):
            atoms = [int(x) for x in atom_id]
        elif isinstance(atom_id, str):
            # Prefer robust utility if available
            try:
                atoms = get_list_from_string_range(atom_id)
            except Exception:
                # Fallback simple parser: split on commas/spaces, handle ranges like '1-5' or '1:5'
                tokens = re.split(r"[,\s]+", atom_id.strip())
                atoms = []
                for tok in tokens:
                    if not tok:
                        continue
                    if "-" in tok:
                        a, b = tok.split("-", 1)
                        atoms.extend(range(int(a), int(b) + 1))
                    elif ":" in tok:
                        a, b = tok.split(":", 1)
                        atoms.extend(range(int(a), int(b) + 1))
                    else:
                        atoms.append(int(tok))
        else:
            # Try to iterate and coerce
            try:
                atoms = [int(x) for x in atom_id]
            except Exception:
                raise TypeError(
                    "Unsupported atom specification type for partition string"
                )

        # make unique sorted list of ints
        atoms = sorted(set(int(x) for x in atoms))
        if not atoms:
            return None

        # compress contiguous sequences into start:end or single integer
        ranges = []
        start = prev = atoms[0]
        for a in atoms[1:]:
            if a == prev + 1:
                prev = a
                continue
            # break in sequence
            if start == prev:
                ranges.append(str(start))
            else:
                ranges.append(f"{start}:{prev}")
            start = prev = a
        # finalize last segment
        if start == prev:
            ranges.append(str(start))
        else:
            ranges.append(f"{start}:{prev}")

        return " ".join(ranges)

    def _write_qmmm_block(self):
        """
        Generate complete %qmmm block for ORCA multiscale calculations.

        Constructs the full %qmmm block containing all necessary parameters
        for QM/MM calculations including atom partitions, charges, force fields,
        embedding options, and crystal-specific parameters.

        Returns:
            str: Complete %qmmm block ready for ORCA input file

        Example output:
            %qmmm
            QMAtoms {16:33 68:82} end
            QM2Atoms {0:12 83:104} end
            Charge_Medium 0
            ORCAFFFilename "system.prms"
            end
        """
        full_qm_block = "%qmmm\n"

        # Add atom partition specifications
        partition_block = self._get_partition_string()
        if partition_block:
            full_qm_block += partition_block

        # Add charge/multiplicity for medium or total system
        charge_str, mult_str = self._get_charge_and_multiplicity()
        # Only include lines that do not contain the string 'None'
        if charge_str and "None" not in str(charge_str):
            full_qm_block += f"{charge_str}\n"
        if mult_str and "None" not in str(mult_str):
            full_qm_block += f"{mult_str}\n"

        # Add intermediate-level solvation if specified
        if self.intermediate_solv_scheme is not None:
            full_qm_block += f"solv_scheme {self.intermediate_solv_scheme}\n"

        # Add active atoms specification
        active_fmt = self._get_formatted_partition_strings(self.active_atoms)
        if active_fmt is not None:
            full_qm_block += f"ActiveAtoms {{{active_fmt}}} end\n"

        # Add intermediate-level method/basis specifications
        if (
            self.intermediate_level_method is not None
            and self.intermediate_level_method.strip()
        ):
            # Only write QM2CustomFile if method is NOT in the built-in list
            if (
                self.intermediate_level_method.upper()
                not in self.INTERMEDIATE_LEVEL_BUILT_IN_METHODS
            ):
                # intermediate-level method provided via external file
                full_qm_block += (
                    f'QM2CustomFile "{self.intermediate_level_method}" end\n'
                )
        else:
            # If custom intermediate-level functional/basis provided, include them
            if (
                self.intermediate_level_functional is not None
                and self.intermediate_level_functional.strip()
            ):
                full_qm_block += f'QM2CUSTOMMETHOD "{self.intermediate_level_functional}" end\n'
            if (
                self.intermediate_level_basis is not None
                and self.intermediate_level_basis.strip()
            ):
                full_qm_block += (
                    f'QM2CUSTOMBASIS "{self.intermediate_level_basis}" end\n'
                )

        # Add force field for specific job types
        # assert (
        #         self.low_level_force_field is not None
        # ), f"Force field file missing for {self.jobtype} job!"
        if self.jobtype and self.jobtype.upper() in [
            "QMMM" "QM/MM",
            "QM/QM2/MM",
            "IONIC-CRYSTAL-QMMM",
        ]:
            assert (
                self.low_level_force_field is not None
            ), f"Force field file missing for {self.jobtype} job!"
            full_qm_block += f'ORCAFFFilename "{self.low_level_force_field}"\n'

        # Fixed atoms options
        if self.use_active_info_from_pbc:
            full_qm_block += "Use_Active_InfoFromPDB true\n"
        elif self.optregion_fixed_atoms:
            fixed_fmt = self._get_formatted_partition_strings(
                self.optregion_fixed_atoms
            )
            if fixed_fmt is not None:
                full_qm_block += f"OptRegion_FixedAtoms {{{fixed_fmt}}} end\n"

        # Add custom H-bond lengths
        if self.high_level_h_bond_length is not None:
            h_block = self._get_h_bond_length()
            if h_block:
                # _get_h_bond_length may return a multi-line string
                full_qm_block += (
                    f"{h_block}\n" if not h_block.endswith("\n") else h_block
                )

        # Add double-counting removal options
        if self.delete_la_double_counting is True:
            full_qm_block += "Delete_LA_Double_Counting true\n"
        if self.delete_la_bond_double_counting_atoms:
            full_qm_block += "DeleteLABondDoubleCounting true\n"
        if self.embedding_type is not None:
            full_qm_block += f"{self._get_embedding_type()}\n"

        # Add crystal QM/MM parameters if any
        crystal_sub = self._write_crystal_qmmm_subblock()
        if crystal_sub is not None:
            full_qm_block += crystal_sub

        full_qm_block += "end\n"
        return full_qm_block

    def _write_crystal_qmmm_subblock(self):
        """
        Generate crystal-specific parameters for %qmmm block.

        Creates parameter block for MOL-CRYSTAL-QMMM and IONIC-CRYSTAL-QMMM
        calculations including charge convergence, ECP settings, and scaling factors.

        Returns:
            str: Crystal QM/MM parameter block
        """
        crystal_qmmm_subblock = ""
        if not self.conv_charges:
            crystal_qmmm_subblock += "Conv_Charges False \n"
            crystal_qmmm_subblock += (
                f'ORCAFFFilename "{self.low_level_force_field}" \n'
            )
        if self.conv_charges_max_n_cycles is not None:
            crystal_qmmm_subblock += (
                f"Conv_Charges_MaxNCycles {self.conv_charges_max_n_cycles} \n"
            )
        if self.conv_charges_conv_thresh is not None:
            crystal_qmmm_subblock += (
                f"Conv_Charges_ConvThresh {self.conv_charges_conv_thresh} \n"
            )
        if self.scale_formal_charge_mm_atom is not None:
            crystal_qmmm_subblock += f"Scale_FormalCharge_MMAtom {self.scale_formal_charge_mm_atom} \n"
        if self.n_unit_cell_atoms is not None:
            crystal_qmmm_subblock += (
                f"NumUnitCellAtoms {self.n_unit_cell_atoms} \n"
            )
        if self.ecp_layer_ecp is not None:
            crystal_qmmm_subblock += f"cECPs {self.ecp_layer_ecp} \n"
        if self.ecp_layer is not None:
            crystal_qmmm_subblock += f"ECPLayers {self.ecp_layer} \n"
        if self.scale_formal_charge_ecp_atom is not None:
            crystal_qmmm_subblock += f"Scale_FormalCharge_ECPAtom {self.scale_formal_charge_ecp_atom} \n"
        return crystal_qmmm_subblock

    def __eq__(self, other):
        """
        Compare two ORCAQMMMJobSettings objects for equality.

        Compares all attributes between two ORCA QMMM settings objects, including the
        QMMM-specific attributes that are not present in the parent class.

        Args:
            other (ORCAQMMMJobSettings): Settings object to compare with.

        Returns:
            bool or NotImplemented: True if equal, False if different,
                NotImplemented if types don't match.
        """
        if type(self) is not type(other):
            return NotImplemented

        # Get dictionaries of both objects
        self_dict = self.__dict__.copy()
        other_dict = other.__dict__.copy()

        # Exclude append_additional_info from the comparison (inherited behavior)
        self_dict.pop("append_additional_info", None)
        other_dict.pop("append_additional_info", None)

        return self_dict == other_dict
