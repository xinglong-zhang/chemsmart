"""
Gaussian job configuration and settings management.

This module provides the GaussianJobSettings class for configuring
Gaussian computational chemistry calculations. It handles all aspects
of job setup including calculation methods, basis sets, solvent models,
optimization parameters, and advanced options like GenECP definitions.

The settings class supports merging configurations, validation of
parameter combinations, and conversion to/from various input formats.
"""

import copy
import logging
import os
import re

from chemsmart.io.gaussian import GAUSSIAN_SOLVATION_MODELS
from chemsmart.io.gaussian.gengenecp import GenGenECPSection
from chemsmart.jobs.settings import MolecularJobSettings
from chemsmart.utils.periodictable import PeriodicTable
from chemsmart.utils.repattern import (
    gaussian_freq_keywords_pattern,
    gaussian_opt_keywords_pattern,
    multiple_spaces_pattern,
)

pt = PeriodicTable()


logger = logging.getLogger(__name__)


class GaussianJobSettings(MolecularJobSettings):
    """
    Configuration settings for Gaussian computational chemistry jobs.

    Manages all parameters needed to configure Gaussian calculations
    including quantum chemistry methods, basis sets, solvation models,
    optimization options, and advanced features like pseudopotentials.
    Provides validation and merging capabilities for job configurations.

    Inherits common calculation fields from `MolecularJobSettings`, such as
    `ab_initio`, `functional`, `basis`, `semiempirical`, `charge`,
    `multiplicity`, `jobtype`, `title`, `freq`, `numfreq`, `solvent_model`,
    `solvent_id`, `additional_solvent_options`, `additional_opt_options_in_route`,
    `append_additional_info`, `gen_genecp_file`, `heavy_elements`,
    `heavy_elements_basis`, `light_elements_basis`, `custom_solvent`, `forces`,
    and `input_string`.

    Attributes:
        chk (bool): Whether to use checkpoint files.
        dieze_tag (str): Calculation level specification.
        additional_solvent_options (str): Extra solvent parameters.
        additional_opt_options_in_route (str): Extra optimization options.
        append_additional_info (str): Additional input file content.
        gen_genecp_file (str): Path to GenECP definition file.
    """

    def __init__(
        self,
        ab_initio=None,
        functional=None,
        basis=None,
        semiempirical=None,
        charge=None,
        multiplicity=None,
        chk=True,
        jobtype=None,
        title=None,
        freq=False,
        numfreq=False,
        dieze_tag=None,
        solvent_model=None,
        solvent_id=None,
        additional_solvent_options=None,
        additional_opt_options_in_route=None,
        additional_route_parameters=None,
        route_to_be_written=None,
        modred=None,
        gen_genecp_file=None,
        heavy_elements=None,
        heavy_elements_basis=None,
        light_elements_basis=None,
        custom_solvent=None,
        append_additional_info=None,
        forces=False,
        input_string=None,
        **kwargs,
    ):
        """
        Initialize Gaussian job settings with calculation parameters.

        Sets up all configuration options for Gaussian calculations
        including validation of parameter combinations and path expansion.

        Args:
            ab_initio (str, optional): Ab initio method (e.g., 'HF', 'MP2').
            functional (str, optional): DFT functional (e.g., 'B3LYP').
            basis (str, optional): Basis set (e.g., '6-31G*').
            semiempirical (str, optional): Semi-empirical method.
            charge (int, optional): Molecular charge.
            multiplicity (int, optional): Spin multiplicity.
            chk (bool): Use checkpoint files (default True).
            jobtype (str, optional): Calculation type (e.g., 'opt', 'freq').
            title (str, optional): Job title.
            freq (bool): Include frequency calculations.
            numfreq (bool): Use numerical frequencies.
            dieze_tag (str, optional): Calculation level tag.
            solvent_model (str, optional): Solvation model.
            solvent_id (str, optional): Solvent identifier.
            additional_solvent_options (str, optional): Extra solvent options.
            additional_opt_options_in_route (str, optional): Extra opt options.
            additional_route_parameters (str, optional): Extra route params.
            route_to_be_written (str, optional): Custom route string.
            modred (dict, optional): Modified redundant coordinates.
            gen_genecp_file (str, optional): GenECP definition file path.
            heavy_elements (list, optional): Heavy element symbols.
            heavy_elements_basis (str, optional): Heavy element basis set.
            light_elements_basis (str, optional): Light element basis set.
            custom_solvent (dict, optional): Custom solvent parameters.
            append_additional_info (str, optional): Additional input content.
            forces (bool): Calculate forces.
            input_string (str, optional): Custom input string.
            **kwargs: Additional keyword arguments.

        Raises:
            ValueError: If incompatible options are specified (freq + forces).
        """
        super().__init__(
            ab_initio=ab_initio,
            functional=functional,
            basis=basis,
            semiempirical=semiempirical,
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
        self.chk = chk
        self.dieze_tag = dieze_tag
        self.additional_solvent_options = additional_solvent_options
        self.additional_opt_options_in_route = additional_opt_options_in_route
        self.append_additional_info = append_additional_info
        self._route_string = None

        if gen_genecp_file is not None and "~" in gen_genecp_file:
            gen_genecp_file = os.path.expanduser(gen_genecp_file)
            logger.debug(f"Expanded GenECP file path: {gen_genecp_file}")
        self.gen_genecp_file = gen_genecp_file

        # Validate that frequency and force calculations are not both enabled
        if forces is True and freq is True:
            logger.error(
                "Attempted to enable both frequency and force calculations"
            )
            raise ValueError(
                "Frequency and Force calculations cannot be performed by "
                "Gaussian at the same time!\n"
                'Such an input file will give "Illegal IType or MSType '
                'generated by parse." error.'
            )

    @property
    def genecp(self):
        """
        Check if GenECP pseudopotentials are configured.

        Determines whether the job uses GenECP pseudopotential
        definitions either from a file or heavy element specifications.

        Returns:
            bool: True if GenECP is configured, False otherwise.
        """
        return (
            self.gen_genecp_file is not None or self.heavy_elements is not None
        )

    def merge(
        self,
        other,
        keywords=("charge", "multiplicity", "title"),
        merge_all=False,
    ):
        """
        Merge current settings with another settings object.

        Combines settings from two sources, with the other settings
        taking precedence. Supports selective merging of specific
        keywords or complete merging of all parameters.

        Args:
            keywords (list): Specific list of keywords to merge.
                Defaults to charge and multiplicity.
                If None, all settings will be merged (Caution: may cause issue if e.g.,
                genecp log file used to prepare input without genecp).
            other (JobSettings, dict): Settings to merge. Can also take the form of a dictionary
            merge_all (bool): If True, merge all settings.
            If False, only merge the settings specified in keywords.
        """
        other_dict = other if isinstance(other, dict) else other.__dict__

        if merge_all:
            # Update self with other for all
            logger.debug(
                f"Merging all settings from {type(other).__name__} "
                f"into {type(self).__name__}"
            )
            merged_dict = self.__dict__.copy()
            merged_dict.update(other_dict)
            return type(self)(**merged_dict)

        if keywords is not None:
            logger.debug(f"Merging specific keywords: {keywords}")
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

        Returns a completely independent copy of the settings
        that can be modified without affecting the original.

        Returns:
            GaussianJobSettings: Deep copy of this settings object.
        """
        return copy.deepcopy(self)

    def __getitem__(self, key):
        """
        Get setting value by key using dictionary-style access.

        Allows settings to be accessed like a dictionary for
        compatibility with existing code patterns.

        Args:
            key (str): Setting name to retrieve.

        Returns:
            Any: Value of the requested setting.
        """
        return self.__dict__[key]

    def __eq__(self, other):
        """
        Compare two settings objects for equality.

        Checks if all attributes are equal between two settings objects,
        excluding the append_additional_info field from comparison.

        Args:
            other (GaussianJobSettings): Settings object to compare with.

        Returns:
            bool or NotImplemented: True if equal, False if different,
                NotImplemented if types don't match.
        """
        if type(self) is not type(other):
            return NotImplemented

        # Exclude append_additional_info from the comparison
        self_dict = self.__dict__
        self_dict.pop("append_additional_info")

        other_dict = other.__dict__
        other_dict.pop("append_additional_info")

        is_equal = self_dict == other_dict
        if not is_equal:
            import dictdiffer

            logger.info("Gaussian job settings are not equal.")
            for diff in list(dictdiffer.diff(self_dict, other_dict)):
                logger.info(f"Difference: {diff}")
        return self_dict == other_dict

    @classmethod
    def from_comfile(cls, filename):
        """
        Create Gaussian settings from a Gaussian .com input file.

        Parses a Gaussian input file to extract calculation settings
        and returns a configured GaussianJobSettings object.

        Args:
            filename (str): Path to the .com file to be parsed.

        Returns:
            GaussianJobSettings: Settings object extracted from input file.
        """
        from chemsmart.io.gaussian.input import Gaussian16Input

        com_path = os.path.abspath(filename)
        gaussian_settings_from_comfile = Gaussian16Input(
            filename=com_path
        ).read_settings()
        return gaussian_settings_from_comfile

    @classmethod
    def from_inpfile(cls, filename):
        """
        Create Gaussian settings from an ORCA .inp input file.

        Converts ORCA input file settings to equivalent Gaussian
        settings by parsing the file and merging with default values.

        Args:
            filename (str): Path to the .inp file to be parsed.

        Returns:
            GaussianJobSettings: Gaussian-compatible settings object.
        """
        from chemsmart.io.orca.input import ORCAInput

        inp_path = os.path.abspath(filename)
        logger.info(f"Return Settings object from inp file: {inp_path}")
        orca_settings_from_inpfile = ORCAInput(
            filename=inp_path
        ).read_settings()
        gaussian_default_settings = cls.default()
        gaussian_settings_from_inpfile = gaussian_default_settings.merge(
            orca_settings_from_inpfile, merge_all=True
        )
        logger.info(
            f"with settings: {gaussian_settings_from_inpfile.__dict__}"
        )
        return gaussian_settings_from_inpfile

    @classmethod
    def from_logfile(cls, filename):
        """
        Create Gaussian settings from a Gaussian .log output file.

        Extracts calculation settings from a completed Gaussian
        log file, including support for periodic boundary conditions.

        Args:
            filename (str): Path to the .log file to be parsed.

        Returns:
            GaussianJobSettings: Settings object extracted from log file.
        """
        log_path = os.path.abspath(filename)
        from chemsmart.io.gaussian.output import (
            Gaussian16Output,
            Gaussian16OutputWithPBC,
        )

        logger.info(f"Return Settings object from logfile: {log_path}")
        try:
            settings = Gaussian16Output(filename=log_path).read_settings()
        except ValueError:
            settings = Gaussian16OutputWithPBC(
                filename=log_path
            ).read_settings()

        return settings

    @classmethod
    def from_outfile(cls, filename):
        """
        Create Gaussian settings from an ORCA .out output file.

        Converts ORCA output file settings to equivalent Gaussian
        settings by parsing the file and merging with defaults.

        Args:
            filename (str): Path to the ORCA .out file to be parsed.

        Returns:
            GaussianJobSettings: Gaussian-compatible settings object.
        """
        from chemsmart.io.orca.output import ORCAOutput

        out_path = os.path.abspath(filename)
        logger.info(
            f"Return Settings object from ORCA .out filename: {out_path}"
        )
        orca_settings_from_outfile = ORCAOutput(
            filename=out_path
        ).read_settings()
        gaussian_default_settings = cls.default()
        gaussian_settings_from_outfile = gaussian_default_settings.merge(
            orca_settings_from_outfile, merge_all=True
        )
        logger.info(
            f"with settings: {gaussian_settings_from_outfile.__dict__}"
        )
        return gaussian_settings_from_outfile

    @classmethod
    def default(cls):
        """
        Create default Gaussian job settings.

        Returns a GaussianJobSettings object with sensible default
        values for common calculations including frequency analysis.

        Returns:
            GaussianJobSettings: Default settings configuration.
        """
        return cls(
            ab_initio=None,
            functional=None,
            basis=None,
            semiempirical=None,
            charge=None,
            multiplicity=None,
            chk=True,
            jobtype=None,
            title="Gaussian job with default settings",
            freq=True,
            numfreq=False,
            dieze_tag=None,
            solvent_model=None,
            solvent_id=None,
            additional_solvent_options=None,
            additional_opt_options_in_route=None,
            additional_route_parameters=None,
            route_to_be_written=None,
            modred=None,
            gen_genecp_file=None,
            heavy_elements=None,
            heavy_elements_basis=None,
            light_elements_basis=None,
            custom_solvent=None,
            append_additional_info=None,
            forces=False,
            input_string=None,
        )

    @classmethod
    def from_filepath(cls, filepath, **kwargs):
        """
        Create settings from a file based on its extension.

        Automatically determines the file type and uses the appropriate
        parser to create settings from the input file.

        Args:
            filepath (str): Path to input file (.com, .gjf, .inp, .log).
            **kwargs: Additional keyword arguments for file parsing.

        Returns:
            GaussianJobSettings: Settings object from the input file.

        Raises:
            ValueError: If file extension is not supported.
        """
        if filepath.endswith((".com", ".gjf")):
            logger.debug(
                f"Loading settings from Gaussian input file: {filepath}"
            )
            return cls.from_comfile(filepath)
        if filepath.endswith(".inp"):
            logger.debug(f"Loading settings from ORCA input file: {filepath}")
            return cls.from_inpfile(filepath)
        if filepath.endswith(".log"):
            logger.debug(
                f"Loading settings from Gaussian log file: {filepath}"
            )
            return cls.from_logfile(filepath)
        raise ValueError(f"Could not create {cls} from {filepath}")

    @property
    def route_string(self):
        """
        Generate the Gaussian route string for the calculation.

        Constructs the complete route line (#-line) for the Gaussian
        input file based on the configured job settings.

        Returns:
            str: Complete route string for Gaussian input file.
        """
        if self.route_to_be_written is not None:
            route_string = self._get_route_string_from_user_input()
        else:
            route_string = self._get_route_string_from_jobtype()
        logger.debug(f"Route for settings {self}: {route_string}")
        return route_string

    @route_string.setter
    def route_string(self, value):
        self._route_string = value

    def _get_route_string_from_user_input(self):
        """
        Generate route string from user-provided route specification.

        Processes user-defined route string and adds appropriate
        prefix (#-tag) if not already present.

        Returns:
            str: Formatted route string with proper prefix.
        """
        route_string = self.route_to_be_written
        if not route_string.startswith("#"):
            route_string = (
                f"#{self.dieze_tag} {route_string}"
                if self.dieze_tag is not None
                else f"# {route_string}"
            )
        return route_string

    def get_light_elements(self, molecule):
        """
        Identify light elements in molecule for GenECP calculations.

        Determines which elements in the molecular structure are
        considered light elements (not in heavy_elements list) for
        mixed basis set calculations with pseudopotentials.

        Args:
            molecule: Molecule object containing atomic information.

        Returns:
            list or None: Sorted list of light element symbols,
                or None if no heavy elements are specified.
        """
        if self.heavy_elements is None:
            return None

        unique_atoms = set(molecule.chemical_symbols)
        light_elements_set = unique_atoms - set(self.heavy_elements)
        light_elements_list = list(light_elements_set)

        sorted_light_elements_list = pt.sorted_periodic_table_list(
            light_elements_list
        )
        logger.info(
            f"Light elements in structure: {sorted_light_elements_list}"
        )
        return sorted_light_elements_list

    def _get_route_string_from_jobtype(self):
        """
        Generate route string based on job type and settings.

        Constructs the Gaussian route string by analyzing the job type
        and incorporating all relevant calculation parameters including
        optimization options, frequency calculations, and solvation.

        Returns:
            str: Complete route string for the specified job type.

        Raises:
            ValueError: If required parameters are missing or incompatible
                options are specified.
        """
        route_string = ""

        # Get dieze tag with job route and freq string
        dieze_tag = self._get_dieze_tag()
        route_string += dieze_tag

        level_of_theory_string = self._get_level_of_theory_string()
        route_string += level_of_theory_string

        return route_string

    def _get_dieze_tag(self):
        """Get dieze tag from job type."""
        route_string = ""
        # Add #-tag prefix for calculation level specification
        if self.dieze_tag is not None:
            route_string += (
                f"#{self.dieze_tag}"  # e.g. dieze_tag='p' to get '#p'
            )
            logger.debug(f"Added dieze tag: {self.dieze_tag}")
        else:
            route_string += "#"

        # Write optimization keywords with additional options
        # e.g., maxstep, calcall etc
        if self.additional_opt_options_in_route is not None:
            logger.debug(
                f"Adding additional opt options: "
                f"{self.additional_opt_options_in_route}"
            )
            if self.jobtype == "opt":
                route_string += (
                    f" opt=({self.additional_opt_options_in_route})"
                )
            elif self.jobtype == "ts":
                if "calcall" not in self.additional_opt_options_in_route:
                    route_string += f" opt=(ts,calcfc,noeigentest,{self.additional_opt_options_in_route})"
                else:
                    route_string += f" opt=(ts,noeigentest,{self.additional_opt_options_in_route})"
            elif self.jobtype == "modred":
                route_string += f" opt=(modredundant,{self.additional_opt_options_in_route})"
                self.freq = True
            elif self.jobtype == "scan":
                route_string += f" opt=(modredundant,{self.additional_opt_options_in_route})"
                self.freq = False
            elif self.jobtype == "sp":
                route_string += ""
                self.freq = False  # turn off freq calculation for sp job
        elif self.additional_opt_options_in_route is None:
            if self.jobtype == "opt":
                route_string += " opt"
            elif self.jobtype == "ts":
                route_string += " opt=(ts,calcfc,noeigentest)"
            elif self.jobtype == "modred":
                route_string += " opt=modredundant"
                self.freq = True
            elif self.jobtype == "scan":
                route_string += " opt=modredundant"
                self.freq = False
            elif self.jobtype == "sp":
                route_string += ""
                self.freq = False  # turn off freq calculation for sp job

        # Write frequency calculation keywords
        if self.freq and not self.numfreq:
            route_string += " freq"
            logger.debug("Added frequency calculation")
        elif not self.freq and self.numfreq:
            route_string += " freq=numer"
            logger.debug("Added numerical frequency calculation")
        elif self.freq and self.numfreq:
            raise ValueError(
                "Both freq and numfreq cannot be True at the same time!"
            )
        return route_string

    def _get_level_of_theory_string(self):
        """Get level of theory string for route."""
        route_string = ""

        # Determine computational method and add to route string
        if self.semiempirical is not None:
            # Semiempirical methods do not require a basis set
            if self.basis is not None:
                logger.warning(
                    "Basis set provided but not required for semiempirical "
                    "methods."
                )
            route_string += f" {self.semiempirical}"
            logger.debug(f"Added semiempirical method: {self.semiempirical}")

        elif self.ab_initio is not None and self.functional is None:
            # Ab initio method requires a basis set
            if self.basis is None:
                logger.error("Basis set required for ab initio methods")
                raise ValueError(
                    "Error: Basis set is required for ab initio methods."
                )
            route_string += f" {self.ab_initio} {self.basis}"
            logger.debug(
                f"Added ab initio method: {self.ab_initio} with basis: "
                f"{self.basis}"
            )

        elif self.functional is not None and self.ab_initio is None:
            # DFT method requires a basis set
            if self.basis is None:
                logger.error("Basis set required for DFT methods")
                raise ValueError(
                    "Error: Basis set is required for DFT methods."
                )
            route_string += f" {self.functional} {self.basis}"
            logger.debug(
                f"Added DFT functional: {self.functional} with basis: "
                f"{self.basis}"
            )

        elif self.ab_initio is not None and self.functional is not None:
            logger.error(
                "Both ab initio and DFT functional provided - invalid"
            )
            raise ValueError(
                "Error: Both ab initio and DFT functional provided.\n"
                "Specify only one."
            )

        else:
            raise ValueError("Error: No computational method provided.")

        # Write forces calculation keyword
        if self.forces:
            route_string += " force"
            logger.debug("Added force calculation")

        # Handle solvation model configuration
        if self.custom_solvent is not None:
            logger.debug("Using custom solvent parameters")
            if self.solvent_model is None and self.solvent_id is None:
                route_string += " scrf=(pcm,read"  # using pcm model as default
                logger.debug("Using default PCM model for custom solvent")
            else:
                # Set default values if any of solvent_model or solvent_id
                # are None
                solvent_model = self.solvent_model or "pcm"
                solvent_id = self.solvent_id or "generic,read"
                route_string += f" scrf=({solvent_model},solvent={solvent_id}"
                logger.debug(
                    f"Using custom solvent with model: {solvent_model}, "
                    f"ID: {solvent_id}"
                )
        elif (
            self.solvent_model is not None and self.solvent_id is not None
        ):  # solvation is turned on
            route_string += (
                f" scrf=({self.solvent_model},solvent={self.solvent_id}"
            )
            logger.debug(
                f"Added solvation: {self.solvent_model} with solvent "
                f"{self.solvent_id}"
            )
            if self.additional_solvent_options is not None:
                route_string += f",{self.additional_solvent_options})"
                logger.debug(
                    f"Added additional solvent options: "
                    f"{self.additional_solvent_options}"
                )
            else:
                route_string += ")"
        elif (self.solvent_model is not None and self.solvent_id is None) or (
            self.solvent_model is None and self.solvent_id is not None
        ):  # if one is provided but the other not
            logger.error(
                f"Incomplete solvation specification - solvent_model: "
                f"{self.solvent_model}, solvent_id: {self.solvent_id}"
            )
            raise ValueError(
                f"Both solvent model (via -sm or --solvent-model)"
                f"and solvent ID (via -si or --solvent-id) "
                f"need to be specified.\n"
                f"Currently, solvent model is {self.solvent_model} and "
                f"solvent id is {self.solvent_id}!"
            )

        # Write additional parameters for route
        if self.additional_route_parameters is not None:
            route_string += f" {self.additional_route_parameters}"
            logger.debug(
                f"Added additional route parameters: "
                f"{self.additional_route_parameters}"
            )

        # Write job type specific route keywords
        if self.jobtype == "nci":
            route_string += " output=wfn"  # output wavefunction file for NCI
            logger.debug("Added NCI-specific output=wfn keyword")
        elif self.jobtype == "wbi":
            route_string += " pop=nboread"  # write bond order matrix
            logger.debug("Added WBI-specific pop=nboread keyword")
        return route_string

    @property
    def _genecp_elements_specified(self):
        """
        Check if GenECP elements are properly specified.

        Verifies that both heavy elements and their basis set are
        defined for GenECP pseudopotential calculations.

        Returns:
            bool: True if GenECP elements are properly specified.
        """
        return (
            self.heavy_elements is not None
            and self.heavy_elements_basis is not None
        )

    @property
    def _genecp_file_specified(self):
        """
        Check if GenECP file is properly specified and exists.

        Verifies that a GenECP file path is provided and the file
        exists on the filesystem for pseudopotential definitions.

        Returns:
            bool: True if GenECP file is specified and exists.
        """
        return self.gen_genecp_file is not None and os.path.exists(
            self.gen_genecp_file
        )

    def get_genecp_section(self, molecule):
        """
        Generate GenECP section for mixed basis set calculations.

        Creates the GenECP pseudopotential section either from
        explicitly specified elements and basis sets or from a
        pre-defined GenECP file.

        Args:
            molecule: Molecule object for element analysis.

        Returns:
            GenGenECPSection: GenECP section for Gaussian input.

        Raises:
            ValueError: If neither GenECP elements nor file are specified.
        """
        if self._genecp_elements_specified:
            logger.info(
                f"GENECP elements specified:\n"
                f"Heavy elements: {self.heavy_elements}\n"
                f"Heavy elements basis: {self.heavy_elements_basis}\n"
                f"Light elements basis: {self.light_elements_basis}\n"
            )
            # Method 1 for getting genecp
            # Need to supply self.heavy_elements, self.heavy_elements_basis
            # and self.light_elements_basis
            heavy_elements_in_structure = self.prune_heavy_elements(molecule)

            genecp_section = GenGenECPSection.from_bse_api(
                light_elements=self.get_light_elements(molecule),
                light_elements_basis=self.light_elements_basis,
                heavy_elements=heavy_elements_in_structure,
                heavy_elements_basis=self.heavy_elements_basis,
            )

        elif self._genecp_file_specified:
            logger.info(f"GENECP file specified: {self.gen_genecp_file}")
            # Method 2 for getting genecp:
            # Supplied path to genecp file
            genecp_section = GenGenECPSection.from_genecp_path(
                genecp_path=self.gen_genecp_file
            )
        else:
            raise ValueError("Could not get GenECPSection")
        return genecp_section

    def prune_heavy_elements(self, molecule):
        """
        Filter heavy elements list to those present in the molecule.

        Removes heavy elements from the settings list that are not
        actually present in the current molecular structure, ensuring
        only relevant pseudopotentials are included.

        Args:
            molecule: Molecule object containing atomic information.

        Returns:
            list or None: Heavy elements present in the molecule,
                or None if no heavy elements are specified.
        """
        # heavy atoms list supplied from settings contains all heavy atoms
        # needed for heavy_atom_basis but in each structure, some heavy atoms
        # supplied from settings may not appear in the structure
        if self.heavy_elements is None:
            return None
        return list(
            set(molecule.chemical_symbols).intersection(self.heavy_elements)
        )

    def determine_basis_keyword(self, molecule):
        """
        Determine the appropriate basis keyword (gen or genecp) for the molecule.

        Based on the heavy elements present in the molecule, determines whether
        to use 'gen' (for elements that don't require ECPs) or 'genecp' (for
        elements that require ECPs). Elements with atomic number > 36 require
        ECPs and thus need 'genecp', while elements <= 36 only need 'gen'.

        If no heavy elements are present in the molecule, returns the light
        elements basis keyword (with hyphens removed).

        Args:
            molecule: Molecule object containing atomic information.

        Returns:
            str: 'gen' if all heavy elements have atomic number <= 36,
                 'genecp' if any heavy element has atomic number > 36,
                 light_elements_basis (formatted) if no heavy elements present,
                 or the original basis keyword if not using gen/genecp.
        """
        # Only applies when basis is gen or genecp
        if self.basis not in ["gen", "genecp"]:
            return self.basis

        # Get heavy elements actually present in the molecule
        heavy_elements_in_structure = self.prune_heavy_elements(molecule)

        # If no heavy elements specified or present, use light elements basis
        if (
            heavy_elements_in_structure is None
            or len(heavy_elements_in_structure) == 0
        ):
            # Return light elements basis if available,
            # otherwise return original basis
            if self.light_elements_basis is not None:
                # Remove hyphens for Gaussian compatibility
                # (def2-SVP -> def2svp)
                return self.light_elements_basis.replace("-", "").lower()
            return self.basis

        # Check if any heavy element requires ECP (atomic number > 36)
        for element in heavy_elements_in_structure:
            if pt.requires_ecp(element):
                logger.debug(
                    f"Element {element} requires ECP (Z > 36), using 'genecp'"
                )
                return "genecp"

        # All heavy elements have atomic number <= 36, use 'gen'
        logger.debug(
            f"All heavy elements {heavy_elements_in_structure} have Z <= 36, "
            f"using 'gen'"
        )
        return "gen"

    def _check_solvent(self, solvent_model):
        """
        Validate that the specified solvent model is supported.

        Checks if the provided solvent model is in the list of
        supported Gaussian solvation models.

        Args:
            solvent_model (str): Solvent model to validate.

        Raises:
            ValueError: If solvent model is not supported.
        """
        if solvent_model.lower() not in GAUSSIAN_SOLVATION_MODELS:
            raise ValueError(
                f"The specified solvent model {solvent_model} is not in \n"
                f"the available solvent models: {GAUSSIAN_SOLVATION_MODELS}"
            )


class GaussianIRCJobSettings(GaussianJobSettings):
    """
    Specialized settings for Gaussian IRC (Intrinsic Reaction Coordinate) jobs.

    Extends GaussianJobSettings with IRC-specific parameters for tracking
    reaction pathways from transition states to reactants and products.
    Includes options for predictors, recorrection, and flat IRC calculations.

    Attributes:
        predictor (str): IRC predictor method (e.g., 'LQA', 'Euler').
        recorrect (str): Recorrection strategy (e.g., 'never', 'always').
        direction (str): IRC direction ('forward' or 'reverse').
        recalc_step (int): Interval between energy recalculations.
        maxpoints (int): Maximum number of IRC points to follow.
        maxcycles (int): Max optimization cycles per IRC point.
        stepsize (int): IRC integration step size.
    """

    def __init__(
        self,
        predictor=None,
        recorrect=None,
        recalc_step=6,
        direction=None,
        maxpoints=512,
        maxcycles=128,
        stepsize=20,
        flat_irc=False,
        **kwargs,
    ):
        """
        Initialize IRC-specific Gaussian job settings.

        Configures parameters for IRC calculations including predictor
        methods, recorrection strategies, and integration limits.

        Args:
            predictor (str, optional): Predictor method for IRC integration.
            recorrect (str, optional): Recorrection strategy.
            recalc_step (int): Steps between energy recalculations.
            direction (str, optional): IRC direction ('forward'/'reverse'). If None, run both directions.
            maxpoints (int): Maximum number of IRC points.
            maxcycles (int): Maximum optimization cycles per point.
            stepsize (int): IRC integration step size.
            flat_irc (bool): Enable flat IRC calculations.
            **kwargs: Additional arguments for parent class.
        """
        super().__init__(**kwargs)
        self.predictor = predictor
        self.recorrect = recorrect
        self.recalc_step = recalc_step
        self.direction = direction
        self.maxpoints = maxpoints
        self.maxcycles = maxcycles
        self.stepsize = stepsize
        self.flat_irc = flat_irc
        self.freq = False  # turn off freq calc for IRC jobs
        self.forces = False  # turn off forces calculations
        self.route_to_be_written = None

    def _get_route_string_from_jobtype(self):
        """
        Generate IRC-specific route string for Gaussian calculations.

        Constructs the route string with IRC-specific keywords including
        predictor methods, recorrection settings, and integration
        parameters. Handles both regular and flat IRC calculations.

        Returns:
            str: Complete IRC route string.

        Raises:
            ValueError: If predictor and recorrect are inconsistently
                specified.
        """
        route_string = super()._get_route_string_from_jobtype()

        # Configure flat IRC parameters with appropriate defaults
        if self.flat_irc:
            logger.debug("Configuring flat IRC calculation parameters")
            # default route for flat IRC run if no options are specified
            self.predictor = (
                "LQA" if self.predictor is None else self.predictor
            )
            self.recorrect = (
                "never" if self.recorrect is None else self.recorrect
            )
            self.recalc_step = (
                -5 if self.recalc_step == 6 else self.recalc_step
            )
            logger.debug(
                f"Flat IRC settings - predictor: {self.predictor}, "
                f"recorrect: {self.recorrect}, recalc_step: {self.recalc_step}"
            )

        # Write job type specific route for IRC direction
        if self.jobtype == "ircf":
            self.direction = "forward"
            logger.debug("Set IRC direction to forward")
        elif self.jobtype == "ircr":
            self.direction = "reverse"
            logger.debug("Set IRC direction to reverse")

        if self.predictor is not None and self.recorrect is not None:
            route_string += (
                f" irc({self.predictor},calcfc,recorrect={self.recorrect},"
                f"recalc={self.recalc_step},"
                f"stepsize={self.stepsize},{self.direction},"
                f"maxpoints={self.maxpoints},maxcycle={self.maxcycles})"
            )
            logger.debug(
                f"Added IRC route with predictor {self.predictor} and "
                f"recorrect {self.recorrect}"
            )
        elif self.predictor is None and self.recorrect is None:
            route_string += (
                f" irc(calcfc,recalc={self.recalc_step},{self.direction},"
                f"maxpoints={self.maxpoints},maxcycle={self.maxcycles})"
            )
            logger.debug("Added basic IRC route without predictor/recorrect")
        else:
            logger.error(
                f"Inconsistent IRC predictor/recorrect specification - "
                f"predictor: {self.predictor}, recorrect: {self.recorrect}"
            )
            raise ValueError(
                f"Only one of predictor type and recorrect is specified, "
                f"please check!\n"
                f"Predictor: {self.predictor}; Recorrect: {self.recorrect}"
            )

        if self.additional_route_parameters is not None:
            # Check if the additional parameters are already in the route string
            # to avoid duplication (e.g., scf=qc appearing twice)
            additional_params = self.additional_route_parameters.strip()
            if additional_params not in route_string:
                route_string += f" {additional_params}"
            else:
                logger.debug(
                    f"Additional route parameters '{additional_params}' already present in route string"
                )

        return route_string


class GaussianLinkJobSettings(GaussianJobSettings):
    """
    Specialized settings for Gaussian multi-step link calculations.

    Extends GaussianJobSettings for calculations that require multiple
    linked steps, such as stability analysis followed by optimization.
    Manages checkpoint file usage and guess orbital specifications.

    Attributes:
        link (bool): Whether to use link job functionality.
        link_route (str): Custom route string for the link step.
        stable (str): Stability analysis type ('opt', 'qrhf', etc.).
        guess (str): Initial guess method ('mix', 'read', etc.).
    """

    def __init__(
        self,
        link=True,
        link_route=None,
        stable="opt",
        guess="mix",
        # IRC-specific parameters for IRC link jobs
        predictor=None,
        recorrect=None,
        recalc_step=6,
        direction=None,
        maxpoints=512,
        maxcycles=128,
        stepsize=20,
        flat_irc=False,
        **kwargs,
    ):
        """
        Initialize link job specific Gaussian settings.

        Configures parameters for multi-step calculations that use
        checkpoint files to pass information between calculation steps.

        Args:
            link (bool): Enable link job functionality.
            link_route (str, optional): Custom route for link step.
            stable (str): Stability analysis method.
            guess (str): Initial orbital guess method.
            predictor (str, optional): IRC predictor method for IRC link jobs.
            recorrect (str, optional): IRC recorrection strategy for IRC link jobs.
            recalc_step (int): Steps between energy recalculations for IRC.
            direction (str, optional): IRC direction ('forward'/'reverse').
            maxpoints (int): Maximum number of IRC points.
            maxcycles (int): Maximum optimization cycles per IRC point.
            stepsize (int): IRC integration step size.
            flat_irc (bool): Enable flat IRC calculations.
            **kwargs: Additional arguments for parent class.
        """
        super().__init__(**kwargs)
        self.link = link
        self.link_route = link_route
        self.stable = stable
        self.guess = guess

        # IRC-specific parameters (will be None for non-IRC link jobs)
        self.predictor = predictor
        self.recorrect = recorrect
        self.recalc_step = recalc_step
        self.direction = direction
        self.maxpoints = maxpoints
        self.maxcycles = maxcycles
        self.stepsize = stepsize
        self.flat_irc = flat_irc

    @property
    def link_route_string(self):
        """
        Generate the route string for the link calculation step.

        Creates the route string for the second step in a link job,
        ensuring proper geometry and orbital guess specifications
        for continuation from the previous step.

        Returns:
            str: Route string for the link calculation step.
        """
        if self.link_route is not None:
            link_route_string = self.link_route
            if self.functional not in self.link_route:
                link_route_string += f" {self.functional}"
            if self.basis not in self.link_route:
                link_route_string += f" {self.basis}"
            if "geom=check" not in self.link_route:
                link_route_string += " geom=check"
            if "guess=read" not in self.link_route:
                link_route_string += " guess=read"
            logger.debug(
                f"Link route for settings {self}: {link_route_string}"
            )
            return link_route_string

        link_route_string = self._get_link_route_string_from_jobtype()
        logger.debug(f"Link route for settings {self}: {link_route_string}")
        return link_route_string

    def _get_route_string_from_jobtype(self):
        """
        Generate route string for the initial stability analysis step.

        Creates the first route string in a link job by removing
        optimization and frequency keywords and adding stability
        analysis and guess method specifications.

        Returns:
            str: Route string for stability analysis step.
        """
        route_string = super()._get_route_string_from_jobtype()
        # Remove opt keywords
        route_string_final = re.sub(
            gaussian_opt_keywords_pattern,
            " ",
            route_string,
            flags=re.IGNORECASE,
        )
        # Remove freq keywords
        route_string_final = re.sub(
            gaussian_freq_keywords_pattern,
            " ",
            route_string_final,
            flags=re.IGNORECASE,
        )
        # Clean up multiple spaces
        route_string_final = re.sub(
            multiple_spaces_pattern, " ", route_string_final
        ).strip()

        if self.stable:
            logger.debug(f"Stable: {self.stable}")
            route_string_final += f" stable={self.stable}"
        if self.guess:
            logger.debug(f"Guess: {self.guess}")
            route_string_final += f" guess={self.guess}"

        return route_string_final

    def _get_link_route_string_from_jobtype(self):
        """
        Generate route string for the optimization/IRC step in link jobs.

        Creates the second route string that continues from the
        stability analysis by adding checkpoint geometry and orbital
        reading specifications. Handles IRC jobs by using GaussianIRCJobSettings.

        Returns:
            str: Route string for the optimization/IRC step.
        """
        # Special handling for IRC jobs - use GaussianIRCJobSettings
        if self.jobtype in ["ircf", "ircr"]:
            # Create a temporary GaussianIRCJobSettings instance with current settings
            irc_settings = GaussianIRCJobSettings(**self.__dict__)
            # Get the IRC route string from the specialized class
            route_string = irc_settings._get_route_string_from_jobtype()
        else:
            # For non-IRC jobs, use the existing logic
            route_string = super()._get_route_string_from_jobtype()

        if "geom=check" not in route_string:
            route_string += " geom=check"
        if "guess=read" not in route_string:
            route_string += " guess=read"
        return route_string


class GaussianTDDFTJobSettings(GaussianJobSettings):
    """
    Specialized settings for Gaussian TD-DFT excited state calculations.

    Extends GaussianJobSettings with parameters specific to time-dependent
    density functional theory calculations for excited states, including
    state selection, solvation effects, and transition analysis.

    Attributes:
        states (str): Type of excited states ('singlets', 'triplets').
        root (int): Specific excited state root to analyze.
        nstates (int): Number of excited states to calculate.
        eqsolv (str): Equilibrium solvation treatment.
    """

    def __init__(
        self, states="singlets", root=1, nstates=3, eqsolv=None, **kwargs
    ):
        """
        Initialize TD-DFT specific Gaussian job settings.

        Configures parameters for excited state calculations including
        the number of states, state types, and solvation effects.

        Args:
            states (str): Type of excited states to calculate.
            root (int): Specific excited state root for analysis.
            nstates (int): Number of excited states to compute.
            eqsolv (str, optional): Equilibrium solvation option
                ('eqsolv' or 'noneqsolv').
            **kwargs: Additional arguments for parent class.
        """
        super().__init__(**kwargs)
        self.states = states
        self.root = root
        self.nstates = nstates
        self.eqsolv = eqsolv

    def _get_route_string_from_jobtype(self):
        """
        Generate TD-DFT specific route string for excited state calculations.

        Constructs the route string with TD-DFT keywords including
        state types, number of states, root selection, and solvation
        equilibrium options.

        Returns:
            str: Complete TD-DFT route string.

        Raises:
            AssertionError: If eqsolv option is not valid.
        """
        route_string = super()._get_route_string_from_jobtype()

        if self.eqsolv is None:
            eqsolv = ""
        else:
            eqsolv_options = ["eqsolv", "noneqsolv"]
            assert (
                self.eqsolv.lower() in eqsolv_options
            ), f"Possible equilibrium solvation options are: {eqsolv_options}!"
            eqsolv = f",{self.eqsolv}"

        route_string += f" TD({self.states},nstates={self.nstates},root={self.root}{eqsolv})"

        return route_string


class GaussianQMMMJobSettings(GaussianJobSettings):
    """
    Configuration settings for Gaussian QM/MM calculations using ONIOM methodology.

    This class manages all parameters needed to set up multi-layer ONIOM calculations
    in Gaussian, which partition molecular systems into different regions treated with
    varying levels of theory. The ONIOM approach enables accurate quantum mechanical
    treatment of chemically active regions while efficiently handling large molecular
    environments with molecular mechanics.

    The class supports 2-layer and 3-layer ONIOM calculations:

    **2-Layer ONIOM**: High(QM):Low(MM)
        - High level: Quantum mechanics (DFT, ab initio, etc.)
        - Low level: Molecular mechanics force fields

    **3-Layer ONIOM**: High(QM):Medium(QM):Low(MM)
        - High level: High-accuracy quantum mechanics
        - Medium level: Lower-cost quantum mechanics
        - Low level: Molecular mechanics force fields

    Key Features:
    - Flexible layer definition with atom selection
    - Support for mixed QM/MM and QM/QM/MM schemes
    - Automatic link atom handling for covalent boundaries
    - Customizable scale factors for link atom placement
    - Integration with popular force fields (AMBER, UFF, etc.)
    - Multiple charge/multiplicity specifications per layer

    Attributes:
        jobtype (str, optional): Type of ONIOM calculation ('sp', 'opt', 'freq', 'ts', 'irc',
            'modred', 'scan'). When using CLI commands, this is automatically inferred from
            the parent command (e.g., 'chemsmart sub gaussian opt qmmm' sets jobtype='opt').

        Level-specific theory parameters:
            high_level_functional (str): DFT functional for high layer (e.g., 'B3LYP', 'M06-2X')
            high_level_basis (str): Basis set for high layer (e.g., '6-31G*', 'def2-TZVP')
            high_level_force_field (str): Force field for high layer (if MM)
            medium_level_functional (str): DFT functional for medium layer
            medium_level_basis (str): Basis set for medium layer
            medium_level_force_field (str): Force field for medium layer
            low_level_functional (str): DFT functional for low layer (if QM)
            low_level_basis (str): Basis set for low layer (if QM)
            low_level_force_field (str): Force field for low layer (usually MM)

        Charge and multiplicity specifications:
            charge_total/mult_total (int): Full system properties (legacy: charge_total/real_multiplicity)
            charge_intermediate/mult_intermediate (int): Intermediate system properties (legacy: int_charge/int_multiplicity)
            charge_high/mult_high (int): Model/high-layer properties (legacy: model_charge/model_multiplicity)

        Atom partitioning:
            high_level_atoms (list/str): Atoms in high layer (1-indexed)
            medium_level_atoms (list/str): Atoms in medium layer (1-indexed)
            low_level_atoms (list/str): Atoms in low layer (1-indexed)
            bonded_atoms (list): Covalent bonds crossing layer boundaries
            scale_factors (dict): Custom scale factors for link atom placement

    Examples:
        2-layer enzyme active site calculation:
        >>> settings = GaussianQMMMJobSettings(
        ...     jobtype='opt',
        ...     high_level_functional='B3LYP',
        ...     high_level_basis='6-31G*',
        ...     low_level_force_field='AMBER=HardFirst',
        ...     charge_total=0,
        ...     real_multiplicity=1,
        ...     high_level_atoms=[1, 2, 3, 4, 5],  # Active site residues
        ...     bonded_atoms=[(5, 6)],  # QM-MM boundary bond
        ... )

        3-layer organometallic catalyst:
        >>> settings = GaussianQMMMJobSettings(
        ...     jobtype='freq',
        ...     high_level_functional='M06-2X',
        ...     high_level_basis='def2-TZVP',
        ...     medium_level_functional='B3LYP',
        ...     medium_level_basis='6-31G*',
        ...     low_level_force_field='UFF',
        ...     charge_total=-1,
        ...     real_multiplicity=2,
        ...     high_level_atoms='1-10',     # Metal center and ligands
        ...     medium_level_atoms='11-50',  # Extended coordination sphere
        ...     bonded_atoms=[(10,11), (50,51)],
        ...     scale_factors={(10,11): [0.709, 0.709, 0.709]}
        ... )

    Note:
        The real/intermediate/model charge and multiplicity specifications follow
        ONIOM conventions where:
        - Real system: Complete molecular system
        - Intermediate system: High+Medium layers (3-layer only)
        - Model system: High layer only

        Scale factors control link atom placement and default to covalent radii
        ratios if not specified. The format is {(atom1,atom2): [low,medium,high]}.

    See Also:
        - GaussianQMMMJob: Job execution class for QM/MM calculations
        - GaussianJobSettings: Base class for Gaussian job configuration
        - Molecule: Molecular structure with QM/MM partitioning information
    """

    def __init__(
        self,
        jobtype=None,
        parent_jobtype=None,
        high_level_functional=None,
        high_level_basis=None,
        high_level_force_field=None,
        medium_level_functional=None,
        medium_level_basis=None,
        medium_level_force_field=None,
        low_level_functional=None,
        low_level_basis=None,
        low_level_force_field=None,
        charge_total=None,
        mult_total=None,
        charge_intermediate=None,
        mult_intermediate=None,
        charge_high=None,
        mult_high=None,
        real_charge=None,
        real_multiplicity=None,
        int_charge=None,
        int_multiplicity=None,
        model_charge=None,
        model_multiplicity=None,
        high_level_atoms=None,
        medium_level_atoms=None,
        low_level_atoms=None,
        bonded_atoms=None,
        scale_factors=None,
        **kwargs,
    ):
        """
        Initialize Gaussian QM/MM job settings for ONIOM calculations.

        Args:
            jobtype (str, optional): Type of ONIOM calculation to perform.
                Options: 'sp' (single-point), 'opt' (optimization), 'freq' (frequency),
                'ts' (transition state), 'irc' (intrinsic reaction coordinate),
                'modred' (redundant coordinates), 'scan' (coordinate scan)

                Note: When using the CLI with `chemsmart sub gaussian <jobtype> qmmm`,
                the jobtype is automatically inferred from the parent command and does
                not need to be specified manually.

            Theory level parameters:
                high_level_functional (str): DFT functional for high QM layer
                    (e.g., 'B3LYP', 'M06-2X', 'wB97X-D', 'PBE0')
                high_level_basis (str): Basis set for high QM layer
                    (e.g., '6-31G*', '6-311++G(d,p)', 'def2-TZVP', 'cc-pVTZ')
                high_level_force_field (str): Force field for high MM layer (rare)
                medium_level_functional (str): DFT functional for medium QM layer
                medium_level_basis (str): Basis set for medium QM layer
                medium_level_force_field (str): Force field for medium MM layer
                low_level_functional (str): DFT functional for low QM layer (uncommon)
                low_level_basis (str): Basis set for low QM layer (uncommon)
                low_level_force_field (str): Force field for low MM layer
                    (e.g., 'AMBER=HardFirst', 'UFF', 'DREIDING', 'MM3')

            Charge and multiplicity:
                charge_total (int): Total charge of complete molecular system (legacy: charge_total)
                mult_total (int): Spin multiplicity of complete system (legacy: real_multiplicity)
                charge_intermediate (int): Charge of high+medium layers (legacy: int_charge)
                mult_intermediate (int): Multiplicity of high+medium layers (legacy: int_multiplicity)
                charge_high (int): Charge of high layer only (legacy: model_charge)
                mult_high (int): Multiplicity of high layer only (legacy: model_multiplicity)

            Atom partitioning:
                high_level_atoms (list/str): Atoms in high layer. Can be:
                    - List of integers: [1, 2, 3, 5, 8]
                    - Range string: "1-10,15,20-25"
                    - List of ranges: ["1-10", "15", "20-25"]
                medium_level_atoms (list/str): Atoms in medium layer (3-layer only)
                low_level_atoms (list/str): Atoms in low layer (usually auto-assigned)
                bonded_atoms (list): Covalent bonds crossing layer boundaries.
                    Format: [(atom1, atom2), (atom3, atom4)] where atoms are 1-indexed
                scale_factors (dict): Custom link atom scale factors. Format:
                    {(atom1, atom2): [scale_low, scale_medium, scale_high]}
                    Values typically range from 0.5-1.0, common value is 0.709

            **kwargs: Additional keyword arguments passed to parent GaussianJobSettings

        Note:
            - All atom indices are 1-based following Gaussian conventions
            - If only 2 layers specified, use high_level_* and low_level_* parameters
            - Force fields must be available in your Gaussian installation
            - Link atoms are automatically placed for bonded_atoms specifications
            - The parent class 'functional' and 'basis' attributes are set to high_level values
        """
        super().__init__(**kwargs)
        self.jobtype = jobtype
        self.parent_jobtype = parent_jobtype
        self.high_level_functional = high_level_functional
        self.high_level_basis = high_level_basis
        self.high_level_force_field = high_level_force_field
        self.medium_level_functional = medium_level_functional
        self.medium_level_basis = medium_level_basis
        self.medium_level_force_field = medium_level_force_field
        self.low_level_functional = low_level_functional
        self.low_level_basis = low_level_basis
        self.low_level_force_field = low_level_force_field
        # Canonical charge/multiplicity names (aligned with ORCA) with legacy fallbacks
        self.charge_total = charge_total
        if self.charge_total is None:
            self.charge_total = real_charge
        self.mult_total = mult_total
        if self.mult_total is None:
            self.mult_total = real_multiplicity

        self.charge_intermediate = charge_intermediate
        if self.charge_intermediate is None:
            self.charge_intermediate = int_charge
        self.mult_intermediate = mult_intermediate
        if self.mult_intermediate is None:
            self.mult_intermediate = int_multiplicity

        self.charge_high = charge_high
        if self.charge_high is None:
            self.charge_high = model_charge
        self.mult_high = mult_high
        if self.mult_high is None:
            self.mult_high = model_multiplicity

        # Maintain legacy attribute names for backward compatibility
        self.real_charge = self.charge_total
        self.real_multiplicity = self.mult_total
        self.int_charge = self.charge_intermediate
        self.int_multiplicity = self.mult_intermediate
        self.model_charge = self.charge_high
        self.model_multiplicity = self.mult_high
        self.high_level_atoms = high_level_atoms
        self.medium_level_atoms = medium_level_atoms
        self.low_level_atoms = low_level_atoms
        self.bonded_atoms = bonded_atoms
        self.scale_factors = scale_factors

        self.title = "Gaussian QM/MM job"

        if self.charge_total is not None and self.mult_total is not None:
            # the charge and multiplicity of the real system equal to
            # that of the low_level_charge and low_level_multiplicity
            self.charge = self.charge_total
            self.multiplicity = self.mult_total

    @property
    def charge_and_multiplicity_string(self):
        """Obtain charge and multiplicity string."""
        return self._get_charge_and_multiplicity()

    def validate_and_assign_level(
        self, functional, basis, force_field, level_name
    ):
        """Validates functional and basis set for a given level
        and returns formatted theory string.
        Return level of theory if both functional and basis are specified,
        or force field if both are not specified.
        """

        if functional and basis and force_field:
            raise ValueError(
                f"For {level_name} level of theory, one should specify only functional/basis or force field!"
            )

        if force_field:
            assert functional is None and basis is None, (
                f"Force field is given for {level_name} level of theory, "
                f"thus no functional and basis should be given!"
            )
            level_of_theory = force_field
        else:
            # if force field is not given, then functional and basis can be given,
            # so that level of theory takes functional and basis set
            if functional and basis:
                level_of_theory = f"{functional}/{basis}"
            else:
                # but functional and basis set can also not be given, in which case,
                # all 3 are None and overall level of theory for that layer is None.
                level_of_theory = None

        logger.debug(
            f"Obtained level of theory {level_of_theory} for {level_name} level."
        )

        return level_of_theory

    def _get_route_string_from_jobtype(self):
        """Generate QM/MM route string with job type, freq, and ONIOM specification."""
        route_string = "#"
        if self.dieze_tag:
            route_string += self.dieze_tag

        # Add job type
        parent_jobtype = (
            self.parent_jobtype.lower() if self.parent_jobtype else None
        )
        if parent_jobtype == "opt":
            route_string += " opt"
        elif parent_jobtype == "scan" or parent_jobtype == "modred":
            route_string += " opt=modredundant"
        elif parent_jobtype == "ts":
            route_string += " opt=(ts,calcfc,noeigentest)"
        elif parent_jobtype == "irc":
            route_string += " irc"
        # sp doesn't add any keyword

        # Add freq if enabled (and not already a freq job)
        if (self.freq or self.numfreq) and parent_jobtype != "freq":
            route_string += " freq"

        # Add ONIOM level of theory string
        oniom_string = self._get_oniom_string()
        route_string += oniom_string

        # Add solvation if specified
        if self.solvent_model is not None and self.solvent_id is not None:
            route_string += (
                f" scrf=({self.solvent_model},solvent={self.solvent_id})"
            )

        return route_string

    def _get_oniom_string(self):
        """Get ONIOM level of theory string."""
        high_level_of_theory = self.validate_and_assign_level(
            self.high_level_functional,
            self.high_level_basis,
            self.high_level_force_field,
            level_name="high",
        )

        medium_level_of_theory = self.validate_and_assign_level(
            self.medium_level_functional,
            self.medium_level_basis,
            self.medium_level_force_field,
            level_name="medium",
        )

        low_level_of_theory = self.validate_and_assign_level(
            self.low_level_functional,
            self.low_level_basis,
            self.low_level_force_field,
            level_name="low",
        )

        # Build ONIOM string with proper parentheses handling
        levels = []
        if high_level_of_theory is not None:
            levels.append(high_level_of_theory)
        if medium_level_of_theory is not None:
            levels.append(medium_level_of_theory)
        if low_level_of_theory is not None:
            levels.append(low_level_of_theory)

        if levels:
            oniom_string = f" oniom({':'.join(levels)})"
        else:
            oniom_string = " oniom"

        return oniom_string

    def get_qmmm_level_of_theory_string(self):
        """Get ONIOM level of theory for route string.

        Deprecated: Use _get_route_string_from_jobtype instead.

        Note: jobtype is now inferred from the parent command when using CLI,
        so it should always be set. If not set, defaults to basic ONIOM route.
        """
        if self.jobtype is None:
            logger.warning(
                "Job type not specified for ONIOM job. "
                "Using basic route string without job-specific keywords."
            )
        return self._get_route_string_from_jobtype()

    def _get_charge_and_multiplicity(self):
        """Obtain charge and multiplicity string.
        For two-layer ONIOM jobs, the format for this input line is:

        chrg_real-low spin_real-low [chrg_model-high spin_model-high
                                    [chrg_model-low spin_model-low [chrg_real-high spin_real-high]]]

        Fourth pair applies only to ONIOM=SValue calculations.
        When only a single value pair is specified, all levels will use those values.
        If two pairs of values are included, then third pair defaults to same values as second pair.
        If final pair is omitted for an S-value job, it defaults to values for the real system at low level.
        For such two-layer ONIOM jobs, users are required to specify the charge and multiplicity of high-level
         layer and low-level layer, instead high and medium level.

        For 3-layers ONIOM, the format is:
        cRealL sRealL [cIntM sIntM [cIntL sIntL [cModH sModH [cModM sModM [cModL sModL]]]]]
        Real, Int=Intermediate system, and Mod=Model system, and second character
        is one of: H, M and L for the High, Medium and Low levels).
        """
        assert (
            self.charge_total is not None and self.mult_total is not None
        ), "Charge and multiplicity for the real system must be specified!"
        real_low_charge = self.charge_total
        real_low_multiplicity = self.mult_total
        int_med_charge = self.charge_intermediate
        int_med_multiplicity = self.mult_intermediate
        int_low_charge = self.charge_intermediate
        int_low_multiplicity = self.mult_intermediate
        model_high_charge = self.charge_high
        model_high_multiplicity = self.mult_high
        model_med_charge = self.charge_high
        model_med_multiplicity = self.mult_high
        model_low_charge = self.charge_high
        model_low_multiplicity = self.mult_high

        # two-layer ONIOM model
        if (
            self.validate_and_assign_level(
                self.medium_level_functional,
                self.medium_level_basis,
                self.medium_level_force_field,
                level_name="medium",
            )
            is None
            or self.validate_and_assign_level(
                self.low_level_functional,
                self.low_level_basis,
                self.low_level_force_field,
                level_name="low",
            )
            is None
        ):
            charge_and_multiplicity_list = [
                real_low_charge,
                real_low_multiplicity,
                model_high_charge,
                model_high_multiplicity,
                model_low_charge,
                model_low_multiplicity,
            ]
            if all(var is None for var in charge_and_multiplicity_list[2:]):
                for i in range(2, len(charge_and_multiplicity_list), 2):
                    charge_and_multiplicity_list[i] = real_low_charge
                    charge_and_multiplicity_list[i + 1] = real_low_multiplicity
            elif all(var is None for var in charge_and_multiplicity_list[4:]):
                for i in range(4, len(charge_and_multiplicity_list), 2):
                    charge_and_multiplicity_list[i] = model_high_charge
                    charge_and_multiplicity_list[i + 1] = (
                        model_high_multiplicity
                    )
            elif all(var is not None for var in charge_and_multiplicity_list):
                pass
            else:
                raise ValueError(
                    "The charge and multiplicity of lower level-of-theory cannot override the higher ones!"
                )
            updated_list = []
            for charge_and_multiplicity in charge_and_multiplicity_list:
                updated_list.append(str(charge_and_multiplicity))
            charge_and_multiplicity = " ".join(updated_list)
        else:
            # three-layer ONIOM model
            charge_and_multiplicity_list = [
                real_low_charge,
                real_low_multiplicity,
                int_med_charge,
                int_med_multiplicity,
                int_low_charge,
                int_low_multiplicity,
                model_high_charge,
                model_high_multiplicity,
                model_med_charge,
                model_med_multiplicity,
                model_low_charge,
                model_low_multiplicity,
            ]
            # Defaults for missing charge / spin multiplicity pairs are taken from the next highest
            # calculation level and / or system size.
            if all(var is None for var in charge_and_multiplicity_list[2:]):
                # only charge and multiplicity of real system is specified,
                # so the charge and multiplicity of other systems will be the same as the real system
                for i in range(2, len(charge_and_multiplicity_list), 2):
                    charge_and_multiplicity_list[i] = real_low_charge
                    charge_and_multiplicity_list[i + 1] = real_low_multiplicity
            elif all(var is None for var in charge_and_multiplicity_list[4:]):
                # only charge and multiplicity of real system and that of intermediate layer,
                # medium level-of-theory are specified, the charge and multiplicity of other
                # systems will be the same as the intermediate layer
                for i in range(4, len(charge_and_multiplicity_list), 2):
                    charge_and_multiplicity_list[i] = int_med_charge
                    charge_and_multiplicity_list[i + 1] = int_med_multiplicity
            elif all(var is None for var in charge_and_multiplicity_list[6:]):
                # only charge and multiplicity of real system, intermediate layer, medium level-of-theory
                # and intermediate layer, medium level-of-theory are specified, the charge and multiplicity of other
                # systems will be the same as intermediate layer, medium level-of-theory,...
                for i in range(6, len(charge_and_multiplicity_list), 2):
                    charge_and_multiplicity_list[i] = int_med_charge
                    charge_and_multiplicity_list[i + 1] = int_med_multiplicity
            elif all(var is None for var in charge_and_multiplicity_list[8:]):
                # the rest systems will follow the model system, high level-of-theory
                for i in range(8, len(charge_and_multiplicity_list), 2):
                    charge_and_multiplicity_list[i] = model_high_charge
                    charge_and_multiplicity_list[i + 1] = (
                        model_high_multiplicity
                    )
            elif all(var is None for var in charge_and_multiplicity_list[10:]):
                charge_and_multiplicity_list[-2] = model_med_charge
                charge_and_multiplicity_list[-1] = model_med_multiplicity
            elif all(var is not None for var in charge_and_multiplicity_list):
                pass
            else:
                raise ValueError(
                    "The charge and multiplicity of lower level-of-theory cannot override the higher ones!"
                )
            updated_list = []
            for charge_and_multiplicity in charge_and_multiplicity_list:
                updated_list.append(str(charge_and_multiplicity))
            charge_and_multiplicity = " ".join(updated_list)
        return charge_and_multiplicity

    def __eq__(self, other):
        """
        Compare two GaussianQMMMJobSettings objects for equality.

        Compares all attributes between two QMMM settings objects, including the
        QMMM-specific attributes that are not present in the parent class.

        Args:
            other (GaussianQMMMJobSettings): Settings object to compare with.

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

        is_equal = self_dict == other_dict
        if not is_equal:
            import dictdiffer

            logger.info("Gaussian QMMM job settings are not equal.")
            for diff in list(dictdiffer.diff(self_dict, other_dict)):
                logger.info(f"Difference: {diff}")

        return self_dict == other_dict
