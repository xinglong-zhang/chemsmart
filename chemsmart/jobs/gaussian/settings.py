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
    `multiplicity`, `job_type`, `title`, `freq`, `numfreq`, `solvent_model`,
    `solvent_id`, `custom_solvent`, `forces`, and `input_string`.

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
        job_type=None,
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
            job_type (str, optional): Calculation type (e.g., 'opt', 'freq').
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
            job_type=job_type,
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
            job_type=None,
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
            if self.job_type == "opt":
                route_string += (
                    f" opt=({self.additional_opt_options_in_route})"
                )
            elif self.job_type == "ts":
                if "calcall" not in self.additional_opt_options_in_route:
                    route_string += f" opt=(ts,calcfc,noeigentest,{self.additional_opt_options_in_route})"
                else:
                    route_string += f" opt=(ts,noeigentest,{self.additional_opt_options_in_route})"
            elif self.job_type == "modred":
                route_string += f" opt=(modredundant,{self.additional_opt_options_in_route})"
                self.freq = True
            elif self.job_type == "scan":
                route_string += f" opt=(modredundant,{self.additional_opt_options_in_route})"
                self.freq = False
            elif self.job_type == "sp":
                route_string += ""
                self.freq = False  # turn off freq calculation for sp job
        elif self.additional_opt_options_in_route is None:
            if self.job_type == "opt":
                route_string += " opt"
            elif self.job_type == "ts":
                route_string += " opt=(ts,calcfc,noeigentest)"
            elif self.job_type == "modred":
                route_string += " opt=modredundant"
                self.freq = True
            elif self.job_type == "scan":
                route_string += " opt=modredundant"
                self.freq = False
            elif self.job_type == "sp":
                route_string += ""
                self.freq = False  # turn off freq calculation for sp job

        # Write frequency calculation keywords
        if self.freq and not self.numfreq:
            route_string += " freq"
            logger.debug("Added frequency calculation")
        elif not self.freq and self.numfreq:
            route_string += " freq=numer"
            logger.debug("Added numerical frequency calculation")

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
        elif (self.solvent_model is not None and self.solvent_id is None) or (
            self.solvent_model is None and self.solvent_id is not None
        ):  # if one is provided but the other not
            logger.error(
                f"Incomplete solvation specification - model: "
                f"{self.solvent_model}, ID: {self.solvent_id}"
            )
            raise ValueError(
                f"Both solvent model and solvent ID need to be specified.\n"
                f"Currently, solvent model is {self.solvent_model} and "
                f"solvent id is {self.solvent_id}!"
            )

        # Add additional solvation options if present
        if "scrf" in route_string:
            if self.additional_solvent_options is not None:
                route_string += f",{self.additional_solvent_options})"
                logger.debug(
                    f"Added additional solvent options: "
                    f"{self.additional_solvent_options}"
                )
            else:
                route_string += ")"

        # Write additional parameters for route
        if self.additional_route_parameters is not None:
            route_string += f" {self.additional_route_parameters}"
            logger.debug(
                f"Added additional route parameters: "
                f"{self.additional_route_parameters}"
            )

        # Write job type specific route keywords
        if self.job_type == "nci":
            route_string += " output=wfn"  # output wavefunction file for NCI
            logger.debug("Added NCI-specific output=wfn keyword")
        elif self.job_type == "wbi":
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
            direction (str, optional): IRC direction ('forward'/'reverse').
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
        if self.job_type == "ircf":
            self.direction = "forward"
            logger.debug("Set IRC direction to forward")
        elif self.job_type == "ircr":
            self.direction = "reverse"
            logger.debug("Set IRC direction to reverse")

        # Build IRC route string with predictor and recorrection options
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
        # Remove opt keywords (opt, opt=(...), opt=word)
        pattern = re.compile(
            r"\bopt\s*(=\s*(\([^)]*\)|\w+))?\s*", re.IGNORECASE
        )
        route_string_final = re.sub(pattern, " ", route_string)
        # Remove freq keywords for all job types
        route_string_final = re.sub(
            r"\bfreq\s*(=\s*\w+)?\s*",
            " ",
            route_string_final,
            flags=re.IGNORECASE,
        )
        # Clean up multiple spaces
        route_string_final = re.sub(r"\s+", " ", route_string_final).strip()

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
        if self.job_type in ["irc", "ircf", "ircr"]:
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
