import copy
import logging
import os

from chemsmart.io.gaussian import GAUSSIAN_SOLVATION_MODELS
from chemsmart.io.gaussian.gengenecp import GenGenECPSection
from chemsmart.jobs.settings import MolecularJobSettings
from chemsmart.utils.periodictable import PeriodicTable

pt = PeriodicTable()


logger = logging.getLogger(__name__)


class XTBJobSettings(MolecularJobSettings):
    def __init__(
        self,
        xtb_version="gfn2",  # use gfn2 by default
        optimization_level=None,
        charge=None,
        uhf=None,
        job_type=None,
        title=None,
        freq=False,
        solvent_model=None,
        solvent_id=None,
        **kwargs,
    ):
        """Initialize XTB job settings.
        uhf(int): Number of unpaired electrons.
        """
        super().__init__(
            charge=charge,
            freq=freq,
            job_type=job_type,
            title=title,
            solvent_model=solvent_model,
            solvent_id=solvent_id,
            **kwargs,
        )
        self.xtb_version = xtb_version
        self.optimization_level = optimization_level
        self.uhf = uhf

    def copy(self):
        return copy.deepcopy(self)

    def __getitem__(self, key):
        return self.__dict__[key]

    def __eq__(self, other):
        """Two settings objects are equal if all their attributes are equal."""
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
        """Return Gaussian settings object from a given gaussian.com file.

        Args:
            filename (str): file path of the .com file string to be supplied.
        """
        from chemsmart.io.gaussian.input import Gaussian16Input

        com_path = os.path.abspath(filename)
        gaussian_settings_from_comfile = Gaussian16Input(
            filename=com_path
        ).read_settings()
        return gaussian_settings_from_comfile

    @classmethod
    def from_inpfile(cls, filename):
        """Return Gaussian settings object from a given orca.inp file.

        Args:
            filename (str): file path of the .inp file string to be supplied.
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
        """Return Gaussian settings object from a given gaussian.log file.

        Args:
            filename (str): file path of the .log file to be supplied.
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
        """Return Gaussian job settings from ORCA output file."""
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
        return cls(
            xtb_version="gfn2",
            optimization_level="vtight",
            charge=None,
            uhf=None,
            job_type=None,
            title="xtb job with default settings",
            freq=False,
            solvent_model=None,
            solvent_id=None,
        )

    @classmethod
    def from_filepath(cls, filepath, **kwargs):
        if filepath.endswith((".com", ".gjf")):
            return cls.from_comfile(filepath)
        if filepath.endswith(".inp"):
            return cls.from_inpfile(filepath)
        if filepath.endswith(".log"):
            return cls.from_logfile(filepath)
        if filepath.endswith(".out"):
            return cls.from_outfile(filepath)
        if filepath.endswith(".xyz"):
            return cls.default()
        raise ValueError(f"Could not create {cls} from {filepath}")

    @property
    def route_string(self):
        if self.route_to_be_written is not None:
            route_string = self._get_route_string_from_user_input()
        else:
            route_string = self._get_route_string_from_jobtype()
        logger.debug(f"Route for settings {self}: {route_string}")
        return route_string

    def _get_route_string_from_user_input(self):
        route_string = self.route_to_be_written
        if not route_string.startswith("#"):
            route_string = (
                f"#{self.dieze_tag} {route_string}"
                if self.dieze_tag is not None
                else f"# {route_string}"
            )
        return route_string

    def get_light_elements(self, molecule):
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
        route_string = ""
        if self.dieze_tag is not None:
            route_string += (
                f"#{self.dieze_tag}"  # e.g. dieze_tag='p' to get '#p'
            )
        else:
            route_string += "#"

        # write opt with additional options e.g., maxstep, calcall etc
        if self.additional_opt_options_in_route is not None:
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

        # write frequency
        if self.freq and not self.numfreq:
            route_string += " freq"
        elif not self.freq and self.numfreq:
            route_string += " freq=numer"

        # write functional and basis
        if self.basis is None:
            raise ValueError("Warning: Basis is missing!")
        if self.ab_initio is not None and self.functional is None:
            method = self.ab_initio
        elif self.ab_initio is None and self.functional is not None:
            method = self.functional
        elif self.ab_initio is not None and self.functional is not None:
            raise ValueError(
                "Warning: Both ab initio and DFT functional are provided!"
            )
        else:
            raise ValueError(
                "Warning: Both ab initio and DFT functional are missing!"
            )

        # write basis set
        route_string += f" {method} {self.basis}"

        # write forces calculation
        if self.forces:
            route_string += " force"

        if self.custom_solvent is not None:
            if self.solvent_model is None and self.solvent_id is None:
                route_string += (
                    " scrf=(pcm,read)"  # using pcm model as default
                )
            else:
                # Set default values if any of solvent_model or solvent_id are None
                solvent_model = self.solvent_model or "pcm"
                solvent_id = self.solvent_id or "generic,read"
                route_string += f" scrf=({solvent_model},solvent={solvent_id})"
        elif (
            self.solvent_model is not None and self.solvent_id is not None
        ):  # solvation is turned on
            route_string += (
                f" scrf=({self.solvent_model},solvent={self.solvent_id})"
            )
        elif (self.solvent_model is not None and self.solvent_id is None) or (
            self.solvent_model is None and self.solvent_id is not None
        ):  # if one is provided but the other not
            raise ValueError(
                f"Both solvent model and solvent ID need to be specified.\n"
                f"Currently, solvent model is {self.solvent_model} and solvent id is {self.solvent_id}!"
            )

        # write additional parameters for route
        if self.additional_route_parameters is not None:
            route_string += f" {self.additional_route_parameters}"

        # write job type specific route
        if self.job_type == "nci":
            route_string += " output=wfn"  # output wavefunction file for NCI
        elif self.job_type == "wbi":
            route_string += " pop=nboread"  # write bond order matrix
        return route_string

    @property
    def _genecp_elements_specified(self):
        return (
            self.heavy_elements is not None
            and self.heavy_elements_basis is not None
        )

    @property
    def _genecp_file_specified(self):
        return self.gen_genecp_file is not None and os.path.exists(
            self.gen_genecp_file
        )

    def get_genecp_section(self, molecule):
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
        # heavy atoms list supplied from settings contains all heavy atoms needed for
        # heavy_atom_basis but in each structure, some heave atoms supplied from settings
        # may not appear in the structure
        if self.heavy_elements is None:
            return None
        return list(
            set(molecule.chemical_symbols).intersection(self.heavy_elements)
        )

    def _check_solvent(self, solvent_model):
        if solvent_model.lower() not in GAUSSIAN_SOLVATION_MODELS:
            raise ValueError(
                f"The specified solvent model {solvent_model} is not in \n"
                f"the available solvent models: {GAUSSIAN_SOLVATION_MODELS}"
            )
