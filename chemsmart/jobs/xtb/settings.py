"""
xTB job settings implementation.

This module contains the XTBJobSettings class for configuring xTB jobs.
It manages parameters such as GFN version, optimization level, and solvent models.
"""

import copy
import logging
import os

from chemsmart.io.xtb import (
    XTB_ALL_METHODS,
    XTB_ALL_OPT_LEVELS,
    XTB_ALL_SOLVENT_IDS,
    XTB_ALL_SOLVENT_MODELS,
)

logger = logging.getLogger(__name__)


class XTBJobSettings:
    """
    Settings for xTB jobs.

    This class manages xTB-specific parameters including GFN versions,
    optimization levels, and solvation models.

    Key xTB-specific parameters:
    - gfn_version: GFN0, GFN1, GFN2, GFN-FF
    - optimization_level: crude to extreme (xTB-specific convergence criteria)
    - solvent_model: GBSA, ALPB, COSMO, etc.
    - job_type: sp, opt, hess, md, ohess, omd, path, modef

    Attributes:
        gfn_version (str): GFN-xTB version ('gfn0', 'gfn1', 'gfn2', 'gfnff').
        optimization_level (str): Optimization convergence level.
        charge (int): Molecular charge.
        multiplicity (int): Spin multiplicity (uhf parameter in xTB).
        job_type (str): Type of calculation to perform.
        title (str): Job title/description.
        opt (bool): Whether to perform a prior optimization.
        grad (bool): Whether to calculate gradient (forces).
        solvent_model (str): Implicit solvation model.
        solvent_id (str): Solvent identifier/name.
        electronic_temperature (float): Electronic temperature (K).
        accuracy (float): Numerical accuracy for xTB calculations.
        iterations (int): Maximum SCC iterations.
        input_string (str): Custom input file content.
    """

    def __init__(
        self,
        gfn_version="gfn2",
        optimization_level="vtight",
        charge=0,
        multiplicity=1,
        job_type=None,
        title=None,
        opt=False,
        grad=False,
        solvent_model=None,
        solvent_id=None,
        electronic_temperature=None,
        accuracy=1.0,
        iterations=250,
        input_string=None,
        **kwargs,
    ):
        """
        Initialize xTB job settings.

        Args:
            gfn_version (str): GFN method version. Options: 'gfn0', 'gfn1', 'gfn2', 'gfnff'.
                Defaults to 'gfn2' (most accurate and balanced).
            optimization_level (str): Optimization convergence level.
                Options: 'crude', 'sloppy', 'loose', 'lax', 'normal', 'tight', 'vtight', 'extreme'.
                Defaults to 'vtight'.
            charge (int): Molecular charge. Defaults to 0.
            multiplicity (int): Spin multiplicity (2S+1). Defaults to 1 (singlet).
            job_type (str): Calculation type. Options: 'sp', 'opt', 'hess', 'md', 'path', 'modef'.
            title (str): Job title/description.
            opt (bool): Whether to perform a prior optimization. Defaults to False.
            grad (bool): Calculate gradient (forces). Defaults to False.
            solvent_model (str): Implicit solvent model.
                Options: 'gbsa', 'alpb', 'cosmo', 'tmcosmo', 'cpcmx'.
            solvent_id (str): Solvent name (e.g., 'water', 'acetone', 'dmso').
            electronic_temperature (float): Electronic temperature (K). Defaults to None.
            accuracy (float): Numerical accuracy. Defaults to 1.0.
            iterations (int): Maximum SCC iterations. Defaults to 250.
            input_string (str): Custom xcontrol input content.
            **kwargs: Additional arguments (for compatibility).
        """
        # Validate xTB-specific parameters
        if gfn_version.lower() not in XTB_ALL_METHODS:
            logger.warning(
                f"GFN version '{gfn_version}' not in standard xTB methods: {XTB_ALL_METHODS}. "
                f"Proceeding anyway."
            )

        if (
            optimization_level is not None
            and optimization_level.lower() not in XTB_ALL_OPT_LEVELS
        ):
            logger.warning(
                f"Optimization level '{optimization_level}' not in standard xTB levels: "
                f"{XTB_ALL_OPT_LEVELS}. Proceeding anyway."
            )

        if (
            solvent_model is not None
            and solvent_model.lower() not in XTB_ALL_SOLVENT_MODELS
        ):
            logger.warning(
                f"Solvent model '{solvent_model}' not in standard xTB models: "
                f"{XTB_ALL_SOLVENT_MODELS}. Proceeding anyway."
            )

        # Core xTB parameters
        self.gfn_version = gfn_version
        self.optimization_level = optimization_level
        self.charge = charge
        self.multiplicity = multiplicity
        self.job_type = job_type
        self.title = title if title is not None else "xTB calculation"
        self.opt = opt
        self.grad = grad

        # Solvation parameters
        self.solvent_model = solvent_model
        self.solvent_id = solvent_id

        # Advanced xTB parameters
        self.electronic_temperature = electronic_temperature
        self.accuracy = accuracy
        self.iterations = iterations
        self.input_string = input_string

    def remove_solvent(self):
        """Remove implicit solvation from the calculation."""
        self.solvent_model = None
        self.solvent_id = None
        logger.debug("Removed solvent settings")

    def update_solvent(self, solvent_model=None, solvent_id=None):
        """
        Update implicit solvation settings.

        Args:
            solvent_model (str, optional): Solvent model (e.g., 'gbsa', 'alpb', 'cosmo').
            solvent_id (str, optional): Solvent identifier (e.g., 'water', 'acetone').
        """
        if solvent_model is not None:
            if solvent_model.lower() not in XTB_ALL_SOLVENT_MODELS:
                logger.warning(
                    f"Solvent model '{solvent_model}' not in standard xTB models: "
                    f"{XTB_ALL_SOLVENT_MODELS}"
                )
            self.solvent_model = solvent_model
            logger.debug(f"Updated solvent model to: {solvent_model}")

        if solvent_id is not None:
            if solvent_id.lower() not in XTB_ALL_SOLVENT_IDS:
                logger.warning(
                    f"Solvent ID '{solvent_id}' not in standard xTB solvents: "
                    f"{XTB_ALL_SOLVENT_IDS}"
                )
            self.solvent_id = solvent_id
            logger.debug(f"Updated solvent ID to: {solvent_id}")

    def modify_solvent(self, remove_solvent=False, **kwargs):
        """
        Modify solvation settings.

        Args:
            remove_solvent (bool): If True, remove solvation. Otherwise update.
            **kwargs: Solvent parameters to update (solvent_model, solvent_id).
        """
        if not remove_solvent:
            self.update_solvent(**kwargs)
        else:
            self.remove_solvent()

    @classmethod
    def from_dict(cls, settings_dict):
        return cls(**settings_dict)

    def merge(
        self,
        other,
        keywords=("charge", "multiplicity", "title"),
        merge_all=False,
    ):
        """Overwrite self settings with other settings."""

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
            XTBJobSettings: A new instance with copied settings object.
        """
        return copy.deepcopy(self)

    def __getitem__(self, key):
        return self.__dict__[key]

    def __eq__(self, other):
        """Two settings objects are equal if all their attributes are equal."""
        if type(self) is not type(other):
            return NotImplemented

        # Exclude append_additional_info from the comparison
        self_dict = self.__dict__.copy()
        self_dict.pop("append_additional_info", None)

        other_dict = other.__dict__.copy()
        other_dict.pop("append_additional_info", None)

        is_equal = self_dict == other_dict
        if not is_equal:
            import dictdiffer

            logger.info("xTB job settings are not equal.")
            for diff in list(dictdiffer.diff(self_dict, other_dict)):
                logger.info(f"Difference: {diff}")
        return is_equal

    @classmethod
    def from_comfile(cls, filename):
        """Return xTB job settings from Gaussian .com file.

        Args:
            filename (str): Path to the Gaussian .com file.

        Returns:
            XTBJobSettings: Settings object from com file.
        """
        from chemsmart.io.gaussian.input import Gaussian16Input

        com_path = os.path.abspath(filename)
        logger.info(
            f"Reading molecular info from Gaussian com file: {com_path}"
        )
        gaussian_settings_from_comfile = Gaussian16Input(
            filename=com_path
        ).read_settings()
        xtb_default_settings = cls.default()
        basic_properties = ("charge", "multiplicity", "title")
        return xtb_default_settings.merge(
            gaussian_settings_from_comfile, keywords=basic_properties
        )

    @classmethod
    def from_inpfile(cls, filename):
        """Return xTB job settings from ORCA .inp file.

        Args:
            filename (str): Path to the ORCA .inp file.

        Returns:
            XTBJobSettings: Settings object from inp file.
        """
        from chemsmart.io.orca.input import ORCAInput

        inp_path = os.path.abspath(filename)
        logger.info(f"Reading molecular info from ORCA inp file: {inp_path}")
        orca_settings_from_inpfile = ORCAInput(
            filename=inp_path
        ).read_settings()
        xtb_default_settings = cls.default()
        basic_properties = ("charge", "multiplicity", "title")
        return xtb_default_settings.merge(
            orca_settings_from_inpfile, keywords=basic_properties
        )

    @classmethod
    def from_logfile(cls, filename):
        """Return xTB job settings from Gaussian .log file.

        Args:
            filename (str): Path to the Gaussian .log file.

        Returns:
            XTBJobSettings: Settings object from log file.
        """
        log_path = os.path.abspath(filename)
        from chemsmart.io.gaussian.output import (
            Gaussian16Output,
            Gaussian16OutputWithPBC,
        )

        log_path = os.path.abspath(filename)
        logger.info(
            f"Reading molecular info from Gaussian log file: {log_path}"
        )
        try:
            gaussian_settings_from_logfile = Gaussian16Output(
                filename=log_path
            ).read_settings()
        except ValueError:
            gaussian_settings_from_logfile = Gaussian16OutputWithPBC(
                filename=log_path
            ).read_settings()

        xtb_default_settings = cls.default()
        basic_properties = ("charge", "multiplicity", "title")
        return xtb_default_settings.merge(
            gaussian_settings_from_logfile, keywords=basic_properties
        )

    @classmethod
    def from_outfile(cls, filename):
        """
        Return xTB job settings from ORCA .out file.

        Args:
            filename (str): Path to the ORCA .out file.

        Returns:
            XTBJobSettings: Settings object from out file.
        """
        from chemsmart.io.orca.output import ORCAOutput

        out_path = os.path.abspath(filename)
        logger.info(f"Reading molecular info from ORCA out file: {out_path}")
        orca_settings_from_outfile = ORCAOutput(
            filename=out_path
        ).read_settings()
        xtb_default_settings = cls.default()
        basic_properties = ("charge", "multiplicity", "title")
        return xtb_default_settings.merge(
            orca_settings_from_outfile, keywords=basic_properties
        )

    @classmethod
    def default(cls):
        """
        Get default xTB job settings.

        Returns:
            XTBJobSettings: Default settings object.
        """
        return cls(
            gfn_version="gfn2",
            optimization_level="vtight",
            charge=0,
            multiplicity=1,
            job_type=None,
            title="xTB calculation with default settings",
            opt=False,
            grad=False,
            solvent_model=None,
            solvent_id=None,
            electronic_temperature=None,
            accuracy=1.0,
            iterations=250,
            input_string=None,
        )

    @classmethod
    def from_filepath(cls, filepath, **kwargs):
        """
        Create settings from any supported file type.

        Args:
            filepath (str): Path to the input file (.com, .gjf, .inp, .log, .out, .xyz).
            **kwargs: Additional arguments.

        Returns:
            XTBJobSettings: The created settings object.

        Raises:
            ValueError: If the file extension is not supported.
        """
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
