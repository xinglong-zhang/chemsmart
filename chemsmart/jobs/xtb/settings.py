"""
xTB job settings implementation.

This module contains the XTBJobSettings class for configuring xTB jobs.
It manages parameters such as GFN version, optimization level, and solvent models.
"""

import copy
import logging
import os

from chemsmart.jobs.settings import MolecularJobSettings
from chemsmart.utils.periodictable import PeriodicTable

pt = PeriodicTable()


logger = logging.getLogger(__name__)


class XTBJobSettings(MolecularJobSettings):
    """
    Settings for xTB jobs.

    Manages configuration parameters for xTB calculations, including
    GFN version, optimization level, charge, and solvent settings.

    Attributes:
        gfn_version (str): GFN-xTB version (e.g., 'gfn2').
        optimization_level (str): Optimization level (e.g., 'vtight').
    """

    def __init__(
        self,
        gfn_version="gfn2",  # use gfn2 by default
        optimization_level=None,
        charge=None,
        multiplicity=None,
        job_type=None,
        title=None,
        freq=False,
        solvent_model=None,
        solvent_id=None,
        **kwargs,
    ):
        """
        Initialize xTB job settings.

        Args:
            gfn_version (str, optional): GFN-xTB version. Defaults to "gfn2".
            optimization_level (str, optional): Optimization level.
            charge (int, optional): Molecular charge.
            multiplicity (int, optional): Spin multiplicity.
            job_type (str, optional): Type of job.
            title (str, optional): Job title.
            freq (bool, optional): Whether to calculate frequencies. Defaults to False.
            solvent_model (str, optional): Solvent model to use.
            solvent_id (str, optional): Solvent identifier.
            **kwargs: Additional arguments.
        """
        super().__init__(
            charge=charge,
            multiplicity=multiplicity,
            freq=freq,
            job_type=job_type,
            title=title,
            solvent_model=solvent_model,
            solvent_id=solvent_id,
            **kwargs,
        )
        self.gfn_version = gfn_version
        self.optimization_level = optimization_level

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
        self_dict = self.__dict__
        self_dict.pop("append_additional_info")

        other_dict = other.__dict__
        other_dict.pop("append_additional_info")

        is_equal = self_dict == other_dict
        if not is_equal:
            import dictdiffer

            logger.info("xTB job settings are not equal.")
            for diff in list(dictdiffer.diff(self_dict, other_dict)):
                logger.info(f"Difference: {diff}")
        return self_dict == other_dict

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
        gaussian_settings_from_comfile = Gaussian16Input(
            filename=com_path
        ).read_settings()
        xtb_default_settings = cls.default()
        return xtb_default_settings.merge(
            gaussian_settings_from_comfile, merge_all=True
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
        logger.info(f"Return Settings object from inp file: {inp_path}")
        orca_settings_from_inpfile = ORCAInput(
            filename=inp_path
        ).read_settings()
        xtb_default_settings = cls.default()
        return xtb_default_settings.merge(
            orca_settings_from_inpfile, merge_all=True
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

        logger.info(f"Return Settings object from logfile: {log_path}")
        try:
            settings = Gaussian16Output(filename=log_path).read_settings()
        except ValueError:
            settings = Gaussian16OutputWithPBC(
                filename=log_path
            ).read_settings()

        xtb_default_settings = cls.default()
        return xtb_default_settings.merge(settings, merge_all=True)

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
        logger.info(
            f"Return Settings object from ORCA .out filename: {out_path}"
        )
        orca_settings_from_outfile = ORCAOutput(
            filename=out_path
        ).read_settings()
        xtb_default_settings = cls.default()
        return xtb_default_settings.merge(
            orca_settings_from_outfile, merge_all=True
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
            charge=None,
            multiplicity=None,
            job_type=None,
            title="xtb job with default settings",
            freq=False,
            solvent_model=None,
            solvent_id=None,
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
