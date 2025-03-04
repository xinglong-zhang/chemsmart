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
        gfn_version="gfn2",  # use gfn2 by default
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
        self.gfn_version = gfn_version
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
        xtb_default_settings = cls.default()
        return xtb_default_settings.merge(gaussian_settings_from_comfile, merge_all=True)

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
        xtb_default_settings = cls.default()
        return xtb_default_settings.merge(orca_settings_from_inpfile, merge_all=True)

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

        xtb_default_settings = cls.default()
        return xtb_default_settings.merge(settings, merge_all=True)

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
        xtb_default_settings = cls.default()
        return xtb_default_settings.merge(
            orca_settings_from_outfile, merge_all=True
        )

    @classmethod
    def default(cls):
        return cls(
            gfn_version="gfn2",
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
