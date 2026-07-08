import copy
import logging
import os

from chemsmart.io.xtb import (
    XTB_ALL_JOB_TYPES,
    XTB_ALL_METHODS,
    XTB_ALL_OPT_LEVELS,
    XTB_ALL_SOLVENT_IDS,
    XTB_ALL_SOLVENT_MODELS,
)

logger = logging.getLogger(__name__)


class XTBJobSettings:
    """Settings for xTB command-line jobs."""

    def __init__(
        self,
        gfn_version="gfn2",
        optimization_level="vtight",
        charge=0,
        multiplicity=1,
        jobtype=None,
        title=None,
        grad=False,
        solvent_model=None,
        solvent_id=None,
        input_string=None,
        **kwargs,
    ):
        if gfn_version is not None:
            gfn_version = gfn_version.lower()
        if optimization_level is not None:
            optimization_level = optimization_level.lower()
        if jobtype is not None:
            jobtype = jobtype.lower()
        if solvent_model is not None:
            solvent_model = solvent_model.lower()
        if solvent_id is not None:
            solvent_id = solvent_id.lower()

        self._warn_if_unknown(gfn_version, XTB_ALL_METHODS, "GFN version")
        self._warn_if_unknown(
            optimization_level, XTB_ALL_OPT_LEVELS, "optimization level"
        )
        self._warn_if_unknown(jobtype, XTB_ALL_JOB_TYPES, "job type")
        self._warn_if_unknown(
            solvent_model, XTB_ALL_SOLVENT_MODELS, "solvent model"
        )
        self._warn_if_unknown(solvent_id, XTB_ALL_SOLVENT_IDS, "solvent id")

        self.gfn_version = gfn_version
        self.optimization_level = optimization_level
        self.charge = charge
        self.multiplicity = multiplicity
        self.jobtype = jobtype
        self.title = title
        self.grad = grad
        self.solvent_model = solvent_model
        self.solvent_id = solvent_id
        self.input_string = input_string

    @staticmethod
    def _warn_if_unknown(value, known_values, label):
        if value is not None and value not in known_values:
            logger.warning(
                f"{label} {value!r} is not in the known xTB values: "
                f"{known_values}"
            )

    @classmethod
    def default(cls):
        return cls()

    @classmethod
    def from_dict(cls, settings_dict):
        return cls(**settings_dict)

    def copy(self):
        return copy.deepcopy(self)

    def merge(
        self,
        other,
        keywords=("charge", "multiplicity", "title"),
        merge_all=False,
    ):
        other_dict = other if isinstance(other, dict) else other.__dict__
        if merge_all:
            merged_dict = self.__dict__.copy()
            merged_dict.update(other_dict)
            logger.debug(f"Merged all xTB settings: {merged_dict}")
            return type(self)(**merged_dict)

        if keywords is not None:
            other_dict = {
                key: other_dict[key]
                for key in keywords
                if key in other_dict and other_dict[key] is not None
            }
        merged_dict = self.__dict__.copy()
        merged_dict.update(other_dict)
        logger.debug(
            f"Merged xTB settings with keywords {keywords}: {merged_dict}"
        )
        return type(self)(**merged_dict)

    def remove_solvent(self):
        self.solvent_model = None
        self.solvent_id = None

    def update_solvent(self, solvent_model=None, solvent_id=None):
        if solvent_model is not None:
            self._warn_if_unknown(
                solvent_model.lower(), XTB_ALL_SOLVENT_MODELS, "solvent model"
            )
            self.solvent_model = solvent_model.lower()
        if solvent_id is not None:
            self._warn_if_unknown(
                solvent_id.lower(), XTB_ALL_SOLVENT_IDS, "solvent id"
            )
            self.solvent_id = solvent_id.lower()

    def modify_solvent(self, remove_solvent=False, **kwargs):
        if remove_solvent:
            self.remove_solvent()
        else:
            self.update_solvent(**kwargs)

    @classmethod
    def from_filepath(cls, filepath, **kwargs):
        filepath = os.path.abspath(filepath)
        if filepath.endswith((".com", ".gjf")):
            from chemsmart.io.gaussian.input import Gaussian16Input

            return cls.default().merge(
                Gaussian16Input(filename=filepath).read_settings(),
                keywords=("charge", "multiplicity", "title"),
            )
        if filepath.endswith(".log"):
            from chemsmart.io.gaussian.output import Gaussian16Output

            return cls.default().merge(
                Gaussian16Output(filename=filepath).read_settings(),
                keywords=("charge", "multiplicity", "title"),
            )
        if filepath.endswith(".inp"):
            from chemsmart.io.orca.input import ORCAInput

            return cls.default().merge(
                ORCAInput(filename=filepath).read_settings(),
                keywords=("charge", "multiplicity", "title"),
            )
        if filepath.endswith(".out"):
            from chemsmart.io.orca.output import ORCAOutput

            return cls.default().merge(
                ORCAOutput(filename=filepath).read_settings(),
                keywords=("charge", "multiplicity", "title"),
            )
        return cls.default()

    def __eq__(self, other):
        if type(self) is not type(other):
            return NotImplemented
        return self.__dict__ == other.__dict__
