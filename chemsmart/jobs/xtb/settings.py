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
        additional_route_parameters=None,
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
        self.additional_route_parameters = additional_route_parameters

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

    @classmethod
    def from_database(
        cls,
        filepath,
        record_index=None,
        record_id=None,
        structure_index="-1",
        structure_id=None,
    ):
        """Create job settings from a chemsmart database file.

        With record selectors (record_index/record_id), this fills
        charge/multiplicity from the selected structure and, when present in
        record metadata, xTB-relevant fields (gfn_version, solvent).
        With a global structure selector (structure_id), this uses defaults
        and fills only charge/multiplicity from the selected structure.
        """
        from chemsmart.database.database import Database
        from chemsmart.database.utils import resolve_record
        from chemsmart.utils.utils import string2index_1based

        if not os.path.isfile(filepath):
            raise FileNotFoundError(f"Database file not found: {filepath}")

        db = Database(filepath)
        record_selected = record_index is not None or record_id is not None
        if structure_id is not None and record_selected:
            raise ValueError(
                "Use either structure_id or record_index/record_id, not both."
            )

        settings = cls.default()

        if structure_id is not None:
            full_sid = db.get_structure_by_partial_id(structure_id)
            structure = db.get_structure(full_sid)
            if structure is None:
                raise ValueError(
                    f"No structure found with ID '{structure_id}'."
                )
            settings.charge = structure.get("charge")
            settings.multiplicity = structure.get("multiplicity")
            settings.title = (
                "Job prepared from chemsmart database "
                f"{os.path.basename(filepath)}"
            )
            logger.info(
                "Created JobSettings from database: "
                f"charge={settings.charge}, "
                f"multiplicity={settings.multiplicity}"
            )
            return settings

        record = resolve_record(
            db,
            record_index=record_index,
            record_id=record_id,
            return_list=False,
        )
        if record is None:
            return None

        meta = record.get("meta", {})
        molecules = record.get("molecules", [])
        selected_index = string2index_1based(str(structure_index))
        if isinstance(selected_index, slice):
            raise ValueError(
                "Database-aware jobs support one structure at a time."
            )
        try:
            structure = molecules[selected_index] if molecules else {}
        except IndexError as exc:
            raise ValueError(
                f"Structure index {structure_index} out of range for "
                "selected record."
            ) from exc
        settings.charge = structure.get("charge")
        settings.multiplicity = structure.get("multiplicity")
        method = meta.get("method")
        if method is not None:
            method_lower = str(method).lower()
            if method_lower in XTB_ALL_METHODS:
                settings.gfn_version = method_lower
        settings.solvent_model = meta.get("solvent_model")
        settings.solvent_id = meta.get("solvent_id")
        settings.title = (
            "Job prepared from chemsmart database "
            f"{os.path.basename(filepath)}"
        )
        logger.info(
            "Created JobSettings from database: "
            f"charge={settings.charge}, "
            f"multiplicity={settings.multiplicity}"
        )
        return settings

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
    def from_outfile(cls, filepath):
        """Create xTB settings from an xTB main output file.

        Inherits all xTB-specific settings: GFN version, optimization level,
        charge, multiplicity, jobtype, gradient, and solvent configuration.
        This is same-program inheritance, analogous to Gaussian from_logfile
        or ORCA from_outfile.
        """
        from chemsmart.io.xtb.file import XTBMainOut

        output = XTBMainOut(filename=os.path.abspath(filepath))

        return cls(
            gfn_version=output.method,
            optimization_level=output.optimization_level,
            charge=output.net_charge,
            multiplicity=output.multiplicity,
            jobtype=output.jobtype,
            title=(
                f"Job prepared from xTB file " f"{os.path.basename(filepath)}"
            ),
            grad=output.route_object.grad,
            solvent_model=output.solvent_model,
            solvent_id=output.solvent_id,
        )

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
            from chemsmart.utils.io import get_program_type_from_file

            program = get_program_type_from_file(filepath)
            if program == "xtb":
                return cls.from_outfile(filepath)
            if program == "gaussian":
                from chemsmart.io.gaussian.output import Gaussian16Output

                return cls.default().merge(
                    Gaussian16Output(filename=filepath).read_settings(),
                    keywords=("charge", "multiplicity", "title"),
                )
            if program == "orca":
                from chemsmart.io.orca.output import ORCAOutput

                return cls.default().merge(
                    ORCAOutput(filename=filepath).read_settings(),
                    keywords=("charge", "multiplicity", "title"),
                )
            raise ValueError(
                f"Unsupported .out file program type: {program}. "
                "Only Gaussian, ORCA, and xTB outputs are supported."
            )
        return cls.default()

    def __eq__(self, other):
        if type(self) is not type(other):
            return NotImplemented
        return self.__dict__ == other.__dict__
