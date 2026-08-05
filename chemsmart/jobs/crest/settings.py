import copy
import logging
import os

from chemsmart.io.crest import (
    CREST_ALL_JOB_TYPES,
    CREST_ALL_METHODS,
    CREST_ALL_OPT_LEVELS,
    CREST_ALL_QUICK_MODES,
    CREST_ALL_SOLVENT_IDS,
    CREST_ALL_SOLVENT_MODELS,
)

logger = logging.getLogger(__name__)


class CRESTJobSettings:
    """Settings for CREST conformational search jobs."""

    def __init__(
        self,
        gfn_version="gfn2",
        charge=None,
        multiplicity=None,
        jobtype=None,
        solvent_model=None,
        solvent_id=None,
        energy_window=None,
        rmsd_threshold=None,
        energy_threshold=None,
        bconst_threshold=None,
        population_threshold=None,
        temperature=None,
        md_timestep=None,
        md_length=None,
        md_dump_step=None,
        vbias_dump_interval=None,
        additional_md_temperature=None,
        no_topology_check=False,
        no_reference_topology_check=False,
        optimization_level=None,
        quick_mode=None,
        nci=False,
        constraints=None,
        force_constant=None,
        nprocs=None,
        additional_flags=None,
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
        if quick_mode is not None:
            quick_mode = quick_mode.lower()

        self._warn_if_unknown(gfn_version, CREST_ALL_METHODS, "GFN version")
        self._warn_if_unknown(
            optimization_level, CREST_ALL_OPT_LEVELS, "optimization level"
        )
        self._warn_if_unknown(jobtype, CREST_ALL_JOB_TYPES, "job type")
        self._warn_if_unknown(
            solvent_model, CREST_ALL_SOLVENT_MODELS, "solvent model"
        )
        self._warn_if_unknown(solvent_id, CREST_ALL_SOLVENT_IDS, "solvent id")
        self._warn_if_unknown(quick_mode, CREST_ALL_QUICK_MODES, "quick mode")

        self.gfn_version = gfn_version
        self.charge = charge
        self.multiplicity = multiplicity
        self.jobtype = jobtype
        self.solvent_model = solvent_model
        self.solvent_id = solvent_id
        self.energy_window = energy_window
        self.rmsd_threshold = rmsd_threshold
        self.energy_threshold = energy_threshold
        self.bconst_threshold = bconst_threshold
        self.population_threshold = population_threshold
        self.temperature = temperature
        self.md_timestep = md_timestep
        self.md_length = md_length
        self.md_dump_step = md_dump_step
        self.vbias_dump_interval = vbias_dump_interval
        self.additional_md_temperature = additional_md_temperature
        self.no_topology_check = no_topology_check
        self.no_reference_topology_check = no_reference_topology_check
        self.optimization_level = optimization_level
        self.quick_mode = quick_mode
        self.nci = nci
        self.constraints = constraints
        self.force_constant = force_constant
        self.nprocs = nprocs
        self.additional_flags = additional_flags

    @staticmethod
    def _warn_if_unknown(value, known_values, label):
        if value is not None and value not in known_values:
            logger.warning(
                f"{label} {value!r} is not in the known CREST values: "
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
        record metadata, CREST-relevant fields (gfn_version, solvent).
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
            if method_lower in CREST_ALL_METHODS:
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
        keywords=("charge", "multiplicity"),
        merge_all=False,
    ):
        other_dict = other if isinstance(other, dict) else other.__dict__
        if merge_all:
            merged_dict = self.__dict__.copy()
            merged_dict.update(other_dict)
            logger.debug(f"Merged all CREST settings: {merged_dict}")
            return type(self)(**merged_dict)

        if keywords is not None:
            other_dict = {
                key: other_dict[key] for key in keywords if key in other_dict
            }
        merged_dict = self.__dict__.copy()
        merged_dict.update(other_dict)
        logger.debug(
            f"Merged CREST settings with keywords {keywords}: {merged_dict}"
        )
        return type(self)(**merged_dict)

    def remove_solvent(self):
        self.solvent_model = None
        self.solvent_id = None

    def update_solvent(self, solvent_model=None, solvent_id=None):
        if solvent_model is not None:
            self._warn_if_unknown(
                solvent_model.lower(),
                CREST_ALL_SOLVENT_MODELS,
                "solvent model",
            )
            self.solvent_model = solvent_model.lower()
        if solvent_id is not None:
            self._warn_if_unknown(
                solvent_id.lower(), CREST_ALL_SOLVENT_IDS, "solvent id"
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
                keywords=("charge", "multiplicity"),
            )
        if filepath.endswith(".log"):
            from chemsmart.io.gaussian.output import Gaussian16Output

            return cls.default().merge(
                Gaussian16Output(filename=filepath).read_settings(),
                keywords=("charge", "multiplicity"),
            )
        if filepath.endswith(".inp"):
            from chemsmart.io.orca.input import ORCAInput

            return cls.default().merge(
                ORCAInput(filename=filepath).read_settings(),
                keywords=("charge", "multiplicity"),
            )
        if filepath.endswith(".out"):
            from chemsmart.utils.io import get_program_type_from_file

            program = get_program_type_from_file(filepath)
            if program == "gaussian":
                from chemsmart.io.gaussian.output import Gaussian16Output

                return cls.default().merge(
                    Gaussian16Output(filename=filepath).read_settings(),
                    keywords=("charge", "multiplicity"),
                )
            if program == "orca":
                from chemsmart.io.orca.output import ORCAOutput

                return cls.default().merge(
                    ORCAOutput(filename=filepath).read_settings(),
                    keywords=("charge", "multiplicity"),
                )
            if program == "xtb":
                from chemsmart.jobs.xtb.settings import XTBJobSettings

                return cls.default().merge(
                    XTBJobSettings.from_outfile(filepath),
                    keywords=("charge", "multiplicity"),
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
