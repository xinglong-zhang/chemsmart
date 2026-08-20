"""Strict settings contract for ChemSmart-managed xTB jobs."""

import copy
import logging
import math
import os
import unicodedata
from collections.abc import Mapping
from dataclasses import asdict, dataclass

from chemsmart.io.xtb import (
    XTB_ALL_METHODS,
    XTB_ALL_OPT_LEVELS,
    XTB_ALL_SOLVENT_IDS,
    XTB_ALL_SOLVENT_MODELS,
)

logger = logging.getLogger(__name__)

_XTB_671_NAMED_SOLVENTS = frozenset(XTB_ALL_SOLVENT_IDS)
_XTB_671_ALPB_SOLVENTS = _XTB_671_NAMED_SOLVENTS - {
    "h2o",
    "n-hexane",
}
_XTB_671_GBSA_GFN1_SOLVENTS = frozenset(
    {
        "acetone",
        "acetonitrile",
        "benzene",
        "ch2cl2",
        "chcl3",
        "cs2",
        "dmso",
        "ether",
        "h2o",
        "methanol",
        "thf",
        "toluene",
        "water",
    }
)
_XTB_671_GBSA_GFN2_SOLVENTS = frozenset(
    {
        "acetone",
        "acetonitrile",
        "ch2cl2",
        "chcl3",
        "cs2",
        "dmf",
        "dmso",
        "ether",
        "h2o",
        "methanol",
        "n-hexane",
        "thf",
        "toluene",
        "water",
    }
)


@dataclass(frozen=True)
class XTBSolventCapabilityV1:
    """Pinned-native support observation for one xTB solvent pair.

    ``unknown`` is deliberately different from ``confirmed_unsupported``.
    Unknown, syntactically safe identifiers may be retained in a plan, but a
    real execution preflight must resolve them before launching xTB.
    """

    schema_version: str
    required_version: str
    method: str
    model: str | None
    identifier: str | None
    status: str
    rule_id: str
    evidence: str

    def as_dict(self):
        return asdict(self)


class XTBJobSettings:
    """Validated settings for the executable xTB ``sp``, ``opt`` and ``hess`` surface.

    This class is deliberately narrower than the xTB output parser.  The parser
    understands additional native xTB task families, while the execution plane
    only exposes the three job kinds that have a ChemSmart CLI contract.
    """

    FIELDS = frozenset(
        {
            "gfn_version",
            "optimization_level",
            "charge",
            "multiplicity",
            "jobtype",
            "grad",
            "solvent_model",
            "solvent_id",
        }
    )
    JOBTYPES = ("sp", "opt", "hess")
    GFN_VERSIONS = tuple(XTB_ALL_METHODS)
    OPTIMIZATION_LEVELS = tuple(XTB_ALL_OPT_LEVELS)
    SOLVENT_MODELS = tuple(XTB_ALL_SOLVENT_MODELS)

    def __init__(
        self,
        gfn_version="gfn2",
        optimization_level="vtight",
        charge=0,
        multiplicity=1,
        jobtype=None,
        grad=False,
        solvent_model=None,
        solvent_id=None,
    ):
        self.gfn_version = self._normalized_text(gfn_version)
        self.optimization_level = self._normalized_text(optimization_level)
        self.charge = charge
        self.multiplicity = multiplicity
        self.jobtype = self._normalized_text(jobtype)
        self.grad = grad
        self.solvent_model = self._normalized_text(solvent_model)
        self.solvent_id = self._normalized_solvent_id(solvent_id)
        self.validate()

    @staticmethod
    def _normalized_text(value):
        if value is None:
            return None
        if not isinstance(value, str):
            raise TypeError(f"Expected a string, got {type(value).__name__}.")
        return value.strip().lower()

    @staticmethod
    def _normalized_solvent_id(value):
        if value is None:
            return None
        if isinstance(value, bool) or not isinstance(value, (str, int, float)):
            raise TypeError(
                "xTB solvent_id must be a non-empty string or numeric "
                "dielectric identifier."
            )
        raw_value = str(value)
        if any(
            unicodedata.category(char).startswith("C") for char in raw_value
        ):
            raise ValueError(
                "xTB solvent_id must not contain control characters."
            )
        normalized = raw_value.strip().lower()
        if not normalized:
            raise ValueError("xTB solvent_id must not be empty.")
        if normalized.startswith("-"):
            raise ValueError(
                "xTB solvent_id must be data, not an option-like token."
            )
        return normalized

    @classmethod
    def default(cls):
        """Return ChemSmart's explicit GFN2/vtight default settings."""
        return cls()

    @classmethod
    def from_dict(cls, settings_dict):
        if not isinstance(settings_dict, Mapping):
            raise TypeError(
                "xTB settings must be a mapping, got "
                f"{type(settings_dict).__name__}."
            )
        unknown = sorted(set(settings_dict) - cls.FIELDS)
        if unknown:
            raise ValueError(
                "Unknown xTB setting key(s): " + ", ".join(unknown)
            )
        return cls(**dict(settings_dict))

    def validate(self, expected_jobtype=None):
        """Fail closed on values that cannot be rendered unambiguously."""
        if self.gfn_version not in self.GFN_VERSIONS:
            raise ValueError(
                f"Unsupported xTB gfn_version {self.gfn_version!r}; expected "
                f"one of {self.GFN_VERSIONS}."
            )
        if self.optimization_level not in self.OPTIMIZATION_LEVELS:
            raise ValueError(
                "Unsupported xTB optimization_level "
                f"{self.optimization_level!r}; expected one of "
                f"{self.OPTIMIZATION_LEVELS}."
            )
        if self.jobtype is not None and self.jobtype not in self.JOBTYPES:
            raise ValueError(
                f"Unsupported xTB jobtype {self.jobtype!r}; expected one of "
                f"{self.JOBTYPES}."
            )
        if expected_jobtype is not None:
            expected = self._normalized_text(expected_jobtype)
            if expected not in self.JOBTYPES:
                raise ValueError(
                    f"Unknown expected xTB jobtype: {expected!r}."
                )
            if self.jobtype != expected:
                raise ValueError(
                    "xTB settings/job mismatch: "
                    f"settings declare {self.jobtype!r}, job requires {expected!r}."
                )
        if isinstance(self.charge, bool) or not isinstance(self.charge, int):
            raise TypeError(
                "xTB charge must be an integer (bool is not valid)."
            )
        if (
            isinstance(self.multiplicity, bool)
            or not isinstance(self.multiplicity, int)
            or self.multiplicity < 1
        ):
            raise ValueError("xTB multiplicity must be a positive integer.")
        if not isinstance(self.grad, bool):
            raise TypeError("xTB grad must be a boolean.")
        if self.grad and self.jobtype not in (None, "sp"):
            raise ValueError(
                "xTB --grad is supported only for the sp job kind; it must "
                "not be combined with opt or hess."
            )

        has_model = self.solvent_model is not None
        has_identifier = self.solvent_id is not None
        if has_model != has_identifier:
            raise ValueError(
                "xTB solvation requires solvent_model and solvent_id together."
            )
        if has_model and self.solvent_model not in self.SOLVENT_MODELS:
            raise ValueError(
                f"Unsupported xTB solvent_model {self.solvent_model!r}; "
                f"expected one of {self.SOLVENT_MODELS}."
            )
        if has_model:
            self._validate_solvent_pair()
        return self

    def _validate_solvent_pair(self):
        """Reject malformed or confirmed-unsupported xTB solvent pairs.

        A well-formed identifier absent from the pinned 6.7.1 vocabulary is
        not evidence of non-support.  Such a pair remains ``unknown`` for a
        later environment-bound preflight rather than being falsely rejected
        while the workflow is still being planned.
        """

        capability = self.solvent_capability()
        if capability.status == "confirmed_unsupported":
            raise ValueError(
                "Confirmed unsupported xTB 6.7.1 solvent pair: "
                f"{self.solvent_model}={self.solvent_id!r} for "
                f"{self.gfn_version} ({capability.rule_id})."
            )

    def solvent_capability(self):
        """Return a tri-state pinned-native capability observation."""

        if self.solvent_model is None and self.solvent_id is None:
            return XTBSolventCapabilityV1(
                schema_version="chemsmart.xtb-solvent-capability.v1",
                required_version="6.7.1",
                method=self.gfn_version,
                model=None,
                identifier=None,
                status="not_applicable",
                rule_id="xtb.solvent.not_requested",
                evidence="settings:no_solvent",
            )

        identifier = self.solvent_id
        try:
            numeric_identifier = float(identifier)
            identifier_is_numeric = True
        except (TypeError, ValueError):
            numeric_identifier = None
            identifier_is_numeric = False

        if identifier_is_numeric:
            if self.solvent_model not in ("cosmo", "tmcosmo"):
                return XTBSolventCapabilityV1(
                    schema_version="chemsmart.xtb-solvent-capability.v1",
                    required_version="6.7.1",
                    method=self.gfn_version,
                    model=self.solvent_model,
                    identifier=identifier,
                    status="confirmed_unsupported",
                    rule_id="xtb.solvent.numeric_model_unsupported",
                    evidence="xTB-6.7.1:native-cli-contract",
                )
            if (
                not math.isfinite(numeric_identifier)
                or numeric_identifier <= 0
            ):
                raise ValueError(
                    "An xTB COSMO dielectric identifier must be finite and positive."
                )
            return XTBSolventCapabilityV1(
                schema_version="chemsmart.xtb-solvent-capability.v1",
                required_version="6.7.1",
                method=self.gfn_version,
                model=self.solvent_model,
                identifier=identifier,
                status="confirmed_supported",
                rule_id="xtb.solvent.numeric_dielectric_supported",
                evidence="xTB-6.7.1:native-cli-contract",
            )

        if self.solvent_model == "gbsa":
            allowed = {
                "gfn1": _XTB_671_GBSA_GFN1_SOLVENTS,
                "gfn2": _XTB_671_GBSA_GFN2_SOLVENTS,
            }.get(self.gfn_version)
            if allowed is None:
                return XTBSolventCapabilityV1(
                    schema_version="chemsmart.xtb-solvent-capability.v1",
                    required_version="6.7.1",
                    method=self.gfn_version,
                    model=self.solvent_model,
                    identifier=identifier,
                    status="confirmed_unsupported",
                    rule_id="xtb.solvent.gbsa_method_unsupported",
                    evidence="xTB-6.7.1:pinned-native-registry",
                )
        elif self.solvent_model in ("alpb", "cosmo", "tmcosmo"):
            allowed = _XTB_671_ALPB_SOLVENTS
        elif self.solvent_model == "cpcmx":
            # The execution plane intentionally exposes only its pinned,
            # reviewable vocabulary rather than accepting arbitrary database
            # labels that cannot be validated offline.
            allowed = _XTB_671_NAMED_SOLVENTS
        else:  # guarded by SOLVENT_MODELS above
            allowed = frozenset()

        status = "confirmed_supported" if identifier in allowed else "unknown"
        rule_id = (
            "xtb.solvent.pinned_pair_supported"
            if status == "confirmed_supported"
            else "xtb.solvent.pair_unverified"
        )
        evidence = (
            "xTB-6.7.1:pinned-native-registry"
            if status == "confirmed_supported"
            else "xTB-6.7.1:no-pinned-native-observation"
        )
        return XTBSolventCapabilityV1(
            schema_version="chemsmart.xtb-solvent-capability.v1",
            required_version="6.7.1",
            method=self.gfn_version,
            model=self.solvent_model,
            identifier=identifier,
            status=status,
            rule_id=rule_id,
            evidence=evidence,
        )

    @property
    def requires_solvent_clarification(self):
        return self.solvent_capability().status == "unknown"

    def validate_for_molecule(self, molecule):
        """Validate the electron-count parity implied by charge/multiplicity.

        ChemSmart maps multiplicity to xTB's ``--uhf`` as ``multiplicity - 1``.
        That number must not exceed the electron count and must have the same
        parity as the total number of electrons.
        """
        from ase.data import atomic_numbers

        self.validate()
        try:
            nuclear_charge = sum(
                atomic_numbers[str(symbol)] for symbol in molecule.symbols
            )
        except (AttributeError, KeyError) as exc:
            raise ValueError(
                "Cannot validate xTB electronic state because the molecule "
                "contains an unknown atomic symbol."
            ) from exc

        electron_count = nuclear_charge - self.charge
        unpaired_electrons = self.multiplicity - 1
        if electron_count < 0:
            raise ValueError(
                f"xTB charge {self.charge} implies a negative electron count."
            )
        if unpaired_electrons > electron_count:
            raise ValueError(
                "xTB multiplicity implies more unpaired electrons than the "
                "molecule contains."
            )
        if (electron_count - unpaired_electrons) % 2:
            raise ValueError(
                "xTB charge/multiplicity parity is inconsistent with the "
                f"molecular electron count ({electron_count} electrons, "
                f"multiplicity {self.multiplicity})."
            )
        return self

    def copy(self):
        return copy.deepcopy(self)

    def merge(
        self,
        other,
        keywords=("charge", "multiplicity"),
        merge_all=False,
    ):
        if isinstance(other, Mapping):
            other_dict = dict(other)
        else:
            other_dict = vars(other).copy()

        unknown = sorted(set(other_dict) - self.FIELDS)
        if merge_all and unknown:
            raise ValueError(
                "Unknown xTB setting key(s) in merge: " + ", ".join(unknown)
            )
        if merge_all:
            selected = other_dict
        else:
            requested = tuple(keywords or ())
            unknown_keywords = sorted(set(requested) - self.FIELDS)
            if unknown_keywords:
                raise ValueError(
                    "Unknown xTB merge key(s): " + ", ".join(unknown_keywords)
                )
            selected = {
                key: other_dict[key]
                for key in requested
                if key in self.FIELDS
                and key in other_dict
                and other_dict[key] is not None
            }

        merged = vars(self).copy()
        merged.update(selected)
        logger.debug(f"Merged xTB settings: {merged}")
        return type(self).from_dict(merged)

    def remove_solvent(self):
        self.solvent_model = None
        self.solvent_id = None
        return self

    def update_solvent(self, solvent_model=None, solvent_id=None):
        model = (
            self.solvent_model
            if solvent_model is None
            else self._normalized_text(solvent_model)
        )
        identifier = (
            self.solvent_id
            if solvent_id is None
            else self._normalized_solvent_id(solvent_id)
        )
        candidate = self.copy()
        candidate.solvent_model = model
        candidate.solvent_id = identifier
        candidate.validate()
        self.solvent_model = model
        self.solvent_id = identifier
        return self

    def modify_solvent(self, remove_solvent=False, **kwargs):
        if remove_solvent:
            return self.remove_solvent()
        return self.update_solvent(**kwargs)

    @classmethod
    def from_filepath(cls, filepath, **kwargs):
        """Read only charge and multiplicity from a supported input/output."""
        filepath = os.path.abspath(filepath)
        suffix = os.path.splitext(filepath)[1].lower()
        settings = None

        if suffix in (".com", ".gjf"):
            from chemsmart.io.gaussian.input import Gaussian16Input

            settings = Gaussian16Input(filename=filepath).read_settings()
        elif suffix == ".log":
            from chemsmart.io.gaussian.output import Gaussian16Output

            settings = Gaussian16Output(filename=filepath).read_settings()
        elif suffix == ".inp":
            from chemsmart.io.orca.input import ORCAInput

            settings = ORCAInput(filename=filepath).read_settings()
        elif suffix == ".out":
            from chemsmart.utils.io import get_program_type_from_file

            program = get_program_type_from_file(filepath)
            if program == "orca":
                from chemsmart.io.orca.output import ORCAOutput

                settings = ORCAOutput(filename=filepath).read_settings()
            elif program == "gaussian":
                from chemsmart.io.gaussian.output import Gaussian16Output

                settings = Gaussian16Output(filename=filepath).read_settings()
            elif program == "xtb":
                from chemsmart.io.xtb.output import XTBOutput

                output = XTBOutput(folder=os.path.dirname(filepath))
                settings = {
                    "charge": output.charge,
                    "multiplicity": output.multiplicity,
                }
            else:
                raise ValueError(
                    f"Cannot derive xTB settings from unknown .out file: {filepath}"
                )
        else:
            return cls.default()

        return cls.default().merge(
            settings,
            keywords=("charge", "multiplicity"),
        )

    def __eq__(self, other):
        if type(self) is not type(other):
            return NotImplemented
        return vars(self) == vars(other)
