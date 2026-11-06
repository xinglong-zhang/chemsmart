"""Utility functions and row abstractions for tabular dataset processing."""

from __future__ import annotations

import logging
import os
import re

import numpy as np

logger = logging.getLogger(__name__)


def normalize_table_cell(value):
    """Convert blank/NaN table cells to None for stable downstream handling."""
    if value is None:
        return None
    if isinstance(value, float) and np.isnan(value):
        return None
    try:
        import pandas as pd

        if pd.isna(value):
            return None
    except (TypeError, ImportError):
        pass
    if isinstance(value, str):
        stripped = value.strip()
        return stripped if stripped else None
    return value


class PKaTableEntry:
    """Generic table-row abstraction for pKa job-submission tables.

    This refactor keeps backward compatibility for dict-like access while
    avoiding dynamic attribute lookup.
    """

    _ALIASES = {
        "filepath": ["filepath", "file_path", "path", "ha_file"],
        "proton_index": ["proton_index", "pi"],
        "charge": ["charge", "q"],
        "multiplicity": ["multiplicity", "mult", "m"],
    }

    def __init__(self, *args, row_number: int = None, **kwargs):
        # Explicitly initialized attributes
        self.row_number = row_number
        self.filepath = None
        self.proton_index = None
        self.charge = None
        self.multiplicity = None

        # Canonical + non-canonical storage (no dynamic attrs on self)
        self._data = {}
        self._extra_data = {}

        data = {}
        if len(args) == 4:
            data.update(
                {
                    "filepath": args[0],
                    "proton_index": args[1],
                    "charge": args[2],
                    "multiplicity": args[3],
                }
            )
        elif len(args) == 1 and isinstance(args[0], dict):
            data.update(args[0])
        elif len(args) != 0:
            raise TypeError(
                "PKaTableEntry accepts either (filepath, proton_index, charge, multiplicity), "
                "or a single dict, plus keyword fields."
            )

        data.update(kwargs)
        for key, value in data.items():
            if key == "row_number":
                continue
            self._set_field(key, value)

    @staticmethod
    def normalize_header(header):
        """Normalize a table header into snake_case."""
        header = str(header).strip().lower()
        header = re.sub(r"[^0-9a-zA-Z]+", "_", header)
        header = re.sub(r"_+", "_", header).strip("_")
        return header

    @staticmethod
    def resolve_column(columns, candidates, required=True):
        """Resolve a physical column name from logical candidates."""
        normalized = {PKaTableEntry.normalize_header(c): c for c in columns}
        for candidate in candidates:
            key = PKaTableEntry.normalize_header(candidate)
            if key in normalized:
                return normalized[key]
        if required:
            raise ValueError(
                "Could not resolve required column. Tried: "
                + ", ".join(candidates)
            )
        return None

    @staticmethod
    def parse_table(table_path, delimiter=None, comment="#"):
        """Backward-compatible shim for generic table parsing."""
        from chemsmart.io.datasets import TabularDataset

        return TabularDataset.parse_table(
            table_path=table_path,
            delimiter=delimiter,
            comment=comment,
        )

    @classmethod
    def from_headers_and_row(cls, headers, row, row_number: int = None):
        if len(headers) != len(row):
            raise ValueError(
                f"Header/value length mismatch: {len(headers)} headers vs {len(row)} values"
            )
        row_dict = {
            str(h).strip(): v for h, v in zip(headers, row, strict=False)
        }
        return cls(row_dict, row_number=row_number)

    def _canonical_key(self, key):
        k = str(key)
        nk = self.normalize_header(k)
        for canonical, aliases in self._ALIASES.items():
            if nk == self.normalize_header(canonical):
                return canonical
            for alias in aliases:
                if nk == self.normalize_header(alias):
                    return canonical
        return None

    def _set_field(self, key, value):
        value = normalize_table_cell(value)
        canonical = self._canonical_key(key)
        if canonical == "filepath":
            self.filepath = value
            self._data["filepath"] = value
        elif canonical == "proton_index":
            self.proton_index = value
            self._data["proton_index"] = value
        elif canonical == "charge":
            self.charge = value
            self._data["charge"] = value
        elif canonical == "multiplicity":
            self.multiplicity = value
            self._data["multiplicity"] = value
        else:
            key_str = str(key)
            self._extra_data[key_str] = value
            self._data[key_str] = value

    def __repr__(self):
        return (
            f"PKaTableEntry(row_number={self.row_number}, data={self._data!r})"
        )

    def __getitem__(self, key):
        if key in self._data:
            return self._data[key]
        canonical = self._canonical_key(key)
        if canonical is not None and canonical in self._data:
            return self._data[canonical]
        raise KeyError(key)

    def __setitem__(self, key, value):
        self._set_field(key, value)

    def __contains__(self, key):
        if key in self._data:
            return True
        canonical = self._canonical_key(key)
        if canonical is None:
            return False
        return canonical in self._data

    def get(self, key, default=None):
        if key in self._data:
            return self._data.get(key, default)
        canonical = self._canonical_key(key)
        if canonical is None:
            return default
        return self._data.get(canonical, default)

    def get_canonical(self, canonical, default=None):
        if canonical == "filepath":
            return self.filepath if self.filepath is not None else default
        if canonical == "proton_index":
            return (
                self.proton_index if self.proton_index is not None else default
            )
        if canonical == "charge":
            return self.charge if self.charge is not None else default
        if canonical == "multiplicity":
            return (
                self.multiplicity if self.multiplicity is not None else default
            )
        return self.get(canonical, default)

    def keys(self):
        return self._data.keys()

    def items(self):
        return self._data.items()

    def values(self):
        return self._data.values()

    def to_dict(self):
        out = dict(self._extra_data)
        out["filepath"] = self.filepath
        out["proton_index"] = self.proton_index
        out["charge"] = self.charge
        out["multiplicity"] = self.multiplicity
        return out

    def to_kwargs(self, rename_map=None, drop_none=False):
        out = self.to_dict()
        if rename_map:
            out = {rename_map.get(k, k): v for k, v in out.items()}
        if drop_none:
            out = {k: v for k, v in out.items() if v is not None}
        return out

    @staticmethod
    def _is_cdxml_filepath(filepath):
        return bool(filepath) and str(filepath).lower().endswith(
            (".cdx", ".cdxml")
        )

    def validate(self):
        errors = []
        row_info = f" (row {self.row_number})" if self.row_number else ""

        filepath = self.filepath
        proton_index = self.proton_index
        charge = self.charge
        multiplicity = self.multiplicity
        is_cdxml = self._is_cdxml_filepath(filepath)

        if not filepath:
            errors.append(f"Empty filepath{row_info}")
        elif not os.path.exists(str(filepath)):
            errors.append(f"File not found: {filepath}{row_info}")

        if proton_index is None:
            if not is_cdxml:
                errors.append(f"Missing proton_index{row_info}")
        else:
            try:
                proton_index = int(proton_index)
                if proton_index < 1:
                    errors.append(
                        f"proton_index must be >= 1, got {proton_index}{row_info}"
                    )
            except (TypeError, ValueError):
                errors.append(
                    f"Invalid proton_index: {proton_index!r}{row_info}"
                )

        if charge is None:
            errors.append(f"Missing charge{row_info}")
        else:
            try:
                int(charge)
            except (TypeError, ValueError):
                errors.append(f"Invalid charge: {charge!r}{row_info}")

        if multiplicity is None:
            errors.append(f"Missing multiplicity{row_info}")
        else:
            try:
                multiplicity = int(multiplicity)
                if multiplicity < 1:
                    errors.append(
                        f"multiplicity must be >= 1, got {multiplicity}{row_info}"
                    )
            except (TypeError, ValueError):
                errors.append(
                    f"Invalid multiplicity: {multiplicity!r}{row_info}"
                )

        if errors:
            raise ValueError("; ".join(errors))

    @staticmethod
    def is_submission_table(table_path) -> bool:
        """Return True when *table_path* has pKa submission-table columns."""
        from chemsmart.io.datasets import TabularDataset

        if not table_path:
            return False
        if str(table_path).lower().endswith((".cdx", ".cdxml")):
            return False
        try:
            dataset = TabularDataset.parse_table(
                table_path=table_path,
                comment="#",
            )
            for aliases in PKaTableEntry._ALIASES.values():
                TabularDataset.resolve_column(dataset.columns, aliases)
            return True
        except (ValueError, FileNotFoundError, OSError):
            return False

    @staticmethod
    def parse_pka_table(
        table_path: str,
        delimiter: str = None,
        skip_header: bool = True,
    ) -> list:
        """Thin pKa adapter on top of the generic tabular parser layer."""
        from chemsmart.io.datasets import TabularDataset

        dataset = TabularDataset.parse_table(
            table_path=table_path,
            delimiter=delimiter,
            comment="#",
        )

        try:
            file_col = TabularDataset.resolve_column(
                dataset.columns,
                PKaTableEntry._ALIASES["filepath"],
            )
            proton_col = TabularDataset.resolve_column(
                dataset.columns,
                PKaTableEntry._ALIASES["proton_index"],
            )
            charge_col = TabularDataset.resolve_column(
                dataset.columns,
                PKaTableEntry._ALIASES["charge"],
            )
            mult_col = TabularDataset.resolve_column(
                dataset.columns,
                PKaTableEntry._ALIASES["multiplicity"],
            )
        except ValueError:
            raise ValueError(
                "Invalid table format: expected 4 columns "
                "(filepath, proton_index, charge, multiplicity)."
            )

        # Canonicalize selected columns for generic validation/entry materialization.
        canonical_df = dataset.dataframe.rename(
            columns={
                file_col: "filepath",
                proton_col: "proton_index",
                charge_col: "charge",
                mult_col: "multiplicity",
            }
        )[["filepath", "proton_index", "charge", "multiplicity"]]

        canonical_dataset = TabularDataset(
            canonical_df, source_path=table_path
        )
        canonical_dataset.validate(
            required_columns=[
                "filepath",
                "proton_index",
                "charge",
                "multiplicity",
            ],
            integer_columns=[],
            positive_integer_columns=[],
            path_columns=[],
            check_file_exists=False,
        )

        entries = canonical_dataset.to_entries(
            entry_cls=PKaTableEntry, row_offset=2
        )

        # Coerce numeric fields to preserve historical parse_pka_table behavior.
        for entry in entries:
            line_num = entry.row_number if entry.row_number is not None else 0
            is_cdxml = PKaTableEntry._is_cdxml_filepath(entry.filepath)

            if entry.proton_index is None:
                if not is_cdxml:
                    raise ValueError(
                        f"Missing proton_index at line {line_num} for "
                        f"{entry.filepath!r}"
                    )
            else:
                try:
                    entry["proton_index"] = int(entry["proton_index"])
                except (TypeError, ValueError):
                    raise ValueError(
                        f"Invalid proton_index at line {line_num}: "
                        f"{entry['proton_index']!r} is not an integer"
                    )

            try:
                entry["charge"] = int(entry["charge"])
            except (TypeError, ValueError):
                raise ValueError(
                    f"Invalid charge at line {line_num}: "
                    f"{entry['charge']!r} is not an integer"
                )

            try:
                entry["multiplicity"] = int(entry["multiplicity"])
            except (TypeError, ValueError):
                raise ValueError(
                    f"Invalid multiplicity at line {line_num}: "
                    f"{entry['multiplicity']!r} is not an integer"
                )

        if skip_header and len(entries) >= 0:
            pass

        if not entries:
            raise ValueError(
                f"No valid entries found in pKa table: {table_path}"
            )

        return entries


class PKaOutputTableEntry:
    """Row abstraction for a pKa output table.

    Uses explicit attributes (initialized to None) instead of dynamic
    attribute lookup.
    """

    _ALIASES = {
        "basename": ["basename", "name", "label", "system"],
        "ha_gas": [
            "ha_gas",
            "ha_opt",
            "ha_gas_file",
            "ha_optimization_output",
        ],
        "a_gas": ["a_gas", "a_opt", "a_gas_file", "a_optimization_output"],
        "href_gas": [
            "href_gas",
            "href_opt",
            "href_gas_file",
            "href_optimization_output",
        ],
        "ref_gas": [
            "ref_gas",
            "ref_opt",
            "ref_gas_file",
            "ref_optimization_output",
        ],
        "ha_sp": [
            "ha_sp",
            "ha_solv",
            "ha_solv_file",
            "ha_single_point_output",
        ],
        "a_sp": ["a_sp", "a_solv", "a_solv_file", "a_single_point_output"],
        "href_sp": [
            "href_sp",
            "href_solv",
            "href_solv_file",
            "href_single_point_output",
        ],
        "ref_sp": [
            "ref_sp",
            "ref_solv",
            "ref_solv_file",
            "ref_single_point_output",
        ],
        "pka_ref": ["pka_ref", "pka_reference", "reference_pka", "ref_pka"],
    }

    _REFERENCE_COLUMNS = (
        "href_gas",
        "ref_gas",
        "href_sp",
        "ref_sp",
        "pka_ref",
    )

    _OUTPUT_PATH_FIELDS = (
        "ha_gas",
        "a_gas",
        "ha_sp",
        "a_sp",
        "href_gas",
        "ref_gas",
        "href_sp",
        "ref_sp",
    )

    # Filename suffix patterns for CHEMSMART pKa job outputs (not CSV aliases).
    _OUTPUT_SUFFIX_CANDIDATES = {
        "ha_gas": ["_pka_HA_opt", "_pka_HA", "_pka"],
        "a_gas": ["_pka_A_opt", "_pka_A", "_pka_cb"],
        "ha_sp": ["_pka_HA_sp", "_pka_sp"],
        "a_sp": ["_pka_A_sp", "_pka_cb_sp"],
        "href_gas": ["_pka_HRef_opt", "_HRef_opt", "_pka_HRef", "_HRef"],
        "ref_gas": [
            "_pka_Ref_opt",
            "_Ref_opt",
            "_pka_Ref",
            "_Ref",
            "_pka_cb",
            "_cb",
        ],
        "href_sp": ["_pka_HRef_sp", "_HRef_sp"],
        "ref_sp": ["_pka_Ref_sp", "_Ref_sp", "_pka_cb_sp", "_cb_sp"],
    }

    _TARGET_AUTO_DISCOVER_FIELDS = ("ha_gas", "a_gas", "ha_sp", "a_sp")

    TARGET_SUFFIX_HELP = (
        "  <basename>_pka_A_opt.<ext>   (conjugate base gas-phase)\n"
        "  <basename>_pka_HA_sp.<ext>   (HA solvent single-point)\n"
        "  <basename>_pka_A_sp.<ext>    (conjugate base solvent SP)"
    )

    REFERENCE_SUFFIX_HELP = (
        "  <basename>_pka_Ref_opt.<ext>  (reference conjugate base)\n"
        "  <basename>_pka_HRef_sp.<ext>  (reference acid solvent SP)\n"
        "  <basename>_pka_Ref_sp.<ext>   (reference conjugate base solvent SP)"
    )

    def __init__(self, data: dict, row_number: int = None):
        # Explicitly initialized attributes
        self.row_number = row_number
        self.basename = None
        self.ha_gas = None
        self.a_gas = None
        self.href_gas = None
        self.ref_gas = None
        self.ha_sp = None
        self.a_sp = None
        self.href_sp = None
        self.ref_sp = None
        self.pka_ref = None

        # Explicit helper attributes requested for pKa workflow clarity
        self.ha_basename = None
        self.hb_basename = None
        self.ha_conjugated_base_label = None
        self.hb_conjugated_base_label = None
        self.opt_path = None
        self.sp_path = None

        self._data = {}
        self._extra_data = {}

        for key, value in data.items():
            if key == "row_number":
                continue
            self._set_field(key, value)

        self._derive_helper_fields()

    @staticmethod
    def _normalize_header(header):
        header = str(header).strip().lower()
        header = re.sub(r"[^0-9a-zA-Z]+", "_", header)
        header = re.sub(r"_+", "_", header).strip("_")
        return header

    def _canonical_key(self, key):
        nk = self._normalize_header(key)
        for canonical, aliases in self._ALIASES.items():
            if nk == self._normalize_header(canonical):
                return canonical
            for alias in aliases:
                if nk == self._normalize_header(alias):
                    return canonical
        return None

    def _set_field(self, key, value):
        value = normalize_table_cell(value)
        canonical = self._canonical_key(key)
        if canonical == "basename":
            self.basename = value
            self._data["basename"] = value
        elif canonical == "ha_gas":
            self.ha_gas = value
            self._data["ha_gas"] = value
        elif canonical == "a_gas":
            self.a_gas = value
            self._data["a_gas"] = value
        elif canonical == "href_gas":
            self.href_gas = value
            self._data["href_gas"] = value
        elif canonical == "ref_gas":
            self.ref_gas = value
            self._data["ref_gas"] = value
        elif canonical == "ha_sp":
            self.ha_sp = value
            self._data["ha_sp"] = value
        elif canonical == "a_sp":
            self.a_sp = value
            self._data["a_sp"] = value
        elif canonical == "href_sp":
            self.href_sp = value
            self._data["href_sp"] = value
        elif canonical == "ref_sp":
            self.ref_sp = value
            self._data["ref_sp"] = value
        elif canonical == "pka_ref":
            self.pka_ref = value
            self._data["pka_ref"] = value
        else:
            key_str = str(key)
            self._extra_data[key_str] = value
            self._data[key_str] = value

    def _derive_helper_fields(self):
        # HA basename defaults to basename column
        if self.basename is not None:
            self.ha_basename = str(self.basename)

        # HB basename derived from href_gas file stem when available
        if self.href_gas is not None:
            href_name = os.path.basename(str(self.href_gas))
            self.hb_basename = os.path.splitext(href_name)[0]

        # Conjugated base labels are made explicit
        if self.ha_basename is not None:
            self.ha_conjugated_base_label = (
                f"{self.ha_basename}_conjugated_base"
            )
        if self.hb_basename is not None:
            self.hb_conjugated_base_label = (
                f"{self.hb_basename}_conjugated_base"
            )

        # Convenience paths for target acid outputs
        self.opt_path = self.ha_gas
        self.sp_path = self.ha_sp

    def __getitem__(self, key):
        if key in self._data:
            return self._data[key]
        canonical = self._canonical_key(key)
        if canonical is not None and canonical in self._data:
            return self._data[canonical]
        raise KeyError(key)

    def __setitem__(self, key, value):
        self._set_field(key, value)
        self._derive_helper_fields()

    def __contains__(self, key):
        if key in self._data:
            return True
        canonical = self._canonical_key(key)
        if canonical is None:
            return False
        return canonical in self._data

    def __repr__(self):
        return (
            f"PKaOutputTableEntry(row={self.row_number}, "
            f"basename={self.basename if self.basename is not None else '?'})"
        )

    def _existing_output_paths(self):
        """Return existing output-file paths already listed in this row."""
        paths = []
        for field in self._OUTPUT_PATH_FIELDS:
            val = normalize_table_cell(self.get(field))
            if val and os.path.isfile(str(val)):
                paths.append(str(val))
        return paths

    def _basename_output_paths(self, suffix_candidates):
        """Return existing basename-pattern output files for this row."""
        paths = []
        for suffixes in suffix_candidates.values():
            for suffix in suffixes:
                for ext in (".log", ".out"):
                    candidate = f"{self.basename}{suffix}{ext}"
                    if os.path.isfile(candidate):
                        paths.append(candidate)
        return paths

    def _detect_output_program(self, suffix_candidates):
        """Detect Gaussian/ORCA from basename outputs or explicit row paths."""
        from chemsmart.utils.io import get_program_type_from_file

        allowed_programs = {"gaussian", "orca"}

        for path_group in (
            self._basename_output_paths(suffix_candidates),
            self._existing_output_paths(),
        ):
            if not path_group:
                continue

            programs = set()
            for fp in path_group:
                detected = get_program_type_from_file(fp)
                if detected != "unknown":
                    programs.add(detected)

            if len(programs) != 1:
                continue

            program = next(iter(programs))
            if program not in allowed_programs:
                continue
            return program

        return None

    def _resolve_filenames(self):
        """Auto-discover output files based on basename if not explicitly provided."""
        from chemsmart.utils.io import get_program_output_extensions

        suffix_candidates = {
            key: self._OUTPUT_SUFFIX_CANDIDATES[key]
            for key in self._TARGET_AUTO_DISCOVER_FIELDS
        }

        program = self._detect_output_program(suffix_candidates)
        extensions = get_program_output_extensions(program)
        default_ext = extensions[0]

        for field, suffixes in suffix_candidates.items():
            if normalize_table_cell(self.get(field)) is not None:
                continue

            found = False
            for suffix in suffixes:
                for ext in extensions:
                    candidate = f"{self.basename}{suffix}{ext}"
                    if os.path.exists(candidate):
                        self._set_field(field, candidate)
                        found = True
                        break
                if found:
                    break

            if not found:
                self._set_field(
                    field, f"{self.basename}{suffixes[0]}{default_ext}"
                )

        self._derive_helper_fields()

    def get(self, key, default=None):
        if key in self._data:
            return self._data.get(key, default)
        canonical = self._canonical_key(key)
        if canonical is None:
            return default
        return self._data.get(canonical, default)

    def keys(self):
        return self._data.keys()

    def items(self):
        return self._data.items()

    def to_dict(self):
        out = dict(self._extra_data)
        out.update(
            {
                "basename": self.basename,
                "ha_gas": self.ha_gas,
                "a_gas": self.a_gas,
                "href_gas": self.href_gas,
                "ref_gas": self.ref_gas,
                "ha_sp": self.ha_sp,
                "a_sp": self.a_sp,
                "href_sp": self.href_sp,
                "ref_sp": self.ref_sp,
                "pka_ref": self.pka_ref,
            }
        )
        return out

    def validate(self, check_file_exists=True, scheme="proton exchange"):
        """Validate that all required file paths are present and non-empty."""
        errors = []
        row_info = f" (row {self.row_number})" if self.row_number else ""

        if self.basename is None or self.basename.strip() == "":
            errors.append(f"Missing basename{row_info}")

        # Auto-discover missing files if basename is present
        if self.basename:
            self._resolve_filenames()

        if scheme == "direct":
            required_files = [
                ("ha_gas", self.ha_gas),
                ("a_gas", self.a_gas),
                ("ha_sp", self.ha_sp),
                ("a_sp", self.a_sp),
            ]
        elif scheme == "proton exchange":
            required_files = [
                ("ha_gas", self.ha_gas),
                ("a_gas", self.a_gas),
                ("href_gas", self.href_gas),
                ("ref_gas", self.ref_gas),
                ("ha_sp", self.ha_sp),
                ("a_sp", self.a_sp),
                ("href_sp", self.href_sp),
                ("ref_sp", self.ref_sp),
            ]
        else:
            raise ValueError(f"Unsupported pKa analysis scheme: {scheme!r}")
        for col, val in required_files:
            if val is None or (isinstance(val, float) and np.isnan(val)):
                errors.append(f"Missing {col}{row_info}")
            elif check_file_exists and not os.path.exists(str(val)):
                errors.append(f"File not found for {col}: {val}{row_info}")

        if scheme == "proton exchange":
            if self.pka_ref is None or (
                isinstance(self.pka_ref, float) and np.isnan(self.pka_ref)
            ):
                errors.append(f"Missing pka_ref{row_info}")
            else:
                try:
                    float(self.pka_ref)
                except (TypeError, ValueError):
                    errors.append(
                        f"Invalid pka_ref: {self.pka_ref!r}{row_info}"
                    )

        if errors:
            raise ValueError("; ".join(errors))


class PKaOutputTable:
    """Parsed pKa output table with reference resolution and pKa execution.

    The table object keeps the parsed rows together so callers can pass it
    between workflow classes instead of repeatedly passing raw table paths and
    entry lists.
    """

    def __init__(self, entries: list, source_path: str = None):
        self.entries = list(entries)
        self.source_path = source_path
        self.results = None

    @staticmethod
    def parse_pka_output_table(table_path: str, delimiter: str = None) -> list:
        """Parse an output-table file into :class:`PKaOutputTableEntry` rows."""
        import pandas as pd

        from chemsmart.io.datasets import TabularDataset

        dataset = TabularDataset.parse_table(
            table_path=table_path,
            delimiter=delimiter,
            comment="#",
        )

        rename_map = {}
        for canonical, aliases in PKaOutputTableEntry._ALIASES.items():
            resolved = TabularDataset.resolve_column(
                dataset.columns,
                aliases,
                required=False,
            )
            if resolved is not None:
                rename_map[resolved] = canonical

        if "basename" not in rename_map.values():
            raise ValueError(
                "Output table is missing required column: basename. "
                f"Accepted aliases: {PKaOutputTableEntry._ALIASES['basename']}"
            )

        canonical_df = dataset.dataframe.rename(columns=rename_map)
        canonical_cols = [
            c
            for c in PKaOutputTableEntry._ALIASES
            if c in canonical_df.columns
        ]
        canonical_df = canonical_df[canonical_cols]
        canonical_df = canonical_df.where(pd.notnull(canonical_df), None)

        canonical_dataset = TabularDataset(
            canonical_df, source_path=table_path
        )
        entries = canonical_dataset.to_entries(
            entry_cls=PKaOutputTableEntry,
            row_offset=2,
        )

        if not entries:
            raise ValueError(
                f"No valid entries found in pKa output table: {table_path}"
            )

        return entries

    @classmethod
    def from_file(cls, table_path: str, delimiter: str = None):
        """Parse an output table file into a ``PKaOutputTable`` object."""
        entries = cls.parse_pka_output_table(table_path, delimiter=delimiter)
        return cls(entries=entries, source_path=table_path)

    @staticmethod
    def resolve_pka_output_references(entries: list) -> list:
        """Fill blank reference-acid columns by carrying forward earlier rows."""
        ref_cols = PKaOutputTableEntry._REFERENCE_COLUMNS
        last_ref = {}

        for entry in entries:
            for col in ref_cols:
                val = entry.get(col)
                if val is not None and not (
                    isinstance(val, float) and np.isnan(val)
                ):
                    last_ref[col] = val
                else:
                    if col not in last_ref:
                        raise ValueError(
                            f"Reference column '{col}' is blank in row "
                            f"{entry.row_number} and no previous value exists "
                            f"to carry forward."
                        )
                    entry[col] = last_ref[col]

        return entries

    @staticmethod
    def compute_pka_from_output_table(
        entries: list,
        output_cls,
        temperature: float = 298.15,
        concentration: float = 1.0,
        pressure: float = 1.0,
        cutoff_entropy_grimme: float = 100.0,
        cutoff_enthalpy: float = 100.0,
        entropy_method: str = "grimme",
        scheme: str = "proton exchange",
        delta_G_proton: float = None,
    ) -> list:
        """Compute pKa for every row using the supplied output class."""
        if scheme == "direct" and delta_G_proton is None:
            raise ValueError(
                "delta_G_proton is required when scheme='direct'."
            )

        results = []
        for entry in entries:
            pka_kwargs = dict(
                ha_gas_file=entry["ha_gas"],
                a_gas_file=entry["a_gas"],
                ha_solv_file=entry["ha_sp"],
                a_solv_file=entry["a_sp"],
                temperature=temperature,
                concentration=concentration,
                pressure=pressure,
                cutoff_entropy_grimme=cutoff_entropy_grimme,
                cutoff_enthalpy=cutoff_enthalpy,
                entropy_method=entropy_method,
                scheme=scheme,
            )
            if scheme == "direct":
                pka_kwargs["delta_G_proton"] = delta_G_proton
            elif scheme == "proton exchange":
                pka_kwargs.update(
                    href_gas_file=entry["href_gas"],
                    ref_gas_file=entry["ref_gas"],
                    href_solv_file=entry["href_sp"],
                    ref_solv_file=entry["ref_sp"],
                    pka_reference=float(entry["pka_ref"]),
                )
            else:
                raise ValueError(
                    f"Unsupported pKa analysis scheme: {scheme!r}"
                )
            compute_fn = (
                output_cls.compute_pka
                if isinstance(output_cls, type)
                else output_cls
            )
            pka_result = compute_fn(**pka_kwargs)
            pka_result["basename"] = entry["basename"]
            results.append(pka_result)
        return results

    @staticmethod
    def pka_scheme_delta_g_key(scheme):
        """Return the result-dict key for a scheme-specific ΔG value."""
        keys = {
            "direct": "delta_G_diss_kcal_mol",
            "proton exchange": "delta_G_soln_kcal_mol",
        }
        if scheme is None:
            return None
        return keys.get(scheme)

    @staticmethod
    def pka_scheme_delta_g_value(result, scheme=None):
        """Return the ΔG value appropriate for the analysis scheme."""
        resolved_scheme = scheme or result.get("scheme")
        key = PKaOutputTable.pka_scheme_delta_g_key(resolved_scheme)
        if key is not None:
            return result[key]
        return result.get(
            "delta_G_diss_kcal_mol", result["delta_G_soln_kcal_mol"]
        )

    @staticmethod
    def format_pka_batch_results_table(
        entries,
        results,
        temperature,
        pressure,
        scheme=None,
    ):
        """Return the formatted batch pKa summary table shown on stdout."""
        display_scheme = scheme
        if display_scheme is None and results:
            display_scheme = results[0].get("scheme")

        header = PKaOutputTable._scheme_batch_header(display_scheme)
        dg_label = PKaOutputTable._scheme_delta_g_label(display_scheme)

        lines = [
            "=" * 78,
            header,
            "=" * 78,
            f"Temperature: {temperature} K",
            f"Pressure: {pressure} atm",
            f"{'basename':<30} {'pKa':>10} {dg_label:>20}",
            "-" * 78,
        ]

        for entry, result in zip(entries, results):
            dg_value = PKaOutputTable.pka_scheme_delta_g_value(result, scheme)
            lines.append(
                f"{entry['basename']:<30} "
                f"{result['pKa']:>10.2f} "
                f"{dg_value:>20.4f}"
            )
        lines.append("=" * 78)
        return "\n".join(lines)

    def echo_pka_output_table_results(
        self,
        results,
        output_results,
        temperature,
        pressure,
        scheme=None,
    ):
        table_text = self.format_pka_batch_results_table(
            self.entries,
            results,
            temperature,
            pressure,
            scheme=scheme,
        )
        if output_results is not None:
            self.export_results(
                output_results,
                results=results,
                scheme=scheme,
                temperature=temperature,
                pressure=pressure,
            )
        return table_text

    @staticmethod
    def _scheme_batch_header(scheme):
        if scheme is None:
            return "Batch pKa Results"
        headers = {
            "direct": "Batch pKa Results (Direct Dissociation)",
            "proton exchange": "Batch pKa Results (Dual-level Proton Exchange)",
        }
        return headers.get(scheme, f"Batch pKa Results ({scheme})")

    @staticmethod
    def _scheme_delta_g_label(scheme):
        labels = {
            "direct": "ΔG_diss (kcal/mol)",
            "proton exchange": "ΔG_soln (kcal/mol)",
        }
        return labels.get(scheme, "ΔG (kcal/mol)")

    @staticmethod
    def export_pka_results_table(
        entries: list,
        results: list,
        output_path: str,
        scheme: str = None,
        temperature: float = 298.15,
        pressure: float = 1.0,
    ) -> None:
        """Write the formatted batch pKa summary table to *output_path*."""
        table_text = PKaOutputTable.format_pka_batch_results_table(
            entries,
            results,
            temperature,
            pressure,
            scheme=scheme,
        )
        with open(output_path, "w", encoding="utf-8") as fh:
            fh.write(table_text + "\n")

        logger.info(f"pKa results written to {output_path}")

    def __len__(self):
        return len(self.entries)

    def __iter__(self):
        return iter(self.entries)

    def resolve_references(self):
        """Fill blank reference-acid columns by carrying previous values."""
        self.resolve_pka_output_references(self.entries)
        return self

    def validate(
        self, check_file_exists: bool = True, scheme: str = "proton exchange"
    ):
        """Validate all rows in the table."""
        all_errors = []
        for entry in self.entries:
            try:
                entry.validate(
                    check_file_exists=check_file_exists, scheme=scheme
                )
            except ValueError as exc:
                all_errors.append(str(exc))
        if all_errors:
            raise ValueError(
                "Output table validation failed:\n" + "\n".join(all_errors)
            )
        return self

    def prepare(
        self, check_file_exists: bool = True, scheme: str = "proton exchange"
    ):
        """Resolve shared references and validate the table."""
        if scheme == "direct":
            return self.validate(
                check_file_exists=check_file_exists, scheme=scheme
            )
        return self.resolve_references().validate(
            check_file_exists=check_file_exists, scheme=scheme
        )

    def run_pka(
        self,
        output_cls,
        temperature: float = 298.15,
        concentration: float = 1.0,
        pressure: float = 1.0,
        cutoff_entropy_grimme: float = 100.0,
        cutoff_enthalpy: float = 100.0,
        entropy_method: str = "grimme",
        scheme: str = "proton exchange",
        delta_G_proton: float = None,
    ) -> list:
        """Compute pKa for every row using the supplied output class."""
        self.results = self.compute_pka_from_output_table(
            entries=self.entries,
            output_cls=output_cls,
            temperature=temperature,
            concentration=concentration,
            pressure=pressure,
            cutoff_entropy_grimme=cutoff_entropy_grimme,
            cutoff_enthalpy=cutoff_enthalpy,
            entropy_method=entropy_method,
            scheme=scheme,
            delta_G_proton=delta_G_proton,
        )
        return self.results

    def export_results(
        self,
        output_path: str,
        results: list = None,
        scheme: str = None,
        temperature: float = 298.15,
        pressure: float = 1.0,
    ) -> None:
        """Export the formatted batch pKa summary table for this table."""
        self.export_pka_results_table(
            self.entries,
            self.results if results is None else results,
            output_path,
            scheme=scheme,
            temperature=temperature,
            pressure=pressure,
        )

    @staticmethod
    def validate_pka_table_entries(
        entries: list,
        check_file_exists: bool = True,
    ) -> list:
        """Validate a list of PKaTableEntry objects.

        Args:
            entries: List of PKaTableEntry to validate.
            check_file_exists: If True, verify that each filepath exists.

        Returns:
            list[PKaTableEntry]: The validated entries (same list, for chaining).

        Raises:
            ValueError: If any entry fails validation.
        """
        all_errors = []

        for entry in entries:
            try:
                if check_file_exists:
                    entry.validate()
                else:
                    # Validate without file existence check
                    is_cdxml = PKaTableEntry._is_cdxml_filepath(entry.filepath)
                    if entry.proton_index is None:
                        if not is_cdxml:
                            raise ValueError(
                                f"Missing proton_index for {entry.filepath}"
                            )
                    elif entry.proton_index < 1:
                        raise ValueError(
                            f"Invalid proton_index: {entry.proton_index}"
                        )
                    if entry.charge is None:
                        raise ValueError("Missing charge")
                    if entry.multiplicity is None or entry.multiplicity < 1:
                        raise ValueError(
                            f"Invalid multiplicity: {entry.multiplicity}"
                        )
            except ValueError as e:
                all_errors.append(str(e))

        if all_errors:
            raise ValueError(
                "pKa table validation failed:\n" + "\n".join(all_errors)
            )

        return entries


def pka_output_basename_from_path(filepath, role):
    """Strip a known gas-phase suffix to recover the pKa job basename."""
    stem = os.path.splitext(os.path.basename(str(filepath)))[0]
    for suffix in PKaOutputTableEntry._OUTPUT_SUFFIX_CANDIDATES.get(role, []):
        if stem.endswith(suffix):
            return stem[: -len(suffix)]
    return stem


def discover_pka_output_path(
    basename,
    directory,
    role,
    program=None,
    filepath_hint=None,
):
    """Return the first existing companion output path for *role*."""
    from chemsmart.utils.io import (
        get_program_output_extensions,
        get_program_type_from_file,
    )

    if program is None and filepath_hint is not None:
        program = get_program_type_from_file(filepath_hint)
    extensions = get_program_output_extensions(program)
    suffixes = PKaOutputTableEntry._OUTPUT_SUFFIX_CANDIDATES[role]
    directory = directory or "."
    for suffix in suffixes:
        for ext in extensions:
            candidate = os.path.join(directory, f"{basename}{suffix}{ext}")
            if os.path.isfile(candidate):
                return candidate
    return os.path.join(directory, f"{basename}{suffixes[0]}{extensions[0]}")


def discover_pka_reference_companion_outputs(href_gas_path, program=None):
    """Infer Ref- and reference solvent SP paths from an HRef gas-phase file."""
    from chemsmart.utils.io import get_program_type_from_file

    href_gas_path = str(href_gas_path)
    directory = os.path.dirname(href_gas_path) or "."
    if program is None:
        program = get_program_type_from_file(href_gas_path)
    basename = pka_output_basename_from_path(href_gas_path, "href_gas")
    return {
        "ref": discover_pka_output_path(
            basename, directory, "ref_gas", program=program
        ),
        "href_solv": discover_pka_output_path(
            basename, directory, "href_sp", program=program
        ),
        "ref_solv": discover_pka_output_path(
            basename, directory, "ref_sp", program=program
        ),
    }


PKA_OUTPUT_SUFFIX_CANDIDATES = PKaOutputTableEntry._OUTPUT_SUFFIX_CANDIDATES
PKA_TARGET_SUFFIX_HELP = PKaOutputTableEntry.TARGET_SUFFIX_HELP
PKA_REFERENCE_SUFFIX_HELP = PKaOutputTableEntry.REFERENCE_SUFFIX_HELP

parse_pka_output_table = PKaOutputTable.parse_pka_output_table
resolve_pka_output_references = PKaOutputTable.resolve_pka_output_references
compute_pka_from_output_table = PKaOutputTable.compute_pka_from_output_table
pka_scheme_delta_g_key = PKaOutputTable.pka_scheme_delta_g_key
pka_scheme_delta_g_value = PKaOutputTable.pka_scheme_delta_g_value
export_pka_results_table = PKaOutputTable.export_pka_results_table
