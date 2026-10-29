"""Tabular dataset objects for high-throughput workflow I/O."""

from __future__ import annotations

import logging
import os
import re

logger = logging.getLogger(__name__)


class TabularDataset:
    """Generic DataFrame-backed dataset for high-throughput workflows."""

    def __init__(self, dataframe, source_path=None):
        self.dataframe = dataframe.reset_index(drop=True)
        self.source_path = source_path

    @property
    def columns(self):
        return list(self.dataframe.columns)

    def __len__(self):
        return len(self.dataframe)

    @staticmethod
    def normalize_header(header):
        header = str(header).strip().lower()
        header = re.sub(r"[^0-9a-zA-Z]+", "_", header)
        header = re.sub(r"_+", "_", header).strip("_")
        return header

    @classmethod
    def parse_table(cls, table_path, delimiter=None, comment="#"):
        """Parse .txt/.csv into a generic TabularDataset."""
        import pandas as pd

        if not os.path.exists(table_path):
            raise FileNotFoundError(f"Table file not found: {table_path}")

        sep = delimiter
        if sep is None:
            ext = os.path.splitext(table_path)[1].lower()
            sep = "," if ext == ".csv" else r"\s+"
        try:
            df = pd.read_csv(
                table_path,
                sep=sep,
                engine="python",
                comment=comment,
                skip_blank_lines=True,
            )
        except pd.errors.EmptyDataError:
            raise ValueError(f"No valid entries found in table: {table_path}")
        except pd.errors.ParserError as exc:
            raise ValueError(
                f"Failed to parse submission table '{table_path}': {exc}"
            ) from exc

        if df.empty:
            raise ValueError(f"No valid entries found in table: {table_path}")

        df = df.rename(
            columns={c: cls.normalize_header(c) for c in df.columns}
        )
        return cls(df, source_path=table_path)

    @staticmethod
    def resolve_column(columns, candidates, required=True):
        normalized = {TabularDataset.normalize_header(c): c for c in columns}
        for candidate in candidates:
            key = TabularDataset.normalize_header(candidate)
            if key in normalized:
                return normalized[key]
        if required:
            raise ValueError(
                "Could not resolve required column. Tried: "
                + ", ".join(candidates)
            )
        return None

    def to_entries(self, entry_cls, row_offset=2):
        return [
            entry_cls(row.to_dict(), row_number=idx + row_offset)
            for idx, row in self.dataframe.iterrows()
        ]

    def validate(
        self,
        required_columns=None,
        integer_columns=None,
        positive_integer_columns=None,
        path_columns=None,
        check_file_exists=True,
    ):
        """Generic dataset-level validation without job-specific logic."""
        required_columns = required_columns or []
        integer_columns = integer_columns or []
        positive_integer_columns = positive_integer_columns or []
        path_columns = path_columns or []

        missing = [col for col in required_columns if col not in self.columns]
        if missing:
            raise ValueError(
                "Missing required table columns: " + ", ".join(missing)
            )

        errors = []
        for idx, row in self.dataframe.iterrows():
            row_no = idx + 2

            for col in integer_columns:
                value = row[col]
                if value is None:
                    errors.append(f"Missing {col} (row {row_no})")
                    continue
                try:
                    int(value)
                except (TypeError, ValueError):
                    errors.append(
                        f"Invalid integer for {col} at row {row_no}: {value!r}"
                    )

            for col in positive_integer_columns:
                value = row[col]
                if value is None:
                    errors.append(f"Missing {col} (row {row_no})")
                    continue
                try:
                    parsed = int(value)
                    if parsed < 1:
                        errors.append(
                            f"{col} must be >= 1 at row {row_no}, got {parsed}"
                        )
                except (TypeError, ValueError):
                    errors.append(
                        f"Invalid integer for {col} at row {row_no}: {value!r}"
                    )

            if check_file_exists:
                for col in path_columns:
                    value = row[col]
                    if not value:
                        errors.append(f"Missing {col} (row {row_no})")
                        continue
                    if not os.path.exists(str(value)):
                        errors.append(
                            f"File not found for {col} at row {row_no}: {value}"
                        )

        if errors:
            raise ValueError("Table validation failed:\n" + "\n".join(errors))
        return self
