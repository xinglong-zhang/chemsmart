"""
Direct unit tests for chemsmart.io.datasets.TabularDataset.

No existing test coverage previously existed for this module.
"""

import pandas as pd
import pytest

from chemsmart.io.datasets import TabularDataset


class TestNormalizeHeader:
    @pytest.mark.parametrize(
        "raw,expected",
        [
            ("Total Energy", "total_energy"),
            ("  My-Col!! ", "my_col"),
            ("already_snake", "already_snake"),
            ("Multi   Space", "multi_space"),
            (123, "123"),
        ],
    )
    def test_normalizes_various_headers(self, raw, expected):
        assert TabularDataset.normalize_header(raw) == expected


class TestBasicProperties:
    def test_columns_and_len(self):
        df = pd.DataFrame({"a": [1, 2], "b": [3, 4]})
        dataset = TabularDataset(df)
        assert dataset.columns == ["a", "b"]
        assert len(dataset) == 2

    def test_resets_index(self):
        df = pd.DataFrame({"a": [1, 2]}, index=[5, 9])
        dataset = TabularDataset(df)
        assert list(dataset.dataframe.index) == [0, 1]


class TestParseTable:
    def test_file_not_found_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="Table file not found"):
            TabularDataset.parse_table(str(tmp_path / "missing.csv"))

    def test_parses_csv_with_normalized_headers(self, tmp_path):
        csv_path = tmp_path / "data.csv"
        csv_path.write_text("Total Energy,Job Name\n1.0,job1\n2.0,job2\n")
        dataset = TabularDataset.parse_table(str(csv_path))
        assert dataset.columns == ["total_energy", "job_name"]
        assert len(dataset) == 2
        assert dataset.source_path == str(csv_path)

    def test_parses_whitespace_delimited_txt(self, tmp_path):
        txt_path = tmp_path / "data.txt"
        txt_path.write_text("name  value\nfoo   1\nbar   2\n")
        dataset = TabularDataset.parse_table(str(txt_path))
        assert dataset.columns == ["name", "value"]
        assert len(dataset) == 2

    def test_explicit_delimiter_overrides_extension_default(self, tmp_path):
        txt_path = tmp_path / "data.txt"
        txt_path.write_text("a;b\n1;2\n")
        dataset = TabularDataset.parse_table(str(txt_path), delimiter=";")
        assert dataset.columns == ["a", "b"]

    def test_empty_file_raises_value_error(self, tmp_path):
        csv_path = tmp_path / "empty.csv"
        csv_path.write_text("")
        with pytest.raises(ValueError, match="No valid entries found"):
            TabularDataset.parse_table(str(csv_path))

    def test_header_only_file_raises_value_error(self, tmp_path):
        csv_path = tmp_path / "header_only.csv"
        csv_path.write_text("a,b\n")
        with pytest.raises(ValueError, match="No valid entries found"):
            TabularDataset.parse_table(str(csv_path))

    def test_malformed_file_raises_value_error(self, tmp_path):
        # An unterminated quoted field triggers a ParserError from the
        # python engine ("unexpected end of data").
        csv_path = tmp_path / "bad.csv"
        csv_path.write_text('a,b\n"unterminated,1\n2,3\n')
        with pytest.raises(ValueError, match="Failed to parse"):
            TabularDataset.parse_table(str(csv_path))

    def test_comment_lines_are_skipped(self, tmp_path):
        csv_path = tmp_path / "with_comments.csv"
        csv_path.write_text("# a comment\na,b\n1,2\n# another comment\n3,4\n")
        dataset = TabularDataset.parse_table(str(csv_path))
        assert len(dataset) == 2


class TestResolveColumn:
    def test_resolves_first_matching_candidate(self):
        columns = ["Total Energy", "Job Name"]
        resolved = TabularDataset.resolve_column(
            columns, ["energy", "total_energy"]
        )
        assert resolved == "Total Energy"

    def test_returns_none_when_optional_and_not_found(self):
        resolved = TabularDataset.resolve_column(
            ["a", "b"], ["missing"], required=False
        )
        assert resolved is None

    def test_raises_when_required_and_not_found(self):
        with pytest.raises(ValueError, match="Could not resolve"):
            TabularDataset.resolve_column(["a", "b"], ["missing"])


class SimpleEntry:
    def __init__(self, row_dict, row_number):
        self.row_dict = row_dict
        self.row_number = row_number


class TestToEntries:
    def test_builds_entries_with_row_numbers(self):
        df = pd.DataFrame({"a": [1, 2, 3]})
        dataset = TabularDataset(df)
        entries = dataset.to_entries(SimpleEntry)
        assert [e.row_number for e in entries] == [2, 3, 4]
        assert entries[0].row_dict == {"a": 1}

    def test_custom_row_offset(self):
        df = pd.DataFrame({"a": [1]})
        dataset = TabularDataset(df)
        entries = dataset.to_entries(SimpleEntry, row_offset=10)
        assert entries[0].row_number == 10


class TestValidate:
    def test_missing_required_columns_raises(self):
        df = pd.DataFrame({"a": [1]})
        dataset = TabularDataset(df)
        with pytest.raises(ValueError, match="Missing required table col"):
            dataset.validate(required_columns=["a", "b"])

    def test_passes_with_all_required_columns_present(self):
        df = pd.DataFrame({"a": [1], "b": [2]})
        dataset = TabularDataset(df)
        assert dataset.validate(required_columns=["a", "b"]) is dataset

    def test_integer_column_missing_value(self):
        df = pd.DataFrame({"a": [None]})
        dataset = TabularDataset(df)
        with pytest.raises(ValueError, match="Missing a"):
            dataset.validate(integer_columns=["a"])

    def test_integer_column_invalid_value(self):
        df = pd.DataFrame({"a": ["not-an-int"]})
        dataset = TabularDataset(df)
        with pytest.raises(ValueError, match="Invalid integer for a"):
            dataset.validate(integer_columns=["a"])

    def test_integer_column_valid_value_passes(self):
        df = pd.DataFrame({"a": [5]})
        dataset = TabularDataset(df)
        assert dataset.validate(integer_columns=["a"]) is dataset

    def test_positive_integer_column_missing_value(self):
        df = pd.DataFrame({"a": [None]})
        dataset = TabularDataset(df)
        with pytest.raises(ValueError, match="Missing a"):
            dataset.validate(positive_integer_columns=["a"])

    def test_positive_integer_column_below_one_fails(self):
        df = pd.DataFrame({"a": [0]})
        dataset = TabularDataset(df)
        with pytest.raises(ValueError, match="must be >= 1"):
            dataset.validate(positive_integer_columns=["a"])

    def test_positive_integer_column_invalid_value(self):
        df = pd.DataFrame({"a": ["bad"]})
        dataset = TabularDataset(df)
        with pytest.raises(ValueError, match="Invalid integer for a"):
            dataset.validate(positive_integer_columns=["a"])

    def test_positive_integer_column_valid_value_passes(self):
        df = pd.DataFrame({"a": [3]})
        dataset = TabularDataset(df)
        assert dataset.validate(positive_integer_columns=["a"]) is dataset

    def test_path_column_missing_value(self):
        df = pd.DataFrame({"path": [None]})
        dataset = TabularDataset(df)
        with pytest.raises(ValueError, match="Missing path"):
            dataset.validate(path_columns=["path"])

    def test_path_column_file_not_found(self):
        df = pd.DataFrame({"path": ["/does/not/exist.xyz"]})
        dataset = TabularDataset(df)
        with pytest.raises(ValueError, match="File not found for path"):
            dataset.validate(path_columns=["path"])

    def test_path_column_valid_file_passes(self, tmp_path):
        existing = tmp_path / "exists.xyz"
        existing.write_text("1\ncomment\nH 0 0 0\n")
        df = pd.DataFrame({"path": [str(existing)]})
        dataset = TabularDataset(df)
        assert dataset.validate(path_columns=["path"]) is dataset

    def test_path_column_skipped_when_check_file_exists_false(self):
        df = pd.DataFrame({"path": ["/does/not/exist.xyz"]})
        dataset = TabularDataset(df)
        assert (
            dataset.validate(path_columns=["path"], check_file_exists=False)
            is dataset
        )

    def test_multiple_errors_are_aggregated(self):
        df = pd.DataFrame({"a": ["bad"], "b": [0]})
        dataset = TabularDataset(df)
        with pytest.raises(ValueError) as excinfo:
            dataset.validate(
                integer_columns=["a"], positive_integer_columns=["b"]
            )
        message = str(excinfo.value)
        assert "Invalid integer for a" in message
        assert "must be >= 1" in message
