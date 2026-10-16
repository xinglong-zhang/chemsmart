"""
Direct unit tests filling coverage gaps in ``chemsmart.utils.io`` left by
``tests/test_io_utils.py``: the exception/early-break paths of
``get_program_type_from_file``, ``check_program_availability_in_chemsmart``
(previously untested), ``resolve_output_path`` (previously untested), and
the unreachable-type guard in ``select_items_by_index``.
"""

from unittest.mock import patch

import pytest

from chemsmart.utils.io import (
    check_program_availability_in_chemsmart,
    get_program_type_from_file,
    resolve_output_path,
    select_items_by_index,
)


class TestGetProgramTypeFromFileErrorPaths:
    def test_nonexistent_file_returns_unknown(self):
        assert (
            get_program_type_from_file("/no/such/file/does_not_exist.log")
            == "unknown"
        )

    def test_match_beyond_first_200_lines_is_not_detected(self, tmp_path):
        # The scan stops after 200 lines, so a marker appearing later is
        # never seen and the file is reported as "unknown".
        lines = ["irrelevant text\n"] * 250
        lines[210] = "* O   R   C   A *\n"
        outfile = tmp_path / "late_marker.out"
        outfile.write_text("".join(lines))
        assert get_program_type_from_file(str(outfile)) == "unknown"

    def test_match_within_first_200_lines_is_detected(self, tmp_path):
        lines = ["irrelevant text\n"] * 250
        lines[50] = "* O   R   C   A *\n"
        outfile = tmp_path / "early_marker.out"
        outfile.write_text("".join(lines))
        assert get_program_type_from_file(str(outfile)) == "orca"

    def test_blank_lines_are_skipped(self, tmp_path):
        outfile = tmp_path / "blank_then_marker.out"
        outfile.write_text("\n\n\n* O   R   C   A *\n")
        assert get_program_type_from_file(str(outfile)) == "orca"


class TestCheckProgramAvailability:
    @pytest.mark.parametrize("name", ["gaussian", "Gaussian", "GAUSSIAN"])
    def test_gaussian_case_insensitive_is_accepted(self, name):
        check_program_availability_in_chemsmart(name)

    @pytest.mark.parametrize("name", ["orca", "Orca", "ORCA"])
    def test_orca_case_insensitive_is_accepted(self, name):
        check_program_availability_in_chemsmart(name)

    def test_unsupported_program_raises(self):
        with pytest.raises(ValueError, match="Unsupported program"):
            check_program_availability_in_chemsmart("nwchem")


class TestResolveOutputPath:
    def test_different_paths_returned_unchanged(self, tmp_path):
        in_path = tmp_path / "input.xyz"
        out_path = tmp_path / "output.xyz"
        in_path.write_text("x")
        resolved, changed = resolve_output_path(str(in_path), str(out_path))
        assert resolved == str(out_path)
        assert changed is False

    def test_same_path_gets_numeric_suffix(self, tmp_path):
        path = tmp_path / "same.xyz"
        path.write_text("x")
        resolved, changed = resolve_output_path(str(path), str(path))
        assert resolved == str(tmp_path / "same_1.xyz")
        assert changed is True

    def test_existing_suffix_candidate_is_skipped(self, tmp_path):
        path = tmp_path / "same.xyz"
        path.write_text("x")
        (tmp_path / "same_1.xyz").write_text("taken")
        resolved, changed = resolve_output_path(str(path), str(path))
        assert resolved == str(tmp_path / "same_2.xyz")
        assert changed is True


class TestSelectItemsByIndexUnreachableTypeGuard:
    def test_unexpected_return_type_from_parser_raises(self):
        with patch(
            "chemsmart.utils.utils.parse_index_specification",
            return_value="not-a-valid-type",
        ):
            with pytest.raises(ValueError, match="Unexpected index type"):
                select_items_by_index(["a", "b", "c"], "1,2")
