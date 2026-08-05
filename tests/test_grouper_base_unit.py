"""
Direct unit tests for chemsmart.jobs.grouper.base:
ResultsRecorder, StructureGrouperConfig, and MoleculeGrouper.

These exercise the internal helper methods directly (rather than
via a concrete grouper subclass end-to-end) to close coverage gaps
in error branches, alternate output formats, and shared helpers.
"""

import os

import numpy as np
import pandas as pd
import pytest

from chemsmart.jobs.grouper.base import (
    MoleculeGrouper,
    ResultsRecorder,
    StructureGrouperConfig,
)


class _ConcreteGrouper(MoleculeGrouper):
    """Minimal concrete subclass so the abstract base can be instantiated."""

    def group(self):
        indices = list(range(len(self.molecules)))
        return [self.molecules], [indices]


class TestResultsRecorderConstruction:
    def test_rejects_unsupported_format(self, tmp_path):
        with pytest.raises(ValueError, match="Unsupported output format"):
            ResultsRecorder(output_dir=str(tmp_path), matrix_format="pdf")

    def test_normalizes_format_case_and_dot(self, tmp_path):
        recorder = ResultsRecorder(
            output_dir=str(tmp_path), matrix_format=".XLSX"
        )
        assert recorder.matrix_format == "xlsx"

    def test_creates_output_dir(self, tmp_path):
        target = tmp_path / "nested" / "output"
        ResultsRecorder(output_dir=str(target))
        assert target.is_dir()


class TestResultsRecorderLabels:
    def test_get_labels_uses_conformer_ids_when_length_matches(self, tmp_path):
        recorder = ResultsRecorder(
            output_dir=str(tmp_path), conformer_ids=["a", "b", "c"]
        )
        assert recorder.get_labels(3) == ["a", "b", "c"]

    def test_get_labels_falls_back_to_numeric_when_length_mismatches(
        self, tmp_path
    ):
        recorder = ResultsRecorder(
            output_dir=str(tmp_path), conformer_ids=["a", "b"]
        )
        assert recorder.get_labels(3) == ["1", "2", "3"]

    def test_get_labels_falls_back_when_no_conformer_ids(self, tmp_path):
        recorder = ResultsRecorder(output_dir=str(tmp_path))
        assert recorder.get_labels(2) == ["1", "2"]


class TestResultsRecorderGroupsDataframe:
    def test_build_groups_dataframe_basic(self, tmp_path):
        recorder = ResultsRecorder(output_dir=str(tmp_path))
        df = recorder.build_groups_dataframe(
            index_groups=[[0, 1], [2]], n_molecules=3
        )
        assert list(df["Group"]) == [1, 2]
        assert list(df["Members"]) == ["1, 2", "3"]

    def test_build_groups_dataframe_with_extra_columns(self, tmp_path):
        recorder = ResultsRecorder(output_dir=str(tmp_path))
        df = recorder.build_groups_dataframe(
            index_groups=[[0], [1]],
            n_molecules=2,
            extra_columns={"Score": [0.5]},
        )
        # Second row's extra column value is skipped because the values
        # list is shorter than the number of groups (exercises the
        # `if i < len(values)` false branch).
        assert df.loc[0, "Score"] == 0.5
        assert "Score" not in df.columns or pd.isna(df.loc[1, "Score"])


class TestResultsRecorderFilename:
    def test_get_filename_with_label_and_suffix(self, tmp_path):
        recorder = ResultsRecorder(
            output_dir=str(tmp_path), label="mysys", matrix_format="csv"
        )
        filename = recorder.get_filename("RMSDGrouper", suffix="T0.5")
        assert filename == os.path.join(
            str(tmp_path), "mysys_RMSDGrouper_T0.5.csv"
        )

    def test_get_filename_without_label_or_suffix(self, tmp_path):
        recorder = ResultsRecorder(output_dir=str(tmp_path))
        filename = recorder.get_filename("RMSDGrouper")
        assert filename == os.path.join(str(tmp_path), "RMSDGrouper.xlsx")


class TestResultsRecorderRecordResults:
    def _make_common_args(self):
        header_info = [
            ("Title", None),
            ("", "Some Title Line"),
            ("Key", ""),
            ("Method", "RMSD"),
        ]
        sheets_data = {
            "Groups": pd.DataFrame({"Group": [1], "Members": ["1"]})
        }
        matrix = np.array([[0.0, np.inf], [np.inf, 0.0]])
        matrix_data = ("Matrix", matrix, ["1", "2"])
        return header_info, sheets_data, matrix_data

    def test_record_results_xlsx_with_matrix(self, tmp_path):
        recorder = ResultsRecorder(
            output_dir=str(tmp_path), matrix_format="xlsx"
        )
        header_info, sheets_data, matrix_data = self._make_common_args()
        # Include a sheet sharing the matrix's name, to hit the "skip if
        # already written as matrix" continue branch in _write_xlsx.
        sheets_data[matrix_data[0]] = pd.DataFrame({"A": [1]})
        path = recorder.record_results(
            grouper_name="TestGrouper",
            header_info=header_info,
            sheets_data=sheets_data,
            matrix_data=matrix_data,
        )
        assert os.path.isfile(path)

    def test_record_results_xlsx_without_matrix(self, tmp_path):
        recorder = ResultsRecorder(
            output_dir=str(tmp_path), matrix_format="xlsx"
        )
        header_info, sheets_data, _ = self._make_common_args()
        sheets_data["Other"] = pd.DataFrame({"A": [1, 2]})
        path = recorder.record_results(
            grouper_name="TestGrouper",
            header_info=header_info,
            sheets_data=sheets_data,
            startrow=3,
        )
        assert os.path.isfile(path)

    def test_record_results_csv(self, tmp_path):
        recorder = ResultsRecorder(
            output_dir=str(tmp_path), matrix_format="csv"
        )
        header_info, sheets_data, matrix_data = self._make_common_args()
        sheets_data["Extra"] = pd.DataFrame({"A": [1]})
        # Also include a sheet with the same name as the matrix sheet, to
        # exercise the "skip if already written as matrix" continue branch.
        sheets_data[matrix_data[0]] = pd.DataFrame({"A": [1]})
        path = recorder.record_results(
            grouper_name="TestGrouper",
            header_info=header_info,
            sheets_data=sheets_data,
            matrix_data=matrix_data,
        )
        assert os.path.isfile(path)
        extra_file = path[:-4] + "_Extra.csv"
        assert os.path.isfile(extra_file)

    def test_record_results_csv_without_matrix(self, tmp_path):
        recorder = ResultsRecorder(
            output_dir=str(tmp_path), matrix_format="csv"
        )
        header_info, sheets_data, _ = self._make_common_args()
        path = recorder.record_results(
            grouper_name="TestGrouper",
            header_info=header_info,
            sheets_data=sheets_data,
        )
        assert os.path.isfile(path)

    def test_record_results_txt(self, tmp_path):
        recorder = ResultsRecorder(
            output_dir=str(tmp_path), matrix_format="txt"
        )
        header_info, sheets_data, matrix_data = self._make_common_args()
        sheets_data["Extra"] = pd.DataFrame({"A": [1]})
        sheets_data[matrix_data[0]] = pd.DataFrame({"A": [1]})
        path = recorder.record_results(
            grouper_name="TestGrouper",
            header_info=header_info,
            sheets_data=sheets_data,
            matrix_data=matrix_data,
        )
        assert os.path.isfile(path)
        content = open(path).read()
        assert "Matrix" in content
        assert "inf" in content

    def test_record_results_txt_without_matrix(self, tmp_path):
        recorder = ResultsRecorder(
            output_dir=str(tmp_path), matrix_format="txt"
        )
        header_info, sheets_data, _ = self._make_common_args()
        path = recorder.record_results(
            grouper_name="TestGrouper",
            header_info=header_info,
            sheets_data=sheets_data,
        )
        assert os.path.isfile(path)
        content = open(path).read()
        assert "Groups" in content
        assert "Matrix" not in content

    def test_record_results_invalid_format_raises(self, tmp_path):
        recorder = ResultsRecorder(
            output_dir=str(tmp_path), matrix_format="csv"
        )
        # Force an unsupported format after construction to hit the
        # otherwise-unreachable else branch in record_results.
        recorder.matrix_format = "bogus"
        with pytest.raises(ValueError, match="Unsupported output format"):
            recorder.record_results(
                grouper_name="TestGrouper",
                header_info=[],
                sheets_data={},
            )


class TestStructureGrouperConfig:
    def test_defaults(self):
        config = StructureGrouperConfig()
        assert config.ltol == 0.1
        assert config.stol == 0.18
        assert config.angle_tol == 1

    def test_custom_values(self):
        config = StructureGrouperConfig(ltol=0.2, stol=0.3, angle_tol=5)
        assert config.ltol == 0.2
        assert config.stol == 0.3
        assert config.angle_tol == 5


class TestMoleculeGrouperValidation:
    # Note: the `not isinstance(self.molecules, Iterable)` branch in
    # _validate_inputs is unreachable through the public API: by the
    # time it runs, __init__ has already done `self.molecules =
    # list(molecules)`, which either raised TypeError itself (for a
    # genuinely non-iterable argument) or produced a real `list`
    # (always Iterable). So only the "non-Molecule items" branch below
    # is reachable in practice.

    def test_rejects_non_molecule_items(self):
        with pytest.raises(TypeError, match="Molecule instances"):
            _ConcreteGrouper(molecules=[1, 2, 3])

    def test_accepts_valid_molecules(self, water_molecule):
        grouper = _ConcreteGrouper(molecules=[water_molecule])
        assert grouper.molecules == [water_molecule]
        assert grouper.num_procs == 1

    def test_num_procs_floored_at_one(self, water_molecule):
        grouper = _ConcreteGrouper(molecules=[water_molecule], num_procs=0)
        assert grouper.num_procs == 1


class TestMoleculeGrouperOutputDir:
    def test_get_output_dir_with_label(self, water_molecule, tmp_path):
        grouper = _ConcreteGrouper(
            molecules=[water_molecule],
            label="sys1",
            output_dir=str(tmp_path),
        )
        assert grouper._get_output_dir() == os.path.join(
            str(tmp_path), "sys1_group_result"
        )

    def test_get_output_dir_without_label_or_output_dir(self, water_molecule):
        grouper = _ConcreteGrouper(molecules=[water_molecule])
        assert grouper._get_output_dir() == "group_result"

    def test_get_results_recorder(self, water_molecule, tmp_path):
        grouper = _ConcreteGrouper(
            molecules=[water_molecule],
            label="sys1",
            output_dir=str(tmp_path),
        )
        recorder = grouper._get_results_recorder()
        assert isinstance(recorder, ResultsRecorder)
        assert recorder.label == "sys1"


class TestMoleculeGrouperRecordDelegation:
    def test_record_delegates_to_record_results(self, water_molecule):
        grouper = _ConcreteGrouper(molecules=[water_molecule])
        with pytest.raises(NotImplementedError, match="_record_results"):
            grouper.record()


class TestMoleculeGrouperHeaderHelpers:
    @pytest.mark.parametrize("energy_type", ["QHH", "QHG", "SP_QHG", "qhh"])
    def test_append_thermo_header_included_for_thermo_energy_types(
        self, water_molecule, energy_type
    ):
        grouper = _ConcreteGrouper(
            molecules=[water_molecule],
            energy_type=energy_type,
            thermo_parameters="s=1,T=298.15",
        )
        header_info = []
        grouper._append_thermo_header(header_info)
        assert header_info == [("Thermochemistry Parameters", "s=1,T=298.15")]

    def test_append_thermo_header_skipped_for_non_thermo_energy_type(
        self, water_molecule
    ):
        grouper = _ConcreteGrouper(
            molecules=[water_molecule],
            energy_type="E",
            thermo_parameters="s=1,T=298.15",
        )
        header_info = []
        grouper._append_thermo_header(header_info)
        assert header_info == []

    def test_append_thermo_header_skipped_without_thermo_parameters(
        self, water_molecule
    ):
        grouper = _ConcreteGrouper(
            molecules=[water_molecule], energy_type="QHH"
        )
        header_info = []
        grouper._append_thermo_header(header_info)
        assert header_info == []

    def test_append_input_usage_header_with_skipped_ids(self, water_molecule):
        grouper = _ConcreteGrouper(
            molecules=[water_molecule], skipped_ids=["c3", "c4"]
        )
        header_info = []
        grouper._append_input_usage_header(header_info)
        assert ("Used Molecules", 1) in header_info
        assert ("Skipped Molecules", 2) in header_info
        assert ("Skipped IDs", "c3, c4") in header_info

    def test_append_input_usage_header_without_skipped_ids(
        self, water_molecule
    ):
        grouper = _ConcreteGrouper(molecules=[water_molecule])
        header_info = []
        grouper._append_input_usage_header(header_info)
        assert ("Skipped IDs", "None") in header_info


class TestMoleculeGrouperUnique:
    def test_unique_writes_xyz_and_returns_lowest_energy(
        self, water_molecule, methane_molecule, tmp_path
    ):
        mol_low = water_molecule
        mol_low.energy = -1.0
        mol_high = methane_molecule
        mol_high.energy = 5.0

        grouper = _ConcreteGrouper(
            molecules=[mol_low, mol_high], output_dir=str(tmp_path)
        )
        unique = grouper.unique(output_dir=str(tmp_path))

        assert len(unique) == 1
        assert unique[0] is mol_low

        group_file = os.path.join(str(tmp_path), "group_result", "group_1.xyz")
        assert os.path.isfile(group_file)
        content = open(group_file).read()
        assert "Group 1 Member 1" in content
        assert "Group 1 Member 2" in content

    def test_unique_uses_cached_groups_on_second_call(
        self, water_molecule, tmp_path
    ):
        water_molecule.energy = None
        grouper = _ConcreteGrouper(
            molecules=[water_molecule],
            label="cached",
            output_dir=str(tmp_path),
        )
        first = grouper.unique(output_dir=str(tmp_path))
        assert grouper._cached_groups is not None

        # Second call should hit the cached-groups branch without
        # recomputation (verified indirectly: result is consistent).
        second = grouper.unique(output_dir=str(tmp_path))
        assert first[0] is second[0]

    def test_unique_with_conformer_ids_labels_original_index(
        self, water_molecule, tmp_path
    ):
        water_molecule.energy = None
        grouper = _ConcreteGrouper(
            molecules=[water_molecule],
            conformer_ids=["confA"],
            output_dir=str(tmp_path),
        )
        grouper.unique(output_dir=str(tmp_path))
        group_file = os.path.join(str(tmp_path), "group_result", "group_1.xyz")
        content = open(group_file).read()
        assert "Original_Index: confA" in content
        assert "N/A" in content
