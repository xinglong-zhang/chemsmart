"""
Direct unit tests for chemsmart.jobs.orca.runner.FakeORCA, the fake
ORCA execution simulator used by FakeORCAJobRunner. No prior test
coverage existed for this class at all.
"""

import os
from shutil import copyfile

import pytest

from chemsmart.jobs.orca.runner import FakeORCA


@pytest.fixture()
def water_sp_inp_copy(tmp_path, water_sp_input_path):
    target = tmp_path / "water_sp.inp"
    copyfile(water_sp_input_path, target)
    return str(target)


class TestFakeORCAConstruction:
    def test_raises_when_file_missing(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="not found"):
            FakeORCA(str(tmp_path / "missing.inp"))

    def test_stores_absolute_path(self, water_sp_inp_copy):
        fake = FakeORCA(water_sp_inp_copy)
        assert fake.file_to_run == os.path.abspath(water_sp_inp_copy)


class TestFakeORCAProperties:
    def test_file_folder_and_filename(self, water_sp_inp_copy):
        fake = FakeORCA(water_sp_inp_copy)
        assert fake.file_folder == os.path.dirname(water_sp_inp_copy)
        assert fake.filename == "water_sp.inp"

    def test_input_filepath(self, water_sp_inp_copy):
        fake = FakeORCA(water_sp_inp_copy)
        assert fake.input_filepath == os.path.abspath(water_sp_inp_copy)

    def test_output_filepath_uses_out_extension(self, water_sp_inp_copy):
        fake = FakeORCA(water_sp_inp_copy)
        assert fake.output_filepath == os.path.join(
            fake.file_folder, "water_sp.out"
        )

    def test_input_contents_matches_file(self, water_sp_inp_copy):
        fake = FakeORCA(water_sp_inp_copy)
        assert "!HF DEF2-SVP" in fake.input_contents
        assert any("MAXITER 500" in line for line in fake.input_contents)

    def test_molecule_charge_and_multiplicity(self, water_sp_inp_copy):
        fake = FakeORCA(water_sp_inp_copy)
        assert fake.charge == 0
        assert fake.multiplicity == 1
        assert fake.molecule is not None

    def test_spin_restricted_for_singlet(self, water_sp_inp_copy):
        fake = FakeORCA(water_sp_inp_copy)
        assert fake.spin == "R"

    def test_num_atoms_and_symbols(self, water_sp_inp_copy):
        fake = FakeORCA(water_sp_inp_copy)
        assert fake.num_atoms == 3
        assert fake.atomic_symbols == ["O", "H", "H"]

    def test_atomic_numbers(self, water_sp_inp_copy):
        fake = FakeORCA(water_sp_inp_copy)
        assert fake.atomic_numbers == [8, 1, 1]

    def test_atomic_coordinates_shape(self, water_sp_inp_copy):
        fake = FakeORCA(water_sp_inp_copy)
        assert len(fake.atomic_coordinates) == 3

    def test_empirical_formula(self, water_sp_inp_copy):
        fake = FakeORCA(water_sp_inp_copy)
        assert fake.empirical_formula == "H2O"


class TestFakeORCASpinUnrestricted:
    def test_spin_unrestricted_for_non_singlet(self, tmp_path):
        inp_path = tmp_path / "radical.inp"
        inp_path.write_text(
            "!HF DEF2-SVP\n* xyz 0 2\nH   0.0000   0.0000   0.0000\n*\n"
        )
        fake = FakeORCA(str(inp_path))
        assert fake.spin == "U"


class TestFakeORCARun:
    def test_run_writes_fake_output_file(self, water_sp_inp_copy):
        fake = FakeORCA(water_sp_inp_copy)
        returncode = fake.run()
        assert returncode is None or returncode == 0
        assert os.path.exists(fake.output_filepath)
        with open(fake.output_filepath) as f:
            content = f.read()
        assert "O   R   C   A" in content
        assert "ORCA TERMINATED NORMALLY" in content
        assert "Number of atoms" in content
        assert "3" in content
