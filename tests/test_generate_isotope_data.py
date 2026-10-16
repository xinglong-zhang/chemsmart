import importlib.util
from pathlib import Path

import pytest

import chemsmart
from chemsmart.scripts.generate_isotope_data import (
    generate_isotope_data_file,
    parse_isotope_file,
)

ISOTOPES_TXT = Path(chemsmart.__file__).parent / "utils" / "isotopes.txt"


class TestParseIsotopeFile:
    def test_parses_hydrogen_isotopes(self):
        isotopes = parse_isotope_file(ISOTOPES_TXT)
        assert 1 in isotopes
        assert "most_abundant" in isotopes[1]
        assert isotopes[1]["most_abundant"]["mass_number"] == 1


class TestGenerateIsotopeDataFile:
    def test_writes_importable_module_with_expected_data(self, tmp_path):
        output_file = tmp_path / "isotopes_data.py"

        generate_isotope_data_file(ISOTOPES_TXT, output_file)

        assert output_file.exists()
        content = output_file.read_text()
        assert "isotopes = {" in content
        assert "weighted_atomic_mass" in content

        spec = importlib.util.spec_from_file_location(
            "generated_isotopes_data", output_file
        )
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        assert 1 in module.isotopes
        hydrogen = module.isotopes[1]
        assert "most_abundant" in hydrogen
        assert "weighted_atomic_mass" in hydrogen
        assert hydrogen["weighted_atomic_mass"] == pytest.approx(
            1.008, abs=0.05
        )
