import importlib.util
from pathlib import Path
from unittest.mock import patch

import pytest

import chemsmart
from chemsmart.scripts.generate_isotope_data import (
    generate_isotope_data_file,
    most_stable_mass_numbers,
    parse_isotope_file,
)

ISOTOPES_TXT = Path(chemsmart.__file__).parent / "utils" / "isotopes.txt"


def _write_isotope_entry(
    path, atomic_number, mass_number, mass, composition=""
):
    with open(path, "a") as f:
        f.write(
            f"Atomic Number = {atomic_number}\n"
            f"Mass Number = {mass_number}\n"
            f"Relative Atomic Mass = {mass}\n"
            f"Isotopic Composition = {composition}\n"
            "\n"
        )


class TestParseIsotopeFile:
    def test_parses_hydrogen_isotopes(self):
        isotopes = parse_isotope_file(ISOTOPES_TXT)
        assert 1 in isotopes
        assert "most_abundant" in isotopes[1]
        assert isotopes[1]["most_abundant"]["mass_number"] == 1

    def test_weighted_mass_from_abundances(self, tmp_path):
        # Two isotopes with known abundances -> most_abundant and
        # weighted_atomic_mass should follow the abundance-weighted branch.
        input_file = tmp_path / "custom.txt"
        _write_isotope_entry(input_file, 6, 12, 12.0, "0.9")
        _write_isotope_entry(input_file, 6, 13, 13.0, "0.1")

        isotopes = parse_isotope_file(input_file)

        assert isotopes[6]["most_abundant"]["mass_number"] == 12
        assert isotopes[6]["weighted_atomic_mass"] == pytest.approx(
            12.0 * 0.9 + 13.0 * 0.1
        )

    def test_blank_composition_treated_as_zero_abundance(self, tmp_path):
        # A blank "Isotopic Composition =" line (as NIST uses for some
        # minor isotopes, e.g. tritium) has no regex match and falls
        # back to 0.0 abundance rather than raising.
        input_file = tmp_path / "custom.txt"
        _write_isotope_entry(input_file, 6, 12, 12.0, "0.9")
        _write_isotope_entry(input_file, 6, 14, 14.0, composition="")

        isotopes = parse_isotope_file(input_file)

        assert isotopes[6][14]["abundance"] == 0.0

    def test_fully_radioactive_element_falls_back_to_most_stable(
        self, tmp_path
    ):
        # An element with no isotope abundances at all (e.g. Tc, Z=43)
        # cannot pick a most-abundant isotope by abundance, so both the
        # most_abundant and weighted_atomic_mass calculations fall back
        # to the hardcoded most_stable_mass_numbers entry.
        atomic_number = 43
        stable_mass_number = most_stable_mass_numbers[atomic_number]
        input_file = tmp_path / "custom.txt"
        _write_isotope_entry(
            input_file, atomic_number, stable_mass_number, 97.9, ""
        )
        _write_isotope_entry(input_file, atomic_number, 99, 98.9, "")

        isotopes = parse_isotope_file(input_file)

        isotope_data = isotopes[atomic_number]
        assert (
            isotope_data["most_abundant"]["mass_number"] == stable_mass_number
        )
        assert isotope_data["weighted_atomic_mass"] == pytest.approx(97.9)


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

    def test_omits_optional_fields_when_absent_from_parsed_data(
        self, tmp_path
    ):
        # The writer defensively checks for "mass"/"abundance" per
        # isotope and "most_abundant"/"weighted_atomic_mass" per element
        # before emitting them, even though the current parser always
        # supplies all of them together. Exercise those omission
        # branches directly against a hand-built isotopes dict.
        output_file = tmp_path / "isotopes_data.py"
        fake_isotopes = {
            1: {
                1: {"mass": 1.0},  # missing "abundance" -> entry skipped
            },
        }

        with patch(
            "chemsmart.scripts.generate_isotope_data.parse_isotope_file",
            return_value=fake_isotopes,
        ):
            generate_isotope_data_file("unused", output_file)

        content = output_file.read_text()
        assert "1: {" in content
        assert "most_abundant" not in content
        assert "weighted_atomic_mass" not in content
