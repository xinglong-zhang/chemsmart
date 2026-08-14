"""
Tests for the larger standalone analysis CLI scripts under
``chemsmart/scripts`` (fmo, fukui, plot_dias, structure_filter,
generate_isotope_data, get_thermochemistry).

Underlying analysis classes are mocked so these tests exercise only the
CLI argument-parsing/dispatch logic, without real quantum-chemistry
output files or real computation.
"""

import textwrap
from unittest.mock import MagicMock, patch

import pytest
from click.testing import CliRunner


class TestFmoScript:
    def test_closed_shell_system(self):
        from chemsmart.scripts.fmo import entry_point

        runner = CliRunner()
        with (
            patch(
                "chemsmart.scripts.fmo.get_program_type_from_file",
                return_value="gaussian",
            ),
            patch("chemsmart.scripts.fmo.Gaussian16Output") as mock_cls,
        ):
            mock_output = MagicMock()
            mock_output.multiplicity = 1
            mock_output.homo_energy = -5.0
            mock_output.lumo_energy = -1.0
            mock_output.fmo_gap = 4.0
            mock_cls.return_value = mock_output
            result = runner.invoke(
                entry_point,
                ["-f", "some.log"],
                catch_exceptions=False,
            )

        assert result.exit_code == 0, result.output

    def test_closed_shell_missing_energies_returns_early(self):
        from chemsmart.scripts.fmo import entry_point

        runner = CliRunner()
        with (
            patch(
                "chemsmart.scripts.fmo.get_program_type_from_file",
                return_value="gaussian",
            ),
            patch("chemsmart.scripts.fmo.Gaussian16Output") as mock_cls,
        ):
            mock_output = MagicMock()
            mock_output.multiplicity = 1
            mock_output.homo_energy = None
            mock_output.lumo_energy = None
            mock_cls.return_value = mock_output
            result = runner.invoke(
                entry_point,
                ["-f", "some.log"],
                catch_exceptions=False,
            )

        assert result.exit_code == 0, result.output

    def test_open_shell_system_kcal_mol_unit(self):
        from chemsmart.scripts.fmo import entry_point

        runner = CliRunner()
        with (
            patch(
                "chemsmart.scripts.fmo.get_program_type_from_file",
                return_value="orca",
            ),
            patch("chemsmart.scripts.fmo.ORCAOutput") as mock_cls,
        ):
            mock_output = MagicMock()
            mock_output.multiplicity = 2
            mock_output.alpha_homo_energy = -5.0
            mock_output.alpha_lumo_energy = -1.0
            mock_output.beta_homo_energy = -6.0
            mock_output.beta_lumo_energy = -0.5
            mock_output.somo_energies = [-3.0]
            mock_output.highest_somo_energy = -3.0
            mock_output.lowest_somo_energy = -3.0
            mock_output.num_unpaired_electrons = 1
            mock_output.alpha_fmo_gap = 4.0
            mock_output.beta_fmo_gap = 5.5
            mock_output.fmo_gap = 2.0
            mock_cls.return_value = mock_output
            result = runner.invoke(
                entry_point,
                ["-f", "some.out", "-u", "kcal/mol"],
                catch_exceptions=False,
            )

        assert result.exit_code == 0, result.output

    def test_unknown_program_raises(self):
        from chemsmart.scripts.fmo import entry_point

        runner = CliRunner()
        with patch(
            "chemsmart.scripts.fmo.get_program_type_from_file",
            return_value="unknown",
        ):
            with pytest.raises(TypeError, match="unknown filetype"):
                runner.invoke(
                    entry_point,
                    ["-f", "some.txt"],
                    catch_exceptions=False,
                )


class TestFukuiScript:
    def test_requires_at_least_one_radical_file(self):
        from chemsmart.scripts.fukui import entry_point

        runner = CliRunner()
        with pytest.raises(ValueError, match="At least one"):
            runner.invoke(
                entry_point,
                ["-n", "neutral.log"],
                catch_exceptions=False,
            )

    def test_mulliken_mode_with_cation_and_anion(self):
        from chemsmart.scripts.fukui import entry_point

        runner = CliRunner()

        def make_output(shift):
            output = MagicMock()
            output.energies = [-100.0 + shift]
            output.mulliken_atomic_charges = {"1C": 0.1 + shift}
            return output

        with (
            patch(
                "chemsmart.scripts.fukui.get_program_type_from_file",
                return_value="gaussian",
            ),
            patch("chemsmart.scripts.fukui.Gaussian16WBIOutput") as mock_cls,
        ):
            mock_cls.side_effect = [
                make_output(0.0),
                make_output(0.1),
                make_output(-0.1),
            ]
            result = runner.invoke(
                entry_point,
                [
                    "-n",
                    "neutral.log",
                    "-c",
                    "cation.log",
                    "-a",
                    "anion.log",
                ],
                catch_exceptions=False,
            )

        assert result.exit_code == 0, result.output

    def test_unknown_mode_rejected_by_cli(self):
        """``-m`` is restricted to a fixed choice set at the CLI level."""
        from chemsmart.scripts.fukui import entry_point

        runner = CliRunner()
        result = runner.invoke(
            entry_point,
            [
                "-n",
                "neutral.log",
                "-c",
                "cation.log",
                "-m",
                "bogus",
            ],
        )
        assert result.exit_code != 0

    def test_unknown_program_raises(self):
        from chemsmart.scripts.fukui import entry_point

        runner = CliRunner()
        with patch(
            "chemsmart.scripts.fukui.get_program_type_from_file",
            return_value="unknown",
        ):
            with pytest.raises(TypeError, match="unknown filetype"):
                runner.invoke(
                    entry_point,
                    ["-n", "neutral.txt", "-c", "cation.txt"],
                    catch_exceptions=False,
                )


class TestPlotDiasScript:
    def test_gaussian_dispatch(self):
        from chemsmart.scripts.plot_dias import entry_point

        runner = CliRunner()
        with patch(
            "chemsmart.scripts.plot_dias.GaussianDIASLogFolder"
        ) as mock_cls:
            mock_folder = MagicMock()
            mock_cls.return_value = mock_folder
            result = runner.invoke(
                entry_point,
                ["-p", "gaussian", "-a", "1", "-b", "2"],
                catch_exceptions=False,
            )

        assert result.exit_code == 0, result.output
        mock_folder.write_data.assert_called_once()
        mock_folder.plot_dias.assert_called_once()

    def test_orca_dispatch(self):
        from chemsmart.scripts.plot_dias import entry_point

        runner = CliRunner()
        with patch(
            "chemsmart.scripts.plot_dias.ORCADIASOutFolder"
        ) as mock_cls:
            mock_folder = MagicMock()
            mock_cls.return_value = mock_folder
            result = runner.invoke(
                entry_point,
                ["-p", "orca", "-a", "1", "-b", "2"],
                catch_exceptions=False,
            )

        assert result.exit_code == 0, result.output
        mock_folder.write_data.assert_called_once()

    def test_unknown_program_raises(self):
        from chemsmart.scripts.plot_dias import entry_point

        runner = CliRunner()
        with pytest.raises(TypeError, match="Unknown program"):
            runner.invoke(
                entry_point,
                ["-p", "xtb", "-a", "1", "-b", "2"],
                catch_exceptions=False,
            )

    def test_requires_atom_numbers(self):
        from chemsmart.scripts.plot_dias import entry_point

        runner = CliRunner()
        result = runner.invoke(entry_point, ["-p", "gaussian"])
        assert result.exit_code != 0


class TestStructureFilterScript:
    def test_basic_grouping_run(self):
        from chemsmart.scripts.structure_filter import entry_point

        runner = CliRunner()
        with runner.isolated_filesystem():
            with (
                patch(
                    "chemsmart.scripts.structure_filter.BaseFolder"
                ) as mock_folder_cls,
                patch(
                    "chemsmart.scripts.structure_filter.Molecule.from_filepath"
                ) as mock_from_filepath,
                patch(
                    "chemsmart.scripts.structure_filter.StructureGrouperFactory"
                ) as mock_factory,
            ):
                mock_folder = MagicMock()
                mock_folder.get_all_files_in_current_folder_by_suffix.return_value = [
                    "mol_c1.log",
                    "mol_c2.log",
                ]
                mock_folder_cls.return_value = mock_folder
                mock_from_filepath.return_value = MagicMock()

                mock_grouper = MagicMock()
                mock_grouper.group.return_value = ([[0], [1]], [[0], [1]])
                mock_grouper.unique.return_value = [MagicMock(), MagicMock()]
                mock_factory.create.return_value = mock_grouper

                result = runner.invoke(
                    entry_point,
                    ["-d", ".", "-t", "log"],
                    catch_exceptions=False,
                )

            assert result.exit_code == 0, result.output
            mock_factory.create.assert_called_once()

    def test_grouper_creation_error_propagates(self):
        from chemsmart.scripts.structure_filter import entry_point

        runner = CliRunner()
        with runner.isolated_filesystem():
            with (
                patch(
                    "chemsmart.scripts.structure_filter.BaseFolder"
                ) as mock_folder_cls,
                patch(
                    "chemsmart.scripts.structure_filter.Molecule.from_filepath"
                ) as mock_from_filepath,
                patch(
                    "chemsmart.scripts.structure_filter.StructureGrouperFactory"
                ) as mock_factory,
            ):
                mock_folder = MagicMock()
                mock_folder.get_all_files_in_current_folder_by_suffix.return_value = [
                    "mol_c1.log",
                ]
                mock_folder_cls.return_value = mock_folder
                mock_from_filepath.return_value = MagicMock()
                mock_factory.create.side_effect = RuntimeError("boom")

                with pytest.raises(RuntimeError, match="boom"):
                    runner.invoke(
                        entry_point,
                        ["-d", ".", "-t", "log"],
                        catch_exceptions=False,
                    )


class TestGenerateIsotopeData:
    def test_parse_isotope_file_with_natural_abundance(self, tmp_path):
        from chemsmart.scripts.generate_isotope_data import (
            parse_isotope_file,
        )

        isotope_txt = tmp_path / "isotopes.txt"
        isotope_txt.write_text(textwrap.dedent("""\
                Atomic Number = 1
                Mass Number = 1
                Relative Atomic Mass = 1.00782503(1)
                Isotopic Composition = 0.999885(70)
                Atomic Number = 1
                Mass Number = 2
                Relative Atomic Mass = 2.01410178(1)
                Isotopic Composition = 0.000115(70)
                """))

        isotopes = parse_isotope_file(str(isotope_txt))

        assert 1 in isotopes
        assert isotopes[1]["most_abundant"]["mass_number"] == 1
        assert isotopes[1]["weighted_atomic_mass"] == pytest.approx(
            1.00782503, rel=1e-3
        )

    def test_parse_isotope_file_radioactive_element_fallback(self, tmp_path):
        from chemsmart.scripts.generate_isotope_data import (
            parse_isotope_file,
        )

        # Technetium (Z=43): no natural abundance, falls back to most
        # stable mass number 98 per `most_stable_mass_numbers`.
        isotope_txt = tmp_path / "isotopes_tc.txt"
        isotope_txt.write_text(textwrap.dedent("""\
                Atomic Number = 43
                Mass Number = 98
                Relative Atomic Mass = 97.9072124(1)
                Isotopic Composition =
                """))

        isotopes = parse_isotope_file(str(isotope_txt))

        assert isotopes[43]["most_abundant"]["mass_number"] == 98
        assert isotopes[43]["weighted_atomic_mass"] == pytest.approx(
            97.9072124, rel=1e-6
        )


class TestGetThermochemistryScript:
    def _make_thermo_mock(self):
        thermo = MagicMock()
        thermo.electronic_energy = -100.0
        thermo.zero_point_energy = 0.05
        thermo.enthalpy = -99.9
        thermo.qrrho_enthalpy = -99.9
        thermo.entropy_times_temperature = -0.02
        thermo.qrrho_entropy_times_temperature = -0.02
        thermo.gibbs_free_energy = -99.92
        thermo.qrrho_gibbs_free_energy = -99.92
        thermo.qrrho_gibbs_free_energy_qh = -99.92
        thermo.qrrho_gibbs_free_energy_qs = -99.92
        thermo.imaginary_frequencies = []
        thermo.jobtype = "opt"
        thermo.vibrational_frequencies = [100.0]
        return thermo

    def test_default_run_writes_results(self, tmp_path):
        from chemsmart.scripts.get_thermochemistry import entry_point

        log_file = tmp_path / "mol.log"
        log_file.write_text("dummy")

        runner = CliRunner()
        with runner.isolated_filesystem():
            with patch(
                "chemsmart.scripts.get_thermochemistry.Thermochemistry"
            ) as mock_cls:
                mock_cls.return_value = self._make_thermo_mock()
                result = runner.invoke(
                    entry_point,
                    ["-f", str(log_file)],
                    catch_exceptions=False,
                )

            assert result.exit_code == 0, result.output
            mock_cls.assert_called_once()

    def test_quasi_rrho_flag_run(self, tmp_path):
        from chemsmart.scripts.get_thermochemistry import entry_point

        log_file = tmp_path / "mol.log"
        log_file.write_text("dummy")

        runner = CliRunner()
        with runner.isolated_filesystem():
            with patch(
                "chemsmart.scripts.get_thermochemistry.Thermochemistry"
            ) as mock_cls:
                mock_cls.return_value = self._make_thermo_mock()
                result = runner.invoke(
                    entry_point,
                    ["-f", str(log_file), "-q"],
                    catch_exceptions=False,
                )

            assert result.exit_code == 0, result.output

    def test_directory_without_filetype_asserts(self, tmp_path):
        from chemsmart.scripts.get_thermochemistry import entry_point

        runner = CliRunner()
        with runner.isolated_filesystem():
            with pytest.raises(AssertionError):
                runner.invoke(
                    entry_point,
                    ["-d", str(tmp_path)],
                    catch_exceptions=False,
                )

    def test_unsupported_file_extension_logs_error(self, tmp_path):
        from chemsmart.scripts.get_thermochemistry import entry_point

        bad_file = tmp_path / "mol.txt"
        bad_file.write_text("dummy")

        runner = CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(
                entry_point,
                ["-f", str(bad_file)],
                catch_exceptions=False,
            )

        assert result.exit_code == 0, result.output
