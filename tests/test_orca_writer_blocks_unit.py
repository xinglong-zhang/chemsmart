"""
Direct unit tests for private block-writing methods of
:class:`ORCAInputWriter` in ``chemsmart.jobs.orca.writer`` that are not
exercised by the higher-level, file-comparison-based tests in
``test_ORCAWriter.py``.

Covers ``_write_mdci_block``, ``_write_elprop_block``,
``_write_qmmm_block``, ``_write_modred``/``_write_modred_if_list``/
``_write_modred_if_dict``, ``_write_hessian_block_for_ts``,
``_write_irc_block_for_irc``, the ``neb_block`` property and
``_write_neb_block_for_neb``, ``_write_constrained_atoms``, and
``_write_charge_and_multiplicity``.

Each writer method is called directly against an ``io.StringIO`` target
using a lightweight fake job (only ``settings``/``molecule``/``jobrunner``
are needed), avoiding the full CLI/project-YAML/file-comparison machinery
used elsewhere.
"""

import io
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest

from chemsmart.jobs.orca.settings import (
    ORCAIRCJobSettings,
    ORCAJobSettings,
    ORCANEBJobSettings,
    ORCAQMMMJobSettings,
    ORCATSJobSettings,
)
from chemsmart.jobs.orca.writer import ORCAInputWriter


def _writer(settings, molecule=None):
    job = SimpleNamespace(
        settings=settings,
        molecule=(
            molecule if molecule is not None else MagicMock(frozen_atoms=None)
        ),
        jobrunner=MagicMock(num_cores=4, mem_gb=8),
        label="test_job",
        folder="/tmp",
    )
    return ORCAInputWriter(job=job)


class TestMdciBlock:
    def test_no_cutoff_writes_nothing(self):
        writer = _writer(ORCAJobSettings())
        f = io.StringIO()
        writer._write_mdci_block(f)
        assert f.getvalue() == ""

    def test_invalid_cutoff_raises(self):
        writer = _writer(ORCAJobSettings(mdci_cutoff="ultra"))
        with pytest.raises(AssertionError):
            writer._write_mdci_block(io.StringIO())

    @pytest.mark.parametrize(
        "cutoff,expected_pairs",
        [
            ("loose", "TCutPairs 1e-3"),
            ("normal", "TCutPairs 1e-4"),
            ("tight", "TCutPairs 1e-5"),
        ],
    )
    def test_cutoff_levels(self, cutoff, expected_pairs):
        writer = _writer(ORCAJobSettings(mdci_cutoff=cutoff))
        f = io.StringIO()
        writer._write_mdci_block(f)
        out = f.getvalue()
        assert out.startswith("%mdci\n")
        assert expected_pairs in out
        assert out.endswith("end\n")

    def test_invalid_density_raises(self):
        writer = _writer(
            ORCAJobSettings(mdci_cutoff="normal", mdci_density="bogus")
        )
        with pytest.raises(AssertionError):
            writer._write_mdci_block(io.StringIO())

    @pytest.mark.parametrize(
        "density,expected",
        [
            ("none", "Density None"),
            ("unrelaxed", "Density Unrelaxed"),
            ("relaxed", "Density Relaxed"),
        ],
    )
    def test_density_levels(self, density, expected):
        writer = _writer(
            ORCAJobSettings(mdci_cutoff="normal", mdci_density=density)
        )
        f = io.StringIO()
        writer._write_mdci_block(f)
        assert expected in f.getvalue()


class TestElpropBlock:
    def test_neither_dipole_nor_quadrupole_writes_nothing(self):
        writer = _writer(ORCAJobSettings())
        f = io.StringIO()
        writer._write_elprop_block(f)
        assert f.getvalue() == ""

    def test_dipole_true_quadrupole_false(self):
        writer = _writer(ORCAJobSettings(dipole=True))
        f = io.StringIO()
        writer._write_elprop_block(f)
        out = f.getvalue()
        assert "Dipole True" in out
        assert "Quadrupole False" in out

    def test_quadrupole_true_dipole_false(self):
        writer = _writer(ORCAJobSettings(quadrupole=True))
        f = io.StringIO()
        writer._write_elprop_block(f)
        out = f.getvalue()
        assert "Dipole False" in out
        assert "Quadrupole True" in out


class TestQmmmBlockDispatch:
    def test_non_qmmm_settings_writes_nothing(self):
        writer = _writer(ORCAJobSettings())
        f = io.StringIO()
        writer._write_qmmm_block(f)
        assert f.getvalue() == ""

    def test_qmmm_settings_delegates_to_settings_qmmm_block(self):
        settings = ORCAQMMMJobSettings(
            jobtype="QMMM",
            low_level_method="ff.prms",
            charge_total=0,
            mult_total=1,
        )
        writer = _writer(settings)
        f = io.StringIO()
        writer._write_qmmm_block(f)
        assert f.getvalue() == f"{settings.qmmm_block}\n"


class TestModredBlock:
    def test_list_format_writes_constraints(self):
        writer = _writer(ORCAJobSettings())
        f = io.StringIO()
        writer._write_modred_if_list(f, [[1, 2], [3, 4, 5]])
        out = f.getvalue()
        assert out.startswith("  Constraints\n")
        assert "{B 0 1 C}" in out
        assert "{A 2 3 4 C}" in out
        assert out.rstrip().endswith("end")

    def test_dict_format_scan_only(self):
        writer = _writer(ORCAJobSettings())
        f = io.StringIO()
        modred = {
            "coords": [[1, 2]],
            "num_steps": [10],
            "dist_start": [1.5],
            "dist_end": [3.5],
        }
        writer._write_modred_if_dict(f, modred)
        out = f.getvalue()
        assert "Constraints" not in out
        assert "  Scan\n" in out
        assert "B 0 1 = 1.5, 3.5, 10" in out
        assert "bond distance" in out

    def test_dict_format_with_constrained_coordinates(self):
        writer = _writer(ORCAJobSettings())
        f = io.StringIO()
        modred = {
            "constrained_coordinates": [[5, 6]],
            "coords": [[1, 2, 3]],
            "num_steps": [5],
            "dist_start": [90.0],
            "dist_end": [120.0],
        }
        writer._write_modred_if_dict(f, modred)
        out = f.getvalue()
        assert "Constraints" in out
        assert "{B 4 5 C}" in out
        assert "A 0 1 2 = 90.0, 120.0, 5" in out
        assert "angle" in out

    def test_dihedral_scan_variable_name(self):
        writer = _writer(ORCAJobSettings())
        f = io.StringIO()
        modred = {
            "coords": [[1, 2, 3, 4]],
            "num_steps": [8],
            "dist_start": [0.0],
            "dist_end": [180.0],
        }
        writer._write_modred_if_dict(f, modred)
        out = f.getvalue()
        assert "dihedral" in out

    def test_write_modred_dispatches_by_type(self):
        writer = _writer(ORCAJobSettings())
        f_list = io.StringIO()
        writer._write_modred(f_list, [1, 2])
        assert "Constraints" in f_list.getvalue()

        f_dict = io.StringIO()
        writer._write_modred(
            f_dict,
            {
                "coords": [[1, 2]],
                "num_steps": [4],
                "dist_start": [1.0],
                "dist_end": [2.0],
            },
        )
        assert "Scan" in f_dict.getvalue()

    def test_modred_block_wraps_in_geom_section(self):
        writer = _writer(ORCAJobSettings(modred=[1, 2]))
        f = io.StringIO()
        writer._write_modred_block(f)
        out = f.getvalue()
        assert out.startswith("%geom\n")
        assert out.rstrip().endswith("end")

    def test_no_modred_writes_nothing(self):
        writer = _writer(ORCAJobSettings())
        f = io.StringIO()
        writer._write_modred_block(f)
        assert f.getvalue() == ""


class TestHessianBlockForTS:
    def test_minimal_ts_settings(self):
        settings = ORCATSJobSettings()
        writer = _writer(settings)
        f = io.StringIO()
        writer._write_hessian_block_for_ts(f)
        out = f.getvalue()
        assert out.startswith("%geom\n")
        assert "Calc_Hess True" in out
        assert f"NumHess {settings.numhess}" in out
        assert f"Recalc_Hess {settings.recalc_hess}" in out
        assert out.rstrip().endswith("end")

    def test_inhess_requires_existing_file(self, tmp_path):
        hess_file = tmp_path / "guess.hess"
        hess_file.write_text("fake hessian")
        settings = ORCATSJobSettings(
            inhess=True, inhess_filename=str(hess_file)
        )
        writer = _writer(settings)
        f = io.StringIO()
        writer._write_hessian_block_for_ts(f)
        out = f.getvalue()
        assert "InHess Read" in out
        assert f'InHessName "{hess_file}"' in out

    def test_inhess_without_filename_raises(self):
        settings = ORCATSJobSettings(inhess=True)
        writer = _writer(settings)
        with pytest.raises(AssertionError, match="No Hessian file"):
            writer._write_hessian_block_for_ts(io.StringIO())

    def test_inhess_missing_file_raises(self):
        settings = ORCATSJobSettings(
            inhess=True, inhess_filename="/no/such/file.hess"
        )
        writer = _writer(settings)
        with pytest.raises(AssertionError, match="is not found"):
            writer._write_hessian_block_for_ts(io.StringIO())

    def test_hybrid_hess_converts_to_0_indexed(self):
        settings = ORCATSJobSettings(
            hybrid_hess=True, hybrid_hess_atoms=[1, 2, 5]
        )
        writer = _writer(settings)
        f = io.StringIO()
        writer._write_hessian_block_for_ts(f)
        assert "Hybrid_Hess {0 1 4} end" in f.getvalue()

    def test_hybrid_hess_without_atoms_raises(self):
        settings = ORCATSJobSettings(hybrid_hess=True)
        writer = _writer(settings)
        with pytest.raises(AssertionError, match="No atoms"):
            writer._write_hessian_block_for_ts(io.StringIO())

    def test_trust_radius_negative_is_fixed(self):
        settings = ORCATSJobSettings(trust_radius=-0.3)
        writer = _writer(settings)
        f = io.StringIO()
        writer._write_hessian_block_for_ts(f)
        out = f.getvalue()
        assert "Trust -0.3" in out
        assert "fixed trust radius" in out

    def test_trust_radius_positive_is_update(self):
        settings = ORCATSJobSettings(trust_radius=0.3)
        writer = _writer(settings)
        f = io.StringIO()
        writer._write_hessian_block_for_ts(f)
        out = f.getvalue()
        assert "Trust 0.3" in out
        assert "trust radius update" in out

    def test_scants_requires_modred(self):
        settings = ORCATSJobSettings(tssearch_type="scants")
        writer = _writer(settings)
        with pytest.raises(AssertionError, match="No modred"):
            writer._write_hessian_block_for_ts(io.StringIO())

    def test_scants_writes_scan_block(self):
        settings = ORCATSJobSettings(
            tssearch_type="scants",
            scants_modred={
                "coords": [[1, 2]],
                "num_steps": [4],
                "dist_start": [1.0],
                "dist_end": [2.0],
            },
        )
        writer = _writer(settings)
        f = io.StringIO()
        writer._write_hessian_block_for_ts(f)
        assert "Scan" in f.getvalue()

    def test_full_scan_writes_flag(self):
        settings = ORCATSJobSettings(full_scan=True)
        writer = _writer(settings)
        f = io.StringIO()
        writer._write_hessian_block_for_ts(f)
        assert "fullScan True" in f.getvalue()


class TestHessianBlockDispatch:
    def test_non_ts_settings_no_hessian_block(self):
        writer = _writer(ORCAJobSettings())
        f = io.StringIO()
        writer._write_hessian_block(f)
        assert f.getvalue() == ""

    def test_ts_settings_dispatches(self):
        writer = _writer(ORCATSJobSettings())
        f = io.StringIO()
        writer._write_hessian_block(f)
        assert f.getvalue().startswith("%geom\n")


class TestIrcBlockForIrc:
    def test_no_irc_specific_options_writes_nothing(self):
        settings = ORCAIRCJobSettings()
        writer = _writer(settings)
        f = io.StringIO()
        writer._write_irc_block_for_irc(f)
        assert f.getvalue() == ""

    def test_maxiter_option_written(self):
        settings = ORCAIRCJobSettings(maxiter=20)
        writer = _writer(settings)
        f = io.StringIO()
        writer._write_irc_block_for_irc(f)
        out = f.getvalue()
        assert out.startswith("%irc\n")
        assert "maxiter 20" in out
        assert out.rstrip().endswith("end")

    def test_inithess_read_requires_hess_filename(self, tmp_path):
        hess_file = tmp_path / "start.hess"
        hess_file.write_text("fake hessian")
        settings = ORCAIRCJobSettings(
            inithess="read", hess_filename=str(hess_file)
        )
        writer = _writer(settings)
        f = io.StringIO()
        writer._write_irc_block_for_irc(f)
        out = f.getvalue()
        assert "inithess read" in out
        assert f'Hess_Filename "{hess_file}"' in out

    def test_inithess_read_without_hess_filename_raises(self):
        settings = ORCAIRCJobSettings(inithess="read")
        writer = _writer(settings)
        with pytest.raises(AssertionError, match="No Hessian file"):
            writer._write_irc_block_for_irc(io.StringIO())

    def test_inithess_calc_anfreq_no_filename_needed(self):
        settings = ORCAIRCJobSettings(inithess="calc_anfreq")
        writer = _writer(settings)
        f = io.StringIO()
        writer._write_irc_block_for_irc(f)
        assert "inithess calc_anfreq" in f.getvalue()

    def test_monitor_internals_requires_internal_modred(self):
        settings = ORCAIRCJobSettings(monitor_internals=True)
        writer = _writer(settings)
        with pytest.raises(AssertionError, match="internal modred"):
            writer._write_irc_block_for_irc(io.StringIO())

    def test_monitor_internals_writes_coordinates(self):
        settings = ORCAIRCJobSettings(
            monitor_internals=True, internal_modred=[[1, 2], [5, 6]]
        )
        writer = _writer(settings)
        f = io.StringIO()
        writer._write_irc_block_for_irc(f)
        out = f.getvalue()
        assert "True" in out
        assert "{ B 0 1 }" in out
        assert "{ B 4 5 }" in out

    @pytest.mark.parametrize(
        "key,expected_text",
        [
            ("adapt_scale_displ", "Adapt_Scale_Displ True"),
            ("sd_parabolicfit", "SD_ParabolicFit True"),
            ("interpolate_only", "Interpolate_only True"),
            ("do_sd_corr", "Do_SD_Corr True"),
            ("sd_corr_parabolicfit", "SD_Corr_ParabolicFit True"),
        ],
    )
    def test_boolean_flag_options(self, key, expected_text):
        settings = ORCAIRCJobSettings(**{key: True})
        writer = _writer(settings)
        f = io.StringIO()
        writer._write_irc_block_for_irc(f)
        assert expected_text in f.getvalue()

    @pytest.mark.parametrize(
        "key",
        [
            "adapt_scale_displ",
            "sd_parabolicfit",
            "interpolate_only",
            "do_sd_corr",
            "sd_corr_parabolicfit",
        ],
    )
    def test_boolean_flag_options_false_omitted(self, key):
        settings = ORCAIRCJobSettings(**{key: False})
        writer = _writer(settings)
        f = io.StringIO()
        writer._write_irc_block_for_irc(f)
        # False is not None, so the key is "set", but the specific
        # branches only write when the value is True.
        assert f.getvalue() == "%irc\nend\n"


class TestIrcBlockDispatch:
    def test_non_irc_settings_writes_nothing(self):
        writer = _writer(ORCAJobSettings())
        f = io.StringIO()
        writer._write_irc_block(f)
        assert f.getvalue() == ""

    def test_irc_settings_dispatches(self):
        writer = _writer(ORCAIRCJobSettings(maxiter=5))
        f = io.StringIO()
        writer._write_irc_block(f)
        assert "maxiter 5" in f.getvalue()


class TestNebBlock:
    def test_missing_nimages_raises(self):
        settings = ORCANEBJobSettings(ending_xyzfile="product.xyz")
        writer = _writer(settings)
        with pytest.raises(AssertionError, match="number of images"):
            _ = writer.neb_block

    def test_missing_geometry_raises(self):
        settings = ORCANEBJobSettings(nimages=8)
        writer = _writer(settings)
        with pytest.raises(AssertionError, match="valid input geometry"):
            _ = writer.neb_block

    def test_ending_xyzfile_basic_block(self):
        settings = ORCANEBJobSettings(
            nimages=8, ending_xyzfile="/some/path/product.xyz"
        )
        writer = _writer(settings)
        block = writer.neb_block
        assert block.startswith("%NEB")
        assert 'NEB_END_XYZFILE "product.xyz"' in block
        assert "NImages 8" in block
        assert "PREOPT_ENDS False" in block
        assert block.endswith("end")

    def test_preopt_ends_true(self):
        settings = ORCANEBJobSettings(
            nimages=8,
            ending_xyzfile="product.xyz",
            preopt_ends=True,
        )
        writer = _writer(settings)
        assert "PREOPT_ENDS True" in writer.neb_block

    def test_intermediate_xyzfile_included(self):
        settings = ORCANEBJobSettings(
            nimages=8,
            ending_xyzfile="product.xyz",
            intermediate_xyzfile="/some/ts_guess.xyz",
        )
        writer = _writer(settings)
        assert 'NEB_TS_XYZFILE "ts_guess.xyz"' in writer.neb_block

    def test_restart_file_takes_precedence_and_skips_preopt(self):
        settings = ORCANEBJobSettings(
            nimages=8,
            ending_xyzfile="product.xyz",
            restarting_xyzfile="/some/restart.allxyz",
        )
        writer = _writer(settings)
        block = writer.neb_block
        assert 'Restart_ALLXYZFile "restart.allxyz"' in block
        assert "NEB_END_XYZFILE" not in block
        assert "PREOPT_ENDS" not in block

    def test_write_neb_block_for_neb_matches_property(self):
        settings = ORCANEBJobSettings(nimages=8, ending_xyzfile="product.xyz")
        writer = _writer(settings)
        f = io.StringIO()
        writer._write_neb_block_for_neb(f)
        written = f.getvalue().rstrip("\n")
        assert written == writer.neb_block


class TestNebBlockDispatch:
    def test_non_neb_settings_writes_nothing(self):
        writer = _writer(ORCAJobSettings())
        f = io.StringIO()
        writer._write_neb_block(f)
        assert f.getvalue() == ""

    def test_neb_settings_dispatches(self):
        settings = ORCANEBJobSettings(nimages=4, ending_xyzfile="product.xyz")
        writer = _writer(settings)
        f = io.StringIO()
        writer._write_neb_block(f)
        assert f.getvalue().startswith("%NEB\n")


class TestConstrainedAtoms:
    def test_no_frozen_atoms_writes_nothing(self):
        molecule = MagicMock(frozen_atoms=None)
        writer = _writer(ORCAJobSettings(), molecule=molecule)
        f = io.StringIO()
        writer._write_constrained_atoms(f)
        assert f.getvalue() == ""

    def test_frozen_atoms_written_0_indexed(self):
        molecule = MagicMock(frozen_atoms=[-1, 0, -1])
        writer = _writer(ORCAJobSettings(), molecule=molecule)
        f = io.StringIO()
        writer._write_constrained_atoms(f)
        out = f.getvalue()
        assert out.startswith("%geom\n")
        assert "{ C 0 C }" in out
        assert "{ C 2 C }" in out
        assert "{ C 1 C }" not in out
        assert out.rstrip().endswith("end")

    def test_invert_constraints_flag(self):
        molecule = MagicMock(frozen_atoms=[-1])
        writer = _writer(
            ORCAJobSettings(invert_constraints=True), molecule=molecule
        )
        f = io.StringIO()
        writer._write_constrained_atoms(f)
        assert "InvertConstraints True" in f.getvalue()


class TestChargeAndMultiplicity:
    def test_direct_charge_and_multiplicity(self):
        writer = _writer(ORCAJobSettings(charge=1, multiplicity=2))
        f = io.StringIO()
        writer._write_charge_and_multiplicity(f)
        assert f.getvalue() == "* xyz 1 2\n"

    def test_missing_both_raises(self):
        writer = _writer(ORCAJobSettings())
        with pytest.raises(AssertionError, match="Charge and multiplicity"):
            writer._write_charge_and_multiplicity(io.StringIO())

    def test_falls_back_to_qmmm_total_fields(self):
        settings = ORCAQMMMJobSettings(
            jobtype="MOL-CRYSTAL-QMMM",
            n_unit_cell_atoms=12,
            low_level_method="ff.prms",
        )
        # charge/multiplicity are None on this settings object; total
        # fields were also not supplied, so patch them in directly to
        # exercise the QMMM fallback path without invoking
        # check_crystal_qmmm (which forbids setting multiplicity here).
        settings.charge_total = 0
        settings.mult_total = 1
        writer = _writer(settings)
        f = io.StringIO()
        writer._write_charge_and_multiplicity(f)
        assert f.getvalue() == "* xyz 0 1\n"
        # settings object is updated in place for downstream use
        assert settings.charge == 0
        assert settings.multiplicity == 1
