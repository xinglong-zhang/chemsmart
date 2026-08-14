"""
Direct unit tests for ``ORCAQMMMJobSettings`` and ``ORCANEBJobSettings``
in ``chemsmart.jobs.orca.settings``.

These two classes are not covered by ``test_orca_settings_unit.py`` (which
focuses on ``ORCAJobSettings``, ``ORCApKaJobSettings``, ``ORCATSJobSettings``,
and ``ORCAIRCJobSettings``). They implement ORCA's multiscale QM/MM and
Nudged Elastic Band calculation configuration respectively, including
charge/multiplicity fallback logic, intermediate-layer validation, atom
partition-string formatting, and full ``%qmmm`` block generation.
"""

import pytest

from chemsmart.jobs.orca.settings import (
    ORCANEBJobSettings,
    ORCAQMMMJobSettings,
)


class TestQMMMChargeMultiplicityFallback:
    def test_high_level_takes_priority(self):
        s = ORCAQMMMJobSettings(
            jobtype="QMMM",
            low_level_method="ff.prms",
            charge_high=0,
            mult_high=1,
            charge_intermediate=1,
            mult_intermediate=2,
            charge_total=2,
            mult_total=3,
        )
        assert s.charge == 0
        assert s.multiplicity == 1

    def test_intermediate_used_when_no_high(self):
        s = ORCAQMMMJobSettings(
            jobtype="QM/QM2",
            intermediate_level_method="XTB2",
            intermediate_level_atoms=[1, 2],
            charge_intermediate=1,
            mult_intermediate=2,
            charge_total=2,
            mult_total=3,
        )
        assert s.charge == 1
        assert s.multiplicity == 2

    def test_total_used_when_no_high_or_intermediate(self):
        s = ORCAQMMMJobSettings(
            jobtype="QMMM",
            low_level_method="ff.prms",
            charge_total=2,
            mult_total=3,
        )
        assert s.charge == 2
        assert s.multiplicity == 3

    def test_none_when_nothing_specified(self):
        s = ORCAQMMMJobSettings(jobtype="QMMM", low_level_method="ff.prms")
        assert s.charge is None
        assert s.multiplicity is None

    def test_functional_and_basis_aliased_from_high_level(self):
        s = ORCAQMMMJobSettings(
            jobtype="QMMM",
            high_level_functional="B3LYP",
            high_level_basis="def2-SVP",
            low_level_method="ff.prms",
            charge_total=0,
            mult_total=1,
        )
        assert s.functional == "B3LYP"
        assert s.basis == "def2-SVP"


class TestQMMMIntermediateValidation:
    def test_qm_qm2_without_intermediate_params_raises(self):
        with pytest.raises(ValueError, match="requires QM2"):
            ORCAQMMMJobSettings(
                jobtype="QM/QM2",
                high_level_functional="B3LYP",
                high_level_basis="def2-SVP",
                charge_high=0,
                mult_high=1,
            )

    def test_qm_qm2_mm_without_intermediate_params_raises(self):
        with pytest.raises(ValueError, match="requires QM2"):
            ORCAQMMMJobSettings(
                jobtype="QM/QM2/MM",
                low_level_method="ff.prms",
                charge_total=0,
                mult_total=1,
            )

    def test_intermediate_method_without_atoms_but_with_low_level_raises(
        self,
    ):
        with pytest.raises(ValueError, match="intermediate_level_atoms"):
            ORCAQMMMJobSettings(
                jobtype="QM/QM2/MM",
                intermediate_level_method="XTB2",
                low_level_method="ff.prms",
                charge_total=0,
                mult_total=1,
            )

    def test_intermediate_method_without_atoms_or_low_level_is_allowed(self):
        # The atoms-required guard is only triggered when a low-level
        # method is also supplied; a pure QM/QM2 job with just an
        # intermediate method and no MM layer sails through.
        s = ORCAQMMMJobSettings(
            jobtype="QM/QM2",
            high_level_functional="B3LYP",
            high_level_basis="def2-SVP",
            intermediate_level_method="XTB2",
            charge_high=0,
            mult_high=1,
        )
        assert s.intermediate_level_atoms is None

    def test_intermediate_params_for_non_qm2_jobtype_with_low_level_raises(
        self,
    ):
        with pytest.raises(ValueError, match="does not require QM2"):
            ORCAQMMMJobSettings(
                jobtype="QMMM",
                intermediate_level_functional="B3LYP",
                intermediate_level_basis="def2-SVP",
                intermediate_level_atoms=[1, 2],
                low_level_method="ff.prms",
                charge_total=0,
                mult_total=1,
            )

    def test_intermediate_params_for_non_qm2_jobtype_without_low_level_is_silently_accepted(
        self,
    ):
        """Documents a validation gap (see BUGS_FOUND.md #11): the
        "does not require QM2 but intermediate params were provided"
        check is only reached when ``low_level_method`` is also set, so a
        plain ``QMMM`` job can carry a fully-specified (and meaningless)
        intermediate QM2 layer without any error being raised."""
        s = ORCAQMMMJobSettings(
            jobtype="QMMM",
            intermediate_level_functional="B3LYP",
            intermediate_level_basis="def2-SVP",
            charge_total=0,
            mult_total=1,
        )
        assert s.intermediate_level_functional == "B3LYP"

    def test_re_init_and_validate_recomputes_charge_and_reruns_validation(
        self,
    ):
        s = ORCAQMMMJobSettings(
            jobtype="QM/QM2",
            high_level_functional="B3LYP",
            high_level_basis="def2-SVP",
            intermediate_level_method="XTB2",
            intermediate_level_atoms=[1, 2],
            charge_high=0,
            mult_high=1,
        )
        s.charge_high = 1
        s.mult_high = 3
        s.high_level_functional = "PBE0"
        s.re_init_and_validate()
        assert s.charge == 1
        assert s.multiplicity == 3
        assert s.functional == "PBE0"

        s.intermediate_level_method = None
        s.intermediate_level_functional = None
        s.intermediate_level_basis = None
        with pytest.raises(ValueError, match="requires QM2"):
            s.re_init_and_validate()


class TestQMMMLevelOfTheoryString:
    def test_additive_qmmm_route_string(self):
        s = ORCAQMMMJobSettings(
            jobtype="QMMM",
            parent_jobtype="opt",
            high_level_functional="B3LYP",
            high_level_basis="def2-SVP",
            low_level_method="ff.prms",
            charge_total=0,
            mult_total=1,
        )
        route = s._get_level_of_theory_string()
        assert route.startswith("! Opt QMMM")
        assert "B3LYP def2-SVP" in route

    def test_qm_qm2_route_string(self):
        s = ORCAQMMMJobSettings(
            jobtype="QM/QM2",
            parent_jobtype="sp",
            high_level_functional="B3LYP",
            high_level_basis="def2-SVP",
            intermediate_level_method="XTB2",
            intermediate_level_atoms=[1, 2],
            charge_high=0,
            mult_high=1,
        )
        route = s._get_level_of_theory_string()
        assert "QM/XTB2" in route
        assert "MM" not in route.replace("QM/XTB2", "")

    def test_qm_qm2_mm_route_string_includes_low_level(self):
        s = ORCAQMMMJobSettings(
            jobtype="QM/QM2/MM",
            parent_jobtype="ts",
            high_level_functional="B3LYP",
            high_level_basis="def2-SVP",
            intermediate_level_method="XTB2",
            intermediate_level_atoms=[1, 2],
            low_level_method="ff.prms",
            charge_total=0,
            mult_total=1,
        )
        route = s._get_level_of_theory_string()
        assert route.startswith("! OptTS")
        assert "QM/XTB2/MM" in route

    def test_irc_parent_jobtype(self):
        s = ORCAQMMMJobSettings(
            jobtype="QMMM",
            parent_jobtype="irc",
            high_level_functional="B3LYP",
            high_level_basis="def2-SVP",
            low_level_method="ff.prms",
            freq=True,
            charge_total=0,
            mult_total=1,
        )
        route = s._get_level_of_theory_string()
        assert route.startswith("! IRC Freq")

    def test_crystal_jobtype_route_string_is_just_the_jobtype(self):
        s = ORCAQMMMJobSettings(
            jobtype="MOL-CRYSTAL-QMMM",
            n_unit_cell_atoms=12,
            low_level_method="ff.prms",
        )
        assert s._get_level_of_theory_string() == "! MOL-CRYSTAL-QMMM"

    def test_qmmm_route_string_property_delegates(self):
        s = ORCAQMMMJobSettings(
            jobtype="MOL-CRYSTAL-QMMM",
            n_unit_cell_atoms=12,
            low_level_method="ff.prms",
        )
        assert s.qmmm_route_string == s._get_level_of_theory_string()


class TestValidateAndAssignLevel:
    def test_functional_and_basis_and_builtin_together_raises(self):
        s = ORCAQMMMJobSettings(jobtype="QMMM", low_level_method="ff.prms")
        with pytest.raises(ValueError, match="not both"):
            s.validate_and_assign_level(
                "B3LYP", "def2-SVP", "XTB2", level_name="high_level"
            )

    def test_builtin_method_for_intermediate_level(self):
        s = ORCAQMMMJobSettings(jobtype="QMMM", low_level_method="ff.prms")
        result = s.validate_and_assign_level(
            None, None, "XTB2", level_name="intermediate_level"
        )
        assert result == "XTB2"

    def test_unknown_builtin_method_for_intermediate_level_returns_empty(
        self,
    ):
        s = ORCAQMMMJobSettings(jobtype="QMMM", low_level_method="ff.prms")
        result = s.validate_and_assign_level(
            None, None, "NOT-A-METHOD", level_name="intermediate_level"
        )
        assert result == ""

    def test_functional_and_basis_for_intermediate_level_is_qm2(self):
        s = ORCAQMMMJobSettings(jobtype="QMMM", low_level_method="ff.prms")
        result = s.validate_and_assign_level(
            "B3LYP", "def2-SVP", None, level_name="intermediate_level"
        )
        assert result == "QM2"

    def test_functional_and_basis_for_high_level(self):
        s = ORCAQMMMJobSettings(jobtype="QMMM", low_level_method="ff.prms")
        result = s.validate_and_assign_level(
            "B3LYP", "def2-SVP", None, level_name="high_level"
        )
        assert result == "B3LYP def2-SVP"

    def test_low_level_always_mm(self):
        s = ORCAQMMMJobSettings(jobtype="QMMM", low_level_method="ff.prms")
        result = s.validate_and_assign_level(
            None, None, None, level_name="low_level"
        )
        assert result == "MM"


class TestCheckCrystalQmmm:
    def test_mol_crystal_qmmm_requires_no_multiplicity(self):
        s = ORCAQMMMJobSettings(
            jobtype="MOL-CRYSTAL-QMMM",
            n_unit_cell_atoms=12,
            low_level_method="ff.prms",
            mult_total=1,
        )
        with pytest.raises(AssertionError, match="Multiplicity"):
            s.check_crystal_qmmm()

    def test_mol_crystal_qmmm_sets_multiplicity_to_zero(self):
        s = ORCAQMMMJobSettings(
            jobtype="MOL-CRYSTAL-QMMM",
            n_unit_cell_atoms=12,
            low_level_method="ff.prms",
        )
        s.check_crystal_qmmm()
        assert s.multiplicity == 0

    def test_mol_crystal_qmmm_requires_n_unit_cell_atoms(self):
        s = ORCAQMMMJobSettings(
            jobtype="MOL-CRYSTAL-QMMM", low_level_method="ff.prms"
        )
        with pytest.raises(AssertionError, match="number of atoms"):
            s.check_crystal_qmmm()

    def test_ionic_crystal_qmmm_requires_ecp_layer_ecp(self):
        s = ORCAQMMMJobSettings(
            jobtype="IONIC-CRYSTAL-QMMM", low_level_method="ff.prms"
        )
        with pytest.raises(AssertionError, match="cECPs"):
            s.check_crystal_qmmm()

    def test_ionic_crystal_qmmm_rejects_n_unit_cell_atoms(self):
        s = ORCAQMMMJobSettings(
            jobtype="IONIC-CRYSTAL-QMMM",
            ecp_layer_ecp="def2-ECP",
            n_unit_cell_atoms=5,
            low_level_method="ff.prms",
        )
        with pytest.raises(AssertionError, match="only applicable"):
            s.check_crystal_qmmm()

    def test_non_crystal_jobtype_is_a_no_op(self):
        s = ORCAQMMMJobSettings(
            jobtype="QMMM",
            low_level_method="ff.prms",
            charge_total=0,
            mult_total=1,
        )
        s.check_crystal_qmmm()
        assert s.multiplicity == 1


class TestFormattedPartitionStrings:
    def _settings(self):
        return ORCAQMMMJobSettings(jobtype="QMMM", low_level_method="ff.prms")

    def test_none_returns_none(self):
        s = self._settings()
        assert s._get_formatted_partition_strings(None) is None

    def test_empty_list_returns_none(self):
        s = self._settings()
        assert s._get_formatted_partition_strings([]) is None

    def test_list_of_ints_compresses_ranges_and_shifts_to_zero_index(self):
        s = self._settings()
        result = s._get_formatted_partition_strings([1, 2, 3, 5, 7, 8, 9])
        assert result == "0:2 4 6:8"

    def test_string_range_is_parsed_and_shifted(self):
        s = self._settings()
        result = s._get_formatted_partition_strings("1-15,37,39")
        assert result == "0:14 36 38"

    def test_duplicates_are_deduplicated(self):
        s = self._settings()
        result = s._get_formatted_partition_strings([1, 1, 2, 2, 3])
        assert result == "0:2"

    def test_unsupported_type_raises_type_error(self):
        s = self._settings()
        with pytest.raises(TypeError):
            s._get_formatted_partition_strings(object())


class TestQmmmBlockGeneration:
    def test_additive_qmmm_full_block(self):
        s = ORCAQMMMJobSettings(
            jobtype="QMMM",
            high_level_atoms=[1, 2, 3],
            low_level_method="system.prms",
            charge_total=0,
            mult_total=1,
            embedding_type="electronic",
            delete_la_double_counting=True,
        )
        block = s._write_qmmm_block()
        assert block.startswith("%qmmm\n")
        assert block.endswith("end\n")
        assert "QMAtoms {0:2} end" in block
        assert "Charge_Total 0" in block
        assert "Mult_Total 1" in block
        assert 'ORCAFFFilename "system.prms"' in block
        assert "Delete_LA_Double_Counting true" in block
        assert "Embedding Electronic" in block

    def test_qmmm_block_property_delegates(self):
        s = ORCAQMMMJobSettings(
            jobtype="QMMM",
            low_level_method="system.prms",
            charge_total=0,
            mult_total=1,
        )
        assert s.qmmm_block == s._write_qmmm_block()

    def test_qmmm_requires_force_field_file(self):
        s = ORCAQMMMJobSettings(jobtype="QMMM", charge_total=0, mult_total=1)
        with pytest.raises(AssertionError, match="Force field file"):
            s._write_qmmm_block()

    def test_intermediate_atoms_uses_intermediate_charge_line(self):
        s = ORCAQMMMJobSettings(
            jobtype="QM/QM2",
            intermediate_level_method="XTB2",
            intermediate_level_atoms=[1, 2],
            charge_intermediate=1,
            mult_intermediate=2,
        )
        block = s._write_qmmm_block()
        assert "Charge_Intermediate 1" in block
        assert "Mult_Intermediate 2" in block

    def test_none_charge_lines_are_omitted(self):
        s = ORCAQMMMJobSettings(jobtype="QMMM", low_level_method="ff.prms")
        block = s._write_qmmm_block()
        assert "Charge_Total None" not in block
        assert "None" not in block

    def test_custom_intermediate_method_file_written_when_not_builtin(self):
        s = ORCAQMMMJobSettings(
            jobtype="QM/QM2",
            intermediate_level_method="my_custom_method",
            intermediate_level_atoms=[1, 2],
            charge_intermediate=0,
            mult_intermediate=1,
        )
        block = s._write_qmmm_block()
        assert 'QM2CustomFile "my_custom_method" end' in block

    def test_custom_functional_basis_written_when_no_builtin_method(self):
        s = ORCAQMMMJobSettings(
            jobtype="QM/QM2",
            intermediate_level_functional="B3LYP",
            intermediate_level_basis="def2-SVP",
            intermediate_level_atoms=[1, 2],
            charge_intermediate=0,
            mult_intermediate=1,
        )
        block = s._write_qmmm_block()
        assert 'QM2CUSTOMMETHOD "B3LYP" end' in block
        assert 'QM2CUSTOMBASIS "def2-SVP" end' in block

    def test_use_active_info_from_pbc_line(self):
        s = ORCAQMMMJobSettings(
            jobtype="QMMM",
            low_level_method="ff.prms",
            use_active_info_from_pbc=True,
        )
        block = s._write_qmmm_block()
        assert "Use_Active_InfoFromPDB true" in block

    def test_optregion_fixed_atoms_line(self):
        s = ORCAQMMMJobSettings(
            jobtype="QMMM",
            low_level_method="ff.prms",
            optregion_fixed_atoms=[1, 2, 3],
        )
        block = s._write_qmmm_block()
        assert "OptRegion_FixedAtoms {0:2} end" in block

    def test_crystal_subblock_included(self):
        s = ORCAQMMMJobSettings(
            jobtype="MOL-CRYSTAL-QMMM",
            n_unit_cell_atoms=12,
            low_level_method="ff.prms",
            conv_charges=False,
        )
        block = s._write_qmmm_block()
        assert "Conv_Charges False" in block
        assert "NumUnitCellAtoms 12" in block


class TestHBondLength:
    def _settings(self, **kwargs):
        return ORCAQMMMJobSettings(
            jobtype="QMMM", low_level_method="ff.prms", **kwargs
        )

    def test_dict_form(self):
        s = self._settings(
            high_level_h_bond_length={("C", "H"): 1.09, ("O", "H"): 0.98}
        )
        result = s._get_h_bond_length()
        assert "Dist_C_H 1.09" in result
        assert "Dist_O_H 0.98" in result

    def test_string_form_requires_existing_file(self, tmp_path):
        missing = tmp_path / "does_not_exist.txt"
        s = self._settings(high_level_h_bond_length=str(missing))
        with pytest.raises(AssertionError, match="does not exist"):
            s._get_h_bond_length()

    def test_string_form_with_existing_file(self, tmp_path):
        present = tmp_path / "dist.txt"
        present.write_text("Dist_C_HLA 1.09\n")
        s = self._settings(high_level_h_bond_length=str(present))
        result = s._get_h_bond_length()
        assert result == f'H_Dist_FileName "{present}"'


class TestEmbeddingType:
    def test_electronic(self):
        s = ORCAQMMMJobSettings(
            jobtype="QMMM",
            low_level_method="ff.prms",
            embedding_type="electronic",
        )
        assert s._get_embedding_type() == "Embedding Electronic"

    def test_mechanical(self):
        s = ORCAQMMMJobSettings(
            jobtype="QMMM",
            low_level_method="ff.prms",
            embedding_type="Mechanical",
        )
        assert s._get_embedding_type() == "Embedding Mechanical"

    def test_invalid_raises(self):
        s = ORCAQMMMJobSettings(
            jobtype="QMMM",
            low_level_method="ff.prms",
            embedding_type="bogus",
        )
        with pytest.raises(ValueError, match="Invalid embedding type"):
            s._get_embedding_type()


class TestQMMMEquality:
    def test_equal_settings_are_equal(self):
        kwargs = dict(
            jobtype="QMMM",
            low_level_method="ff.prms",
            charge_total=0,
            mult_total=1,
        )
        assert ORCAQMMMJobSettings(**kwargs) == ORCAQMMMJobSettings(**kwargs)

    def test_different_settings_are_not_equal(self):
        s1 = ORCAQMMMJobSettings(
            jobtype="QMMM",
            low_level_method="ff.prms",
            charge_total=0,
            mult_total=1,
        )
        s2 = ORCAQMMMJobSettings(
            jobtype="QMMM",
            low_level_method="ff.prms",
            charge_total=1,
            mult_total=2,
        )
        assert s1 != s2

    def test_different_type_returns_not_implemented(self):
        s = ORCAQMMMJobSettings(
            jobtype="QMMM",
            low_level_method="ff.prms",
            charge_total=0,
            mult_total=1,
        )
        assert s.__eq__(5) is NotImplemented
        assert s != 5


class TestORCANEBJobSettings:
    def test_basic_construction(self):
        s = ORCANEBJobSettings(
            joboption="NEB-TS",
            nimages=8,
            ending_xyzfile="product.xyz",
        )
        assert s.joboption == "NEB-TS"
        assert s.nimages == 8
        assert s.ending_xyzfile == "product.xyz"
        assert s.preopt_ends is False

    def test_route_string_with_dft(self):
        s = ORCANEBJobSettings(
            joboption="NEB-CI",
            functional="B3LYP",
            basis="def2-SVP",
        )
        route = s.route_string
        assert "B3LYP def2-SVP NEB-CI" in route

    def test_route_string_with_semiempirical(self):
        s = ORCANEBJobSettings(joboption="NEB-TS", semiempirical="XTB2")
        route = s.route_string
        assert "XTB2 NEB-TS" in route

    def test_route_string_with_freq(self):
        s = ORCANEBJobSettings(
            joboption="NEB",
            functional="B3LYP",
            basis="def2-SVP",
            freq=True,
        )
        assert "Freq" in s.route_string

    def test_route_string_freq_and_numfreq_conflict_raises(self):
        s = ORCANEBJobSettings(
            joboption="NEB",
            functional="B3LYP",
            basis="def2-SVP",
        )
        s.freq = True
        s.numfreq = True
        with pytest.raises(ValueError, match="Cannot specify both"):
            s.route_string

    def test_route_string_appends_solvent_with_id(self):
        s = ORCANEBJobSettings(
            joboption="NEB",
            functional="B3LYP",
            basis="def2-SVP",
            solvent_model="CPCM",
            solvent_id="water",
        )
        assert "CPCM(water)" in s.route_string

    def test_equal_settings_are_equal(self):
        kwargs = dict(joboption="NEB-TS", nimages=8)
        assert ORCANEBJobSettings(**kwargs) == ORCANEBJobSettings(**kwargs)

    def test_different_type_returns_not_implemented(self):
        s = ORCANEBJobSettings(joboption="NEB-TS", nimages=8)
        assert s.__eq__("not a settings object") is NotImplemented
        assert s != "not a settings object"
