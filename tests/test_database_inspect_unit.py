"""
Direct unit tests for chemsmart.database.inspect.DatabaseInspector.

Complements tests/test_database.py::TestDatabaseInspect (which exercises
overview/record/structure/molecule views through a real assembled
database) by mocking the Database layer and feeding hand-built
record/molecule/structure dicts to reach formatting branches (open-shell
electronic structure, custom basis/solvent detail, mixed electronic
states, physical properties, Mulliken analysis, vibrational data, etc.)
that are hard to set up with real quantum-chemistry fixtures.
"""

from unittest.mock import MagicMock, patch

import pytest

from chemsmart.database.inspect import DatabaseInspector


def _make_inspector(**kwargs):
    with patch("chemsmart.database.inspect.Database") as mock_db_cls:
        mock_db_cls.return_value = MagicMock()
        inspector = DatabaseInspector("fake.db", **kwargs)
    return inspector


class TestResolveId:
    def test_index_not_found_raises(self):
        inspector = _make_inspector(index=3)
        inspector.db.get_record.return_value = None
        with pytest.raises(ValueError, match="No record found with index"):
            inspector.resolve_id()

    def test_index_found_returns_record_id(self):
        inspector = _make_inspector(index=3)
        inspector.db.get_record.return_value = {"record_id": "rid123"}
        assert inspector.resolve_id() == "rid123"

    def test_record_id_path_delegates_to_partial_id_lookup(self):
        inspector = _make_inspector(record_id="abc")
        inspector.db.get_record_by_partial_id.return_value = "abcfull"
        assert inspector.resolve_id() == "abcfull"

    def test_neither_index_nor_record_id_raises(self):
        inspector = _make_inspector()
        with pytest.raises(ValueError, match="Either index or record_id"):
            inspector.resolve_id()


class TestRecordDetail:
    def test_record_not_found_raises(self):
        inspector = _make_inspector(index=1)
        inspector.db.get_record.side_effect = [
            {"record_id": "rid1"},  # resolve_id's lookup
            None,  # record_detail's own lookup
        ]
        with pytest.raises(ValueError, match="Record not found"):
            inspector.record_detail()


class TestFormatOverview:
    def _base_stats(self, **overrides):
        stats = {
            "db_file": "my.db",
            "db_size": 2048,
            "num_records": 2,
            "num_molecules": 1,
            "num_structures": 3,
            "programs": [("gaussian", 2)],
            "methods": [],
            "basis_sets": [],
            "assembled_first": "2024-01-01T00:00:00",
            "assembled_last": "2024-01-02T00:00:00",
            "gas_phase_records": 2,
            "solvated_records": 0,
            "solvents": [],
            "jobtypes": [("opt", 2)],
            "energy_count": 2,
            "forces_count": 0,
            "thermo_count": 1,
        }
        stats.update(overrides)
        return stats

    def test_empty_methods_and_basis_sets_show_null(self):
        inspector = _make_inspector()
        with patch.object(
            inspector, "overview", return_value=self._base_stats()
        ):
            text = inspector.format_overview()
        assert "Methods" in text
        assert "NULL" in text
        assert "Basis Sets" in text

    def test_methods_basis_and_solvents_present(self):
        inspector = _make_inspector()
        stats = self._base_stats(
            methods=[("b3lyp", 2)],
            basis_sets=[("sto-3g", 2)],
            solvents=[("water", 1)],
            solvated_records=1,
        )
        with patch.object(inspector, "overview", return_value=stats):
            text = inspector.format_overview()
        assert "b3lyp (2)" in text
        assert "sto-3g (2)" in text
        assert "water (1)" in text


class TestFormatRecordDetail:
    def _base_record(self, **overrides):
        record = {
            "record_index": 1,
            "record_id": "rid123456789",
            "meta": {
                "method": "b3lyp",
                "basis": "sto-3g",
                "spin": "restricted",
                "jobtype": "opt",
                "solvent_on": False,
                "route_string": "#p b3lyp/sto-3g opt",
            },
            "results": {
                "total_energy": -1.0,
                "num_unpaired_electrons": 0,
                "homo_energy": -5.0,
                "lumo_energy": 1.0,
                "fmo_gap": 6.0,
            },
            "provenance": {
                "source": "/tmp/job.log",
                "program": "gaussian",
                "program_version": "16",
                "parser": "cclib",
                "chemsmart_version": "1.0",
                "source_file_hash": "abc",
                "source_size": 1024,
                "source_date": "2024-01-01",
                "assembled_at": "2024-01-01T00:00:00",
            },
            "molecules": [],
        }
        record["meta"].update(overrides.pop("meta", {}))
        record["results"].update(overrides.pop("results", {}))
        record["provenance"].update(overrides.pop("provenance", {}))
        record.update(overrides)
        return record

    def _run(self, inspector, record):
        with patch.object(inspector, "record_detail", return_value=record):
            return inspector.format_record_detail()

    def test_closed_shell_basic_fields(self):
        inspector = _make_inspector()
        text = self._run(inspector, self._base_record())
        assert "HOMO Energy (eV)" in text
        assert "LUMO Energy (eV)" in text
        assert "α-HOMO" not in text

    def test_open_shell_shows_alpha_beta_and_somo(self):
        inspector = _make_inspector()
        record = self._base_record(
            results={
                "num_unpaired_electrons": 1,
                "alpha_homo_energy": -5.0,
                "beta_homo_energy": -4.5,
                "alpha_lumo_energy": 1.0,
                "beta_lumo_energy": 1.5,
                "alpha_fmo_gap": 6.0,
                "beta_fmo_gap": 6.0,
                "somo_energies": [-4.8],
            }
        )
        text = self._run(inspector, record)
        assert "α-HOMO Energy (eV)" in text
        assert "β-HOMO Energy (eV)" in text
        assert "SOMO Energies (eV)" in text

    def test_open_shell_without_somo_energies(self):
        inspector = _make_inspector()
        record = self._base_record(
            results={
                "num_unpaired_electrons": 1,
                "alpha_homo_energy": -5.0,
                "beta_homo_energy": -4.5,
                "alpha_lumo_energy": 1.0,
                "beta_lumo_energy": 1.5,
                "alpha_fmo_gap": 6.0,
                "beta_fmo_gap": 6.0,
            }
        )
        text = self._run(inspector, record)
        assert "SOMO Energies (eV)" not in text

    def test_customized_basis_without_detail_shows_unavailable_note(self):
        inspector = _make_inspector()
        record = self._base_record(meta={"basis": "customized_basis"})
        text = self._run(inspector, record)
        assert "unavailable" in text

    def test_customized_basis_with_full_detail(self):
        inspector = _make_inspector()
        record = self._base_record(
            meta={
                "basis": "customized_basis",
                "custom_basis": {
                    "light_elements": ["H", "C"],
                    "light_elements_basis": "6-31G",
                    "heavy_elements": ["Fe"],
                    "heavy_elements_basis": {
                        "Fe": [{"shell": "S"}, {"shell": "P"}]
                    },
                    "heavy_elements_ecp": {"Fe": "def2-ecp"},
                },
            }
        )
        text = self._run(inspector, record)
        assert "Light Elements" in text
        assert "[6-31G]" in text
        assert "Heavy Elements" in text
        assert "Fe Shells" in text
        assert "S, P" in text
        assert "Heavy Elements ECP" in text

    def test_customized_basis_light_elements_only_no_light_basis(self):
        inspector = _make_inspector()
        record = self._base_record(
            meta={
                "basis": "customized_basis",
                "custom_basis": {
                    "light_elements": ["H", "C"],
                    "light_elements_basis": None,
                },
            }
        )
        text = self._run(inspector, record)
        assert "Light Elements" in text
        assert "[" not in text.split("Light Elements")[1].split("\n")[0]
        assert "Heavy Elements" not in text

    def test_customized_basis_heavy_elements_only_no_light_elements(self):
        inspector = _make_inspector()
        record = self._base_record(
            meta={
                "basis": "customized_basis",
                "custom_basis": {
                    "light_elements": None,
                    "heavy_elements": ["Fe"],
                },
            }
        )
        text = self._run(inspector, record)
        assert "Light Elements" not in text
        assert "Heavy Elements" in text

    def test_solvent_on_with_custom_solvent_detail(self):
        inspector = _make_inspector()
        record = self._base_record(
            meta={
                "solvent_on": True,
                "solvent_model": "PCM",
                "solvent_id": "water",
                "custom_solvent": {
                    "SolventName": "MySolvent",
                    "Eps": 78.0,
                    "EpsInf": 1.8,
                    "HbondAcidity": 0.1,
                    "HbondBasicity": 0.2,
                    "SurfaceTensionAtInterface": 0.3,
                    "CarbonAromaticity": 0.4,
                    "ElectronegativeHalogenicity": 0.5,
                },
            }
        )
        text = self._run(inspector, record)
        assert "Solvent Model" in text
        assert "Custom Solvent Name" in text
        assert "Eps" in text

    def test_solvent_on_custom_solvent_generic_name_is_skipped(self):
        inspector = _make_inspector()
        record = self._base_record(
            meta={
                "solvent_on": True,
                "solvent_model": "PCM",
                "solvent_id": "water",
                "custom_solvent": {"SolventName": "Generic"},
            }
        )
        text = self._run(inspector, record)
        assert "Custom Solvent Name" not in text

    def test_solvent_on_without_custom_solvent_dict(self):
        inspector = _make_inspector()
        record = self._base_record(
            meta={
                "solvent_on": True,
                "solvent_model": "PCM",
                "solvent_id": "water",
                "custom_solvent": None,
            }
        )
        text = self._run(inspector, record)
        assert "Solvent Model" in text
        assert "Custom Solvent Name" not in text

    def test_performance_section_shown_when_core_hours_present(self):
        inspector = _make_inspector()
        record = self._base_record(
            results={
                "total_core_hours": 12.5,
                "total_elapsed_walltime": 3.2,
            }
        )
        text = self._run(inspector, record)
        assert "Core Hours" in text
        assert "Elapsed Walltime" in text

    def test_thermochemistry_section_shown_when_gibbs_present(self):
        inspector = _make_inspector()
        record = self._base_record(
            meta={"temperature_in_K": 298.15, "pressure_in_atm": 1.0},
            results={
                "gibbs_free_energy": -1.5,
                "zero_point_energy": 0.01,
                "internal_energy": -1.4,
                "enthalpy": -1.3,
                "entropy": 0.002,
            },
        )
        text = self._run(inspector, record)
        assert "Thermochemistry" in text
        assert "Gibbs Free Energy" in text

    def test_normal_termination_shown_when_present(self):
        inspector = _make_inspector()
        record = self._base_record(provenance={"normal_termination": True})
        text = self._run(inspector, record)
        assert "Normal Termination" in text

    def test_uniform_structures_table(self):
        inspector = _make_inspector()
        record = self._base_record(
            molecules=[
                {
                    "index": 1,
                    "structure_id": "sid1",
                    "chemical_formula": "H2O",
                    "charge": 0,
                    "multiplicity": 1,
                    "number_of_atoms": 3,
                    "energy": -1.0,
                    "is_optimized_structure": True,
                },
                {
                    "index": 2,
                    "structure_id": "sid2",
                    "chemical_formula": "H2O",
                    "charge": 0,
                    "multiplicity": 1,
                    "number_of_atoms": 3,
                    "energy": None,
                    "is_optimized_structure": None,
                },
            ]
        )
        text = self._run(inspector, record)
        assert "Formula" in text
        assert "H2O" in text

    def test_non_uniform_structures_table(self):
        inspector = _make_inspector()
        record = self._base_record(
            molecules=[
                {
                    "index": 1,
                    "structure_id": "sid1",
                    "chemical_formula": "H2O",
                    "charge": 0,
                    "multiplicity": 1,
                    "number_of_atoms": 3,
                    "energy": -1.0,
                },
                {
                    "index": 2,
                    "structure_id": "sid2",
                    "chemical_formula": "H2O+",
                    "charge": 1,
                    "multiplicity": 2,
                    "number_of_atoms": 3,
                    "energy": None,
                },
            ]
        )
        text = self._run(inspector, record)
        assert "Idx" in text
        assert "H2O+" in text


class TestFormatStructureDetail:
    def _base_record_and_struct(self, **struct_overrides):
        struct = {
            "structure_index_in_file": 5,
            "energy": -1.0,
            "is_optimized_structure": True,
            "structure_id": "sid123",
            "charge": 0,
            "multiplicity": 1,
            "molecule_id": "mid123",
            "chemical_formula": "H2O",
            "mass": 18.0153,
            "smiles": "O",
            "chemical_symbols": ["O", "H", "H"],
            "positions": [
                [0.0, 0.0, 0.0],
                [0.0, 0.76, 0.59],
                [0.0, -0.76, 0.59],
            ],
        }
        struct.update(struct_overrides)
        record = {
            "record_index": 1,
            "record_id": "rid123",
            "provenance": {"source": "/tmp/job.log"},
        }
        return record, struct

    def _run(self, inspector, record, struct):
        with patch.object(
            inspector, "structure_detail", return_value=(record, struct)
        ):
            return inspector.format_structure_detail()

    def test_full_detail_with_all_optional_sections(self):
        inspector = _make_inspector()
        record, struct = self._base_record_and_struct(
            frozen_atoms=[1, 2],
            center_of_mass=[0.0, 0.0, 0.1],
            moments_of_inertia=[1.0, 2.0, 3.0],
            rotational_symmetry_number=2,
            rotational_constants=[1.0e9, 2.0e9, 3.0e9],
            point_group="C2v",
            dipole_moment=[0.0, 0.0, 1.5],
            dipole_moment_magnitude=1.5,
            mulliken_atomic_charges={"O1": -0.5, "H2": 0.25, "H3": 0.25},
            mulliken_spin_densities={"O1": 0.0, "H2": 0.0, "H3": 0.0},
            num_vibrational_modes=3,
            vibrational_frequencies=[100.0, 200.0, 300.0],
        )
        text = self._run(inspector, record, struct)
        assert "Optimized" in text
        assert "Coordinates" in text
        assert "Frozen Atoms" in text
        assert "Point Group" in text
        assert "Rotational Symmetry Number" in text
        assert "Center of Mass" in text
        assert "Moments of Inertia" in text
        assert "Rotational Constants" in text
        assert "Dipole Moment" in text
        assert "Mulliken Population Analysis" in text
        assert "Vibrational Frequencies" in text

    def test_physical_properties_point_group_only(self):
        inspector = _make_inspector()
        record, struct = self._base_record_and_struct(point_group="C2v")
        text = self._run(inspector, record, struct)
        assert "Point Group" in text
        assert "Center of Mass" not in text
        assert "Moments of Inertia" not in text
        assert "Rotational Constants" not in text
        assert "Rotational Symmetry Number" not in text

    def test_mulliken_charges_only(self):
        inspector = _make_inspector()
        record, struct = self._base_record_and_struct(
            mulliken_atomic_charges={"O1": -0.5}
        )
        text = self._run(inspector, record, struct)
        assert "Charge" in text

    def test_mulliken_spin_only(self):
        inspector = _make_inspector()
        record, struct = self._base_record_and_struct(
            mulliken_spin_densities={"O1": 0.1}
        )
        text = self._run(inspector, record, struct)
        assert "Spin" in text

    def test_vibrational_modes_without_frequencies_list(self):
        inspector = _make_inspector()
        record, struct = self._base_record_and_struct(
            num_vibrational_modes=0, vibrational_frequencies=[]
        )
        text = self._run(inspector, record, struct)
        assert "Vibrational Frequencies" in text

    def test_minimal_fields_skip_optional_sections(self):
        inspector = _make_inspector()
        record, struct = self._base_record_and_struct(
            is_optimized_structure=None,
            chemical_symbols=[],
            positions=[],
        )
        text = self._run(inspector, record, struct)
        assert "Coordinates" not in text
        assert "Frozen Atoms" not in text
        assert "Physical Properties" not in text
        assert "Dipole Moment" not in text
        assert "Mulliken Population Analysis" not in text


class TestFormatMoleculeDetail:
    def _base_molecule(self, **overrides):
        molecule = {
            "molecule_id": "mid123",
            "chemical_formula": "H2O",
            "mass": 18.0153,
            "is_aromatic": False,
            "is_chiral": False,
            "is_linear": False,
            "is_multicomponent": False,
            "smiles": "O",
            "inchi": "InChI=1S/H2O/h1H2",
        }
        molecule.update(overrides)
        return molecule

    def test_molecule_not_found_raises(self):
        inspector = _make_inspector(molecule_id="mid123")
        inspector.db.get_molecule_by_partial_id.return_value = "mid123full"
        inspector.db.get_molecule.return_value = None
        with pytest.raises(ValueError, match="Molecule not found"):
            inspector.molecule_detail()

    def test_no_structures_shows_none(self):
        inspector = _make_inspector()
        molecule = self._base_molecule()
        with patch.object(
            inspector, "molecule_detail", return_value=(molecule, [], [])
        ):
            text = inspector.format_molecule_detail()
        assert "(none)" in text

    def test_element_counts_shown(self):
        inspector = _make_inspector()
        molecule = self._base_molecule(
            element_counts={"O": 1, "H": 2}, number_of_atoms=3
        )
        with patch.object(
            inspector, "molecule_detail", return_value=(molecule, [], [])
        ):
            text = inspector.format_molecule_detail()
        assert "Composition" in text
        assert "H 2" in text

    def test_all_tags_active(self):
        inspector = _make_inspector()
        molecule = self._base_molecule(
            is_aromatic=True,
            is_chiral=True,
            is_linear=True,
            is_multicomponent=True,
        )
        with patch.object(
            inspector, "molecule_detail", return_value=(molecule, [], [])
        ):
            text = inspector.format_molecule_detail()
        assert "aromatic" in text
        assert "chiral" in text

    def test_uniform_structures_with_primary_mb(self):
        inspector = _make_inspector()
        molecule = self._base_molecule()
        structures = [
            {
                "structure_id": "sid1",
                "charge": 0,
                "multiplicity": 1,
                "primary_method_basis": ("b3lyp", "sto-3g"),
                "primary_energy": -2.0,
            },
            {
                "structure_id": "sid2",
                "charge": 0,
                "multiplicity": 1,
                "primary_method_basis": ("b3lyp", "sto-3g"),
                "primary_energy": -1.0,
            },
        ]
        with (
            patch.object(
                inspector,
                "molecule_detail",
                return_value=(molecule, structures, []),
            ),
            patch(
                "chemsmart.database.inspect.sort_structure_dicts_by_energy",
                return_value=structures,
            ),
        ):
            text = inspector.format_molecule_detail()
        assert "Energy[b3lyp/sto-3g]" in text
        assert "ΔE (kcal/mol)" in text

    def test_structures_without_primary_mb_uses_default_header(self):
        inspector = _make_inspector()
        molecule = self._base_molecule()
        structures = [
            {
                "structure_id": "sid1",
                "charge": 0,
                "multiplicity": 1,
                "primary_method_basis": None,
                "primary_energy": None,
            }
        ]
        with (
            patch.object(
                inspector,
                "molecule_detail",
                return_value=(molecule, structures, []),
            ),
            patch(
                "chemsmart.database.inspect.sort_structure_dicts_by_energy",
                return_value=structures,
            ),
        ):
            text = inspector.format_molecule_detail()
        assert "Energy (Eh)" in text
        assert "N/A" in text

    def test_structures_present_but_sorted_structures_empty(self):
        # sort_structure_dicts_by_energy could in principle return an
        # empty list even though `structures` itself is non-empty; the
        # `if sorted_structures:` guard exists specifically for that.
        inspector = _make_inspector()
        molecule = self._base_molecule()
        structures = [
            {
                "structure_id": "sid1",
                "charge": 0,
                "multiplicity": 1,
                "primary_method_basis": None,
                "primary_energy": None,
            }
        ]
        with (
            patch.object(
                inspector,
                "molecule_detail",
                return_value=(molecule, structures, []),
            ),
            patch(
                "chemsmart.database.inspect.sort_structure_dicts_by_energy",
                return_value=[],
            ),
        ):
            text = inspector.format_molecule_detail()
        assert "Energy (Eh)" in text

    def test_mixed_electronic_states_note(self):
        inspector = _make_inspector()
        molecule = self._base_molecule()
        structures = [
            {
                "structure_id": "sid1",
                "charge": 0,
                "multiplicity": 1,
                "primary_method_basis": None,
                "primary_energy": -1.0,
            },
            {
                "structure_id": "sid2",
                "charge": 1,
                "multiplicity": 2,
                "primary_method_basis": None,
                "primary_energy": -0.5,
            },
        ]
        with (
            patch.object(
                inspector,
                "molecule_detail",
                return_value=(molecule, structures, []),
            ),
            patch(
                "chemsmart.database.inspect.sort_structure_dicts_by_energy",
                return_value=structures,
            ),
        ):
            text = inspector.format_molecule_detail()
        assert "ΔE (—)" in text
        assert "Mixed electronic states detected" in text

    def test_related_records_table_populated(self):
        inspector = _make_inspector()
        molecule = self._base_molecule()
        records = [
            {
                "record_id": "rid123456789",
                "jobtype": "opt",
                "program": "gaussian",
                "method": "b3lyp",
                "basis": "sto-3g",
                "total_energy": -1.0,
            }
        ]
        with patch.object(
            inspector,
            "molecule_detail",
            return_value=(molecule, [], records),
        ):
            text = inspector.format_molecule_detail()
        assert "rid123456789"[:12] in text
        assert "gaussian" in text


class TestStandaloneStructureDetail:
    def test_structure_not_found_raises(self):
        inspector = _make_inspector(structure_id="sid123")
        inspector.db.get_structure_by_partial_id.return_value = "sid123full"
        inspector.db.get_structure.return_value = None
        with pytest.raises(ValueError, match="Structure not found"):
            inspector.standalone_structure_detail()


class TestFormatStandaloneStructureDetail:
    def _base_struct(self, **overrides):
        struct = {
            "structure_id": "sid123",
            "charge": 0,
            "multiplicity": 1,
            "molecule_id": "mid123",
            "chemical_formula": "H2O",
            "mass": 18.0153,
            "smiles": "O",
            "chemical_symbols": ["O", "H", "H"],
            "positions": [
                [0.0, 0.0, 0.0],
                [0.0, 0.76, 0.59],
                [0.0, -0.76, 0.59],
            ],
        }
        struct.update(overrides)
        return struct

    def _run(self, inspector, struct, records=None):
        with patch.object(
            inspector,
            "standalone_structure_detail",
            return_value=(struct, records or []),
        ):
            return inspector.format_standalone_structure_detail()

    def test_coordinates_and_no_physical_properties(self):
        inspector = _make_inspector()
        text = self._run(inspector, self._base_struct())
        assert "Coordinates" in text
        assert "Physical Properties" not in text

    def test_com_only(self):
        inspector = _make_inspector()
        struct = self._base_struct(center_of_mass=[0.0, 0.0, 0.1])
        text = self._run(inspector, struct)
        assert "Physical Properties" in text
        assert "Center of Mass" in text
        assert "Moments of Inertia" not in text

    def test_moi_only(self):
        inspector = _make_inspector()
        struct = self._base_struct(moments_of_inertia=[1.0, 2.0, 3.0])
        text = self._run(inspector, struct)
        assert "Moments of Inertia" in text
        assert "Center of Mass" not in text

    def test_no_coordinates(self):
        inspector = _make_inspector()
        struct = self._base_struct(chemical_symbols=[], positions=[])
        text = self._run(inspector, struct)
        assert "Coordinates" not in text

    def test_related_records_included(self):
        inspector = _make_inspector()
        records = [
            {
                "record_id": "rid987654321",
                "jobtype": "sp",
                "program": "orca",
                "method": "wb97xd",
                "basis": "def2tzvp",
                "total_energy": -2.0,
            }
        ]
        text = self._run(inspector, self._base_struct(), records)
        assert "orca" in text


class TestFormatRelatedRecordsTable:
    def test_empty_records_shows_none(self):
        inspector = _make_inspector()
        lines = inspector._format_related_records_table([])
        assert any("(none)" in line for line in lines)

    def test_records_with_missing_optional_fields(self):
        inspector = _make_inspector()
        lines = inspector._format_related_records_table(
            [{"record_id": "rid1"}]
        )
        text = "\n".join(lines)
        assert "rid1" in text
