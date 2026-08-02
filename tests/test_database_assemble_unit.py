"""
Direct unit tests for chemsmart/database/assemble.py's less-exercised
branches: BaseAssembler's failed-calculation/include_failed handling,
empty-molecules-list short-circuit, custom-basis/custom-solvent/
mulliken-spin-density/modredundant meta-data branches, and the
SingleFileAssembler/SingleFolderAssembler dispatch and
exception-handling paths.

tests/test_database.py already exercises the assembler pipeline
end-to-end against real Gaussian/ORCA output files; these tests
instead construct BaseAssembler/GaussianAssembler directly (bypassing
__init__ and the real OUTPUT_CLASS parsing) with a controlled fake
``.output`` object, to reach branches that would otherwise require
hard-to-obtain test fixtures (a failed calculation, a custom basis
set, a custom solvent, non-null Mulliken spin densities, etc.).
"""

from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest

from chemsmart.database.assemble import (
    BaseAssembler,
    GaussianAssembler,
    SingleFileAssembler,
    SingleFolderAssembler,
)


def _make_base_assembler(output, include_failed=False):
    assembler = BaseAssembler.__new__(BaseAssembler)
    assembler.filename = "some_file.log"
    assembler.folder = None
    assembler.index = ":"
    assembler.include_failed = include_failed
    assembler.target = assembler.filename
    assembler.output = output
    return assembler


class TestAssembleFailedTermination:
    def test_failed_termination_without_include_failed_returns_none(self):
        output = MagicMock()
        output.normal_termination = False
        assembler = _make_base_assembler(output, include_failed=False)
        assert assembler.assemble() is None

    def test_failed_termination_with_include_failed_still_assembles(self):
        """include_failed=True logs a different warning but still
        proceeds to assemble partial data instead of bailing out."""
        output = MagicMock()
        output.normal_termination = False
        output.route_string = "opt"
        output.method = "b3lyp"
        output.basis = "def2svp"
        output.jobtype = "opt"
        output.solvent_model = None
        output.solvent_id = None
        output.solvent_on = False
        output.freq = False
        output.num_basis_functions = 10
        output.spin = "singlet"
        molecule = SimpleNamespace(structure_id="mol-1")
        assembler = _make_base_assembler(output, include_failed=True)
        with (
            patch.object(
                BaseAssembler,
                "molecules_list",
                new_callable=property,
                fget=lambda self: [molecule],
            ),
            patch.object(
                BaseAssembler,
                "get_calculation_results",
                return_value={},
            ),
            patch.object(
                BaseAssembler,
                "get_molecule_info",
                return_value={},
            ),
            patch.object(
                BaseAssembler,
                "build_provenance",
                return_value={"program": "Gaussian"},
            ),
        ):
            record = assembler.assemble()
        assert record is not None


class TestAssembleEmptyMolecules:
    def test_empty_molecules_list_returns_none(self):
        output = MagicMock()
        output.normal_termination = True
        assembler = _make_base_assembler(output)
        with patch.object(
            BaseAssembler,
            "molecules_list",
            new_callable=property,
            fget=lambda self: [],
        ):
            assert assembler.assemble() is None


class TestGetMetaDataCustomBasisAndSolvent:
    def _make_output_for_meta(self, **overrides):
        output = MagicMock()
        output.basis = "6-31g(d)"
        output.method = "b3lyp"
        output.num_basis_functions = 50
        output.spin = "singlet"
        output.jobtype = "opt"
        output.solvent_on = False
        output.solvent_model = None
        output.solvent_id = None
        output.custom_solvent = None
        output.route_string = "opt b3lyp/6-31g(d)"
        output.freq = False
        for key, value in overrides.items():
            setattr(output, key, value)
        return output

    def test_custom_basis_uses_customized_basis_label(self):
        output = self._make_output_for_meta(basis="gen")
        assembler = _make_base_assembler(output)
        with patch(
            "chemsmart.database.assemble.is_custom_basis", return_value=True
        ):
            meta = assembler.get_meta_data()
        assert meta["basis"] == "customized_basis"

    def test_standard_basis_is_standardized_not_customized(self):
        output = self._make_output_for_meta(basis="6-31g(d)")
        assembler = _make_base_assembler(output)
        with patch(
            "chemsmart.database.assemble.is_custom_basis", return_value=False
        ):
            meta = assembler.get_meta_data()
        assert meta["basis"] != "customized_basis"

    def test_custom_solvent_uses_customized_solvent_label(self):
        output = self._make_output_for_meta(
            basis="6-31g(d)",
            solvent_on=True,
            solvent_model="smd",
            solvent_id="mysolvent",
            custom_solvent="eps=5.0",
        )
        assembler = _make_base_assembler(output)
        with (
            patch(
                "chemsmart.database.assemble.is_custom_basis",
                return_value=False,
            ),
            patch(
                "chemsmart.database.assemble.is_custom_solvent",
                return_value=True,
            ),
        ):
            meta = assembler.get_meta_data()
        assert meta["solvent_id"] == "customized_solvent"

    def test_standard_solvent_keeps_solvent_id(self):
        output = self._make_output_for_meta(
            basis="6-31g(d)",
            solvent_on=True,
            solvent_model="smd",
            solvent_id="water",
            custom_solvent=None,
        )
        assembler = _make_base_assembler(output)
        with (
            patch(
                "chemsmart.database.assemble.is_custom_basis",
                return_value=False,
            ),
            patch(
                "chemsmart.database.assemble.is_custom_solvent",
                return_value=False,
            ),
        ):
            meta = assembler.get_meta_data()
        assert meta["solvent_id"] == "water"


class TestGetMoleculeInfoMullikenSpinDensities:
    def test_non_none_spin_densities_are_included(self):
        assembler = _make_base_assembler(MagicMock())
        mol = MagicMock()
        mol.mulliken_atomic_charges = None
        mol.mulliken_spin_densities = [0.1, -0.1]
        mol.rotational_symmetry_number = None
        mol.rotational_constants = None
        mol.point_group = None
        mol.vibrational_frequencies = None
        mol.dipole_moment = None
        mol.dipole_moment_magnitude = None

        info = assembler.get_molecule_info(mol)
        assert info["mulliken_spin_densities"] == [0.1, -0.1]


class TestGaussianAssemblerMetaDataBranches:
    def _make_gaussian_assembler(self, **overrides):
        assembler = GaussianAssembler.__new__(GaussianAssembler)
        output = MagicMock()
        output.basis = "6-31g(d)"
        output.method = "b3lyp"
        output.num_basis_functions = 50
        output.spin = "singlet"
        output.jobtype = "opt"
        output.solvent_on = False
        output.solvent_model = None
        output.solvent_id = None
        output.custom_solvent = None
        output.route_string = "opt b3lyp/6-31g(d)"
        output.freq = False
        output.num_primitive_gaussians = 100
        output.num_cartesian_basis_functions = 60
        output.modredundant_group = None
        for key, value in overrides.items():
            setattr(output, key, value)
        assembler.output = output
        assembler.filename = "some_file.log"
        assembler.target = assembler.filename
        assembler.include_failed = False
        return assembler

    def test_modredundant_group_included_when_present(self):
        assembler = self._make_gaussian_assembler(modredundant_group="1 2 3 F")
        with patch(
            "chemsmart.database.assemble.is_custom_basis", return_value=False
        ):
            meta = assembler.get_meta_data()
        assert meta["modredundant_group"] == "1 2 3 F"

    def test_modredundant_group_absent_when_none(self):
        assembler = self._make_gaussian_assembler(modredundant_group=None)
        with patch(
            "chemsmart.database.assemble.is_custom_basis", return_value=False
        ):
            meta = assembler.get_meta_data()
        assert "modredundant_group" not in meta

    def test_custom_basis_includes_heavy_and_light_element_details(self):
        assembler = self._make_gaussian_assembler(
            basis="gen",
            heavy_elements="Fe",
            heavy_elements_basis="lanl2dz",
            heavy_elements_ecp="lanl2dz",
            light_elements="C,H,O",
            light_elements_basis="6-31g(d)",
        )
        with patch(
            "chemsmart.database.assemble.is_custom_basis", return_value=True
        ):
            meta = assembler.get_meta_data()
        assert meta["custom_basis"]["heavy_elements"] == "Fe"
        assert meta["custom_basis"]["light_elements"] == "C,H,O"


class TestSingleFileAssemblerErrorHandling:
    def test_assemble_data_returns_none_and_logs_on_exception(self):
        single = SingleFileAssembler(filename="broken.log")
        broken_assembler = MagicMock()
        broken_assembler.assemble.side_effect = RuntimeError("parse failed")
        with patch.object(
            SingleFileAssembler,
            "_get_assembler",
            return_value=broken_assembler,
        ):
            assert single.assemble_data is None

    def test_unsupported_file_type_raises(self):
        single = SingleFileAssembler(filename="mystery.txt")
        with patch(
            "chemsmart.database.assemble.get_program_type_from_file",
            return_value="unknown",
        ):
            with pytest.raises(ValueError, match="Unsupported file"):
                single._get_assembler("mystery.txt")


class TestSingleFolderAssemblerErrorHandling:
    def test_assemble_data_returns_none_and_logs_on_exception(self):
        single = SingleFolderAssembler(folder="broken_folder")
        broken_assembler = MagicMock()
        broken_assembler.assemble.side_effect = RuntimeError("parse failed")
        with patch.object(
            SingleFolderAssembler,
            "_get_assembler",
            return_value=broken_assembler,
        ):
            assert single.assemble_data is None

    def test_mixed_program_folder_raises(self):
        single = SingleFolderAssembler(folder="mixed_folder")
        with patch(
            "chemsmart.database.assemble.BaseFolder"
        ) as mock_base_folder_cls:
            mock_base_folder_cls.return_value.get_program_type_from_folder.return_value = (
                "mixed"
            )
            with pytest.raises(ValueError, match="multiple programs"):
                single._get_assembler("mixed_folder")

    def test_unsupported_folder_type_raises(self):
        single = SingleFolderAssembler(folder="gaussian_folder")
        with patch(
            "chemsmart.database.assemble.BaseFolder"
        ) as mock_base_folder_cls:
            mock_base_folder_cls.return_value.get_program_type_from_folder.return_value = (
                "gaussian"
            )
            with pytest.raises(ValueError, match="Unsupported folder"):
                single._get_assembler("gaussian_folder")
