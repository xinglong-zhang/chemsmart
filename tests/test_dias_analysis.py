"""
Direct unit tests for ``chemsmart.analysis.dias``.

The module has three classes: the abstract ``DIASOutputFolder`` base
class, and its two concrete subclasses ``GaussianDIASLogFolder`` and
``ORCADIASOutFolder`` which discover DI-AS calculation output files in
a folder (by filename pattern), extract energies/geometries from them,
and combine everything into strain/interaction energy decompositions.

``Gaussian16Output`` and ``ORCAOutput`` are mocked throughout so these
tests don't need real quantum chemistry output files -- only filenames
that match the expected DI-AS naming convention.
"""

import os
from unittest.mock import MagicMock, patch

import matplotlib

matplotlib.use("Agg", force=True)

import numpy as np
import pytest
from ase import units

from chemsmart.analysis.dias import (
    DIASOutputFolder,
    GaussianDIASLogFolder,
    ORCADIASOutFolder,
)


def _touch(path):
    open(path, "w").close()


def _to_kcal_per_mol(hartree_energy):
    energy = hartree_energy * units.Hartree
    energy /= units.kcal / units.mol
    return energy


class _FakeMolecule:
    def __init__(self, energy, distance):
        self.energy = energy
        self.final_energy = energy
        self.positions = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, distance]])
        self.chemical_symbols = ["C", "H"]


# basename -> (energy in Hartree, atom1-atom2 distance in Angstrom)
#
# NOTE: full-molecule filenames deliberately contain the literal
# substring "dias_p" (e.g. "rxndias_p1.log") rather than the more
# natural "rxn_p1.log". See
# TestGaussianDIASBugFullMoleculeGroupIndex below:
# chemsmart/analysis/dias.py:439 always reads match.group(1) from
# `gaussian_dias_filename_point_without_fragment_without_reactant`,
# but that pattern has two alternatives and only the first (the one
# requiring a literal "dias_p" in the filename) populates group(1).
# Real DI-AS jobs name full-molecule files "{label}_p{i}.log" with
# no "dias" substring unless the user's job label happens to contain
# it, so this is the crashing path in normal use.
GAUSSIAN_FILES = {
    "rxndias_p1.log": (-1.00, 1.0),
    "rxndias_p2.log": (-1.10, 1.2),
    "rxndias_p1_f1.log": (-0.50, 1.0),
    "rxndias_p2_f1.log": (-0.55, 1.2),
    "rxndias_p1_f2.log": (-0.40, 1.0),
    "rxndias_p2_f2.log": (-0.45, 1.2),
    "rxndias_r1.log": (-0.48, 1.0),
    "rxndias_r2.log": (-0.39, 1.0),
}

ORCA_FILES = {
    "rxn_p1.out": (-1.00, 1.0),
    "rxn_p2.out": (-1.10, 1.2),
    "rxn_p1_f1.out": (-0.50, 1.0),
    "rxn_p2_f1.out": (-0.55, 1.2),
    "rxn_p1_f2.out": (-0.40, 1.0),
    "rxn_p2_f2.out": (-0.45, 1.2),
    "rxn_r1.out": (-0.48, 1.0),
    "rxn_r2.out": (-0.39, 1.0),
}


def _output_factory(file_map):
    def factory(filename):
        basename = os.path.basename(filename)
        energy, distance = file_map[basename]
        mock = MagicMock()
        mock.molecule = _FakeMolecule(energy=energy, distance=distance)
        return mock

    return factory


@pytest.fixture()
def gaussian_folder(tmp_path):
    for name in GAUSSIAN_FILES:
        _touch(tmp_path / name)
    return tmp_path


@pytest.fixture()
def orca_folder(tmp_path):
    for name in ORCA_FILES:
        _touch(tmp_path / name)
    return tmp_path


@pytest.fixture()
def patched_gaussian_output():
    with patch(
        "chemsmart.analysis.dias.Gaussian16Output",
        side_effect=_output_factory(GAUSSIAN_FILES),
    ) as mock_cls:
        yield mock_cls


@pytest.fixture()
def patched_orca_output():
    with patch(
        "chemsmart.analysis.dias.ORCAOutput",
        side_effect=_output_factory(ORCA_FILES),
    ) as mock_cls:
        yield mock_cls


@pytest.fixture()
def gaussian_dias_folder(gaussian_folder, patched_gaussian_output):
    return GaussianDIASLogFolder(folder=str(gaussian_folder), atom1=1, atom2=2)


@pytest.fixture()
def orca_dias_folder(orca_folder, patched_orca_output):
    return ORCADIASOutFolder(folder=str(orca_folder), atom1=1, atom2=2)


class TestDIASOutputFolderAbstract:
    def test_base_class_cannot_discover_files(self, tmp_path):
        with pytest.raises(NotImplementedError):
            DIASOutputFolder(folder=str(tmp_path), atom1=1, atom2=2)

    def test_abstract_methods_all_raise_not_implemented(self, tmp_path):
        # __init__ short-circuits on the very first abstract method it
        # calls, so build a bare, un-initialized instance to exercise
        # each of the other stub methods individually.
        folder = object.__new__(DIASOutputFolder)
        folder.folder = str(tmp_path)
        with pytest.raises(NotImplementedError):
            folder._get_all_files_fragment1()
        with pytest.raises(NotImplementedError):
            folder._get_all_files_fragment2()
        with pytest.raises(NotImplementedError):
            folder._get_all_files_reactants()
        with pytest.raises(NotImplementedError):
            folder._get_all_energies([])
        # `_get_ref_file` is an unimplemented no-op stub, not abstract.
        assert folder._get_ref_file() is None


class TestGaussianDIASLogFolderFileDiscovery:
    def test_full_molecule_files_sorted_by_point_index(
        self, gaussian_folder, patched_gaussian_output
    ):
        folder = GaussianDIASLogFolder(
            folder=str(gaussian_folder), atom1=1, atom2=2
        )
        names = [
            os.path.basename(f) for f in folder._get_all_files_full_molecule()
        ]
        assert names == ["rxndias_p1.log", "rxndias_p2.log"]

    def test_fragment1_files_sorted_by_point_index(self, gaussian_dias_folder):
        names = [
            os.path.basename(f)
            for f in gaussian_dias_folder._get_all_files_fragment1()
        ]
        assert names == ["rxndias_p1_f1.log", "rxndias_p2_f1.log"]

    def test_fragment2_files_sorted_by_point_index(self, gaussian_dias_folder):
        names = [
            os.path.basename(f)
            for f in gaussian_dias_folder._get_all_files_fragment2()
        ]
        assert names == ["rxndias_p1_f2.log", "rxndias_p2_f2.log"]

    def test_reactant_files_sorted_by_reactant_index(
        self, gaussian_dias_folder
    ):
        names = [
            os.path.basename(f)
            for f in gaussian_dias_folder._get_all_files_reactants()
        ]
        assert names == ["rxndias_r1.log", "rxndias_r2.log"]

    def test_job_basename_from_fragment1_file(self, gaussian_dias_folder):
        assert gaussian_dias_folder.job_basename == "rxndias"

    def test_job_basename_raises_when_no_fragment1_file(self, tmp_path):
        _touch(tmp_path / "unrelated.log")
        # Bypass __init__ (which would itself fail while trying to
        # discover full-molecule/fragment files) to exercise the
        # job_basename property's own assertion in isolation.
        folder = object.__new__(GaussianDIASLogFolder)
        folder.folder = str(tmp_path)
        with pytest.raises(AssertionError):
            _ = folder.job_basename


class TestGaussianDIASLogFolderEnergies:
    def test_get_all_energies_converts_hartree_to_kcal_per_mol(
        self, gaussian_dias_folder
    ):
        energies = gaussian_dias_folder._get_all_energies(
            ["rxndias_p1.log", "rxndias_p2.log"]
        )
        assert energies == pytest.approx(
            [_to_kcal_per_mol(-1.00), _to_kcal_per_mol(-1.10)]
        )

    def test_reactant_energies_list_order(self, gaussian_dias_folder):
        assert gaussian_dias_folder.reactant_energies_list == pytest.approx(
            [_to_kcal_per_mol(-0.48), _to_kcal_per_mol(-0.39)]
        )

    def test_irc_points_matches_number_of_points(self, gaussian_dias_folder):
        assert gaussian_dias_folder.irc_points == 2

    def test_relative_total_energies(self, gaussian_dias_folder):
        reactant_sum = sum(gaussian_dias_folder.reactant_energies_list)
        expected = [
            _to_kcal_per_mol(-1.00) - reactant_sum,
            _to_kcal_per_mol(-1.10) - reactant_sum,
        ]
        assert gaussian_dias_folder.list_rel_total_energies == pytest.approx(
            expected
        )

    def test_strain_energies(self, gaussian_dias_folder):
        r1, r2 = gaussian_dias_folder.reactant_energies_list
        f1 = gaussian_dias_folder.fragment1_energies_list
        f2 = gaussian_dias_folder.fragment2_energies_list
        expected = [(f1[i] - r1) + (f2[i] - r2) for i in range(2)]
        assert (
            gaussian_dias_folder.list_total_strain_energies
            == pytest.approx(expected)
        )

    def test_strain_energies_for_fragments(self, gaussian_dias_folder):
        r1, r2 = gaussian_dias_folder.reactant_energies_list
        f1 = gaussian_dias_folder.fragment1_energies_list
        f2 = gaussian_dias_folder.fragment2_energies_list
        strain1, strain2 = (
            gaussian_dias_folder.get_strain_energies_for_fragments()
        )
        assert strain1 == pytest.approx([f1[i] - r1 for i in range(2)])
        assert strain2 == pytest.approx([f2[i] - r2 for i in range(2)])

    def test_interaction_energies(self, gaussian_dias_folder):
        full = gaussian_dias_folder.full_molecule_energies_list
        f1 = gaussian_dias_folder.fragment1_energies_list
        f2 = gaussian_dias_folder.fragment2_energies_list
        expected = [full[i] - f1[i] - f2[i] for i in range(2)]
        assert (
            gaussian_dias_folder.list_total_interaction_energies
            == pytest.approx(expected)
        )

    def test_irc_points_mismatch_raises_assertion_error(self, tmp_path):
        # Two full-molecule points, but only one fragment1 point.
        subset = {
            "rxndias_p1.log": (-1.00, 1.0),
            "rxndias_p2.log": (-1.10, 1.2),
            "rxndias_p1_f1.log": (-0.50, 1.0),
            "rxndias_p1_f2.log": (-0.40, 1.0),
            "rxndias_p2_f2.log": (-0.45, 1.2),
            "rxndias_r1.log": (-0.48, 1.0),
            "rxndias_r2.log": (-0.39, 1.0),
        }
        for name in subset:
            _touch(tmp_path / name)
        with patch(
            "chemsmart.analysis.dias.Gaussian16Output",
            side_effect=_output_factory(subset),
        ):
            with pytest.raises(AssertionError):
                GaussianDIASLogFolder(folder=str(tmp_path), atom1=1, atom2=2)


class TestGaussianDIASLogFolderReactionCoordinates:
    def test_reaction_coordinates_use_atom_distance(
        self, gaussian_dias_folder
    ):
        assert gaussian_dias_folder.list_rc == pytest.approx([1.0, 1.2])

    def test_unknown_file_extension_raises_value_error(
        self, gaussian_dias_folder
    ):
        with pytest.raises(ValueError, match="unknown format"):
            gaussian_dias_folder._get_reaction_coordinates(["bad.txt"])

    def test_ref_file_properties_none_by_default(self, gaussian_dias_folder):
        assert gaussian_dias_folder.ref_file_rel_energy is None
        assert gaussian_dias_folder.ref_file_rc is None

    def test_ref_file_properties_with_ref_file(
        self, gaussian_folder, patched_gaussian_output
    ):
        ref_file = str(gaussian_folder / "rxndias_p1.log")
        folder = GaussianDIASLogFolder(
            folder=str(gaussian_folder),
            atom1=1,
            atom2=2,
            ref_file=ref_file,
        )
        reactant_sum = sum(folder.reactant_energies_list)
        assert folder.ref_file_rel_energy == pytest.approx(
            _to_kcal_per_mol(-1.00) - reactant_sum
        )
        assert folder.ref_file_rc == pytest.approx(1.0)


class TestOutputNaming:
    def test_outputname_nonzero_ref_default(
        self, gaussian_folder, patched_gaussian_output
    ):
        folder = GaussianDIASLogFolder(
            folder=str(gaussian_folder), atom1=1, atom2=2
        )
        assert folder.outputname == "dias_nonzero_ref"

    def test_outputname_zero_ref(
        self, gaussian_folder, patched_gaussian_output
    ):
        folder = GaussianDIASLogFolder(
            folder=str(gaussian_folder), atom1=1, atom2=2, zero=True
        )
        assert folder.outputname == "dias_zero_ref"

    def test_outputname_with_ref_file_suffix(
        self, gaussian_folder, patched_gaussian_output
    ):
        ref_file = str(gaussian_folder / "rxndias_p1.log")
        folder = GaussianDIASLogFolder(
            folder=str(gaussian_folder),
            atom1=1,
            atom2=2,
            zero=True,
            ref_file=ref_file,
        )
        assert folder.outputname == "dias_zero_ref_user"


class TestGetData:
    def test_no_zero_reference_returns_raw_relative_energies(
        self, gaussian_dias_folder
    ):
        rel_total, strain, interaction = gaussian_dias_folder.get_data()
        assert rel_total == pytest.approx(
            gaussian_dias_folder.list_rel_total_energies
        )
        assert strain == pytest.approx(
            gaussian_dias_folder.list_total_strain_energies
        )
        assert interaction == pytest.approx(
            gaussian_dias_folder.list_total_interaction_energies
        )

    def test_zero_reference_without_ref_file_uses_min(
        self, gaussian_dias_folder
    ):
        lowest = min(gaussian_dias_folder.list_rel_total_energies)
        gaussian_dias_folder.zero = True
        rel_total, strain, interaction = gaussian_dias_folder.get_data()
        assert rel_total == pytest.approx(
            [e - lowest for e in gaussian_dias_folder.list_rel_total_energies]
        )
        assert min(rel_total) == pytest.approx(0.0)

    def test_zero_reference_with_ref_file_uses_ref_energy(
        self, gaussian_folder, patched_gaussian_output
    ):
        ref_file = str(gaussian_folder / "rxndias_p1.log")
        folder = GaussianDIASLogFolder(
            folder=str(gaussian_folder),
            atom1=1,
            atom2=2,
            zero=True,
            ref_file=ref_file,
        )
        ref_energy = folder.ref_file_rel_energy
        rel_total, _, _ = folder.get_data()
        assert rel_total == pytest.approx(
            [e - ref_energy for e in folder.list_rel_total_energies]
        )


class TestWriteData:
    def test_write_data_creates_expected_file(self, gaussian_dias_folder):
        gaussian_dias_folder.write_data()
        outfile = os.path.join(
            gaussian_dias_folder.folder,
            f"{gaussian_dias_folder.job_basename}_"
            f"{gaussian_dias_folder.outputname}_data.txt",
        )
        assert os.path.isfile(outfile)
        with open(outfile) as f:
            lines = f.readlines()
        assert lines[0].startswith("#")
        # header + one line per IRC point
        assert len(lines) == 1 + gaussian_dias_folder.irc_points


class TestPlotDias:
    def test_plot_dias_with_extrapolation(self, gaussian_dias_folder):
        with patch("chemsmart.analysis.dias.plt.show"):
            gaussian_dias_folder.plot_dias(
                extrapolate=True, new_length=10, k=1
            )
        outfile = os.path.join(
            gaussian_dias_folder.folder,
            f"{gaussian_dias_folder.job_basename}_"
            f"{gaussian_dias_folder.outputname}_plot.pdf",
        )
        assert os.path.isfile(outfile)

    def test_plot_dias_without_extrapolation(self, gaussian_dias_folder):
        with patch("chemsmart.analysis.dias.plt.show"):
            gaussian_dias_folder.plot_dias(extrapolate=False, reversed=False)
        outfile = os.path.join(
            gaussian_dias_folder.folder,
            f"{gaussian_dias_folder.job_basename}_"
            f"{gaussian_dias_folder.outputname}_plot.pdf",
        )
        assert os.path.isfile(outfile)

    def test_plot_dias_also_plots_ref_file_point(
        self, gaussian_folder, patched_gaussian_output
    ):
        ref_file = str(gaussian_folder / "rxndias_p1.log")
        folder = GaussianDIASLogFolder(
            folder=str(gaussian_folder),
            atom1=1,
            atom2=2,
            ref_file=ref_file,
        )
        with patch("chemsmart.analysis.dias.plt.show"):
            folder.plot_dias(extrapolate=False)
        outfile = os.path.join(
            folder.folder,
            f"{folder.job_basename}_{folder.outputname}_plot.pdf",
        )
        assert os.path.isfile(outfile)


class TestGaussianDIASBugFullMoleculeGroupIndex:
    """Regression test documenting a real crash in production code.

    ``gaussian_dias_filename_point_without_fragment_without_reactant``
    (chemsmart/utils/repattern.py) is a two-alternative regex: the
    first alternative requires the literal substring "dias_p" in the
    filename and populates capture group 1 with the point number; the
    second (fallback, generic "_p<digits>") alternative populates
    group 3 instead, leaving group 1 as ``None``.
    ``GaussianDIASLogFolder._get_all_files_full_molecule``
    (chemsmart/analysis/dias.py:439) unconditionally does
    ``int(match.group(1))``, so it crashes with a ``TypeError`` for
    any full-molecule filename that doesn't literally contain
    "dias_p" -- which is exactly how real jobs name these files
    (``f"{label}_p{i}.log"``, see chemsmart/jobs/gaussian/dias.py:316)
    unless the user's job label happens to contain "dias".
    """

    def test_full_molecule_scan_crashes_without_dias_in_filename(
        self, tmp_path
    ):
        file_map = {
            "rxn_p1.log": (-1.00, 1.0),
            "rxn_p1_f1.log": (-0.50, 1.0),
            "rxn_p1_f2.log": (-0.40, 1.0),
            "rxn_r1.log": (-0.48, 1.0),
            "rxn_r2.log": (-0.39, 1.0),
        }
        for name in file_map:
            _touch(tmp_path / name)
        with patch(
            "chemsmart.analysis.dias.Gaussian16Output",
            side_effect=_output_factory(file_map),
        ):
            with pytest.raises(TypeError):
                GaussianDIASLogFolder(folder=str(tmp_path), atom1=1, atom2=2)


class TestORCADIASOutFolderFileDiscovery:
    def test_full_molecule_files_sorted_by_point_index(self, orca_dias_folder):
        names = [
            os.path.basename(f)
            for f in orca_dias_folder._get_all_files_full_molecule()
        ]
        assert names == ["rxn_p1.out", "rxn_p2.out"]

    def test_fragment1_files_sorted_by_point_index(self, orca_dias_folder):
        names = [
            os.path.basename(f)
            for f in orca_dias_folder._get_all_files_fragment1()
        ]
        assert names == ["rxn_p1_f1.out", "rxn_p2_f1.out"]

    def test_fragment2_files_sorted_by_point_index(self, orca_dias_folder):
        names = [
            os.path.basename(f)
            for f in orca_dias_folder._get_all_files_fragment2()
        ]
        assert names == ["rxn_p1_f2.out", "rxn_p2_f2.out"]

    def test_reactant_files_sorted_by_reactant_index(self, orca_dias_folder):
        names = [
            os.path.basename(f)
            for f in orca_dias_folder._get_all_files_reactants()
        ]
        assert names == ["rxn_r1.out", "rxn_r2.out"]

    def test_job_basename_from_fragment1_file(self, orca_dias_folder):
        assert orca_dias_folder.job_basename == "rxn"


class TestORCADIASOutFolderEnergies:
    def test_get_all_energies_uses_final_energy(self, orca_dias_folder):
        energies = orca_dias_folder._get_all_energies(
            ["rxn_p1.out", "rxn_p2.out"]
        )
        assert energies == pytest.approx(
            [_to_kcal_per_mol(-1.00), _to_kcal_per_mol(-1.10)]
        )

    def test_interaction_energies(self, orca_dias_folder):
        full = orca_dias_folder.full_molecule_energies_list
        f1 = orca_dias_folder.fragment1_energies_list
        f2 = orca_dias_folder.fragment2_energies_list
        expected = [full[i] - f1[i] - f2[i] for i in range(2)]
        assert (
            orca_dias_folder.list_total_interaction_energies
            == pytest.approx(expected)
        )
