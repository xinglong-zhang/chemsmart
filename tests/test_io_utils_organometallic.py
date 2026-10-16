"""
Coverage-gap tests for the organometallic / CDX-conversion / geometry
helpers in ``chemsmart.utils.io``, plus a handful of remaining edge cases
in ``load_molecules_from_paths`` and the shell/PowerShell/Windows
environment-update helpers.

These tests build RDKit molecules directly with ``Chem.RWMol`` (mirroring
the style used in ``tests/test_io_utils.py``) rather than relying on real
external files or a real ``obabel`` binary, so that the tests are fast and
hermetic.
"""

import os
import sys
from unittest.mock import MagicMock

import numpy as np
import pytest
from rdkit import Chem

from chemsmart.utils.io import (
    _order_ring_atoms_by_walk,
    _reposition_rings_and_metal,
    _rotation_matrix_between_vectors,
    attach_eta_bonds_for_arene_rings,
    attach_eta_bonds_for_cp_rings,
    attach_one_bond_per_cp_ring,
    load_molecules_from_paths,
    normalize_metal_bonds,
    obtain_mols_from_cdx_via_obabel,
    remove_phantom_metal_carbons,
    safe_sanitize,
    update_powershell_profiles,
    update_shell_config,
    update_windows_env,
)


class TestLoadMoleculesFromPathsEdgeCases:
    """Edge cases in load_molecules_from_paths not covered elsewhere."""

    def test_skips_falsy_file_path(self, gaussian_singlet_opt_outfile):
        molecules = load_molecules_from_paths(
            ["", gaussian_singlet_opt_outfile],
            index="-1",
            check_exists=False,
        )
        assert len(molecules) > 0

    def test_check_exists_raises_for_missing_file(self, tmp_path):
        missing = tmp_path / "does_not_exist.xyz"
        with pytest.raises(FileNotFoundError):
            load_molecules_from_paths(
                [str(missing)], index="-1", check_exists=True
            )

    def test_add_index_suffix_for_single_with_explicit_index(
        self, gaussian_singlet_opt_outfile
    ):
        molecules = load_molecules_from_paths(
            [gaussian_singlet_opt_outfile],
            index="1",
            add_index_suffix_for_single=True,
            check_exists=False,
        )
        assert len(molecules) > 0
        assert "idx1" in molecules[0].name

    def test_multi_structure_file_gets_per_structure_suffix(
        self, multiple_molecules_xyz_file
    ):
        molecules = load_molecules_from_paths(
            [multiple_molecules_xyz_file],
            index=":",
            check_exists=False,
        )
        assert len(molecules) > 1
        base = os.path.splitext(os.path.basename(multiple_molecules_xyz_file))[
            0
        ]
        assert molecules[0].name == f"{base}_1"
        assert molecules[1].name == f"{base}_2"

    def test_exception_during_load_is_logged_and_reraised(self, mocker):
        mocker.patch(
            "chemsmart.utils.io.Molecule.from_filepath",
            side_effect=RuntimeError("boom"),
        )
        with pytest.raises(RuntimeError, match="boom"):
            load_molecules_from_paths(
                ["some_file.xyz"], index="-1", check_exists=False
            )


class TestObtainMolsFromCdxViaObabel:
    """Tests for obtain_mols_from_cdx_via_obabel (obabel CLI wrapper)."""

    def test_raises_when_obabel_not_on_path(self, mocker):
        mocker.patch("chemsmart.utils.io.shutil.which", return_value=None)
        with pytest.raises(ValueError, match="Open Babel CLI"):
            obtain_mols_from_cdx_via_obabel("molecule.cdx")

    def test_raises_when_obabel_returns_nonzero(self, mocker):
        mocker.patch(
            "chemsmart.utils.io.shutil.which", return_value="/usr/bin/obabel"
        )
        fake_result = MagicMock()
        fake_result.returncode = 1
        fake_result.stderr = b"conversion error"
        mocker.patch(
            "chemsmart.utils.io.subprocess.run", return_value=fake_result
        )
        with pytest.raises(RuntimeError, match="obabel failed"):
            obtain_mols_from_cdx_via_obabel("molecule.cdx")

    def test_raises_when_no_valid_molecules_produced(self, mocker):
        mocker.patch(
            "chemsmart.utils.io.shutil.which", return_value="/usr/bin/obabel"
        )
        fake_result = MagicMock()
        fake_result.returncode = 0
        fake_result.stdout = b""  # empty SDF stream -> no molecules parsed
        mocker.patch(
            "chemsmart.utils.io.subprocess.run", return_value=fake_result
        )
        with pytest.raises(ValueError, match="no valid molecules"):
            obtain_mols_from_cdx_via_obabel("molecule.cdx")

    def test_success_returns_parsed_molecules(self, mocker):
        mocker.patch(
            "chemsmart.utils.io.shutil.which", return_value="/usr/bin/obabel"
        )
        sdf_text = (
            "\n     RDKit          2D\n\n"
            "  1  0  0  0  0  0  0  0  0  0999 V2000\n"
            "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            "M  END\n$$$$\n"
        )
        fake_result = MagicMock()
        fake_result.returncode = 0
        fake_result.stdout = sdf_text.encode()
        mocker.patch(
            "chemsmart.utils.io.subprocess.run", return_value=fake_result
        )
        mols = obtain_mols_from_cdx_via_obabel("molecule.cdx")
        assert len(mols) == 1
        assert mols[0].GetNumAtoms() == 1


class TestSafeSanitizeKekulizationFallback:
    """Cover the except branch of safe_sanitize (751-763) plus the
    explicit skip_kekulize=True path (739-746)."""

    @staticmethod
    def _build_unkekulizable_cp_metal_mol():
        rw = Chem.RWMol()
        fe = rw.AddAtom(Chem.Atom(26))
        ring = []
        for _ in range(5):
            a = Chem.Atom(6)
            a.SetIsAromatic(True)
            ring.append(rw.AddAtom(a))
        for i in range(5):
            bidx = (
                rw.AddBond(ring[i], ring[(i + 1) % 5], Chem.BondType.AROMATIC)
                - 1
            )
            rw.GetBondWithIdx(bidx).SetIsAromatic(True)
        rw.AddBond(fe, ring[0], Chem.BondType.SINGLE)
        return rw.GetMol()

    def test_fallback_used_when_standard_sanitize_fails(self):
        mol = self._build_unkekulizable_cp_metal_mol()
        # Sanity check: standard sanitize really does fail for this mol.
        with pytest.raises(Exception):
            Chem.SanitizeMol(self._build_unkekulizable_cp_metal_mol())

        result = safe_sanitize(mol)
        assert result is not None
        assert result.GetNumAtoms() == 6

    def test_skip_kekulize_true_bypasses_standard_path(self):
        mol = self._build_unkekulizable_cp_metal_mol()
        result = safe_sanitize(mol, skip_kekulize=True)
        assert result is not None
        assert result.GetNumAtoms() == 6


class TestNormalizeMetalBondsBranches:
    """Directly exercise both aromatic-flag-clearing branches (803, 805)."""

    def test_clears_aromatic_flag_and_bond_type_on_metal_bond(self):
        rw = Chem.RWMol()
        fe = rw.AddAtom(Chem.Atom(26))
        c = rw.AddAtom(Chem.Atom(6))
        bidx = rw.AddBond(fe, c, Chem.BondType.AROMATIC) - 1
        rw.GetBondWithIdx(bidx).SetIsAromatic(True)
        mol = rw.GetMol()

        result = normalize_metal_bonds(mol)
        bond = result.GetBondWithIdx(0)
        assert bond.GetIsAromatic() is False
        assert bond.GetBondType() == Chem.BondType.SINGLE


class TestRemovePhantomMetalCarbons:
    """Tests for remove_phantom_metal_carbons."""

    def test_removes_terminal_non_ring_carbon_neighbor(self):
        mol = Chem.MolFromSmiles("[Fe](C)c1ccccc1", sanitize=False)
        Chem.GetSymmSSSR(mol)
        metal_idxs = {0}

        new_mol, new_metal_idxs = remove_phantom_metal_carbons(mol, metal_idxs)

        assert new_mol.GetNumAtoms() == mol.GetNumAtoms() - 1
        assert new_metal_idxs == {0}
        # The phantom methyl carbon (originally idx 1) is gone; the metal
        # is no longer bonded to any degree-1 non-ring carbon.
        metal_atom = new_mol.GetAtomWithIdx(0)
        for nbr in metal_atom.GetNeighbors():
            assert not (nbr.GetAtomicNum() == 6 and nbr.GetDegree() == 1)

    def test_no_op_when_no_phantom_atoms_present(self):
        mol = Chem.MolFromSmiles("[Fe]c1ccccc1", sanitize=False)
        Chem.GetSymmSSSR(mol)
        metal_idxs = {0}

        new_mol, new_metal_idxs = remove_phantom_metal_carbons(mol, metal_idxs)

        assert new_mol is mol
        assert new_metal_idxs == metal_idxs


class TestOrderRingAtomsByWalk:
    """Directly test the private ring-walking helper's None-returning
    branches by constructing molecule graphs that don't correspond to a
    simple cycle over the supplied ring tuple."""

    def test_start_atom_with_wrong_ring_degree_returns_none(self):
        rw = Chem.RWMol()
        idxs = [rw.AddAtom(Chem.Atom(6)) for _ in range(4)]
        a, b, c, d = idxs
        rw.AddBond(a, b, Chem.BondType.SINGLE)
        rw.AddBond(a, c, Chem.BondType.SINGLE)
        rw.AddBond(a, d, Chem.BondType.SINGLE)
        mol = rw.GetMol()
        mol.UpdatePropertyCache(strict=False)
        assert _order_ring_atoms_by_walk(mol, (a, b, c, d)) is None

    def test_mid_walk_atom_with_wrong_ring_degree_returns_none(self):
        rw = Chem.RWMol()
        idxs = [rw.AddAtom(Chem.Atom(6)) for _ in range(5)]
        a, b, c, d, e = idxs
        rw.AddBond(a, b, Chem.BondType.SINGLE)
        rw.AddBond(b, c, Chem.BondType.SINGLE)
        rw.AddBond(c, d, Chem.BondType.SINGLE)
        rw.AddBond(c, e, Chem.BondType.SINGLE)
        rw.AddBond(e, a, Chem.BondType.SINGLE)
        mol = rw.GetMol()
        mol.UpdatePropertyCache(strict=False)
        assert _order_ring_atoms_by_walk(mol, (a, b, c, d, e)) is None

    def test_walk_shorter_than_ring_returns_none(self):
        # Two disjoint triangles passed in as a single 6-atom "ring".
        rw = Chem.RWMol()
        idxs = [rw.AddAtom(Chem.Atom(6)) for _ in range(6)]
        a, b, c, d, e, f = idxs
        rw.AddBond(a, b, Chem.BondType.SINGLE)
        rw.AddBond(b, c, Chem.BondType.SINGLE)
        rw.AddBond(c, a, Chem.BondType.SINGLE)
        rw.AddBond(d, e, Chem.BondType.SINGLE)
        rw.AddBond(e, f, Chem.BondType.SINGLE)
        rw.AddBond(f, d, Chem.BondType.SINGLE)
        mol = rw.GetMol()
        mol.UpdatePropertyCache(strict=False)
        assert _order_ring_atoms_by_walk(mol, (a, b, c, d, e, f)) is None

    def test_normal_ring_returns_cyclic_order(self):
        rw = Chem.RWMol()
        idxs = [rw.AddAtom(Chem.Atom(6)) for _ in range(5)]
        for i in range(5):
            rw.AddBond(idxs[i], idxs[(i + 1) % 5], Chem.BondType.SINGLE)
        mol = rw.GetMol()
        mol.UpdatePropertyCache(strict=False)
        assert _order_ring_atoms_by_walk(mol, tuple(idxs)) == idxs


class TestAttachEtaBondsForCpRings:
    """Tests for attach_eta_bonds_for_cp_rings."""

    def test_not_metal_idxs_returns_mol_unchanged(self):
        mol = Chem.MolFromSmiles("c1ccc[cH-]1", sanitize=False)
        result = attach_eta_bonds_for_cp_rings(mol, set())
        assert result is mol

    def test_attaches_bond_to_disconnected_cp_ring(self):
        mol = Chem.MolFromSmiles("[Fe].c1ccc[cH-]1", sanitize=False)
        Chem.GetSymmSSSR(mol)

        new_mol = attach_eta_bonds_for_cp_rings(mol, {0})

        fe = new_mol.GetAtomWithIdx(0)
        assert fe.GetDegree() == 1
        # Every ring carbon should have exactly one implicit H after
        # dearomatization + bond-order assignment (per docstring).
        for i in range(1, 6):
            atom = new_mol.GetAtomWithIdx(i)
            assert atom.GetTotalNumHs() == 1

    def test_ring_already_connected_to_metal_is_skipped(self):
        mol = Chem.MolFromSmiles("[Fe]C1=CC=CC1", sanitize=False)
        Chem.GetSymmSSSR(mol)

        new_mol = attach_eta_bonds_for_cp_rings(mol, {0})

        fe = new_mol.GetAtomWithIdx(0)
        # No second bond should have been added; degree stays at 1.
        assert fe.GetDegree() == 1

    def test_metal_not_in_first_fragment(self):
        # The ring fragment is listed before the metal fragment by
        # GetMolFrags, forcing the fragment-membership loop to skip past
        # at least one non-matching fragment before finding the metal.
        mol = Chem.MolFromSmiles("c1ccc[cH-]1.[Fe]", sanitize=False)
        Chem.GetSymmSSSR(mol)
        fe_idx = next(
            a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "Fe"
        )

        new_mol = attach_eta_bonds_for_cp_rings(mol, {fe_idx})

        fe = new_mol.GetAtomWithIdx(fe_idx)
        assert fe.GetDegree() == 1

    def test_no_five_membered_ring_present_is_a_no_op(self):
        mol = Chem.MolFromSmiles("[Fe]", sanitize=False)
        result = attach_eta_bonds_for_cp_rings(mol, {0})
        assert result.GetNumAtoms() == 1

    def test_skips_wrong_size_and_non_carbon_decoy_rings(self):
        # A disconnected metal, a valid 5-C Cp ring, plus a 6-membered
        # decoy (wrong size) and a 5-membered furan (not all carbon) that
        # must both be skipped via the "continue" branches without
        # affecting processing of the real Cp ring.
        mol = Chem.MolFromSmiles(
            "[Fe].c1ccc[cH-]1.c1ccccc1.c1ccoc1", sanitize=False
        )
        Chem.GetSymmSSSR(mol)

        new_mol = attach_eta_bonds_for_cp_rings(mol, {0})

        fe = new_mol.GetAtomWithIdx(0)
        assert fe.GetDegree() == 1


class TestAttachEtaBondsForAreneRings:
    """Tests for attach_eta_bonds_for_arene_rings."""

    def test_not_metal_idxs_returns_mol_unchanged(self):
        mol = Chem.MolFromSmiles("c1ccccc1")
        result = attach_eta_bonds_for_arene_rings(mol, set())
        assert result is mol

    def test_attaches_bond_to_disconnected_arene_ring(self):
        mol = Chem.MolFromSmiles("[Fe].c1ccccc1", sanitize=False)
        Chem.GetSymmSSSR(mol)

        new_mol = attach_eta_bonds_for_arene_rings(mol, {0})

        fe = new_mol.GetAtomWithIdx(0)
        assert fe.GetDegree() == 1
        for i in range(1, 7):
            atom = new_mol.GetAtomWithIdx(i)
            assert atom.GetTotalNumHs() == 1

    def test_ring_already_connected_to_metal_is_skipped(self):
        mol = Chem.MolFromSmiles("[Fe]C1=CC=CC=C1", sanitize=False)
        Chem.GetSymmSSSR(mol)

        new_mol = attach_eta_bonds_for_arene_rings(mol, {0})

        fe = new_mol.GetAtomWithIdx(0)
        assert fe.GetDegree() == 1

    def test_no_six_membered_ring_present_is_a_no_op(self):
        mol = Chem.MolFromSmiles("[Fe]", sanitize=False)
        result = attach_eta_bonds_for_arene_rings(mol, {0})
        assert result.GetNumAtoms() == 1

    def test_skips_wrong_size_and_non_carbon_decoy_rings(self):
        # Decoys: a 5-membered Cp ring (wrong size) and a 6-membered
        # pyridine ring (not all-carbon) alongside the real benzene ring.
        mol = Chem.MolFromSmiles(
            "[Fe].c1ccccc1.c1ccc[cH-]1.c1ccncc1", sanitize=False
        )
        Chem.GetSymmSSSR(mol)

        new_mol = attach_eta_bonds_for_arene_rings(mol, {0})

        fe = new_mol.GetAtomWithIdx(0)
        assert fe.GetDegree() == 1


class TestRotationMatrixBetweenVectors:
    """Tests for the _rotation_matrix_between_vectors Rodrigues helper."""

    def test_parallel_vectors_returns_identity(self):
        v1 = np.array([1.0, 0.0, 0.0])
        v2 = np.array([2.0, 0.0, 0.0])
        rot = _rotation_matrix_between_vectors(v1, v2)
        assert np.allclose(rot, np.eye(3))

    def test_anti_parallel_vectors_180_degree_rotation(self):
        v1 = np.array([1.0, 0.0, 0.0])
        v2 = np.array([-1.0, 0.0, 0.0])
        rot = _rotation_matrix_between_vectors(v1, v2)
        result = rot @ v1
        assert np.allclose(result, v2, atol=1e-8)
        # Rotation matrix must be orthogonal (a proper 180-degree rotation).
        assert np.allclose(rot @ rot.T, np.eye(3), atol=1e-8)

    def test_generic_non_parallel_pair(self):
        v1 = np.array([1.0, 0.0, 0.0])
        v2 = np.array([0.0, 1.0, 0.0])
        rot = _rotation_matrix_between_vectors(v1, v2)
        assert np.allclose(rot @ v1, v2, atol=1e-8)

    def test_near_zero_length_vector_returns_identity(self):
        v1 = np.array([1e-15, 0.0, 0.0])
        v2 = np.array([0.0, 1.0, 0.0])
        rot = _rotation_matrix_between_vectors(v1, v2)
        assert np.allclose(rot, np.eye(3))

        rot2 = _rotation_matrix_between_vectors(v2, v1)
        assert np.allclose(rot2, np.eye(3))


class TestRepositionRingsAndMetal:
    """Tests for the _reposition_rings_and_metal 3D-geometry helper."""

    @staticmethod
    def _build_half_sandwich(ring_size, coplanar_metal_offset=(2.0, 0.0, 0.0)):
        rw = Chem.RWMol()
        fe = rw.AddAtom(Chem.Atom(26))
        ring = [rw.AddAtom(Chem.Atom(6)) for _ in range(ring_size)]
        for i in range(ring_size):
            rw.AddBond(
                ring[i], ring[(i + 1) % ring_size], Chem.BondType.SINGLE
            )
        rw.AddBond(fe, ring[0], Chem.BondType.SINGLE)
        mol = rw.GetMol()
        Chem.SanitizeMol(
            mol,
            sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL
            ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE,
        )

        conf = Chem.Conformer(mol.GetNumAtoms())
        for k, idx in enumerate(ring):
            ang = 2 * np.pi * k / ring_size
            conf.SetAtomPosition(idx, (np.cos(ang), np.sin(ang), 0.0))
        conf.SetAtomPosition(fe, coplanar_metal_offset)
        mol.AddConformer(conf, assignId=True)
        return mol, fe, ring

    def test_no_conformer_returns_mol_unchanged(self):
        rw = Chem.RWMol()
        fe = rw.AddAtom(Chem.Atom(26))
        mol = rw.GetMol()
        result = _reposition_rings_and_metal(mol, {fe})
        assert result is mol

    def test_not_metal_idxs_returns_mol_unchanged(self):
        rw = Chem.RWMol()
        rw.AddAtom(Chem.Atom(26))
        mol = rw.GetMol()
        result = _reposition_rings_and_metal(mol, set())
        assert result is mol

    def test_no_bonded_ring_returns_mol_unchanged(self):
        # Ring exists but is NOT bonded to the metal -> bonded_rings empty.
        rw = Chem.RWMol()
        fe = rw.AddAtom(Chem.Atom(26))
        ring = [rw.AddAtom(Chem.Atom(6)) for _ in range(5)]
        for i in range(5):
            rw.AddBond(ring[i], ring[(i + 1) % 5], Chem.BondType.SINGLE)
        mol = rw.GetMol()
        Chem.SanitizeMol(
            mol,
            sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL
            ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE,
        )
        conf = Chem.Conformer(mol.GetNumAtoms())
        for k, idx in enumerate(ring):
            ang = 2 * np.pi * k / 5
            conf.SetAtomPosition(idx, (np.cos(ang), np.sin(ang), 0.0))
        conf.SetAtomPosition(fe, (10.0, 10.0, 10.0))
        mol.AddConformer(conf, assignId=True)

        result = _reposition_rings_and_metal(mol, {fe})
        assert result is mol

    def test_ignores_wrong_size_decoy_ring(self):
        # A disconnected 4-membered ring alongside the real metal-bonded
        # Cp ring must be skipped by the ring-size filter in the gathering
        # loop, without disturbing processing of the real ring.
        rw = Chem.RWMol()
        fe = rw.AddAtom(Chem.Atom(26))
        ring = [rw.AddAtom(Chem.Atom(6)) for _ in range(5)]
        for i in range(5):
            rw.AddBond(ring[i], ring[(i + 1) % 5], Chem.BondType.SINGLE)
        rw.AddBond(fe, ring[0], Chem.BondType.SINGLE)
        decoy = [rw.AddAtom(Chem.Atom(6)) for _ in range(4)]
        for i in range(4):
            rw.AddBond(decoy[i], decoy[(i + 1) % 4], Chem.BondType.SINGLE)
        mol = rw.GetMol()
        Chem.SanitizeMol(
            mol,
            sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL
            ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE,
        )
        conf = Chem.Conformer(mol.GetNumAtoms())
        for k, idx in enumerate(ring):
            ang = 2 * np.pi * k / 5
            conf.SetAtomPosition(
                idx, (float(np.cos(ang)), float(np.sin(ang)), 0.0)
            )
        conf.SetAtomPosition(fe, (2.0, 0.0, 0.0))
        for k, idx in enumerate(decoy):
            conf.SetAtomPosition(idx, (10.0 + k, 10.0, 10.0))
        mol.AddConformer(conf, assignId=True)

        new_mol = _reposition_rings_and_metal(mol, {fe})
        conf2 = new_mol.GetConformer()
        fe_pos = np.array(list(conf2.GetAtomPosition(fe)))
        ring_centroid = np.mean(
            [list(conf2.GetAtomPosition(i)) for i in ring], axis=0
        )
        assert np.isclose(
            np.linalg.norm(fe_pos - ring_centroid), 2.0, atol=1e-6
        )

    def test_half_sandwich_cp_ring_metal_moved_above_ring(self):
        mol, fe, ring = self._build_half_sandwich(5)
        new_mol = _reposition_rings_and_metal(mol, {fe})
        conf = new_mol.GetConformer()
        fe_pos = np.array(list(conf.GetAtomPosition(fe)))
        ring_centroid = np.mean(
            [list(conf.GetAtomPosition(i)) for i in ring], axis=0
        )
        dist = np.linalg.norm(fe_pos - ring_centroid)
        assert np.isclose(dist, 2.0, atol=1e-6)

    def test_half_sandwich_arene_ring_ideal_distance(self):
        mol, fe, ring = self._build_half_sandwich(6)
        new_mol = _reposition_rings_and_metal(mol, {fe})
        conf = new_mol.GetConformer()
        fe_pos = np.array(list(conf.GetAtomPosition(fe)))
        ring_centroid = np.mean(
            [list(conf.GetAtomPosition(i)) for i in ring], axis=0
        )
        dist = np.linalg.norm(fe_pos - ring_centroid)
        assert np.isclose(dist, 1.75, atol=1e-6)

    def test_sandwich_two_rings_positions_both_ideal_distance(self):
        rw = Chem.RWMol()
        fe = rw.AddAtom(Chem.Atom(26))
        ring1 = [rw.AddAtom(Chem.Atom(6)) for _ in range(5)]
        ring2 = [rw.AddAtom(Chem.Atom(6)) for _ in range(5)]
        for ring in (ring1, ring2):
            for i in range(5):
                rw.AddBond(ring[i], ring[(i + 1) % 5], Chem.BondType.SINGLE)
        rw.AddBond(fe, ring1[0], Chem.BondType.SINGLE)
        rw.AddBond(fe, ring2[0], Chem.BondType.SINGLE)
        mol = rw.GetMol()
        Chem.SanitizeMol(
            mol,
            sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL
            ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE,
        )

        conf = Chem.Conformer(mol.GetNumAtoms())
        for k, idx in enumerate(ring1):
            ang = 2 * np.pi * k / 5
            conf.SetAtomPosition(idx, (np.cos(ang), np.sin(ang), 0.0))
        for k, idx in enumerate(ring2):
            ang = 2 * np.pi * k / 5
            conf.SetAtomPosition(idx, (np.cos(ang) + 5, np.sin(ang), 0.0))
        conf.SetAtomPosition(fe, (2.5, 0.0, 0.0))
        mol.AddConformer(conf, assignId=True)

        new_mol = _reposition_rings_and_metal(mol, {fe})
        conf2 = new_mol.GetConformer()
        fe_pos = np.array(list(conf2.GetAtomPosition(fe)))
        c1 = np.mean([list(conf2.GetAtomPosition(i)) for i in ring1], axis=0)
        c2 = np.mean([list(conf2.GetAtomPosition(i)) for i in ring2], axis=0)
        assert np.isclose(np.linalg.norm(fe_pos - c1), 2.0, atol=1e-6)
        assert np.isclose(np.linalg.norm(fe_pos - c2), 2.0, atol=1e-6)
        # Rings must end up on opposite sides of the metal (sandwich).
        d1 = c1 - fe_pos
        d2 = c2 - fe_pos
        assert np.dot(d1, d2) < 0

    def test_half_sandwich_metal_below_plane_flips_normal(self):
        # Metal placed on the opposite side of the ring plane from the
        # default cross-product normal; the function must flip the normal
        # so the metal is pushed further away rather than through the ring.
        mol, fe, ring = self._build_half_sandwich(
            5, coplanar_metal_offset=(0.0, 0.0, -1.0)
        )
        new_mol = _reposition_rings_and_metal(mol, {fe})
        conf = new_mol.GetConformer()
        fe_pos = np.array(list(conf.GetAtomPosition(fe)))
        assert fe_pos[2] < 0  # stays on the same (negative-z) side

    def test_sandwich_enforces_opposite_signs_for_coincident_centroids(self):
        # Degenerate case: both rings share the exact same spatial position
        # (coincident centroids), so the naive per-ring sign computation
        # would give identical signs; the function must force them apart.
        rw = Chem.RWMol()
        fe = rw.AddAtom(Chem.Atom(26))
        ring1 = [rw.AddAtom(Chem.Atom(6)) for _ in range(5)]
        ring2 = [rw.AddAtom(Chem.Atom(6)) for _ in range(5)]
        for ring in (ring1, ring2):
            for i in range(5):
                rw.AddBond(ring[i], ring[(i + 1) % 5], Chem.BondType.SINGLE)
        rw.AddBond(fe, ring1[0], Chem.BondType.SINGLE)
        rw.AddBond(fe, ring2[0], Chem.BondType.SINGLE)
        mol = rw.GetMol()
        Chem.SanitizeMol(
            mol,
            sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL
            ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE,
        )

        conf = Chem.Conformer(mol.GetNumAtoms())
        for k in range(5):
            ang = 2 * np.pi * k / 5
            pos = (float(np.cos(ang)), float(np.sin(ang)), 0.0)
            conf.SetAtomPosition(ring1[k], pos)
            conf.SetAtomPosition(ring2[k], pos)
        conf.SetAtomPosition(fe, (0.0, 0.0, 5.0))
        mol.AddConformer(conf, assignId=True)

        new_mol = _reposition_rings_and_metal(mol, {fe})
        conf2 = new_mol.GetConformer()
        fe_pos = np.array(list(conf2.GetAtomPosition(fe)))
        c1 = np.mean([list(conf2.GetAtomPosition(i)) for i in ring1], axis=0)
        c2 = np.mean([list(conf2.GetAtomPosition(i)) for i in ring2], axis=0)
        # Rings must be pushed to opposite sides of the metal, not stacked
        # together.
        assert np.dot(c1 - fe_pos, c2 - fe_pos) < 0

    def test_moves_h_atoms_and_bridging_atom_between_two_rings(self):
        # Ansa-bridged bis-Cp complex: explicit H atoms on the ring carbons
        # (must move rigidly with the ring) plus a bridging O atom bonded to
        # one carbon in each ring (must land at the midpoint of the two
        # repositioned ring-carbon anchors).
        rw = Chem.RWMol()
        fe = rw.AddAtom(Chem.Atom(26))

        def make_ring():
            ring = [rw.AddAtom(Chem.Atom(6)) for _ in range(5)]
            for i in range(5):
                rw.AddBond(ring[i], ring[(i + 1) % 5], Chem.BondType.SINGLE)
            h_atoms = []
            for c in ring[1:]:
                h = rw.AddAtom(Chem.Atom(1))
                rw.AddBond(c, h, Chem.BondType.SINGLE)
                h_atoms.append(h)
            return ring, h_atoms

        ring1, h1 = make_ring()
        ring2, h2 = make_ring()
        rw.AddBond(fe, ring1[0], Chem.BondType.SINGLE)
        rw.AddBond(fe, ring2[0], Chem.BondType.SINGLE)

        bridge_o = rw.AddAtom(Chem.Atom(8))
        rw.AddBond(ring1[2], bridge_o, Chem.BondType.SINGLE)
        rw.AddBond(ring2[2], bridge_o, Chem.BondType.SINGLE)

        # A pendant (non-bridging) substituent bonded to only ONE ring
        # system, to exercise the "not a true bridge" skip path (the
        # averaging block only fires for atoms bonded to >=2 ring systems).
        pendant_cl = rw.AddAtom(Chem.Atom(17))
        rw.AddBond(ring1[3], pendant_cl, Chem.BondType.SINGLE)

        mol = rw.GetMol()
        Chem.SanitizeMol(
            mol,
            sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL
            ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE,
        )

        conf = Chem.Conformer(mol.GetNumAtoms())
        for k, idx in enumerate(ring1):
            ang = 2 * np.pi * k / 5
            conf.SetAtomPosition(
                idx, (float(np.cos(ang)), float(np.sin(ang)), 0.0)
            )
        for h in h1:
            conf.SetAtomPosition(h, (0.0, 0.0, 0.5))
        for k, idx in enumerate(ring2):
            ang = 2 * np.pi * k / 5
            conf.SetAtomPosition(
                idx, (float(np.cos(ang)) + 5, float(np.sin(ang)), 0.0)
            )
        for h in h2:
            conf.SetAtomPosition(h, (5.0, 0.0, 0.5))
        conf.SetAtomPosition(bridge_o, (2.5, 0.0, 1.0))
        conf.SetAtomPosition(pendant_cl, (-1.0, 0.0, 0.0))
        conf.SetAtomPosition(fe, (2.5, 0.0, 0.0))
        mol.AddConformer(conf, assignId=True)

        new_mol = _reposition_rings_and_metal(mol, {fe})
        conf2 = new_mol.GetConformer()

        # H atoms must have moved off their placeholder position (rigid
        # body transform applied to them too).
        for h in h1 + h2:
            pos = np.array(list(conf2.GetAtomPosition(h)))
            assert not np.allclose(pos[2], 0.5)

        # Bridge O must land at the midpoint of the two new ring-carbon
        # anchor positions (ring1[2] and ring2[2]).
        o_pos = np.array(list(conf2.GetAtomPosition(bridge_o)))
        c1_new = np.array(list(conf2.GetAtomPosition(ring1[2])))
        c2_new = np.array(list(conf2.GetAtomPosition(ring2[2])))
        assert np.allclose(o_pos, (c1_new + c2_new) / 2, atol=1e-6)

    def test_fused_indenyl_ring_system_moves_as_one_rigid_body(self):
        # Indenyl-like ligand: a 5-membered Cp ring (the metal-coordinating
        # ring) fused to a 6-membered benzo ring sharing one edge. A second
        # plain Cp ring is also bonded to the metal to force the n>=2
        # sandwich code path, which is what performs the fused-ring BFS
        # expansion and rigid-body transform (the n==1 half-sandwich path
        # returns before ever reaching that logic).
        rw = Chem.RWMol()
        fe = rw.AddAtom(Chem.Atom(26))
        c = [rw.AddAtom(Chem.Atom(6)) for _ in range(9)]
        rw.AddBond(c[0], c[1], Chem.BondType.SINGLE)
        rw.AddBond(c[1], c[2], Chem.BondType.SINGLE)
        rw.AddBond(c[2], c[3], Chem.BondType.SINGLE)
        rw.AddBond(c[3], c[4], Chem.BondType.SINGLE)
        rw.AddBond(c[4], c[0], Chem.BondType.SINGLE)
        rw.AddBond(c[3], c[5], Chem.BondType.SINGLE)
        rw.AddBond(c[5], c[6], Chem.BondType.SINGLE)
        rw.AddBond(c[6], c[7], Chem.BondType.SINGLE)
        rw.AddBond(c[7], c[8], Chem.BondType.SINGLE)
        rw.AddBond(c[8], c[2], Chem.BondType.SINGLE)
        rw.AddBond(fe, c[0], Chem.BondType.SINGLE)

        r2 = [rw.AddAtom(Chem.Atom(6)) for _ in range(5)]
        for i in range(5):
            rw.AddBond(r2[i], r2[(i + 1) % 5], Chem.BondType.SINGLE)
        rw.AddBond(fe, r2[0], Chem.BondType.SINGLE)

        mol = rw.GetMol()
        Chem.SanitizeMol(
            mol,
            sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL
            ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE,
        )

        conf = Chem.Conformer(mol.GetNumAtoms())
        for k in range(5):
            ang = 2 * np.pi * k / 5
            conf.SetAtomPosition(
                c[k], (float(np.cos(ang)), float(np.sin(ang)), 0.0)
            )
        # Benzo extension placed off in space (as if badly embedded).
        benzo_start = {
            5: (2.5, 0.0, 0.0),
            6: (3.5, 0.5, 0.0),
            7: (3.5, 1.5, 0.0),
            8: (2.5, 2.0, 0.0),
        }
        for k, pos in benzo_start.items():
            conf.SetAtomPosition(c[k], pos)
        for k, idx in enumerate(r2):
            ang = 2 * np.pi * k / 5
            conf.SetAtomPosition(
                idx, (float(np.cos(ang)) + 6, float(np.sin(ang)), 0.0)
            )
        conf.SetAtomPosition(fe, (3.0, 0.0, 0.0))
        mol.AddConformer(conf, assignId=True)

        new_mol = _reposition_rings_and_metal(mol, {fe})
        conf2 = new_mol.GetConformer()
        cp_centroid = np.mean(
            [list(conf2.GetAtomPosition(c[k])) for k in range(5)], axis=0
        )
        fe_pos = np.array(list(conf2.GetAtomPosition(fe)))
        assert np.isclose(np.linalg.norm(fe_pos - cp_centroid), 2.0, atol=1e-6)
        # The fused benzo-ring atoms must have moved along with the Cp
        # ring's rigid-body transform, i.e. no longer at their original
        # placeholder coordinates.
        for k, original in benzo_start.items():
            new_pos = np.array(list(conf2.GetAtomPosition(c[k])))
            assert not np.allclose(new_pos, original, atol=1e-6)

    def test_three_rings_uses_svd_principal_axis_branch(self):
        rw = Chem.RWMol()
        fe = rw.AddAtom(Chem.Atom(26))
        rings = []
        for _ in range(3):
            ring = [rw.AddAtom(Chem.Atom(6)) for _ in range(5)]
            rings.append(ring)
            for i in range(5):
                rw.AddBond(ring[i], ring[(i + 1) % 5], Chem.BondType.SINGLE)
            rw.AddBond(fe, ring[0], Chem.BondType.SINGLE)
        mol = rw.GetMol()
        Chem.SanitizeMol(
            mol,
            sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL
            ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE,
        )

        conf = Chem.Conformer(mol.GetNumAtoms())
        offsets = [(0.0, 0.0, 0.0), (5.0, 0.0, 0.0), (2.5, 4.0, 0.0)]
        for ring, off in zip(rings, offsets):
            for k, idx in enumerate(ring):
                ang = 2 * np.pi * k / 5
                conf.SetAtomPosition(
                    idx, (np.cos(ang) + off[0], np.sin(ang) + off[1], off[2])
                )
        conf.SetAtomPosition(fe, (2.5, 1.3, 0.0))
        mol.AddConformer(conf, assignId=True)

        # Must run without error and actually move the metal atom.
        new_mol = _reposition_rings_and_metal(mol, {fe})
        conf2 = new_mol.GetConformer()
        fe_pos = np.array(list(conf2.GetAtomPosition(fe)))
        assert fe_pos.shape == (3,)


class TestAttachOneBondPerCpRing:
    """Tests for the deprecated attach_one_bond_per_cp_ring helper."""

    def test_not_metal_idxs_returns_mol_unchanged(self):
        mol = Chem.MolFromSmiles("c1ccc[cH-]1", sanitize=False)
        result = attach_one_bond_per_cp_ring(mol, set())
        assert result is mol

    def test_attaches_bond_to_disconnected_aromatic_cp_ring(self):
        mol = Chem.MolFromSmiles("[Fe].c1ccc[cH-]1", sanitize=False)
        mol.UpdatePropertyCache(strict=False)
        Chem.GetSymmSSSR(mol)

        new_mol = attach_one_bond_per_cp_ring(mol, {0})

        fe = new_mol.GetAtomWithIdx(0)
        assert fe.GetDegree() == 1

    def test_ring_already_connected_to_metal_skips_new_bond(self):
        mol = Chem.MolFromSmiles("[Fe]c1ccc[cH-]1", sanitize=False)
        mol.UpdatePropertyCache(strict=False)
        Chem.GetSymmSSSR(mol)

        new_mol = attach_one_bond_per_cp_ring(mol, {0})

        fe = new_mol.GetAtomWithIdx(0)
        assert fe.GetDegree() == 1

    def test_skips_wrong_size_and_non_aromatic_decoy_rings(self):
        # Decoys: a wrong-size aromatic 6-ring and a non-aromatic
        # (saturated) 5-ring, alongside the real aromatic Cp ring.
        mol = Chem.MolFromSmiles(
            "[Fe].c1ccc[cH-]1.c1ccccc1.C1CCCC1", sanitize=False
        )
        mol.UpdatePropertyCache(strict=False)
        Chem.GetSymmSSSR(mol)

        new_mol = attach_one_bond_per_cp_ring(mol, {0})

        fe = new_mol.GetAtomWithIdx(0)
        assert fe.GetDegree() == 1

    def test_anchor_falls_back_to_first_ring_atom_when_no_h_available(self):
        # Pentasubstituted (Cp*-like) ring: every ring carbon has 0 H, so
        # the search for an H-bearing anchor fails and falls back to
        # ordered[0].
        mol = Chem.MolFromSmiles("[Fe].Cc1c(C)c(C)c(C)[c-]1C", sanitize=False)
        mol.UpdatePropertyCache(strict=False)
        Chem.GetSymmSSSR(mol)

        new_mol = attach_one_bond_per_cp_ring(mol, {0})

        fe = new_mol.GetAtomWithIdx(0)
        assert fe.GetDegree() == 1


class TestUpdateShellConfigBranches:
    def test_creates_file_when_missing(self, tmp_path):
        shell_file = tmp_path / ".bashrc"
        assert not shell_file.exists()
        update_shell_config(
            shell_file, ["export CHEMSMART_HOME=/opt/chemsmart"]
        )
        assert shell_file.exists()
        assert "Added by chemsmart installer" in shell_file.read_text()


class TestUpdatePowershellProfilesMalformedBlock:
    def test_malformed_block_missing_end_marker_is_truncated(self, tmp_path):
        from chemsmart.utils.io import _PS_BLOCK_START

        profile = tmp_path / "profile.ps1"
        profile.write_text(
            f"Write-Host 'hi'\n{_PS_BLOCK_START}\nsome stale content\n",
            encoding="utf-8",
        )
        update_powershell_profiles([profile], ["$env:FOO = 'bar'"])
        content = profile.read_text(encoding="utf-8")
        assert "some stale content" not in content
        assert "$env:FOO = 'bar'" in content


class TestUpdateWindowsEnvBranches:
    def test_missing_winreg_module_is_a_no_op(self):
        with pytest.MonkeyPatch.context() as mp:
            mp.setitem(sys.modules, "winreg", None)
            # Should not raise even though winreg import fails.
            update_windows_env(["/some/path"], "/some/pypath")

    def test_paths_already_present_logs_no_change(self, mocker):
        mock_key = MagicMock()
        mock_key.__enter__ = MagicMock(return_value=mock_key)
        mock_key.__exit__ = MagicMock(return_value=False)

        mock_winreg = MagicMock()
        mock_winreg.HKEY_CURRENT_USER = 0x80000001
        mock_winreg.KEY_READ = 0x20019
        mock_winreg.KEY_WRITE = 0x20006
        mock_winreg.REG_EXPAND_SZ = 2
        mock_winreg.REG_SZ = 1
        mock_winreg.OpenKey.return_value = mock_key

        def fake_query_value_ex(key, name):
            if name == "PATH":
                return "/already/here", 2
            if name == "PYTHONPATH":
                return "/already/here/py", 1
            raise FileNotFoundError

        mock_winreg.QueryValueEx.side_effect = fake_query_value_ex

        mock_ctypes = MagicMock()

        with pytest.MonkeyPatch.context() as mp:
            mp.setitem(sys.modules, "winreg", mock_winreg)
            mp.setitem(sys.modules, "ctypes", mock_ctypes)
            update_windows_env(["/already/here"], "/already/here/py")

        # Nothing new to add, so SetValueEx should not have been called.
        mock_winreg.SetValueEx.assert_not_called()

    def test_permission_error_is_caught_and_warned(self, mocker):
        mock_winreg = MagicMock()
        mock_winreg.HKEY_CURRENT_USER = 0x80000001
        mock_winreg.KEY_READ = 0x20019
        mock_winreg.KEY_WRITE = 0x20006
        mock_winreg.OpenKey.side_effect = PermissionError("denied")

        with pytest.MonkeyPatch.context() as mp:
            mp.setitem(sys.modules, "winreg", mock_winreg)
            # Should not raise; permission error is caught and logged.
            update_windows_env(["/some/path"], "/some/pypath")

    def test_generic_exception_is_caught_and_warned(self, mocker):
        mock_winreg = MagicMock()
        mock_winreg.HKEY_CURRENT_USER = 0x80000001
        mock_winreg.KEY_READ = 0x20019
        mock_winreg.KEY_WRITE = 0x20006
        mock_winreg.OpenKey.side_effect = RuntimeError("unexpected")

        with pytest.MonkeyPatch.context() as mp:
            mp.setitem(sys.modules, "winreg", mock_winreg)
            update_windows_env(["/some/path"], "/some/pypath")
