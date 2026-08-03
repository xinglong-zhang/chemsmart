"""Tests for SOAP descriptor calculations on Molecule objects."""

from unittest.mock import patch

import numpy as np
import pytest

from chemsmart.analysis.soap import calculate_soap
from chemsmart.io.molecules.structure import Molecule

dscribe = pytest.importorskip("dscribe")
from ase import Atoms  # noqa: E402
from dscribe.descriptors import SOAP  # noqa: E402

# Cross-platform numerical tolerances for dense float64 SOAP vectors.
RTOL = 1e-6
ATOL = 1e-8


@pytest.fixture
def water_molecule():
    """Simple water geometry in Å."""
    return Molecule(
        symbols=["O", "H", "H"],
        positions=np.array(
            [
                [0.00000, 0.00000, 0.11779],
                [0.00000, 0.75545, -0.47116],
                [0.00000, -0.75545, -0.47116],
            ],
            dtype=float,
        ),
    )


@pytest.fixture
def methanol_molecule():
    """Methanol geometry with C/H/O species."""
    return Molecule(
        symbols=["C", "O", "H", "H", "H", "H"],
        positions=np.array(
            [
                [-0.6626, 0.0000, 0.0000],
                [0.6626, 0.0000, 0.0000],
                [-1.0280, 1.0280, 0.0000],
                [-1.0280, -0.5140, 0.8900],
                [-1.0280, -0.5140, -0.8900],
                [1.0280, 0.0000, 0.0000],
            ],
            dtype=float,
        ),
    )


def _direct_dscribe_soap(
    molecule,
    *,
    r_cut=6.0,
    n_max=8,
    l_max=6,
    sigma=1.0,
    species=None,
    centers=None,
):
    """Reference SOAP via direct DScribe call (0-based centers)."""
    symbols = list(molecule.chemical_symbols)
    if species is None:
        species = sorted(set(symbols))
    soap = SOAP(
        species=species,
        periodic=False,
        r_cut=r_cut,
        n_max=n_max,
        l_max=l_max,
        sigma=sigma,
        average="off",
        sparse=False,
        dtype="float64",
    )
    atoms = Atoms(symbols=symbols, positions=molecule.positions)
    features = soap.create(atoms, centers=centers)
    features = np.asarray(features, dtype=np.float64)
    if features.ndim == 1:
        features = features.reshape(1, -1)
    return features


def _rotation_matrix(axis, angle_rad):
    """Rodrigues rotation matrix for a unit axis."""
    axis = np.asarray(axis, dtype=float)
    axis = axis / np.linalg.norm(axis)
    x, y, z = axis
    c = np.cos(angle_rad)
    s = np.sin(angle_rad)
    C = 1.0 - c
    return np.array(
        [
            [c + x * x * C, x * y * C - z * s, x * z * C + y * s],
            [y * x * C + z * s, c + y * y * C, y * z * C - x * s],
            [z * x * C - y * s, z * y * C + x * s, c + z * z * C],
        ],
        dtype=float,
    )


class TestCalculateSoapBasics:
    """Basic contract and Molecule delegate tests."""

    def test_shape_dtype_finite(self, water_molecule):
        features = calculate_soap(water_molecule, n_max=4, l_max=2)
        assert features.dtype == np.float64
        assert features.ndim == 2
        assert features.shape[0] == water_molecule.num_atoms
        assert features.shape[1] > 0
        assert np.isfinite(features).all()
        assert np.linalg.norm(features) > 0

    def test_molecule_delegate_matches_module(self, water_molecule):
        via_method = water_molecule.calculate_soap(n_max=4, l_max=2)
        via_module = calculate_soap(water_molecule, n_max=4, l_max=2)
        np.testing.assert_allclose(
            via_method, via_module, rtol=RTOL, atol=ATOL
        )

    def test_parity_with_direct_dscribe(self, water_molecule):
        features = calculate_soap(water_molecule, n_max=4, l_max=2, sigma=0.5)
        expected = _direct_dscribe_soap(
            water_molecule, n_max=4, l_max=2, sigma=0.5
        )
        np.testing.assert_allclose(features, expected, rtol=RTOL, atol=ATOL)

    def test_deterministic(self, water_molecule):
        a = calculate_soap(water_molecule, n_max=4, l_max=2)
        b = calculate_soap(water_molecule, n_max=4, l_max=2)
        np.testing.assert_allclose(a, b, rtol=0.0, atol=0.0)

    def test_single_center_returns_2d(self, water_molecule):
        features = calculate_soap(
            water_molecule, centers=[1], n_max=4, l_max=2
        )
        assert features.shape == (1, features.shape[1])

    def test_centers_1based_order_preserved(self, methanol_molecule):
        features = calculate_soap(
            methanol_molecule, centers=[2, 1], n_max=4, l_max=2
        )
        expected = _direct_dscribe_soap(
            methanol_molecule, centers=[1, 0], n_max=4, l_max=2
        )
        np.testing.assert_allclose(features, expected, rtol=RTOL, atol=ATOL)

    def test_explicit_species_basis(self, water_molecule):
        species = ["H", "C", "O", "N"]
        features = calculate_soap(
            water_molecule, species=species, n_max=4, l_max=2
        )
        expected = _direct_dscribe_soap(
            water_molecule, species=species, n_max=4, l_max=2
        )
        np.testing.assert_allclose(features, expected, rtol=RTOL, atol=ATOL)
        # Shared basis yields a larger feature space than inferred H/O.
        inferred = calculate_soap(water_molecule, n_max=4, l_max=2)
        assert features.shape[1] > inferred.shape[1]


class TestAggregation:
    """Global SOAP aggregation tests."""

    def test_mean_aggregation(self, water_molecule):
        local = calculate_soap(water_molecule, n_max=4, l_max=2)
        global_vec = calculate_soap(
            water_molecule, n_max=4, l_max=2, aggregation="mean"
        )
        assert global_vec.ndim == 1
        np.testing.assert_allclose(
            global_vec, np.mean(local, axis=0), rtol=RTOL, atol=ATOL
        )

    def test_sum_aggregation(self, water_molecule):
        local = calculate_soap(water_molecule, n_max=4, l_max=2)
        global_vec = calculate_soap(
            water_molecule, n_max=4, l_max=2, aggregation="sum"
        )
        assert global_vec.ndim == 1
        np.testing.assert_allclose(
            global_vec, np.sum(local, axis=0), rtol=RTOL, atol=ATOL
        )

    def test_sum_scales_with_center_count(self, water_molecule):
        one = calculate_soap(
            water_molecule, centers=[1], n_max=4, l_max=2, aggregation="sum"
        )
        all_centers = calculate_soap(
            water_molecule, n_max=4, l_max=2, aggregation="sum"
        )
        # Sum over all centers should differ from a single-center sum.
        assert not np.allclose(one, all_centers, rtol=RTOL, atol=ATOL)


class TestInvariances:
    """Rigid-motion and permutation properties."""

    def test_translation_invariance(self, water_molecule):
        shifted = Molecule(
            symbols=water_molecule.chemical_symbols,
            positions=water_molecule.positions + np.array([10.0, -3.5, 2.0]),
        )
        a = calculate_soap(water_molecule, n_max=4, l_max=2)
        b = calculate_soap(shifted, n_max=4, l_max=2)
        np.testing.assert_allclose(a, b, rtol=RTOL, atol=ATOL)

    def test_rotation_invariance(self, water_molecule):
        R = _rotation_matrix([1.0, 2.0, 3.0], np.deg2rad(37.0))
        rotated = Molecule(
            symbols=water_molecule.chemical_symbols,
            positions=water_molecule.positions @ R.T,
        )
        a = calculate_soap(water_molecule, n_max=4, l_max=2)
        b = calculate_soap(rotated, n_max=4, l_max=2)
        np.testing.assert_allclose(a, b, rtol=RTOL, atol=ATOL)

    def test_permutation_local_equivariance(self, methanol_molecule):
        perm = [1, 0, 5, 2, 4, 3]  # 0-based permutation of atom order
        symbols = [methanol_molecule.chemical_symbols[i] for i in perm]
        positions = methanol_molecule.positions[perm]
        permuted = Molecule(symbols=symbols, positions=positions)

        original = calculate_soap(methanol_molecule, n_max=4, l_max=2)
        permuted_features = calculate_soap(permuted, n_max=4, l_max=2)
        # Local SOAP is row-equivariant under atom reordering.
        np.testing.assert_allclose(
            permuted_features, original[perm], rtol=RTOL, atol=ATOL
        )

    def test_permutation_global_invariance(self, methanol_molecule):
        perm = [1, 0, 5, 2, 4, 3]
        symbols = [methanol_molecule.chemical_symbols[i] for i in perm]
        positions = methanol_molecule.positions[perm]
        permuted = Molecule(symbols=symbols, positions=positions)

        a = calculate_soap(
            methanol_molecule, n_max=4, l_max=2, aggregation="mean"
        )
        b = calculate_soap(permuted, n_max=4, l_max=2, aggregation="mean")
        np.testing.assert_allclose(a, b, rtol=RTOL, atol=ATOL)


class TestSensitivity:
    """SOAP should respond to meaningful chemical / geometric changes."""

    def test_geometry_change_alters_descriptor(self, water_molecule):
        stretched = Molecule(
            symbols=water_molecule.chemical_symbols,
            positions=water_molecule.positions * 1.25,
        )
        a = calculate_soap(
            water_molecule, n_max=4, l_max=2, aggregation="mean"
        )
        b = calculate_soap(stretched, n_max=4, l_max=2, aggregation="mean")
        assert not np.allclose(a, b, rtol=RTOL, atol=ATOL)

    def test_species_change_alters_descriptor(self, water_molecule):
        # Replace oxygen with sulfur; use shared species basis.
        sulfur = Molecule(
            symbols=["S", "H", "H"],
            positions=water_molecule.positions.copy(),
        )
        species = ["H", "O", "S"]
        a = calculate_soap(
            water_molecule,
            species=species,
            n_max=4,
            l_max=2,
            aggregation="mean",
        )
        b = calculate_soap(
            sulfur, species=species, n_max=4, l_max=2, aggregation="mean"
        )
        assert not np.allclose(a, b, rtol=RTOL, atol=ATOL)


class TestValidation:
    """Input validation and PBC rejection."""

    def test_rejects_active_pbc(self, water_molecule):
        water_molecule.pbc_conditions = [1, 0, 0]
        with pytest.raises(ValueError, match="Periodic SOAP is not supported"):
            calculate_soap(water_molecule, n_max=4, l_max=2)

    def test_rejects_zero_pbc_flags_as_inactive(self, water_molecule):
        # Explicit all-false PBC should still be treated as finite.
        water_molecule.pbc_conditions = [0, 0, 0]
        features = calculate_soap(water_molecule, n_max=4, l_max=2)
        assert features.shape[0] == 3

    def test_invalid_aggregation(self, water_molecule):
        with pytest.raises(ValueError, match="aggregation"):
            calculate_soap(water_molecule, aggregation="inner")

    def test_invalid_r_cut(self, water_molecule):
        with pytest.raises(ValueError, match="r_cut"):
            calculate_soap(water_molecule, r_cut=0.0)
        with pytest.raises(ValueError, match="r_cut"):
            calculate_soap(water_molecule, r_cut=1.0)
        with pytest.raises(ValueError, match="r_cut"):
            calculate_soap(water_molecule, r_cut=0.5)

    def test_invalid_n_max(self, water_molecule):
        with pytest.raises(ValueError, match="n_max"):
            calculate_soap(water_molecule, n_max=0)
        with pytest.raises(TypeError, match="n_max"):
            calculate_soap(water_molecule, n_max=True)

    def test_invalid_l_max(self, water_molecule):
        with pytest.raises(ValueError, match="l_max"):
            calculate_soap(water_molecule, l_max=-1)
        with pytest.raises(TypeError, match="l_max"):
            calculate_soap(water_molecule, l_max=1.5)

    def test_invalid_sigma(self, water_molecule):
        with pytest.raises(ValueError, match="sigma"):
            calculate_soap(water_molecule, sigma=-0.1)

    def test_centers_out_of_range(self, water_molecule):
        with pytest.raises(ValueError, match="out of range"):
            calculate_soap(water_molecule, centers=[4])

    def test_centers_empty(self, water_molecule):
        with pytest.raises(ValueError, match="non-empty"):
            calculate_soap(water_molecule, centers=[])

    def test_centers_non_integer(self, water_molecule):
        with pytest.raises(TypeError, match="1-based integer"):
            calculate_soap(water_molecule, centers=[1.0])

    def test_species_missing_elements(self, water_molecule):
        with pytest.raises(ValueError, match="not in the SOAP basis"):
            calculate_soap(water_molecule, species=["H", "C"])

    def test_species_empty(self, water_molecule):
        with pytest.raises(ValueError, match="non-empty"):
            calculate_soap(water_molecule, species=[])

    def test_non_molecule_type(self):
        with pytest.raises(TypeError, match="Molecule"):
            calculate_soap("not-a-molecule")

    def test_nonfinite_positions(self, water_molecule):
        bad = Molecule(
            symbols=water_molecule.chemical_symbols,
            positions=np.array(
                [
                    [0.0, 0.0, np.nan],
                    [0.0, 0.7, -0.5],
                    [0.0, -0.7, -0.5],
                ]
            ),
        )
        with pytest.raises(ValueError, match="finite"):
            calculate_soap(bad, n_max=4, l_max=2)

    def test_unrecognized_element_symbol(self):
        mol = Molecule(
            symbols=["O", "H", "D"],
            positions=np.array(
                [
                    [0.0, 0.0, 0.1],
                    [0.0, 0.7, -0.5],
                    [0.0, -0.7, -0.5],
                ]
            ),
        )
        with pytest.raises(ValueError, match="Unrecognized elemental symbol"):
            calculate_soap(mol, n_max=4, l_max=2)

    def test_single_atom_molecule(self):
        mol = Molecule(
            symbols=["C"],
            positions=np.array([[0.0, 0.0, 0.0]]),
        )
        features = calculate_soap(mol, n_max=4, l_max=2)
        assert features.shape[0] == 1
        assert np.isfinite(features).all()

    def test_uses_live_symbols_not_cached(self, water_molecule):
        # Warm any caches that might exist, then mutate symbols in place.
        _ = water_molecule.chemical_symbols
        water_molecule.symbols = ["S", "H", "H"]
        species = ["H", "O", "S"]
        features = calculate_soap(
            water_molecule, species=species, n_max=4, l_max=2
        )
        expected = _direct_dscribe_soap(
            Molecule(
                symbols=["S", "H", "H"],
                positions=water_molecule.positions,
            ),
            species=species,
            n_max=4,
            l_max=2,
        )
        np.testing.assert_allclose(features, expected, rtol=RTOL, atol=ATOL)


class TestDocumentedLimitations:
    """Lock geometry-only SOAP limitations described in the docs."""

    def test_charge_and_multiplicity_invariance(self, water_molecule):
        charged = Molecule(
            symbols=water_molecule.chemical_symbols,
            positions=water_molecule.positions.copy(),
            charge=-1,
            multiplicity=2,
        )
        a = calculate_soap(water_molecule, n_max=4, l_max=2)
        b = calculate_soap(charged, n_max=4, l_max=2)
        np.testing.assert_allclose(a, b, rtol=RTOL, atol=ATOL)

    def test_reflection_invariance(self, methanol_molecule):
        mirrored = Molecule(
            symbols=methanol_molecule.chemical_symbols,
            positions=methanol_molecule.positions * np.array([1.0, 1.0, -1.0]),
        )
        a = calculate_soap(
            methanol_molecule, n_max=4, l_max=2, aggregation="mean"
        )
        b = calculate_soap(mirrored, n_max=4, l_max=2, aggregation="mean")
        np.testing.assert_allclose(a, b, rtol=RTOL, atol=ATOL)


class TestOptionalDependency:
    """DScribe must remain optional at import / call time."""

    def test_import_molecule_without_dscribe(self):
        # Molecule import path must not require DScribe.
        from chemsmart.io.molecules import structure as structure_mod

        assert structure_mod.Molecule is Molecule

    def test_missing_dscribe_raises_actionable_error(self, water_molecule):
        real_import = __import__

        def fake_import(name, *args, **kwargs):
            if name == "dscribe" or name.startswith("dscribe."):
                raise ImportError("No module named 'dscribe'")
            return real_import(name, *args, **kwargs)

        with patch("builtins.__import__", side_effect=fake_import):
            with pytest.raises(ImportError, match=r"chemsmart\[soap\]") as exc:
                calculate_soap(water_molecule, n_max=4, l_max=2)
            assert isinstance(exc.value.__cause__, ImportError)
