"""Tests for SOAP descriptor calculations on Molecule objects."""

from pathlib import Path

import numpy as np
import pytest

from chemsmart.analysis.soap import calculate_soap
from chemsmart.io.molecules.structure import Molecule

# Cross-platform numerical tolerances for dense float64 SOAP vectors.
RTOL = 1e-6
ATOL = 1e-8

_REFERENCE_PATH = (
    Path(__file__).resolve().parent / "data" / "soap_reference.npz"
)


@pytest.fixture(scope="module")
def soap_reference():
    """Golden SOAP vectors generated once from DScribe 2.1.2 GTO SOAP."""
    return np.load(_REFERENCE_PATH)


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

    def test_parity_with_golden_reference(
        self, water_molecule, soap_reference
    ):
        features = calculate_soap(water_molecule, n_max=4, l_max=2, sigma=0.5)
        expected = soap_reference["water_sigma_0_5"]
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

    def test_centers_1based_order_preserved(
        self, methanol_molecule, soap_reference
    ):
        features = calculate_soap(
            methanol_molecule, centers=[2, 1], n_max=4, l_max=2
        )
        expected = soap_reference["methanol_centers_2_1"]
        np.testing.assert_allclose(features, expected, rtol=RTOL, atol=ATOL)

    def test_duplicate_centers_preserved_and_bias_mean(self, water_molecule):
        # Duplicates are kept; mean aggregation overweights that center.
        local = calculate_soap(
            water_molecule, centers=[1, 1], n_max=4, l_max=2
        )
        assert local.shape[0] == 2
        np.testing.assert_allclose(local[0], local[1], rtol=RTOL, atol=ATOL)

        mean_dup = calculate_soap(
            water_molecule,
            centers=[1, 1],
            n_max=4,
            l_max=2,
            aggregation="mean",
        )
        mean_all = calculate_soap(
            water_molecule, n_max=4, l_max=2, aggregation="mean"
        )
        single = calculate_soap(
            water_molecule, centers=[1], n_max=4, l_max=2, aggregation="mean"
        )
        np.testing.assert_allclose(mean_dup, single, rtol=RTOL, atol=ATOL)
        assert not np.allclose(mean_dup, mean_all, rtol=RTOL, atol=ATOL)

    def test_explicit_species_basis(self, water_molecule, soap_reference):
        species = ["H", "C", "O", "N"]
        features = calculate_soap(
            water_molecule, species=species, n_max=4, l_max=2
        )
        expected = soap_reference["water_species_HCON"]
        np.testing.assert_allclose(features, expected, rtol=RTOL, atol=ATOL)
        # Shared basis yields a larger feature space than inferred H/O.
        inferred = calculate_soap(water_molecule, n_max=4, l_max=2)
        assert features.shape[1] > inferred.shape[1]

    def test_parity_public_defaults(self, water_molecule, soap_reference):
        features = calculate_soap(water_molecule)
        expected = soap_reference["water_defaults"]
        np.testing.assert_allclose(features, expected, rtol=RTOL, atol=ATOL)

    def test_parity_lmax4_covers_lpmv_path(
        self, water_molecule, soap_reference
    ):
        features = calculate_soap(water_molecule, n_max=4, l_max=4, sigma=1.0)
        expected = soap_reference["water_lmax4_sigma1"]
        np.testing.assert_allclose(features, expected, rtol=RTOL, atol=ATOL)

    def test_parity_methanol_sigma1(self, methanol_molecule, soap_reference):
        features = calculate_soap(
            methanol_molecule, n_max=4, l_max=2, sigma=1.0
        )
        expected = soap_reference["methanol_sigma1"]
        np.testing.assert_allclose(features, expected, rtol=RTOL, atol=ATOL)

    def test_parity_padding_shell_neighbor(self, soap_reference):
        # He at 8 Å from O sits in DScribe's r_cut+~3.72σ padding shell.
        mol = Molecule(
            symbols=["O", "H", "H", "He"],
            positions=np.array(
                [
                    [0.00000, 0.00000, 0.11779],
                    [0.00000, 0.75545, -0.47116],
                    [0.00000, -0.75545, -0.47116],
                    [8.0, 0.0, 0.11779],
                ],
                dtype=float,
            ),
        )
        features = calculate_soap(
            mol,
            species=["H", "He", "O"],
            n_max=4,
            l_max=2,
            sigma=1.0,
            r_cut=6.0,
        )
        expected = soap_reference["padding_shell_He"]
        np.testing.assert_allclose(features, expected, rtol=RTOL, atol=ATOL)

    def test_padding_shell_changes_descriptor(self, water_molecule):
        # Regression lock: atoms beyond r_cut but inside the padded neighbor
        # shell must contribute (DScribe 2.1.2 GTO convention).
        species = ["H", "He", "O"]
        without = calculate_soap(
            water_molecule,
            species=species,
            n_max=4,
            l_max=2,
            sigma=1.0,
            aggregation="mean",
        )
        with_he = Molecule(
            symbols=["O", "H", "H", "He"],
            positions=np.vstack(
                [
                    water_molecule.positions,
                    [[8.0, 0.0, 0.11779]],
                ]
            ),
        )
        with_ = calculate_soap(
            with_he,
            species=species,
            centers=[1, 2, 3],  # water centers only
            n_max=4,
            l_max=2,
            sigma=1.0,
            aggregation="mean",
        )
        assert not np.allclose(without, with_, rtol=RTOL, atol=ATOL)


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

    def test_rejects_translation_vectors_without_pbc(self, water_molecule):
        # Coordinate-block style: TVs set, pbc_conditions left unset.
        water_molecule.pbc_conditions = None
        water_molecule.translation_vectors = [
            [5.0, 0.0, 0.0],
            [0.0, 5.0, 0.0],
            [0.0, 0.0, 5.0],
        ]
        with pytest.raises(ValueError, match="Periodic SOAP is not supported"):
            calculate_soap(water_molecule, n_max=4, l_max=2)

    def test_empty_translation_vectors_allowed(self, water_molecule):
        water_molecule.translation_vectors = []
        features = calculate_soap(water_molecule, n_max=4, l_max=2)
        assert features.shape[0] == 3

    def test_rejects_translation_vectors_without_len(self, water_molecule):
        # Non-sized TV containers still count as periodic/cell-bearing.
        class _UnsizedTV:
            def __len__(self):
                raise TypeError("unsized")

        water_molecule.pbc_conditions = None
        water_molecule.translation_vectors = _UnsizedTV()
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
        with pytest.raises(ValueError, match="l_max"):
            calculate_soap(water_molecule, l_max=21)
        with pytest.raises(TypeError, match="l_max"):
            calculate_soap(water_molecule, l_max=1.5)

    def test_lmax0_and_lmax1_smoke(self, water_molecule):
        # Exercise solid-harmonics early-return path (l_max < 2).
        for l_max in (0, 1):
            features = calculate_soap(water_molecule, n_max=4, l_max=l_max)
            assert features.shape[0] == water_molecule.num_atoms
            assert np.isfinite(features).all()

    def test_solid_harmonics_survives_cos_theta_roundoff(self):
        # Float noise can push |z/r| slightly outside [-1, 1]; lpmv needs clip.
        from chemsmart.analysis.soap import _solid_harmonics_polynomial

        x = np.array([1e-16, 0.0])
        y = np.array([0.0, 0.0])
        z = np.array([1.0, 0.0])
        harm = _solid_harmonics_polynomial(x, y, z, l_max=4)
        assert harm.shape == (2, 25)
        assert np.isfinite(harm).all()

    def test_solid_harmonics_early_return_empty_or_low_l(self):
        from chemsmart.analysis.soap import _solid_harmonics_polynomial

        empty = _solid_harmonics_polynomial(
            np.array([]), np.array([]), np.array([]), l_max=4
        )
        assert empty.shape == (0, 25)
        low_l = _solid_harmonics_polynomial(
            np.array([1.0]), np.array([0.0]), np.array([0.0]), l_max=1
        )
        assert low_l.shape == (1, 4)
        assert np.all(low_l == 0.0)

    def test_validate_positions_shape_and_count(self, water_molecule):
        from chemsmart.analysis.soap import _validate_positions

        with pytest.raises(ValueError, match=r"shape \(N, 3\)"):
            _validate_positions(np.array([0.0, 0.0, 0.0]), 1)
        with pytest.raises(ValueError, match="does not match"):
            _validate_positions(water_molecule.positions, num_atoms=1)

    def test_rejects_none_or_empty_symbols(self, water_molecule):
        water_molecule.symbols = None
        with pytest.raises(ValueError, match="no atoms"):
            calculate_soap(water_molecule, n_max=4, l_max=2)

        water_molecule.symbols = []
        with pytest.raises(ValueError, match="no atoms"):
            calculate_soap(water_molecule, n_max=4, l_max=2)

    def test_gto_basis_rejects_complex_normalization(self, monkeypatch):
        import chemsmart.analysis.soap as soap_mod

        def _complex_sqrtm(_matrix):
            return np.array([[1.0 + 1.0j]], dtype=np.complex128)

        monkeypatch.setattr(soap_mod, "sqrtm", _complex_sqrtm)
        with pytest.raises(ValueError, match="normalization factors"):
            soap_mod._gto_basis(r_cut=6.0, n_max=1, l_max=0)

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

    def test_species_unrecognized_symbol(self, water_molecule):
        with pytest.raises(ValueError, match="Unrecognized elemental symbol"):
            calculate_soap(water_molecule, species=["H", "O", "Xx"])
        with pytest.raises(ValueError, match="Unrecognized elemental symbol"):
            calculate_soap(water_molecule, species=["H", "O", "Zz"])
        with pytest.raises(ValueError, match="Unrecognized elemental symbol"):
            calculate_soap(water_molecule, species=["H", "O", ""])
        with pytest.raises(ValueError, match="Unrecognized elemental symbol"):
            calculate_soap(water_molecule, species=["H", "O", "---"])

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

    def test_uses_live_symbols_not_cached(
        self, water_molecule, soap_reference
    ):
        # Warm any caches that might exist, then mutate symbols in place.
        _ = water_molecule.chemical_symbols
        water_molecule.symbols = ["S", "H", "H"]
        species = ["H", "O", "S"]
        features = calculate_soap(
            water_molecule, species=species, n_max=4, l_max=2
        )
        expected = soap_reference["sulfur_species_HOS"]
        np.testing.assert_allclose(features, expected, rtol=RTOL, atol=ATOL)

    def test_mixed_case_molecule_symbols_match_canonical(self, water_molecule):
        mixed = Molecule(
            symbols=["o", "h", "h"],
            positions=water_molecule.positions.copy(),
        )
        a = calculate_soap(water_molecule, n_max=4, l_max=2)
        b = calculate_soap(mixed, n_max=4, l_max=2)
        np.testing.assert_allclose(a, b, rtol=RTOL, atol=ATOL)

    def test_mixed_case_species_list_match_canonical(self, water_molecule):
        a = calculate_soap(
            water_molecule, species=["H", "C", "O", "N"], n_max=4, l_max=2
        )
        b = calculate_soap(
            water_molecule, species=["h", "c", "o", "n"], n_max=4, l_max=2
        )
        np.testing.assert_allclose(a, b, rtol=RTOL, atol=ATOL)

    def test_whitespace_padded_symbols_match_canonical(self, water_molecule):
        padded = Molecule(
            symbols=[" O", "H ", "H"],
            positions=water_molecule.positions.copy(),
        )
        a = calculate_soap(water_molecule, n_max=4, l_max=2)
        b = calculate_soap(padded, n_max=4, l_max=2)
        np.testing.assert_allclose(a, b, rtol=RTOL, atol=ATOL)

    def test_lowercase_deuterium_rejected(self):
        mol = Molecule(
            symbols=["O", "H", "d"],
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


class TestBuiltInImplementation:
    """SOAP is implemented in-tree without optional dependencies."""

    def test_import_molecule_without_optional_deps(self):
        from chemsmart.io.molecules import structure as structure_mod

        assert structure_mod.Molecule is Molecule

    def test_dscribe_not_imported_by_soap_module(self):
        import sys

        import chemsmart.analysis.soap as soap_mod

        assert "dscribe" not in sys.modules
        assert soap_mod.__name__ == "chemsmart.analysis.soap"
