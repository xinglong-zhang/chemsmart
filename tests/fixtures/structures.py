"""Structure-file fixtures (xyz, ChemDraw, InChIKey, CXSMILES, canonical geometry, SOAP)."""

import os

import numpy as np
import pytest

from chemsmart.io.molecules.structure import Molecule
from tests.helpers import rdkit_embed_and_optimize, write_xyz_file


# test for structure.py
@pytest.fixture()
def structure_test_directory(test_data_directory):
    return os.path.join(test_data_directory, "StructuresTests")


@pytest.fixture()
def xyz_directory(structure_test_directory):
    return os.path.join(structure_test_directory, "xyz")


@pytest.fixture()
def ts_conformers_log_directory(structure_test_directory):
    """Directory containing TS conformer log files (ch_1c_para_c1.log to c5.log)."""
    return os.path.join(structure_test_directory, "conformers", "log")


@pytest.fixture()
def single_molecule_xyz_file(xyz_directory):
    return os.path.join(xyz_directory, "crest_best.xyz")


@pytest.fixture()
def multiple_molecules_xyz_file(xyz_directory):
    return os.path.join(xyz_directory, "crest_conformers.xyz")


@pytest.fixture()
def two_rotated_molecules_xyz_file(xyz_directory):
    return os.path.join(xyz_directory, "two_rotated_molecules.xyz")


@pytest.fixture()
def xtb_optimized_xyz_file(xyz_directory):
    return os.path.join(xyz_directory, "ts_xtbopt.xyz")


@pytest.fixture()
def chemsmart_generated_xyz_file(xyz_directory):
    return os.path.join(xyz_directory, "frozen_coordinates_opt.xyz")


@pytest.fixture()
def extended_xyz_file(xyz_directory):
    return os.path.join(xyz_directory, "crystal.extxyz")


@pytest.fixture()
def dna_hybrid_visualized_xyz_file(xyz_directory):
    return os.path.join(xyz_directory, "dna_hybrid.xyz")


@pytest.fixture()
def visualized_1_mer_xyz_file(xyz_directory):
    return os.path.join(xyz_directory, "1-mer.xyz")


@pytest.fixture()
def chemdraw_directory(structure_test_directory):
    return os.path.join(structure_test_directory, "chemdraw")


@pytest.fixture()
def single_molecule_cdxml_file_benzene(chemdraw_directory):
    return os.path.join(chemdraw_directory, "benzene.cdxml")


@pytest.fixture()
def single_molecule_cdxml_file_methane(chemdraw_directory):
    return os.path.join(chemdraw_directory, "methane.cdxml")


@pytest.fixture()
def multi_molecule_cdxml_file(chemdraw_directory):
    return os.path.join(chemdraw_directory, "two_molecules.cdxml")


@pytest.fixture()
def single_molecule_cdx_file_imidazole(chemdraw_directory):
    return os.path.join(chemdraw_directory, "imidazole.cdx")


@pytest.fixture()
def complex_molecule_cdxml_file(chemdraw_directory):
    return os.path.join(chemdraw_directory, "complex_molecule.cdxml")


@pytest.fixture()
def chemdraw_expected_directory(chemdraw_directory):
    return os.path.join(chemdraw_directory, "expected")


@pytest.fixture()
def expected_methane_xyz(chemdraw_expected_directory):
    return os.path.join(chemdraw_expected_directory, "methane.xyz")


@pytest.fixture()
def expected_benzene_xyz(chemdraw_expected_directory):
    return os.path.join(chemdraw_expected_directory, "benzene.xyz")


@pytest.fixture()
def expected_imidazole_xyz(chemdraw_expected_directory):
    return os.path.join(chemdraw_expected_directory, "imidazole.xyz")


@pytest.fixture()
def expected_two_molecules_1_xyz(chemdraw_expected_directory):
    return os.path.join(chemdraw_expected_directory, "two_molecules_1.xyz")


@pytest.fixture()
def expected_two_molecules_2_xyz(chemdraw_expected_directory):
    return os.path.join(chemdraw_expected_directory, "two_molecules_2.xyz")


@pytest.fixture()
def expected_benzene_com(chemdraw_expected_directory):
    return os.path.join(chemdraw_expected_directory, "benzene.com")


@pytest.fixture()
def expected_methane_com(chemdraw_expected_directory):
    return os.path.join(chemdraw_expected_directory, "methane.com")


@pytest.fixture()
def expected_two_molecules_1_com(chemdraw_expected_directory):
    return os.path.join(chemdraw_expected_directory, "two_molecules_1.com")


@pytest.fixture()
def expected_two_molecules_2_com(chemdraw_expected_directory):
    return os.path.join(chemdraw_expected_directory, "two_molecules_2.com")


@pytest.fixture()
def expected_complex_molecule_com(chemdraw_expected_directory):
    return os.path.join(chemdraw_expected_directory, "complex_molecule.com")


@pytest.fixture()
def metal_ligand_molecules_cdxml_file(chemdraw_directory):
    return os.path.join(chemdraw_directory, "metal_ligands.cdxml")


@pytest.fixture()
def colored_proton_cdxml_file(chemdraw_directory):
    return os.path.join(chemdraw_directory, "phenol.cdxml")


@pytest.fixture()
def colored_implicit_proton_cdxml_file(chemdraw_directory):
    """Returns the path to a CDXML file with a colored implicit proton, which should be treated as a colored explicit proton."""
    return os.path.join(chemdraw_directory, "phenol_implicit_proton.cdxml")


@pytest.fixture()
def colored_proton_two_molecule_cdxml_file(chemdraw_directory):
    return os.path.join(chemdraw_directory, "phenol_two_molecule.cdxml")


@pytest.fixture()
def pka_scale_cdxml_file(chemdraw_directory):
    """Multi-fragment CDXML with nested groups and coloured acidic protons."""
    return os.path.join(chemdraw_directory, "pka_scale.cdxml")


# ── InChIKey test data ──
@pytest.fixture()
def inchikey_test_directory(structure_test_directory):
    return os.path.join(structure_test_directory, "inchikey")


@pytest.fixture()
def inchikey_normal_file(inchikey_test_directory):
    return os.path.join(
        inchikey_test_directory,
        "normal_testing",
        "inchikey_normal_testing.xyz",
    )


@pytest.fixture()
def inchikey_r_enantiomer_file(inchikey_test_directory):
    return os.path.join(
        inchikey_test_directory,
        "enantiomer_testing",
        "inchikey_r_enantiomer.xyz",
    )


@pytest.fixture()
def inchikey_s_enantiomer_file(inchikey_test_directory):
    return os.path.join(
        inchikey_test_directory,
        "enantiomer_testing",
        "inchikey_s_enantiomer.xyz",
    )


@pytest.fixture()
def inchikey_large_molecule_c3_file(inchikey_test_directory):
    return os.path.join(
        inchikey_test_directory,
        "large_molecule_testing",
        "inchikey_large_molecule_c3.xyz",
    )


@pytest.fixture()
def inchikey_large_molecule_c2_file(inchikey_test_directory):
    return os.path.join(
        inchikey_test_directory,
        "large_molecule_testing",
        "inchikey_large_molecule_c2.xyz",
    )


# ── CXSMILES test data ──
@pytest.fixture()
def cxsmiles_test_directory(structure_test_directory):
    return os.path.join(structure_test_directory, "cxsmiles")


@pytest.fixture()
def cxsmiles_normal_file(cxsmiles_test_directory):
    return os.path.join(
        cxsmiles_test_directory,
        "normal_testing",
        "rotamer_normal_testing.xyz",
    )


@pytest.fixture()
def cxsmiles_r_enantiomer_file(cxsmiles_test_directory):
    return os.path.join(
        cxsmiles_test_directory,
        "enantiomer_testing",
        "rotamer_r_enantiomer.xyz",
    )


@pytest.fixture()
def cxsmiles_s_enantiomer_file(cxsmiles_test_directory):
    return os.path.join(
        cxsmiles_test_directory,
        "enantiomer_testing",
        "rotamer_s_enantiomer.xyz",
    )


@pytest.fixture()
def cxsmiles_r_rotamer_file(cxsmiles_test_directory):
    return os.path.join(
        cxsmiles_test_directory,
        "rotamer_testing",
        "cxsmiles_r_rotamer.xyz",
    )


@pytest.fixture()
def cxsmiles_s_rotamer_file(cxsmiles_test_directory):
    return os.path.join(
        cxsmiles_test_directory,
        "rotamer_testing",
        "cxsmiles_s_rotamer.xyz",
    )


@pytest.fixture()
def cxsmiles_large_molecule_c2_file(cxsmiles_test_directory):
    return os.path.join(
        cxsmiles_test_directory,
        "large_molecule_testing",
        "rotamer_large_molecule_c2.xyz",
    )


@pytest.fixture()
def cxsmiles_large_molecule_c3_file(cxsmiles_test_directory):
    return os.path.join(
        cxsmiles_test_directory,
        "large_molecule_testing",
        "rotamer_large_molecule_c3.xyz",
    )


@pytest.fixture()
def cxsmiles_expected_large_c2_file(cxsmiles_test_directory):
    return os.path.join(
        cxsmiles_test_directory,
        "large_molecule_testing",
        "expected_cxsmiles_c2.txt",
    )


@pytest.fixture()
def cxsmiles_expected_large_c3_file(cxsmiles_test_directory):
    return os.path.join(
        cxsmiles_test_directory,
        "large_molecule_testing",
        "expected_cxsmiles_c3.txt",
    )


############ Canonical Geometry / Structure ID Fixtures ##################
@pytest.fixture()
def canonical_test_directory(structure_test_directory):
    return os.path.join(structure_test_directory, "canonical")


@pytest.fixture()
def canonical_formaldehyde_file(canonical_test_directory):
    """Formaldehyde (CH2O) in its reference orientation — planar, C2v symmetry."""
    return os.path.join(canonical_test_directory, "formaldehyde.xyz")


@pytest.fixture()
def canonical_formaldehyde_trans_rot_file(canonical_test_directory):
    """Formaldehyde translated and rotated from canonical_formaldehyde_file."""
    return os.path.join(canonical_test_directory, "formaldehyde_trans_rot.xyz")


@pytest.fixture()
def canonical_formaldehyde_perturbed_file(canonical_test_directory):
    """Formaldehyde with coordinates perturbed at ~1e-7 Å level from canonical_formaldehyde_file."""
    return os.path.join(canonical_test_directory, "formaldehyde_perturbed.xyz")


@pytest.fixture()
def canonical_methane_file(canonical_test_directory):
    """Methane (CH4) in its reference orientation — tetrahedral, Td symmetry."""
    return os.path.join(canonical_test_directory, "methane.xyz")


@pytest.fixture()
def canonical_methane_trans_rot_file(canonical_test_directory):
    """Methane translated and rotated from canonical_methane_file."""
    return os.path.join(canonical_test_directory, "methane_trans_rot.xyz")


@pytest.fixture()
def canonical_methane_distorted_file(canonical_test_directory):
    """Methane with one C-H bond elongated by ~2e-3 Å from canonical_methane_file."""
    return os.path.join(canonical_test_directory, "methane_distorted.xyz")


@pytest.fixture()
def canonical_3b_file(canonical_test_directory):
    """3b (C17H17NOS) in its reference orientation — large, low-symmetry (C1) molecule."""
    return os.path.join(canonical_test_directory, "3b.xyz")


@pytest.fixture()
def canonical_3b_trans_rot_file(canonical_test_directory):
    """3b translated and rotated from canonical_3b_file."""
    return os.path.join(canonical_test_directory, "3b_trans_rot.xyz")


@pytest.fixture()
def canonical_3b_permuted_file(canonical_test_directory):
    """3b with input atom order permuted — same coordinates as canonical_3b_file."""
    return os.path.join(canonical_test_directory, "3b_permuted.xyz")


@pytest.fixture()
def canonical_r_bromochlorofluoromethane_file(canonical_test_directory):
    """(R)-Bromochlorofluoromethane — one enantiomer of a chiral CHBrClF molecule."""
    return os.path.join(
        canonical_test_directory, "R-Bromochlorofluoromethane.xyz"
    )


@pytest.fixture()
def canonical_s_bromochlorofluoromethane_file(canonical_test_directory):
    """(S)-Bromochlorofluoromethane — non-superimposable mirror image of the R enantiomer."""
    return os.path.join(
        canonical_test_directory, "S-Bromochlorofluoromethane.xyz"
    )


############ SOAP Fixtures ##################
@pytest.fixture()
def soap_test_directory(test_data_directory):
    return os.path.join(test_data_directory, "SOAPTests")


@pytest.fixture()
def soap_reference_filepath(soap_test_directory):
    return os.path.join(soap_test_directory, "soap_reference.npz")


@pytest.fixture()
def soap_reference_data(soap_reference_filepath):
    """Fixture to load the SOAP reference data."""
    return np.load(soap_reference_filepath)


_SOAP_METHANE = (
    ["C", "H", "H", "H", "H"],
    [
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 1.09],
        [0.0, 1.03, -0.36],
        [0.89, -0.51, -0.36],
        [-0.89, -0.51, -0.36],
    ],
)
_SOAP_HALO_CHAIN = (
    ["C"] * 5 + ["F", "Cl"] + ["H"] * 10,
    [
        [-1.0144, 1.0305, 0.9144],
        [-0.5694, 0.4225, -0.4231],
        [-1.2529, -0.9297, -0.664],
        [0.9715, 0.3516, -0.5595],
        [1.6781, -0.5734, 0.4436],
        [1.4033, -1.899, 0.1621],
        [3.441, -0.3363, 0.3383],
        [-0.4859, 1.99, 1.098],
        [-2.1056, 1.2374, 0.8889],
        [-0.8143, 0.3382, 1.7586],
        [-0.9238, 1.1101, -1.2233],
        [-2.3501, -0.7845, -0.7606],
        [-0.8811, -1.3889, -1.6046],
        [-1.071, -1.6266, 0.1808],
        [1.3772, 1.3807, -0.4425],
        [1.2285, 0.0227, -1.5904],
        [1.369, -0.3454, 1.4833],
    ],
)
# Reference key -> ((symbols, positions), calculate_soap keyword arguments).
_SOAP_DSCRIBE_CASES = {
    "methane_rcut5_lmax3": (
        _SOAP_METHANE,
        dict(r_cut=5.0, n_max=4, l_max=3, sigma=1.0),
    ),
    "hcn_linear_sigma_0_5": (
        (
            ["H", "C", "N"],
            [[0.0, 0.0, -1.06], [0.0, 0.0, 0.0], [0.0, 0.0, 1.15]],
        ),
        dict(r_cut=5.0, n_max=4, l_max=3, sigma=0.5),
    ),
    "halo_chain_4species": (
        _SOAP_HALO_CHAIN,
        dict(r_cut=5.0, n_max=4, l_max=3, sigma=1.0),
    ),
    "water_lmax9": (
        (
            ["O", "H", "H"],
            [[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]],
        ),
        dict(r_cut=4.0, n_max=3, l_max=9, sigma=0.7),
    ),
    "methane_lmax0_rcut3": (
        _SOAP_METHANE,
        dict(r_cut=3.0, n_max=2, l_max=0, sigma=0.3),
    ),
    "argon_single_atom": (
        (["Ar"], [[0.0, 0.0, 0.0]]),
        dict(r_cut=5.0, n_max=4, l_max=3, sigma=1.0),
    ),
}


@pytest.fixture(params=sorted(_SOAP_DSCRIBE_CASES))
def soap_dscribe_case(request):
    """One DScribe 2.1.2 reference case: key, symbols, positions, kwargs."""
    (symbols, positions), kwargs = _SOAP_DSCRIBE_CASES[request.param]
    return request.param, symbols, np.array(positions, dtype=float), kwargs


@pytest.fixture()
def clear_pubchem_cache():
    """Clear ``pubchem_search``'s ``lru_cache`` before and after the test so
    cached lookups can't leak between tests. Opt in with
    ``pytestmark = pytest.mark.usefixtures("clear_pubchem_cache")``."""
    from chemsmart.io.molecules.pubchem import pubchem_search

    pubchem_search.cache_clear()
    yield
    pubchem_search.cache_clear()


@pytest.fixture()
def sterimol_methyl_geometry():
    """Symbols and positions of CH4 with H(atom 2)->C(atom 1) along -z."""
    symbols = ["C", "H", "H", "H", "H"]
    positions = np.array(
        [
            [1.0, 2.0, 3.0],
            [1.0, 2.0, 4.09],
            [1.0, 0.97, 2.636],
            [1.892, 1.485, 2.636],
            [0.108, 1.485, 2.636],
        ]
    )
    return symbols, positions


@pytest.fixture()
def sterimol_halo_chain_geometry():
    """Symbols and positions of CH3-CH(CH3)-CH2-CHFCl, with no symmetry.

    Atoms (1-based): C1 and C3 methyl, C2 methine, C4 methylene, C5 CHFCl,
    F6, Cl7. Hydrogens: H8-H10 on C1, H11 on C2, H12-H14 on C3, H15-H16 on C4,
    H17 on C5. The C5-F6 bond is 1.383 A.
    """
    symbols = ["C"] * 5 + ["F", "Cl"] + ["H"] * 10
    positions = np.array(
        [
            [-1.0144, 1.0305, 0.9144],
            [-0.5694, 0.4225, -0.4231],
            [-1.2529, -0.9297, -0.664],
            [0.9715, 0.3516, -0.5595],
            [1.6781, -0.5734, 0.4436],
            [1.4033, -1.899, 0.1621],
            [3.441, -0.3363, 0.3383],
            [-0.4859, 1.99, 1.098],
            [-2.1056, 1.2374, 0.8889],
            [-0.8143, 0.3382, 1.7586],
            [-0.9238, 1.1101, -1.2233],
            [-2.3501, -0.7845, -0.7606],
            [-0.8811, -1.3889, -1.6046],
            [-1.071, -1.6266, 0.1808],
            [1.3772, 1.3807, -0.4425],
            [1.2285, 0.0227, -1.5904],
            [1.369, -0.3454, 1.4833],
        ]
    )
    return symbols, positions


@pytest.fixture()
def sterimol_halo_chain_xyz_file(tmp_path, sterimol_halo_chain_geometry):
    """The halo chain as a single-structure xyz file."""
    symbols, positions = sterimol_halo_chain_geometry
    return write_xyz_file(tmp_path / "chain.xyz", symbols, [positions])


@pytest.fixture()
def sterimol_trajectory_xyz_file(tmp_path, sterimol_halo_chain_geometry):
    """Two frames: the halo chain, then the same with F6 moved 0.05 A in -y.

    F6 sets B5 of the 2->4 substituent, so the two frames give different values.
    """
    symbols, positions = sterimol_halo_chain_geometry
    moved = positions.copy()
    moved[5, 1] -= 0.05
    return write_xyz_file(tmp_path / "traj.xyz", symbols, [positions, moved])


@pytest.fixture()
def sterimol_toluene_molecule():
    """Toluene: C1 methyl, C2-C7 ring, H8-H10 methyl, H11-H15 ring."""
    return Molecule.from_rdkit_mol(rdkit_embed_and_optimize("Cc1ccccc1"))


@pytest.fixture()
def sterimol_toluene_xyz_file(tmp_path, sterimol_toluene_molecule):
    """Toluene as an xyz file."""
    return write_xyz_file(
        tmp_path / "toluene.xyz",
        sterimol_toluene_molecule.chemical_symbols,
        [sterimol_toluene_molecule.positions],
    )
