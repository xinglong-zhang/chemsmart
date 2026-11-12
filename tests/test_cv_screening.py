"""Test molecular geometry measurements and collective_variable screening."""

import numpy as np
import pytest

from chemsmart.analysis.collective_variables import CollectiveVariables
from chemsmart.io.molecules.structure import Molecule

# Geometry measurements for one molecule


def test_distances_include_every_pair_and_update_after_coordinate_edits():
    molecule = Molecule(
        symbols=["H"] * 3,
        positions=[[0, 0, 0], [3, 0, 0], [0, 4, 0]],
    )

    # A 3-4-5 triangle gives independently known pair distances.
    expected = [[0, 3, 4], [3, 0, 5], [4, 5, 0]]
    np.testing.assert_allclose(molecule.distances_matrix, expected)

    # Reading the property again must use the edited coordinates.
    molecule.positions[1] = [0, 0, 4]
    expected = [[0, 4, 4], [4, 0, np.sqrt(32)], [4, np.sqrt(32), 0]]
    np.testing.assert_allclose(molecule.distances_matrix, expected)


@pytest.mark.parametrize(
    "property_name, expected",
    [
        ("bonds_matrix", [[1, 3, 1.5], [1, 4, 1.5], [2, 4, 1.5]]),
        ("angles_matrix", [[1, 4, 2, 90], [3, 1, 4, 90]]),
        ("dihedrals_matrix", [[2, 4, 1, 3, 90]]),
    ],
)
def test_connected_tables_use_sorted_one_based_paths(property_name, expected):
    molecule = Molecule(
        symbols=["C"] * 4,
        positions=[
            [1.5, 0, 0],
            [1.5, 1.5, 1.5],
            [0, 0, 0],
            [1.5, 1.5, 0],
        ],
    )

    # The actual bond detector finds the chain 3-1-4-2.
    # Reversed paths must appear only once, in a reproducible row order.
    np.testing.assert_allclose(getattr(molecule, property_name), expected)


def test_connected_tables_update_when_a_bond_breaks():
    molecule = Molecule(
        symbols=["C"] * 4,
        positions=[[0, 0, 0], [1.5, 0, 0], [1.5, 1.5, 0], [1.5, 1.5, 1.5]],
    )
    assert molecule.bonds_matrix.shape == (3, 3)
    assert molecule.angles_matrix.shape == (2, 4)
    assert molecule.dihedrals_matrix.shape == (1, 5)

    # Move the fourth atom away so the final bond is no longer detected.
    molecule.positions[3][2] = 10.0
    np.testing.assert_allclose(
        molecule.bonds_matrix, [[1, 2, 1.5], [2, 3, 1.5]]
    )
    np.testing.assert_allclose(molecule.angles_matrix, [[1, 2, 3, 90]])
    assert molecule.dihedrals_matrix.shape == (0, 5)


@pytest.mark.parametrize(
    "indices, expected_values",
    [
        ([[3, 4], [1, 3], [4, 3]], [10, np.sqrt(200), 10]),
        ([[3, 2, 1], [1, 3, 2], [1, 2, 3]], [90, 45, 90]),
        ([[4, 3, 2, 1], [1, 2, 3, 4]], [90, 90]),
    ],
)
def test_explicit_tables_preserve_input_order_without_bonds(
    indices, expected_values
):
    molecule = Molecule(
        symbols=["C"] * 4,
        positions=[[0, 0, 0], [10, 0, 0], [10, 10, 0], [10, 10, 10]],
    )
    assert molecule.bonds_matrix.shape == (0, 3)

    # Explicit tuples remain measurable even when their atoms are unbonded.
    table = molecule.geometry_table(indices)
    np.testing.assert_array_equal(table[:, :-1], indices)
    np.testing.assert_allclose(table[:, -1], expected_values)


@pytest.mark.parametrize("height, expected", [(1.5, 90), (-1.5, -90)])
def test_dihedral_sign_and_full_path_reversal(height, expected):
    molecule = Molecule(
        symbols=["C"] * 4,
        positions=[[0, 0, 0], [1.5, 0, 0], [1.5, 1.5, 0], [1.5, 1.5, height]],
    )

    # Reflection changes the sign. Reversing all four labels preserves it.
    table = molecule.geometry_table([[1, 2, 3, 4], [4, 3, 2, 1]])
    np.testing.assert_allclose(table[:, -1], [expected, expected])


def test_undefined_angles_keep_their_rows_as_nan():
    molecule = Molecule(
        symbols=["H"] * 4,
        positions=[[0, 0, 0], [0, 0, 0], [1, 0, 0], [1, 1, 0]],
    )
    indices = [[1, 2, 3], [2, 3, 4]]
    table = molecule.geometry_table(indices)

    # The first angle has a zero-length arm. The second is well defined.
    np.testing.assert_array_equal(table[:, :-1], indices)
    assert np.isnan(table[0, -1])
    assert table[1, -1] == pytest.approx(90)


def test_collinear_dihedral_is_nan_but_straight_angle_is_defined():
    molecule = Molecule(
        symbols=["C"] * 4,
        positions=[[0, 0, 0], [1.5, 0, 0], [3, 0, 0], [4.5, 0, 0]],
    )

    assert molecule.geometry_table([[1, 2, 3]])[0, -1] == pytest.approx(180)
    table = molecule.geometry_table([[1, 2, 3, 4]])
    np.testing.assert_array_equal(table[0, :-1], [1, 2, 3, 4])
    assert np.isnan(table[0, -1])


@pytest.mark.parametrize("order", [2, 3, 4])
def test_empty_explicit_tables_keep_their_column_count(order):
    molecule = Molecule(symbols=["Ar"], positions=[[0, 0, 0]])

    table = molecule.geometry_table(np.empty((0, order), dtype=int))
    assert table.shape == (0, order + 1)


def test_single_atom_has_no_connected_paths():
    molecule = Molecule(symbols=["Ar"], positions=[[0, 0, 0]])

    np.testing.assert_array_equal(molecule.distances_matrix, [[0]])
    assert molecule.bonds_matrix.shape == (0, 3)
    assert molecule.angles_matrix.shape == (0, 4)
    assert molecule.dihedrals_matrix.shape == (0, 5)


@pytest.mark.parametrize(
    "indices, error",
    [
        ([], ValueError),
        ([1, 2], ValueError),
        ([[1]], ValueError),
        ([[1, 2, 3, 4, 5]], ValueError),
        ([[1.0, 2.0]], TypeError),
        ([[True, False]], TypeError),
        ([[0, 2]], IndexError),
        ([[-1, 2]], IndexError),
        ([[1, 5]], IndexError),
        ([[1, 1]], ValueError),
    ],
)
def test_geometry_table_rejects_invalid_indices(indices, error):
    molecule = Molecule(
        symbols=["C"] * 4,
        positions=[[0, 0, 0], [1.5, 0, 0], [3, 0, 0], [4.5, 0, 0]],
    )

    with pytest.raises(error):
        molecule.geometry_table(indices)


@pytest.mark.parametrize(
    "property_name",
    ["distances_matrix", "bonds_matrix", "angles_matrix", "dihedrals_matrix"],
)
@pytest.mark.parametrize("invalid_value", [np.nan, np.inf])
def test_matrix_properties_reject_nonfinite_coordinates(
    property_name, invalid_value
):
    molecule = Molecule(symbols=["C", "C"], positions=[[0, 0, 0], [1.5, 0, 0]])
    molecule.positions[0][0] = invalid_value

    with pytest.raises(ValueError, match="finite N-by-3"):
        getattr(molecule, property_name)


def test_geometry_table_rejects_a_coordinate_count_mismatch():
    molecule = Molecule(symbols=["C", "C"], positions=[[0, 0, 0], [1.5, 0, 0]])
    molecule.positions = [[0, 0, 0]]

    with pytest.raises(ValueError, match="finite N-by-3"):
        molecule.geometry_table([[1, 2]])


# CV profiles across geometries


def test_profiles_align_coordinates_when_a_bond_forms():
    # C4 approaches C3. The C2-C3-C4 angle is measurable in both frames.
    frames = [
        Molecule(
            symbols=["C"] * 4,
            positions=np.array([[0, 0, 0], [1.5, 0, 0], [1.5, 1.5, 0], last]),
        )
        for last in ([1.5, 4.5, -3], [1.5, 1.5, 1.5])
    ]
    analysis = CollectiveVariables(frames)
    angles = analysis.angles
    np.testing.assert_array_equal(
        angles["atom_indices"], [[1, 2, 3], [2, 3, 4]]
    )
    np.testing.assert_array_equal(
        angles["connected_matrix"], [[True, False], [True, True]]
    )
    np.testing.assert_allclose(angles["values_matrix"], [[90, 135], [90, 90]])

    torsions = analysis.dihedrals
    np.testing.assert_array_equal(torsions["atom_indices"], [[1, 2, 3, 4]])
    np.testing.assert_array_equal(
        torsions["connected_matrix"][:, 0], [False, True]
    )
    np.testing.assert_allclose(torsions["values_matrix"][:, 0], [-90, 90])
    np.testing.assert_array_equal(torsions["valid_counts"], [2])


def test_circular_summary_and_refresh_after_a_coordinate_edit():
    # Rotate C4 around the C2-C3 axis. The shortest change across 180° is +2°.
    frames = []
    for degrees in (179, -179, -179):
        angle = np.deg2rad(degrees)
        frames.append(
            Molecule(
                symbols=["C"] * 4,
                positions=np.array(
                    [
                        [0, 1.5, 0],
                        [0, 0, 0],
                        [1.5, 0, 0],
                        [1.5, 1.5 * np.cos(angle), 1.5 * np.sin(angle)],
                    ]
                ),
            )
        )
    analysis = CollectiveVariables(frames)
    data = analysis.dihedrals
    summary = dict(zip(data["summary_columns"], data["summary_matrix"][0]))
    assert summary["signed_change"] == pytest.approx(2)
    assert summary["absolute_change"] == pytest.approx(2)
    assert summary["circular_range"] == pytest.approx(2)

    # A collinear final bond makes the endpoint torsion undefined.
    frames[-1].positions[-1] = [3, 0, 0]
    analysis.refresh()
    data = analysis.dihedrals
    np.testing.assert_allclose(
        data["values_matrix"][:, 0], [179, -179, np.nan], equal_nan=True
    )
    np.testing.assert_array_equal(
        data["connected_matrix"][:, 0], [True, True, True]
    )
    np.testing.assert_array_equal(data["valid_counts"], [2])
    summary = dict(zip(data["summary_columns"], data["summary_matrix"][0]))
    assert np.isnan(summary["end"])
    assert np.isnan(summary["signed_change"])
    assert np.isnan(summary["absolute_change"])
    assert summary["circular_range"] == pytest.approx(2)


def test_xyz_reader_preserves_frame_order_and_source_after_refresh(tmp_path):
    path = tmp_path / "hcl.xyz"
    path.write_text(
        "2\nfirst\nH 0 0 0\nCl 0 0 1.28\n"
        "2\nsecond\nH 0 0 0\nCl 0 0 4.00\n"
        "2\nthird\nH 0 0 0\nCl 0 0 1.20\n",
        encoding="utf-8",
    )
    analysis = CollectiveVariables.from_filepath(path)
    distances = analysis.distances
    np.testing.assert_array_equal(distances["atom_indices"], [[1, 2]])
    np.testing.assert_allclose(
        distances["values_matrix"][:, 0], [1.28, 4, 1.20]
    )
    np.testing.assert_array_equal(
        distances["connected_matrix"][:, 0], [True, False, True]
    )
    np.testing.assert_array_equal(distances["valid_counts"], [3])

    # The intermediate frame sets the range. Only the endpoints set the change.
    summary = dict(
        zip(distances["summary_columns"], distances["summary_matrix"][0])
    )
    assert summary["signed_change"] == pytest.approx(-0.08)
    assert summary["absolute_change"] == pytest.approx(0.08)
    assert summary["range"] == pytest.approx(2.8)

    analysis.molecules[-1].positions[1, 2] = 1.30
    analysis.refresh()
    (selected,) = analysis.select([["cl2"]])
    assert selected["selection_labels"] == ("Cl2",)
    assert selected["distances"]["values_matrix"][-1, 0] == pytest.approx(1.30)
    assert analysis.profile["source_file"] == str(path)


def test_screening_requires_molecules_with_matching_element_order():
    with pytest.raises(ValueError):
        CollectiveVariables([])
    with pytest.raises(TypeError):
        CollectiveVariables([object()])

    positions = [[0, 0, 0], [0, 0, 1.28]]
    first = Molecule(symbols=["H", "Cl"], positions=positions)
    second = Molecule(symbols=["Cl", "H"], positions=positions)
    with pytest.raises(ValueError, match="element sequences"):
        CollectiveVariables([first, second])


# Atom selections and per-atom maxima


def test_select_ordered_triple_and_independent_atom(ethylene_molecule):
    analysis = CollectiveVariables([ethylene_molecule])
    triple, single = analysis.select([["h3", "C1", "C2"], ["H6"]])

    # The middle label defines the angle. All three pair distances survive.
    assert triple["selection_labels"] == ("H3", "C1", "C2")
    np.testing.assert_array_equal(
        triple["distances"]["atom_indices"], [[1, 2], [1, 3], [2, 3]]
    )
    np.testing.assert_array_equal(
        triple["angles"]["atom_indices"], [[2, 1, 3]]
    )
    np.testing.assert_array_equal(
        triple["dihedrals"]["atom_indices"], [[3, 1, 2, 5], [3, 1, 2, 6]]
    )
    assert single["selection_labels"] == ("H6",)
    np.testing.assert_array_equal(
        single["distances"]["atom_indices"],
        [[1, 6], [2, 6], [3, 6], [4, 6], [5, 6]],
    )


def test_select_central_dihedral_bond(ethylene_molecule):
    analysis = CollectiveVariables([ethylene_molecule])
    (outer,) = analysis.select([["H3", "C1"]])
    peripheral, central = analysis.select(
        [["H3", "C1"], ["C2", "C1"]], central_bond_only=True
    )

    assert len(outer["dihedrals"]["atom_indices"]) == 2
    assert peripheral["dihedrals"]["atom_indices"].shape == (0, 4)
    np.testing.assert_array_equal(
        central["dihedrals"]["atom_indices"],
        [[3, 1, 2, 5], [3, 1, 2, 6], [4, 1, 2, 5], [4, 1, 2, 6]],
    )


def test_atom_maxima_keep_ties_and_include_unbonded_pairs_on_request(
    co2_molecule,
):
    expanded = co2_molecule.copy()
    expanded.positions[1:, 2] = [1.26, -1.26]
    analysis = CollectiveVariables([co2_molecule, expanded])

    # Both C–O bonds grow by 0.1 Å. The O–O distance grows by 0.2 Å.
    bonded = analysis.atom_maximum_changes()
    all_pairs = analysis.atom_maximum_changes(connected_distances=False)
    assert bonded["atom_labels"] == ("C1", "O2", "O3")
    np.testing.assert_allclose(bonded["values_matrix"][:, 0], [0.1, 0.1, 0.1])
    np.testing.assert_allclose(
        all_pairs["values_matrix"][:, 0], [0.1, 0.2, 0.2]
    )
    assert bonded["winner_columns"]["distances"][0] == (0, 1)
    assert all_pairs["winner_columns"]["distances"][1] == (2,)
    np.testing.assert_array_equal(bonded["candidate_counts"][:, 2], [0, 0, 0])
    assert np.isnan(bonded["values_matrix"][:, 2]).all()


def test_atom_maxima_report_undefined_linear_dihedral(acetylene_molecule):
    result = CollectiveVariables([acetylene_molecule]).atom_maximum_changes()

    np.testing.assert_array_equal(
        result["candidate_counts"][:, 2], [1, 1, 1, 1]
    )
    np.testing.assert_array_equal(
        result["undefined_counts"][:, 2], [1, 1, 1, 1]
    )
    assert np.isnan(result["values_matrix"][:, 2]).all()
    assert result["winner_columns"]["dihedrals"] == ((), (), (), ())


@pytest.mark.parametrize(
    "groups, central_bond_only",
    [
        (["C1"], False),
        ([["C1", "C1"]], False),
        ([["O1"]], False),
        ([["C7"]], False),
        ([["C0"]], False),
        ([["C1", "C2", "H3"]], True),
    ],
)
def test_select_rejects_invalid_atom_groups(
    ethylene_molecule, groups, central_bond_only
):
    analysis = CollectiveVariables([ethylene_molecule])
    with pytest.raises(ValueError):
        analysis.select(groups, central_bond_only=central_bond_only)
