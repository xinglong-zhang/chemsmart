from chemsmart.utils.periodictable import (
    METALLOIDS,
    NON_METALS_AND_METALLOIDS,
    NONMETALS,
    PeriodicTable,
    is_metal,
    is_metal_atomic_number,
    metal_element_symbols,
)

p = PeriodicTable()


class TestPeriodicTable:
    def test_get_elements(self):
        pass


class TestElementClassifications:
    """Tests for element classification constants."""

    def test_nonmetals_is_frozenset(self):
        """Test that NONMETALS is a frozenset (immutable)."""
        assert isinstance(NONMETALS, frozenset)

    def test_metalloids_is_frozenset(self):
        """Test that METALLOIDS is a frozenset (immutable)."""
        assert isinstance(METALLOIDS, frozenset)

    def test_combined_is_frozenset(self):
        """Test that NON_METALS_AND_METALLOIDS is a frozenset."""
        assert isinstance(NON_METALS_AND_METALLOIDS, frozenset)

    def test_nonmetals_count(self):
        """Test that NONMETALS contains exactly 18 elements."""
        assert len(NONMETALS) == 18

    def test_metalloids_count(self):
        """Test that METALLOIDS contains exactly 7 elements."""
        assert len(METALLOIDS) == 7

    def test_combined_count(self):
        """Test that combined set contains 25 elements (18 + 7)."""
        assert len(NON_METALS_AND_METALLOIDS) == 25

    def test_nonmetals_contains_hydrogen(self):
        """Test that hydrogen (1) is classified as a non-metal."""
        assert 1 in NONMETALS

    def test_nonmetals_contains_noble_gases(self):
        """Test that noble gases are classified as non-metals."""
        noble_gases = {2, 10, 18, 36, 54, 86}  # He, Ne, Ar, Kr, Xe, Rn
        assert noble_gases.issubset(NONMETALS)

    def test_nonmetals_contains_halogens(self):
        """Test that halogens are classified as non-metals."""
        halogens = {9, 17, 35, 53, 85}  # F, Cl, Br, I, At
        assert halogens.issubset(NONMETALS)

    def test_nonmetals_contains_main_group(self):
        """Test that main group non-metals are included."""
        main_group = {6, 7, 8, 15, 16, 34}  # C, N, O, P, S, Se
        assert main_group.issubset(NONMETALS)

    def test_metalloids_contains_expected_elements(self):
        """Test that all expected metalloids are present."""
        expected_metalloids = {
            5,
            14,
            32,
            33,
            51,
            52,
            84,
        }  # B, Si, Ge, As, Sb, Te, Po
        assert METALLOIDS == expected_metalloids

    def test_no_overlap_between_nonmetals_and_metalloids(self):
        """Test that NONMETALS and METALLOIDS don't overlap."""
        assert NONMETALS.isdisjoint(METALLOIDS)

    def test_combined_is_union(self):
        """Test that combined set is the union of nonmetals and metalloids."""
        assert NON_METALS_AND_METALLOIDS == NONMETALS | METALLOIDS

    def test_metal_classification_iron(self):
        """Test that iron (26) is correctly classified as a metal."""
        assert 26 not in NON_METALS_AND_METALLOIDS

    def test_metal_classification_copper(self):
        """Test that copper (29) is correctly classified as a metal."""
        assert 29 not in NON_METALS_AND_METALLOIDS

    def test_metal_classification_gold(self):
        """Test that gold (79) is correctly classified as a metal."""
        assert 79 not in NON_METALS_AND_METALLOIDS

    def test_metal_classification_titanium(self):
        """Test that titanium (22) is correctly classified as a metal."""
        assert 22 not in NON_METALS_AND_METALLOIDS

    def test_nonmetal_classification_carbon(self):
        """Test that carbon (6) is correctly classified as a non-metal."""
        assert 6 in NONMETALS
        assert 6 in NON_METALS_AND_METALLOIDS

    def test_nonmetal_classification_nitrogen(self):
        """Test that nitrogen (7) is correctly classified as a non-metal."""
        assert 7 in NONMETALS
        assert 7 in NON_METALS_AND_METALLOIDS

    def test_nonmetal_classification_oxygen(self):
        """Test that oxygen (8) is correctly classified as a non-metal."""
        assert 8 in NONMETALS
        assert 8 in NON_METALS_AND_METALLOIDS

    def test_metalloid_classification_silicon(self):
        """Test that silicon (14) is correctly classified as a metalloid."""
        assert 14 not in NONMETALS
        assert 14 in METALLOIDS
        assert 14 in NON_METALS_AND_METALLOIDS

    def test_metalloid_classification_arsenic(self):
        """Test that arsenic (33) is correctly classified as a metalloid."""
        assert 33 not in NONMETALS
        assert 33 in METALLOIDS
        assert 33 in NON_METALS_AND_METALLOIDS

    def test_metalloid_classification_antimony(self):
        """Test that antimony (51) is correctly classified as a metalloid."""
        assert 51 not in NONMETALS
        assert 51 in METALLOIDS
        assert 51 in NON_METALS_AND_METALLOIDS

    def test_all_elements_positive(self):
        """Test that all atomic numbers are positive."""
        all_elements = NONMETALS | METALLOIDS
        assert all(element > 0 for element in all_elements)

    def test_all_elements_reasonable_range(self):
        """Test that all atomic numbers are in a reasonable range (1-118)."""
        all_elements = NONMETALS | METALLOIDS
        assert all(1 <= element <= 118 for element in all_elements)

    def test_comprehensive_coverage(self):
        """Test that classification is comprehensive for common elements."""
        # Test a selection of elements across the periodic table
        test_cases = {
            # (atomic_number, expected_classification)
            (1, "nonmetal"),  # H
            (2, "nonmetal"),  # He
            (3, "metal"),  # Li
            (4, "metal"),  # Be
            (5, "metalloid"),  # B
            (6, "nonmetal"),  # C
            (11, "metal"),  # Na
            (12, "metal"),  # Mg
            (13, "metal"),  # Al
            (14, "metalloid"),  # Si
            (15, "nonmetal"),  # P
            (16, "nonmetal"),  # S
            (17, "nonmetal"),  # Cl
            (18, "nonmetal"),  # Ar
            (19, "metal"),  # K
            (20, "metal"),  # Ca
            (26, "metal"),  # Fe
            (29, "metal"),  # Cu
            (30, "metal"),  # Zn
            (34, "nonmetal"),  # Se
            (35, "nonmetal"),  # Br
            (47, "metal"),  # Ag
            (79, "metal"),  # Au
            (82, "metal"),  # Pb
        }

        for atomic_num, expected in test_cases:
            if expected == "nonmetal":
                assert (
                    atomic_num in NONMETALS
                ), f"Element {atomic_num} should be a non-metal"
                assert (
                    atomic_num not in METALLOIDS
                ), f"Element {atomic_num} should not be a metalloid"
            elif expected == "metalloid":
                assert (
                    atomic_num in METALLOIDS
                ), f"Element {atomic_num} should be a metalloid"
                assert (
                    atomic_num not in NONMETALS
                ), f"Element {atomic_num} should not be a non-metal"
            elif expected == "metal":
                assert (
                    atomic_num not in NONMETALS
                ), f"Element {atomic_num} should not be a non-metal"
                assert (
                    atomic_num not in METALLOIDS
                ), f"Element {atomic_num} should not be a metalloid"
                assert (
                    atomic_num not in NON_METALS_AND_METALLOIDS
                ), f"Element {atomic_num} should be a metal"


class TestMetalHelpers:
    def test_is_metal_atomic_number(self):
        assert is_metal_atomic_number(26)
        assert not is_metal_atomic_number(6)
        assert not is_metal_atomic_number(14)

    def test_is_metal_symbol(self):
        assert is_metal("Fe")
        assert is_metal("fe")
        assert is_metal("Mn")
        assert not is_metal("C")
        assert not is_metal("Si")
        assert not is_metal("Xx")
        assert not is_metal("invalid")
        assert not is_metal("Zz")

    def test_periodic_table_is_metal(self):
        assert p.is_metal("Cu")
        assert not p.is_metal("N")

    def test_metal_element_symbols_includes_lanthanides(self):
        symbols = metal_element_symbols()
        assert "La" in symbols
        assert "Mn" in symbols
        assert "C" not in symbols
        assert "Si" not in symbols

    def test_is_metal_empty_string_is_false(self):
        assert is_metal("") is False
        assert is_metal("   ") is False

    def test_is_metal_unparseable_symbol_is_false(self):
        """A string with no alphabetic characters left after cleaning
        makes to_element raise ValueError, which is_metal must catch
        and treat as "not a metal" rather than propagating."""
        assert is_metal("123") is False


class TestPeriodicTableConversions:
    def test_to_element_empty_string_raises(self):
        import pytest

        with pytest.raises(ValueError, match="cannot be empty"):
            p.to_element("")

    def test_to_element_unparseable_raises(self):
        import pytest

        with pytest.raises(ValueError, match="Unable to parse element"):
            p.to_element("123")

    def test_sorted_periodic_table_list(self):
        assert p.sorted_periodic_table_list(["O", "H", "C"]) == [
            "H",
            "C",
            "O",
        ]

    def test_to_weighted_atomic_mass_by_abundance(self):
        mass = p.to_weighted_atomic_mass_by_abundance("C")
        assert isinstance(mass, float)
        assert 11.0 < mass < 13.0

    def test_to_most_abundant_atomic_mass(self):
        mass = p.to_most_abundant_atomic_mass("C")
        assert isinstance(mass, float)
        assert 11.0 < mass < 13.0

    def test_requires_ecp_false_for_light_element(self):
        assert p.requires_ecp("C") is False

    def test_requires_ecp_true_for_heavy_element(self):
        assert p.requires_ecp("Pd") is True

    def test_atomic_masses_property(self):
        masses = p.atomic_masses
        assert masses[p.to_atomic_number("C")] == p.to_atomic_mass("C")

    def test_to_atomic_mass(self):
        mass = p.to_atomic_mass("C")
        assert isinstance(mass, float)
        assert 11.0 < mass < 13.0

    def test_vdw_radius(self):
        radius = p.vdw_radius("C")
        assert isinstance(radius, float)
        assert radius > 0

    def test_covalent_radius(self):
        radius = p.covalent_radius("C")
        assert isinstance(radius, float)
        assert radius > 0
