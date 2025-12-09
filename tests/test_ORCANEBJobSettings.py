"""
Tests for ORCA NEB job settings.
"""

import pytest

from chemsmart.jobs.orca.settings import ORCAJobSettings, ORCANEBJobSettings


class TestORCANEBJobSettings:
    """Test suite for ORCANEBJobSettings class."""

    def test_init_default(self):
        """Test default initialization."""
        settings = ORCANEBJobSettings()

        # Check default values
        assert settings.jobtype is None
        assert settings.nimages is None
        assert settings.ending_xyzfile is None
        assert settings.intermediate_xyzfile is None
        assert settings.starting_xyz is None
        assert settings.restarting_xyzfile is None
        assert settings.preopt_ends is False
        assert settings.semiempirical is None

    def test_init_with_parameters(self):
        """Test initialization with specific parameters."""
        settings = ORCANEBJobSettings(
            jobtype="NEB-TS",
            nimages=8,
            ending_xyzfile="product.xyz",
            intermediate_xyzfile="ts_guess.xyz",
            starting_xyz="reactant.xyz",
            preopt_ends=True,
            semiempirical="XTB2",
        )

        assert settings.jobtype == "NEB-TS"
        assert settings.nimages == 8
        assert settings.ending_xyzfile == "product.xyz"
        assert settings.intermediate_xyzfile == "ts_guess.xyz"
        assert settings.starting_xyz == "reactant.xyz"
        assert settings.preopt_ends is True
        assert settings.semiempirical == "XTB2"

    def test_route_string_with_semiempirical(self):
        """Test route string generation with semiempirical method."""
        settings = ORCANEBJobSettings(jobtype="NEB-CI", semiempirical="XTB2")

        expected = "! XTB2 NEB-CI"
        assert settings.route_string == expected

    def test_route_string_with_dft(self):
        """Test route string generation with DFT method."""
        settings = ORCANEBJobSettings(
            jobtype="NEB-TS", functional="B3LYP", basis="def2-SVP"
        )

        expected = "! B3LYP def2-SVP NEB-TS"
        assert settings.route_string == expected

    def test_route_string_no_semiempirical(self):
        """Test route string generation without semiempirical method."""
        settings = ORCANEBJobSettings(
            jobtype="NEB", functional="PBE0", basis="def2-TZVP"
        )

        expected = "! PBE0 def2-TZVP NEB"
        assert settings.route_string == expected

    def test_neb_block_basic(self):
        """Test basic NEB block generation."""
        settings = ORCANEBJobSettings(
            nimages=5, starting_xyz="start.xyz", ending_xyzfile="end.xyz"
        )

        neb_block = settings.neb_block

        # Check that required components are present
        assert "%neb" in neb_block
        assert "NImages 5" in neb_block
        assert "NEB_END_XYZFile 'end.xyz'" in neb_block
        assert "PREOPT_ENDS False" in neb_block
        assert "END" in neb_block

    def test_neb_block_with_intermediate(self):
        """Test NEB block generation with intermediate geometry."""
        settings = ORCANEBJobSettings(
            nimages=8,
            starting_xyz="reactant.xyz",
            ending_xyzfile="product.xyz",
            intermediate_xyzfile="ts.xyz",
            preopt_ends=True,
        )

        neb_block = settings.neb_block

        assert "NImages 8" in neb_block
        assert "NEB_END_XYZFile 'product.xyz'" in neb_block
        assert "NEB_TS_XYZFILE ts.xyz" in neb_block
        assert "PREOPT_ENDS True" in neb_block

    def test_neb_block_restart(self):
        """Test NEB block generation for restart calculation."""
        settings = ORCANEBJobSettings(
            nimages=6, restarting_xyzfile="restart_path.xyz"
        )

        neb_block = settings.neb_block

        assert "NImages 6" in neb_block
        assert "Restart_ALLXYZFile 'restart_path.xyz'" in neb_block
        assert "END" in neb_block

    def test_neb_block_missing_nimages_raises_error(self):
        """Test that missing nimages raises assertion error."""
        settings = ORCANEBJobSettings(
            starting_xyz="start.xyz", ending_xyzfile="end.xyz"
        )

        with pytest.raises(
            AssertionError, match="The number of images is missing"
        ):
            _ = settings.neb_block

    def test_neb_block_missing_geometries_raises_error(self):
        """Test that missing geometry files raises assertion error."""
        settings = ORCANEBJobSettings(nimages=5)

        with pytest.raises(
            AssertionError, match="No valid input geomertry is given"
        ):
            _ = settings.neb_block

    def test_neb_block_with_restart_only(self):
        """Test that restart file alone is sufficient for geometry."""
        settings = ORCANEBJobSettings(
            nimages=4, restarting_xyzfile="restart.xyz"
        )

        # Should not raise error
        neb_block = settings.neb_block
        assert "Restart_ALLXYZFile 'restart.xyz'" in neb_block

    def test_neb_block_missing_ending_file_raises_error(self):
        """Test that missing ending file (with starting file) raises error."""
        settings = ORCANEBJobSettings(
            nimages=5,
            starting_xyz="start.xyz",
            # missing ending_xyzfile
        )

        with pytest.raises(
            AssertionError, match="No valid input geomertry is given"
        ):
            _ = settings.neb_block

    def test_neb_block_missing_starting_file_raises_error(self):
        """Test that missing starting file (with ending file) raises error."""
        settings = ORCANEBJobSettings(
            nimages=5,
            ending_xyzfile="end.xyz",
            # missing starting_xyz
        )

        with pytest.raises(
            AssertionError, match="No valid input geomertry is given"
        ):
            _ = settings.neb_block

    def test_inheritance_from_orca_job_settings(self):
        """Test that ORCANEBJobSettings properly inherits from ORCAJobSettings."""
        settings = ORCANEBJobSettings(
            functional="B3LYP", basis="def2-SVP", charge=0, multiplicity=1
        )

        # Should have inherited attributes
        assert isinstance(settings, ORCAJobSettings)
        assert settings.functional == "B3LYP"
        assert settings.basis == "def2-SVP"
        assert settings.charge == 0
        assert settings.multiplicity == 1

    def test_supported_job_types(self):
        """Test various supported NEB job types."""
        job_types = [
            "NEB",
            "NEB-CI",
            "NEB-TS",
            "FAST-NEB-TS",
            "TIGHT-NEB-TS",
            "LOOSE-NEB",
            "ZOOM-NEB",
            "ZOOM-NEB-CI",
            "ZOOM-NEB-TS",
            "NEB-IDPP",
        ]

        for job_type in job_types:
            settings = ORCANEBJobSettings(
                jobtype=job_type, semiempirical="XTB2"
            )
            assert job_type in settings.route_string

    def test_neb_block_complete_example(self):
        """Test complete NEB block generation with all options."""
        settings = ORCANEBJobSettings(
            nimages=10,
            starting_xyz="reactant.xyz",
            ending_xyzfile="product.xyz",
            intermediate_xyzfile="ts_initial.xyz",
            preopt_ends=True,
        )

        neb_block = settings.neb_block
        expected_lines = [
            "%neb",
            "NImages 10",
            "NEB_END_XYZFile 'product.xyz'",
            "PREOPT_ENDS True",
            "NEB_TS_XYZFILE ts_initial.xyz",
            "END",
        ]

        for line in expected_lines:
            assert line in neb_block

    def test_merge_functionality(self):
        """Test that merge functionality returns valid NEB settings object."""
        base_settings = ORCANEBJobSettings(
            jobtype="NEB", nimages=5, functional="B3LYP"
        )

        update_dict = {"basis": "def2-TZVP"}

        merged = base_settings.merge(update_dict)

        # Should return a valid ORCANEBJobSettings object
        assert isinstance(merged, ORCANEBJobSettings)

        # Original NEB attributes should be preserved
        assert merged.jobtype == "NEB"
        assert merged.nimages == 5

        # Should have original functional (merge may not update all attributes)
        assert merged.functional == "B3LYP"

    def test_from_dict_creation(self):
        """Test creating ORCANEBJobSettings from dictionary."""
        neb_dict = {
            "jobtype": "NEB-CI",
            "nimages": 8,
            "starting_xyz": "react.xyz",
            "ending_xyzfile": "prod.xyz",
            "semiempirical": "XTB2",
            "preopt_ends": False,
        }

        settings = ORCANEBJobSettings.from_dict(neb_dict)

        assert settings.jobtype == "NEB-CI"
        assert settings.nimages == 8
        assert settings.starting_xyz == "react.xyz"
        assert settings.ending_xyzfile == "prod.xyz"
        assert settings.semiempirical == "XTB2"
        assert settings.preopt_ends is False

    def test_copy_functionality(self):
        """Test that copy functionality works correctly."""
        original = ORCANEBJobSettings(
            jobtype="NEB-TS", nimages=6, ending_xyzfile="end.xyz"
        )

        copied = original.copy()

        # Should be equal but different objects
        assert copied.jobtype == original.jobtype
        assert copied.nimages == original.nimages
        assert copied.ending_xyzfile == original.ending_xyzfile
        assert copied is not original

    def test_route_string_property_consistency(self):
        """Test that route_string property is consistent with internal method."""
        settings = ORCANEBJobSettings(
            jobtype="ZOOM-NEB", functional="M06-2X", basis="def2-TZVP"
        )

        # Both should give the same result
        assert settings.route_string == settings._get_neb_route_string()

    def test_neb_block_property_consistency(self):
        """Test that neb_block property is consistent with internal method."""
        settings = ORCANEBJobSettings(
            nimages=7, starting_xyz="start.xyz", ending_xyzfile="end.xyz"
        )

        # Both should give the same result
        assert settings.neb_block == settings._write_neb_block()

    def test_edge_case_empty_strings(self):
        """Test behavior with empty string values."""
        settings = ORCANEBJobSettings(
            nimages=5,
            starting_xyz="",  # Empty string
            ending_xyzfile="end.xyz",
        )

        # Should treat empty string as falsy in validation
        with pytest.raises(
            AssertionError, match="No valid input geomertry is given"
        ):
            _ = settings.neb_block

    def test_attribute_modification(self):
        """Test that attributes can be modified after initialization."""
        settings = ORCANEBJobSettings()

        # Modify attributes
        settings.jobtype = "NEB-TS"
        settings.nimages = 12
        settings.ending_xyzfile = "modified.xyz"

        assert settings.jobtype == "NEB-TS"
        assert settings.nimages == 12
        assert settings.ending_xyzfile == "modified.xyz"

    def test_inheritance_kwargs_passthrough(self):
        """Test that kwargs are properly passed to parent class."""
        settings = ORCANEBJobSettings(
            jobtype="NEB",
            # These should be passed to ORCAJobSettings
            dispersion="D3BJ",
            defgrid="2",
            forces=True,
        )

        # Should have inherited attributes from kwargs
        assert settings.jobtype == "NEB"  # NEB-specific
        assert settings.dispersion == "D3BJ"  # From parent
        assert settings.defgrid == "2"  # From parent
        assert settings.forces is True  # From parent

    def test_neb_block_validation_combinations(self):
        """Test various valid and invalid combinations for NEB block."""
        # Valid: start + end
        settings1 = ORCANEBJobSettings(
            nimages=3, starting_xyz="start.xyz", ending_xyzfile="end.xyz"
        )
        assert "NEB_END_XYZFile" in settings1.neb_block  # Should work

        # Valid: restart only
        settings2 = ORCANEBJobSettings(
            nimages=3, restarting_xyzfile="restart.xyz"
        )
        assert "Restart_ALLXYZFile" in settings2.neb_block  # Should work

        # Invalid: only start
        settings3 = ORCANEBJobSettings(nimages=3, starting_xyz="start.xyz")
        with pytest.raises(AssertionError):
            _ = settings3.neb_block

        # Invalid: only end
        settings4 = ORCANEBJobSettings(nimages=3, ending_xyzfile="end.xyz")
        with pytest.raises(AssertionError):
            _ = settings4.neb_block

    def test_docstring_completeness(self):
        """Test that the class has proper documentation."""
        assert ORCANEBJobSettings.__doc__ is not None
        assert "NEB" in ORCANEBJobSettings.__doc__
        assert "Nudged Elastic Band" in ORCANEBJobSettings.__doc__
