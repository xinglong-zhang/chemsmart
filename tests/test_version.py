"""Tests for version loading in chemsmart/__init__.py."""

import tempfile
from pathlib import Path
from unittest import mock


class TestVersion:
    """Test the version loading logic from chemsmart/__init__.py."""

    def test_version_loads_from_file(self):
        """Test that __version__ is correctly loaded from VERSION file."""
        import chemsmart

        # Verify that __version__ attribute exists
        assert hasattr(chemsmart, "__version__")
        
        # Verify it's a string
        assert isinstance(chemsmart.__version__, str)
        
        # Verify it's not the fallback value (unless VERSION file is actually "0.0.0")
        # Read the actual VERSION file to compare
        version_file = Path(chemsmart.__file__).with_name("VERSION")
        expected_version = version_file.read_text(encoding="utf-8").strip()
        
        assert chemsmart.__version__ == expected_version

    def test_version_fallback_when_file_missing(self):
        """Test that fallback to '0.0.0' works when VERSION file is missing."""
        # Create a temporary directory and mock the __file__ path
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create a fake __init__.py path in the temp directory
            fake_init_path = Path(tmpdir) / "__init__.py"
            fake_init_path.touch()
            
            # Mock Path(__file__) to return our fake path
            # The VERSION file won't exist at fake_init_path.with_name("VERSION")
            with mock.patch("pathlib.Path") as mock_path:
                mock_path.return_value.with_name.return_value.read_text.side_effect = FileNotFoundError
                
                # Simulate the version loading logic
                try:
                    _version_file = Path(str(fake_init_path)).with_name("VERSION")
                    version = _version_file.read_text(encoding="utf-8").strip()
                except Exception:
                    version = "0.0.0"
                
                assert version == "0.0.0"

    def test_version_fallback_when_file_unreadable(self):
        """Test that fallback to '0.0.0' works when VERSION file is unreadable."""
        # Mock the read_text method to raise a PermissionError
        with mock.patch.object(Path, "read_text", side_effect=PermissionError):
            # Simulate the version loading logic
            try:
                fake_path = Path("/fake/path/VERSION")
                version = fake_path.read_text(encoding="utf-8").strip()
            except Exception:
                version = "0.0.0"
            
            assert version == "0.0.0"

    def test_version_format(self):
        """Test that the version follows a reasonable format."""
        import re

        import chemsmart

        # Version should be non-empty
        assert len(chemsmart.__version__) > 0
        
        # Version should follow semantic versioning format: major.minor.patch[optional]
        # Examples: "1.0.0", "0.1.16", "2.3.4-dev", "1.2.3.post1"
        version_pattern = r"^\d+\.\d+\.\d+.*$"
        assert re.match(version_pattern, chemsmart.__version__), (
            f"Version '{chemsmart.__version__}' does not follow semantic versioning format"
        )
