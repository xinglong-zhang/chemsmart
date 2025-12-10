import tempfile
from pathlib import Path

import pytest
import tomlkit

from chemsmart.cli.update import Updater


class TestUpdater:
    """Tests for the Updater class focusing on dependency normalization."""

    def test_extract_pkg_name_normalization(self):
        """Test that package names are normalized correctly (PEP 503)."""
        updater = Updater()
        
        # Create a temporary requirements file
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as tmp:
            tmp.write("scikit_learn>=1.6.1\n")
            tmp.write("pandas>=2.0.0\n")
            tmp.write("some-package>=1.0.0\n")
            tmp_path = Path(tmp.name)
        
        try:
            # Mock the existing dependencies with dash notation
            existing_deps = {
                "scikit-learn>=1.6.0",  # Different version, dash notation
                "pandas>=2.0.0",  # Same version
                "some_package>=0.9.0",  # Different version, underscore notation
            }
            
            # Temporarily replace the method to return our mock dependencies
            original_method = updater._get_existing_dependencies
            updater._get_existing_dependencies = lambda: existing_deps
            
            # Get missing dependencies
            missing = updater._get_missing_dependencies(tmp_path)
            
            # Restore original method
            updater._get_existing_dependencies = original_method
            
            # Assert no duplicates are found
            # scikit_learn should match scikit-learn
            # pandas should match pandas
            # some-package should match some_package
            assert len(missing) == 0, f"Expected no missing deps, got: {missing}"
            
        finally:
            tmp_path.unlink()
    
    def test_underscore_dash_equivalence(self):
        """Test that underscore and dash package names are treated as equivalent."""
        updater = Updater()
        
        # Create a temporary requirements file with underscore
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as tmp:
            tmp.write("pytest_mock>=3.14.0\n")
            tmp_path = Path(tmp.name)
        
        try:
            # Existing dependency with dash
            existing_deps = {"pytest-mock>=3.0.0"}
            
            original_method = updater._get_existing_dependencies
            updater._get_existing_dependencies = lambda: existing_deps
            
            missing = updater._get_missing_dependencies(tmp_path)
            
            updater._get_existing_dependencies = original_method
            
            # Should not consider pytest_mock as missing since pytest-mock exists
            assert len(missing) == 0, f"Expected no missing deps, got: {missing}"
            
        finally:
            tmp_path.unlink()
    
    def test_dash_underscore_equivalence(self):
        """Test that dash and underscore package names are treated as equivalent."""
        updater = Updater()
        
        # Create a temporary requirements file with dash
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as tmp:
            tmp.write("scikit-learn>=1.6.1\n")
            tmp_path = Path(tmp.name)
        
        try:
            # Existing dependency with underscore
            existing_deps = {"scikit_learn>=1.5.0"}
            
            original_method = updater._get_existing_dependencies
            updater._get_existing_dependencies = lambda: existing_deps
            
            missing = updater._get_missing_dependencies(tmp_path)
            
            updater._get_existing_dependencies = original_method
            
            # Should not consider scikit-learn as missing since scikit_learn exists
            assert len(missing) == 0, f"Expected no missing deps, got: {missing}"
            
        finally:
            tmp_path.unlink()
    
    def test_new_dependency_detected(self):
        """Test that genuinely new dependencies are still detected."""
        updater = Updater()
        
        # Create a temporary requirements file
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as tmp:
            tmp.write("new-package>=1.0.0\n")
            tmp.write("scikit-learn>=1.6.1\n")
            tmp_path = Path(tmp.name)
        
        try:
            # Existing dependencies without new-package
            existing_deps = {"scikit-learn>=1.6.0"}
            
            original_method = updater._get_existing_dependencies
            updater._get_existing_dependencies = lambda: existing_deps
            
            missing = updater._get_missing_dependencies(tmp_path)
            
            updater._get_existing_dependencies = original_method
            
            # Should detect new-package as missing, but not scikit-learn
            assert len(missing) == 1
            assert "new-package>=1.0.0" in missing
            
        finally:
            tmp_path.unlink()
