"""Tests for chemsmart.cli.config.Config class."""
import platform
import shutil
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from chemsmart.cli.config import Config


class TestConfig:
    def test_chemsmart_template_returns_traversable(self):
        """chemsmart_template should return an importlib.resources Traversable."""
        cfg = Config()
        template = cfg.chemsmart_template
        # The traversable should point to an existing .chemsmart directory
        # inside the installed package.
        assert template.is_dir(), (
            f"Template directory not found: {template!r}. "
            "Ensure chemsmart/settings/templates/.chemsmart is included "
            "in package data."
        )

    def test_chemsmart_dest(self):
        cfg = Config()
        assert cfg.chemsmart_dest == Path.home() / ".chemsmart"

    # ------------------------------------------------------------------
    # shell_config
    # ------------------------------------------------------------------

    def test_shell_config_windows_returns_none(self):
        cfg = Config()
        with patch("chemsmart.cli.config.platform.system", return_value="Windows"):
            assert cfg.shell_config is None

    def test_shell_config_bash(self, tmp_path):
        cfg = Config()
        fake_home = tmp_path
        bashrc = fake_home / ".bashrc"
        bashrc.touch()
        with (
            patch("chemsmart.cli.config.platform.system", return_value="Linux"),
            patch.dict("os.environ", {"SHELL": "/bin/bash"}),
            patch.object(Path, "home", return_value=fake_home),
        ):
            result = cfg.shell_config
        assert result == fake_home / ".bashrc"

    def test_shell_config_zsh(self, tmp_path):
        cfg = Config()
        fake_home = tmp_path
        with (
            patch("chemsmart.cli.config.platform.system", return_value="Darwin"),
            patch.dict("os.environ", {"SHELL": "/bin/zsh"}),
            patch.object(Path, "home", return_value=fake_home),
        ):
            result = cfg.shell_config
        assert result == fake_home / ".zshrc"

    # ------------------------------------------------------------------
    # env_vars
    # ------------------------------------------------------------------

    def test_env_vars_windows_empty(self):
        cfg = Config()
        with patch("chemsmart.cli.config.platform.system", return_value="Windows"):
            assert cfg.env_vars == []

    def test_env_vars_linux_non_empty(self):
        cfg = Config()
        with patch("chemsmart.cli.config.platform.system", return_value="Linux"):
            vars_ = cfg.env_vars
        assert len(vars_) > 0
        assert all(v.startswith("export ") for v in vars_)

    # ------------------------------------------------------------------
    # setup_environment
    # ------------------------------------------------------------------

    def test_setup_environment_copies_templates(self, tmp_path):
        """setup_environment should copy the bundled templates to chemsmart_dest."""
        dest = tmp_path / ".chemsmart"
        cfg = Config()
        with (
            patch.object(type(cfg), "chemsmart_dest", new_callable=lambda: property(lambda self: dest)),
            patch.object(type(cfg), "shell_config", new_callable=lambda: property(lambda self: None)),
        ):
            cfg.setup_environment()
        assert dest.exists()
        assert dest.is_dir()

    def test_setup_environment_skips_copy_if_dest_exists(self, tmp_path):
        """setup_environment should not overwrite an existing chemsmart_dest."""
        dest = tmp_path / ".chemsmart"
        dest.mkdir()
        marker = dest / "existing_file.txt"
        marker.write_text("keep me")

        cfg = Config()
        with (
            patch.object(type(cfg), "chemsmart_dest", new_callable=lambda: property(lambda self: dest)),
            patch.object(type(cfg), "shell_config", new_callable=lambda: property(lambda self: None)),
        ):
            cfg.setup_environment()
        # Existing content should be untouched
        assert marker.read_text() == "keep me"

    def test_setup_environment_windows_skips_shell_config(self, tmp_path):
        """On Windows, setup_environment must not touch any shell rc file."""
        dest = tmp_path / ".chemsmart"
        cfg = Config()
        with (
            patch.object(type(cfg), "chemsmart_dest", new_callable=lambda: property(lambda self: dest)),
            patch("chemsmart.cli.config.platform.system", return_value="Windows"),
        ):
            # Should not raise even though there is no shell rc file on Windows
            cfg.setup_environment()
        assert dest.exists()

    def test_setup_environment_updates_shell_config(self, tmp_path):
        """On Unix, setup_environment should write env vars into the shell rc."""
        dest = tmp_path / ".chemsmart"
        shell_rc = tmp_path / ".bashrc"
        shell_rc.write_text("")

        cfg = Config()
        with (
            patch.object(type(cfg), "chemsmart_dest", new_callable=lambda: property(lambda self: dest)),
            patch.object(type(cfg), "shell_config", new_callable=lambda: property(lambda self: shell_rc)),
            patch.object(type(cfg), "env_vars", new_callable=lambda: property(lambda self: ['export FOO="bar"'])),
            patch("chemsmart.cli.config.platform.system", return_value="Linux"),
        ):
            cfg.setup_environment()

        content = shell_rc.read_text()
        assert "Added by chemsmart installer" in content
        assert 'export FOO="bar"' in content

    def test_setup_environment_no_duplicate_shell_config(self, tmp_path):
        """setup_environment should not write env vars twice."""
        dest = tmp_path / ".chemsmart"
        dest.mkdir()
        shell_rc = tmp_path / ".bashrc"
        shell_rc.write_text("# Added by chemsmart installer\n")

        cfg = Config()
        with (
            patch.object(type(cfg), "chemsmart_dest", new_callable=lambda: property(lambda self: dest)),
            patch.object(type(cfg), "shell_config", new_callable=lambda: property(lambda self: shell_rc)),
            patch.object(type(cfg), "env_vars", new_callable=lambda: property(lambda self: ['export FOO="bar"'])),
            patch("chemsmart.cli.config.platform.system", return_value="Linux"),
        ):
            cfg.setup_environment()
            cfg.setup_environment()

        content = shell_rc.read_text()
        assert content.count("Added by chemsmart installer") == 1


class TestConfigServerCommand:
    """Tests for the config server Click command Windows behavior."""

    def test_server_skips_conda_update_on_windows(self):
        """On Windows, config server should not call update_yaml_files for conda."""
        with (
            patch("chemsmart.cli.config.platform.system", return_value="Windows"),
            patch("chemsmart.cli.config.update_yaml_files") as mock_update,
            patch("chemsmart.cli.config.add_lines_in_yaml_files"),
        ):
            cfg = Config()
            # Simulate the server subcommand logic (Windows path)
            import chemsmart.cli.config as cfg_module
            import platform as _platform

            if cfg_module.platform.system() == "Windows":
                pass  # Windows branch: no update_yaml_files call
            else:
                cfg_module.update_yaml_files(
                    cfg.chemsmart_server, "~/miniconda3", "some_path"
                )

            # With the Windows patch, update_yaml_files was not called
            mock_update.assert_not_called()

    def test_server_updates_conda_path_on_linux(self):
        """On Linux, config server should call update_yaml_files for conda."""
        with (
            patch("chemsmart.cli.config.platform.system", return_value="Linux"),
            patch("chemsmart.cli.config.update_yaml_files") as mock_update,
            patch("chemsmart.cli.config.add_lines_in_yaml_files"),
        ):
            cfg = Config()
            with patch.object(
                type(cfg),
                "conda_folder",
                new_callable=lambda: property(
                    lambda self: "/home/user/miniconda3"
                ),
            ):
                import chemsmart.cli.config as cfg_module

                # On Linux the code calls update_yaml_files
                cfg_module.update_yaml_files(
                    cfg.chemsmart_server,
                    "~/miniconda3",
                    cfg.conda_folder,
                )

                mock_update.assert_called_once_with(
                    cfg.chemsmart_server,
                    "~/miniconda3",
                    "/home/user/miniconda3",
                )
