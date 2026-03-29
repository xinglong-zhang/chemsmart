"""Tests for chemsmart.cli.config.Config class."""

from pathlib import Path
from unittest.mock import patch

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
        with patch(
            "chemsmart.cli.config.platform.system", return_value="Windows"
        ):
            assert cfg.shell_config is None

    def test_shell_config_bash(self, tmp_path):
        cfg = Config()
        fake_home = tmp_path
        bashrc = fake_home / ".bashrc"
        bashrc.touch()
        with (
            patch(
                "chemsmart.cli.config.platform.system", return_value="Linux"
            ),
            patch.dict("os.environ", {"SHELL": "/bin/bash"}),
            patch.object(Path, "home", return_value=fake_home),
        ):
            result = cfg.shell_config
        assert result == fake_home / ".bashrc"

    def test_shell_config_zsh(self, tmp_path):
        cfg = Config()
        fake_home = tmp_path
        with (
            patch(
                "chemsmart.cli.config.platform.system", return_value="Darwin"
            ),
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
        with patch(
            "chemsmart.cli.config.platform.system", return_value="Windows"
        ):
            assert cfg.env_vars == []

    def test_env_vars_linux_non_empty(self):
        cfg = Config()
        with patch(
            "chemsmart.cli.config.platform.system", return_value="Linux"
        ):
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
            patch.object(
                type(cfg),
                "chemsmart_dest",
                new_callable=lambda: property(lambda self: dest),
            ),
            patch.object(
                type(cfg),
                "shell_config",
                new_callable=lambda: property(lambda self: None),
            ),
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
            patch.object(
                type(cfg),
                "chemsmart_dest",
                new_callable=lambda: property(lambda self: dest),
            ),
            patch.object(
                type(cfg),
                "shell_config",
                new_callable=lambda: property(lambda self: None),
            ),
        ):
            cfg.setup_environment()
        # Existing content should be untouched
        assert marker.read_text() == "keep me"

    def test_setup_environment_windows_skips_shell_config(self, tmp_path):
        """On Windows, setup_environment must not touch any shell rc file."""
        dest = tmp_path / ".chemsmart"
        cfg = Config()
        with (
            patch.object(
                type(cfg),
                "chemsmart_dest",
                new_callable=lambda: property(lambda self: dest),
            ),
            patch(
                "chemsmart.cli.config.platform.system", return_value="Windows"
            ),
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
            patch.object(
                type(cfg),
                "chemsmart_dest",
                new_callable=lambda: property(lambda self: dest),
            ),
            patch.object(
                type(cfg),
                "shell_config",
                new_callable=lambda: property(lambda self: shell_rc),
            ),
            patch.object(
                type(cfg),
                "env_vars",
                new_callable=lambda: property(
                    lambda self: ['export FOO="bar"']
                ),
            ),
            patch(
                "chemsmart.cli.config.platform.system", return_value="Linux"
            ),
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
            patch.object(
                type(cfg),
                "chemsmart_dest",
                new_callable=lambda: property(lambda self: dest),
            ),
            patch.object(
                type(cfg),
                "shell_config",
                new_callable=lambda: property(lambda self: shell_rc),
            ),
            patch.object(
                type(cfg),
                "env_vars",
                new_callable=lambda: property(
                    lambda self: ['export FOO="bar"']
                ),
            ),
            patch(
                "chemsmart.cli.config.platform.system", return_value="Linux"
            ),
        ):
            cfg.setup_environment()
            cfg.setup_environment()

        content = shell_rc.read_text()
        assert content.count("Added by chemsmart installer") == 1


class TestConfigServerCommand:
    """Tests for the config server Click command behavior."""

    def _invoke_server(self, args=None):
        """Invoke the server subcommand via Click's test runner."""
        from click.testing import CliRunner

        from chemsmart.cli.config import config

        runner = CliRunner()
        return runner.invoke(config, ["server"] + (args or []))

    def test_server_skips_conda_update_on_windows_no_flag(self):
        """On Windows without --conda-path, update_yaml_files is not called."""
        with (
            patch(
                "chemsmart.cli.config.platform.system", return_value="Windows"
            ),
            patch("chemsmart.cli.config.update_yaml_files") as mock_update,
            patch("chemsmart.cli.config.add_lines_in_yaml_files"),
            patch.object(Config, "chemsmart_dest", new_callable=lambda: property(lambda self: Path("/fake/.chemsmart"))),
        ):
            result = self._invoke_server()
            assert result.exit_code == 0
            mock_update.assert_not_called()

    def test_server_uses_conda_path_flag_on_windows(self):
        """On Windows with --conda-path, update_yaml_files is called with provided path."""
        with (
            patch(
                "chemsmart.cli.config.platform.system", return_value="Windows"
            ),
            patch("chemsmart.cli.config.update_yaml_files") as mock_update,
            patch("chemsmart.cli.config.add_lines_in_yaml_files"),
            patch.object(Config, "chemsmart_dest", new_callable=lambda: property(lambda self: Path("/fake/.chemsmart"))),
        ):
            result = self._invoke_server(["--conda-path", "~/miniconda3"])
            assert result.exit_code == 0
            mock_update.assert_called_once()
            call_args = mock_update.call_args[0]
            assert call_args[1] == "~/miniconda3"
            assert call_args[2] == "~/miniconda3"

    def test_server_uses_conda_path_flag_on_linux(self):
        """On Linux with --conda-path, the provided path overrides auto-detection."""
        with (
            patch(
                "chemsmart.cli.config.platform.system", return_value="Linux"
            ),
            patch("chemsmart.cli.config.update_yaml_files") as mock_update,
            patch("chemsmart.cli.config.add_lines_in_yaml_files"),
            patch.object(Config, "chemsmart_dest", new_callable=lambda: property(lambda self: Path("/fake/.chemsmart"))),
        ):
            result = self._invoke_server(["--conda-path", "/opt/conda"])
            assert result.exit_code == 0
            mock_update.assert_called_once()
            call_args = mock_update.call_args[0]
            assert call_args[1] == "~/miniconda3"
            assert call_args[2] == "/opt/conda"

    def test_server_updates_conda_path_on_linux(self):
        """On Linux without --conda-path, update_yaml_files uses auto-detected conda."""
        with (
            patch(
                "chemsmart.cli.config.platform.system", return_value="Linux"
            ),
            patch("chemsmart.cli.config.update_yaml_files") as mock_update,
            patch("chemsmart.cli.config.add_lines_in_yaml_files"),
            patch.object(Config, "chemsmart_dest", new_callable=lambda: property(lambda self: Path("/fake/.chemsmart"))),
            patch.object(
                Config,
                "conda_folder",
                new_callable=lambda: property(
                    lambda self: "/home/user/miniconda3"
                ),
            ),
        ):
            result = self._invoke_server()
            assert result.exit_code == 0
            mock_update.assert_called_once()
            call_args = mock_update.call_args[0]
            assert call_args[1] == "~/miniconda3"
            assert call_args[2] == "/home/user/miniconda3"
