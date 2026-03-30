"""Tests for chemsmart.cli.config.Config class."""

import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

from chemsmart.cli.config import Config
from chemsmart.utils.io import update_powershell_profiles, update_windows_env


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
        """On native Windows (no POSIX shell), shell_config returns None."""
        cfg = Config()
        env_without_shell = {
            k: v for k, v in __import__("os").environ.items() if k != "SHELL"
        }
        with (
            patch(
                "chemsmart.cli.config.platform.system", return_value="Windows"
            ),
            patch.dict("os.environ", env_without_shell, clear=True),
        ):
            assert cfg.shell_config is None

    def test_shell_config_git_bash_returns_bashrc(self, tmp_path):
        """On Windows Git Bash (SHELL is set), shell_config returns ~/.bashrc
        even when the file does not yet exist."""
        cfg = Config()
        fake_home = tmp_path
        # Do NOT pre-create .bashrc — the property should return it regardless.
        with (
            patch(
                "chemsmart.cli.config.platform.system", return_value="Windows"
            ),
            patch.dict("os.environ", {"SHELL": "/usr/bin/bash"}),
            patch.object(Path, "home", return_value=fake_home),
        ):
            result = cfg.shell_config
        assert result == fake_home / ".bashrc"

    def test_shell_config_git_bash_exe_returns_bashrc(self, tmp_path):
        """SHELL=/usr/bin/bash.exe (some Windows Git Bash installs) also maps
        to ~/.bashrc."""
        cfg = Config()
        fake_home = tmp_path
        with (
            patch(
                "chemsmart.cli.config.platform.system", return_value="Windows"
            ),
            patch.dict("os.environ", {"SHELL": "/usr/bin/bash.exe"}),
            patch.object(Path, "home", return_value=fake_home),
        ):
            result = cfg.shell_config
        assert result == fake_home / ".bashrc"

    def test_shell_config_bash(self, tmp_path):
        """On Linux/macOS with bash, shell_config returns ~/.bashrc
        (creates it if absent via update_shell_config)."""
        cfg = Config()
        fake_home = tmp_path
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
        """On native Windows (no POSIX shell), env_vars returns []."""
        cfg = Config()
        env_without_shell = {
            k: v for k, v in __import__("os").environ.items() if k != "SHELL"
        }
        with (
            patch(
                "chemsmart.cli.config.platform.system", return_value="Windows"
            ),
            patch.dict("os.environ", env_without_shell, clear=True),
        ):
            assert cfg.env_vars == []

    def test_env_vars_git_bash_non_empty(self):
        """On Windows Git Bash (SHELL is set), env_vars returns Unix exports."""
        cfg = Config()
        with (
            patch(
                "chemsmart.cli.config.platform.system", return_value="Windows"
            ),
            patch.dict("os.environ", {"SHELL": "/usr/bin/bash"}),
        ):
            vars_ = cfg.env_vars
        assert len(vars_) > 0
        assert all(v.startswith("export ") for v in vars_)

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

    def test_setup_environment_windows_calls_update_env(self, tmp_path):
        """On native Windows (no POSIX shell, no PS), setup_environment
        calls _update_windows_env."""
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
            patch.object(
                type(cfg),
                "powershell_profiles",
                new_callable=lambda: property(lambda self: []),
            ),
            patch.object(cfg, "_update_windows_env") as mock_update_env,
        ):
            cfg.setup_environment()
        assert dest.exists()
        mock_update_env.assert_called_once()

    def test_setup_environment_powershell_writes_ps_profile(self, tmp_path):
        """On Windows PowerShell (PSModulePath set), setup_environment
        writes to PS profiles."""
        dest = tmp_path / ".chemsmart"
        ps_profile = (
            tmp_path
            / "Documents"
            / "WindowsPowerShell"
            / "Microsoft.PowerShell_profile.ps1"
        )
        cfg = Config()
        env_without_shell = {
            k: v for k, v in __import__("os").environ.items() if k != "SHELL"
        }
        env_without_shell["PSModulePath"] = "C:\\Windows\\system32"
        with (
            patch.object(
                type(cfg),
                "chemsmart_dest",
                new_callable=lambda: property(lambda self: dest),
            ),
            patch(
                "chemsmart.cli.config.platform.system", return_value="Windows"
            ),
            patch.dict("os.environ", env_without_shell, clear=True),
            patch.object(
                type(cfg),
                "powershell_profiles",
                new_callable=lambda: property(lambda self: [ps_profile]),
            ),
            patch.object(cfg, "_update_windows_env") as mock_update_env,
        ):
            cfg.setup_environment()
        assert dest.exists()
        mock_update_env.assert_not_called()
        assert ps_profile.exists()
        content = ps_profile.read_text()
        assert "chemsmart initialize" in content
        assert "Set-Alias -Name chemsmart -Value chemsmart.exe" in content

    def test_powershell_profiles_returns_empty_on_linux(self):
        """powershell_profiles returns [] on non-Windows platforms."""
        cfg = Config()
        with patch(
            "chemsmart.cli.config.platform.system", return_value="Linux"
        ):
            assert cfg.powershell_profiles == []

    def test_powershell_profiles_returns_empty_without_psmodulepath(self):
        """powershell_profiles returns [] when PSModulePath is not set."""
        cfg = Config()
        # Override PSModulePath to "" (falsy) rather than clearing the whole
        # environment — clearing is unreliable on Windows CI runners that run
        # inside a PowerShell session.
        with (
            patch(
                "chemsmart.cli.config.platform.system", return_value="Windows"
            ),
            patch.dict("os.environ", {"PSModulePath": ""}, clear=False),
        ):
            assert cfg.powershell_profiles == []

    def test_powershell_profiles_returns_paths_with_psmodulepath(
        self, tmp_path
    ):
        """powershell_profiles returns profile paths when
        PSModulePath is set on Windows."""
        cfg = Config()
        with (
            patch(
                "chemsmart.cli.config.platform.system", return_value="Windows"
            ),
            patch.dict(
                "os.environ", {"PSModulePath": "C:\\Windows"}, clear=False
            ),
            patch.object(Path, "home", return_value=tmp_path),
        ):
            profiles = cfg.powershell_profiles
        assert len(profiles) == 2
        assert any("WindowsPowerShell" in str(p) for p in profiles)
        assert any(
            "PowerShell" in str(p) and "Windows" not in str(p)
            for p in profiles
        )

    def test_update_windows_env_adds_paths(self, tmp_path):
        """_update_windows_env should write chemsmart paths into the registry."""
        cfg = Config()
        fake_pkg_path = str(tmp_path / "chemsmart")

        mock_key = MagicMock()
        # Simulate empty existing PATH and PYTHONPATH
        mock_key.__enter__ = lambda self: self
        mock_key.__exit__ = MagicMock(return_value=False)
        mock_key.QueryValueEx = MagicMock(side_effect=FileNotFoundError)

        captured_set_calls = {}

        def fake_set_value_ex(key, name, reserved, reg_type, value):
            captured_set_calls[name] = value

        mock_winreg = MagicMock()
        mock_winreg.HKEY_CURRENT_USER = 0x80000001
        mock_winreg.KEY_READ = 0x20019
        mock_winreg.KEY_WRITE = 0x20006
        mock_winreg.REG_EXPAND_SZ = 2
        mock_winreg.REG_SZ = 1
        mock_winreg.OpenKey.return_value = mock_key
        mock_winreg.QueryValueEx.side_effect = FileNotFoundError
        mock_winreg.SetValueEx.side_effect = fake_set_value_ex

        mock_ctypes = MagicMock()

        with (
            patch.object(
                type(cfg),
                "chemsmart_package_path",
                new_callable=lambda: property(lambda self: fake_pkg_path),
            ),
            patch.dict(
                sys.modules, {"winreg": mock_winreg, "ctypes": mock_ctypes}
            ),
        ):
            cfg._update_windows_env()

        # SetValueEx should have been called at least for PATH
        mock_winreg.SetValueEx.assert_called()

    def test_update_windows_env_handles_missing_winreg(self):
        """_update_windows_env must not raise when winreg is unavailable."""
        cfg = Config()
        with patch.dict(sys.modules, {"winreg": None}):
            # Should log a warning and return gracefully, not raise
            cfg._update_windows_env()


class TestWindowsUpdateEnvUtil:
    """Tests for the standalone chemsmart.utils.utils.windows_update_env function."""

    def test_update_windows_env_util_adds_paths(self, tmp_path):
        """windows_update_env writes directories to the registry PATH."""
        fake_pkg = str(tmp_path / "pkg")
        paths = [fake_pkg, str(tmp_path / "cli")]

        mock_key = MagicMock()
        mock_key.__enter__ = lambda self: self
        mock_key.__exit__ = MagicMock(return_value=False)

        mock_winreg = MagicMock()
        mock_winreg.HKEY_CURRENT_USER = 0x80000001
        mock_winreg.KEY_READ = 0x20019
        mock_winreg.KEY_WRITE = 0x20006
        mock_winreg.REG_EXPAND_SZ = 2
        mock_winreg.REG_SZ = 1
        mock_winreg.OpenKey.return_value = mock_key
        mock_winreg.QueryValueEx.side_effect = FileNotFoundError

        mock_ctypes = MagicMock()

        with patch.dict(
            sys.modules, {"winreg": mock_winreg, "ctypes": mock_ctypes}
        ):
            update_windows_env(paths, fake_pkg)

        mock_winreg.SetValueEx.assert_called()

    def test_update_windows_env_util_handles_missing_winreg(self):
        """windows_update_env must not raise when winreg is unavailable."""
        with patch.dict(sys.modules, {"winreg": None}):
            update_windows_env(["/some/path"], "/some/path")


class TestUpdatePowershellProfiles:
    """Tests for chemsmart.utils.io.update_powershell_profiles."""

    def test_creates_profile_with_alias(self, tmp_path):
        """Creates a new profile containing the Set-Alias declaration."""
        ps_profile = (
            tmp_path
            / "Documents"
            / "WindowsPowerShell"
            / "Microsoft.PowerShell_profile.ps1"
        )
        update_powershell_profiles(
            [ps_profile],
            ["Set-Alias -Name chemsmart -Value chemsmart.exe"],
        )
        assert ps_profile.exists()
        content = ps_profile.read_text()
        assert "chemsmart initialize" in content
        assert "Set-Alias -Name chemsmart -Value chemsmart.exe" in content

    def test_replaces_existing_block(self, tmp_path):
        """Re-running update_powershell_profiles replaces the old block."""
        ps_profile = tmp_path / "profile.ps1"
        ps_profile.write_text(
            "# >>> chemsmart initialize >>>\n"
            "function chemsmart { python -m chemsmart.cli.main $args }\n"
            "# <<< chemsmart initialize <<<\n",
            encoding="utf-8",
        )
        update_powershell_profiles(
            [ps_profile],
            ["Set-Alias -Name chemsmart -Value chemsmart.exe"],
        )
        content = ps_profile.read_text()
        assert "function chemsmart" not in content
        assert "Set-Alias -Name chemsmart -Value chemsmart.exe" in content
        assert (
            content.count("chemsmart initialize") == 2
        )  # BEGIN + END markers

    def test_migrates_legacy_profile_entries(self, tmp_path):
        """Old-style '# Added by chemsmart installer' blocks are replaced."""
        ps_profile = tmp_path / "profile.ps1"
        ps_profile.write_text(
            "# some existing content\n"
            "\n# Added by chemsmart installer\n"
            '$env:PATH = "C:\\old\\cli;$env:PATH"\n'
            "\n",
            encoding="utf-8",
        )
        update_powershell_profiles(
            [ps_profile],
            ["Set-Alias -Name chemsmart -Value chemsmart.exe"],
        )
        content = ps_profile.read_text()
        # Legacy PATH entry must be removed
        assert "C:\\old\\cli" not in content
        # New alias must be present
        assert "Set-Alias -Name chemsmart -Value chemsmart.exe" in content
        # Original non-chemsmart content must be preserved
        assert "some existing content" in content


class TestConfigSetupEnvironment:
    """Tests for Config.setup_environment Unix shell config update."""

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

    def test_server_skips_conda_update_when_conda_not_found(
        self, invoke_config_server
    ):
        """When conda is not in PATH, update_yaml_files is not called."""
        with (
            patch("chemsmart.cli.config.update_yaml_files") as mock_update,
            patch("chemsmart.cli.config.add_lines_in_yaml_files"),
            patch.object(
                Config,
                "chemsmart_dest",
                new_callable=lambda: property(
                    lambda self: Path("/fake/.chemsmart")
                ),
            ),
            patch.object(
                Config,
                "conda_folder",
                new_callable=lambda: property(
                    lambda self: (_ for _ in ()).throw(
                        FileNotFoundError("Conda not found in PATH.")
                    )
                ),
            ),
        ):
            result = invoke_config_server()
            assert result.exit_code == 0
            mock_update.assert_not_called()

    def test_server_updates_conda_path_on_linux(self, invoke_config_server):
        """On Linux, update_yaml_files uses auto-detected conda."""
        with (
            patch(
                "chemsmart.cli.config.platform.system", return_value="Linux"
            ),
            patch("chemsmart.cli.config.update_yaml_files") as mock_update,
            patch("chemsmart.cli.config.add_lines_in_yaml_files"),
            patch.object(
                Config,
                "chemsmart_dest",
                new_callable=lambda: property(
                    lambda self: Path("/fake/.chemsmart")
                ),
            ),
            patch.object(
                Config,
                "conda_folder",
                new_callable=lambda: property(
                    lambda self: "/home/user/miniconda3"
                ),
            ),
        ):
            result = invoke_config_server()
            assert result.exit_code == 0
            mock_update.assert_called_once()
            call_args = mock_update.call_args[0]
            assert call_args[1] == "~/miniconda3"
            assert call_args[2] == "/home/user/miniconda3"

    def test_server_git_bash_auto_detects_conda(self, invoke_config_server):
        """On Windows Git Bash (SHELL set), conda is auto-detected and YAML files updated."""
        with (
            patch(
                "chemsmart.cli.config.platform.system", return_value="Windows"
            ),
            patch.dict("os.environ", {"SHELL": "/usr/bin/bash"}),
            patch("chemsmart.cli.config.update_yaml_files") as mock_update,
            patch("chemsmart.cli.config.add_lines_in_yaml_files"),
            patch.object(
                Config,
                "chemsmart_dest",
                new_callable=lambda: property(
                    lambda self: Path("/fake/.chemsmart")
                ),
            ),
            patch.object(
                Config,
                "conda_folder",
                new_callable=lambda: property(
                    lambda self: "/c/Users/user/miniconda3"
                ),
            ),
        ):
            result = invoke_config_server()
            assert result.exit_code == 0
            mock_update.assert_called_once()
            call_args = mock_update.call_args[0]
            assert call_args[1] == "~/miniconda3"
            assert call_args[2] == "/c/Users/user/miniconda3"

    def test_server_powershell_auto_detects_conda(self, invoke_config_server):
        """On Windows Conda PowerShell (PSModulePath set), conda is auto-detected."""
        with (
            patch(
                "chemsmart.cli.config.platform.system", return_value="Windows"
            ),
            patch.dict(
                "os.environ",
                {"PSModulePath": "C:\\Windows\\system32"},
                clear=False,
            ),
            patch("chemsmart.cli.config.update_yaml_files") as mock_update,
            patch("chemsmart.cli.config.add_lines_in_yaml_files"),
            patch.object(
                Config,
                "chemsmart_dest",
                new_callable=lambda: property(
                    lambda self: Path("/fake/.chemsmart")
                ),
            ),
            patch.object(
                Config,
                "conda_folder",
                new_callable=lambda: property(
                    lambda self: "C:\\Users\\user\\miniconda3"
                ),
            ),
        ):
            result = invoke_config_server()
            assert result.exit_code == 0
            mock_update.assert_called_once()
            call_args = mock_update.call_args[0]
            assert call_args[1] == "~/miniconda3"
            assert call_args[2] == "C:\\Users\\user\\miniconda3"


class TestConfigureCondaInServerYaml:
    """Tests for Config.configure_conda_in_server_yaml."""

    def test_calls_update_yaml_when_conda_found(self, tmp_path):
        """When conda is found, update_yaml_files is called with the detected folder."""
        cfg = Config()
        with (
            patch.object(
                type(cfg),
                "conda_folder",
                new_callable=lambda: property(
                    lambda self: "/home/user/miniconda3"
                ),
            ),
            patch.object(
                type(cfg),
                "chemsmart_server",
                new_callable=lambda: property(lambda self: tmp_path),
            ),
            patch("chemsmart.cli.config.update_yaml_files") as mock_update,
        ):
            cfg.configure_conda_in_server_yaml()
        mock_update.assert_called_once_with(
            tmp_path, "~/miniconda3", "/home/user/miniconda3"
        )

    def test_logs_message_when_conda_not_found(self, tmp_path):
        """When conda is not in PATH, configure_conda_in_server_yaml logs a message."""
        cfg = Config()
        with (
            patch.object(
                type(cfg),
                "conda_folder",
                new_callable=lambda: property(
                    lambda self: (_ for _ in ()).throw(
                        FileNotFoundError("Conda not found")
                    )
                ),
            ),
            patch("chemsmart.cli.config.update_yaml_files") as mock_update,
        ):
            cfg.configure_conda_in_server_yaml()
        mock_update.assert_not_called()


class TestConfigurePathsInteractively:
    """Tests for Config.configure_paths_interactively."""

    def test_updates_yaml_for_provided_paths(self, tmp_path):
        """When the user provides all three paths, update_yaml_files is called for each."""
        cfg = Config()
        with (
            patch.object(
                type(cfg),
                "chemsmart_server",
                new_callable=lambda: property(lambda self: tmp_path),
            ),
            patch("chemsmart.cli.config.update_yaml_files") as mock_update,
            patch("chemsmart.cli.config.click.prompt") as mock_prompt,
        ):
            mock_prompt.side_effect = [
                "/path/to/g16",
                "/path/to/orca",
                "/path/to/nci",
            ]
            cfg.configure_paths_interactively()

        assert mock_update.call_count == 3
        mock_update.assert_any_call(tmp_path, "~/bin/g16", "/path/to/g16")
        mock_update.assert_any_call(
            tmp_path, "~/bin/orca_6_0_0", "/path/to/orca"
        )
        mock_update.assert_any_call(tmp_path, "~/bin/nciplot", "/path/to/nci")

    def test_skips_when_all_prompts_empty(self, tmp_path):
        """When the user presses Enter for all prompts, no YAML updates are made."""
        cfg = Config()
        with (
            patch.object(
                type(cfg),
                "chemsmart_server",
                new_callable=lambda: property(lambda self: tmp_path),
            ),
            patch("chemsmart.cli.config.update_yaml_files") as mock_update,
            patch("chemsmart.cli.config.click.prompt", return_value=""),
        ):
            cfg.configure_paths_interactively()

        mock_update.assert_not_called()

    def test_partial_paths_updates_only_provided(self, tmp_path):
        """Only paths that the user fills in are updated; skipped ones are ignored."""
        cfg = Config()
        with (
            patch.object(
                type(cfg),
                "chemsmart_server",
                new_callable=lambda: property(lambda self: tmp_path),
            ),
            patch("chemsmart.cli.config.update_yaml_files") as mock_update,
            patch("chemsmart.cli.config.click.prompt") as mock_prompt,
        ):
            # Provide Gaussian, skip ORCA and NCIPLOT
            mock_prompt.side_effect = ["/opt/g16", "", ""]
            cfg.configure_paths_interactively()

        assert mock_update.call_count == 1
        mock_update.assert_called_once_with(tmp_path, "~/bin/g16", "/opt/g16")

    def test_handles_eof_gracefully(self, tmp_path):
        """When stdin raises EOFError (non-interactive), all paths are skipped."""
        cfg = Config()
        with (
            patch.object(
                type(cfg),
                "chemsmart_server",
                new_callable=lambda: property(lambda self: tmp_path),
            ),
            patch("chemsmart.cli.config.update_yaml_files") as mock_update,
            patch("chemsmart.cli.config.click.prompt", side_effect=EOFError),
        ):
            cfg.configure_paths_interactively()

        mock_update.assert_not_called()
