import logging
import os
import platform
import shutil
from importlib import resources
from pathlib import Path

import click
import yaml

from chemsmart.utils.io import (
    update_powershell_profiles,
    update_shell_config,
    update_windows_env,
)
from chemsmart.utils.logger import create_logger

logger = logging.getLogger(__name__)

create_logger(debug=True, stream=True)


class Config:
    """
    Central configuration helper for chemsmart.

    The class is organised into five logical sections:

    1. **Template & destination paths** — where bundled templates live and
       where they are copied to (``~/.chemsmart``).
    2. **Conda detection** — locate the active conda installation via
       ``shutil.which`` (works on Linux, macOS, Git Bash and Conda
       PowerShell).
    3. **POSIX shell management** — detect the active shell (bash / zsh / sh)
       and return the correct rc file path and ``export`` lines.
    4. **Windows PowerShell management** — detect an Anaconda / Miniconda
       PowerShell session and write ``$env:PATH`` entries into the PS profile.
    5. **Windows registry management** — fallback for plain CMD sessions that
       have neither a POSIX shell nor PowerShell; persists PATH via the
       Windows registry.
    """

    # ------------------------------------------------------------------ #
    # 1. Template & destination paths                                       #
    # ------------------------------------------------------------------ #

    @property
    def chemsmart_template(self):
        """
        Return a Traversable pointing to the bundled ``.chemsmart`` template
        directory.

        Uses ``importlib.resources`` so that the path is resolved correctly
        after installation on Windows, macOS, and Linux alike.
        """
        return (
            resources.files("chemsmart.settings") / "templates" / ".chemsmart"
        )

    @property
    def chemsmart_dest(self):
        """Destination path for the user's ``.chemsmart`` configuration."""
        return Path.home() / ".chemsmart"

    @property
    def chemsmart_server(self):
        """Path to the ``server`` sub-directory of the user config."""
        return self.chemsmart_dest / "server"

    @property
    def chemsmart_gaussian(self):
        """Path to the ``gaussian`` sub-directory of the user config."""
        return self.chemsmart_dest / "gaussian"

    @property
    def chemsmart_orca(self):
        """Path to the ``orca`` sub-directory of the user config."""
        return self.chemsmart_dest / "orca"

    @property
    def chemsmart_package_path(self):
        """Absolute path to the root of the chemsmart source/install tree."""
        chemsmart_path = Path(__file__).resolve().parent / ".." / ".."
        chemsmart_path = os.path.abspath(chemsmart_path)
        logger.debug(f"chemsmart package path: {chemsmart_path}")
        return chemsmart_path

    # ------------------------------------------------------------------ #
    # 2. Conda detection                                                    #
    # ------------------------------------------------------------------ #

    @property
    def conda_path(self):
        """
        Locate the ``conda`` executable via ``shutil.which``.

        Works on Linux, macOS, Windows Git Bash, and Anaconda / Miniconda
        PowerShell Prompt (all of which add conda to ``PATH``).

        Raises :class:`FileNotFoundError` if conda is not found.
        """
        conda_path = shutil.which("conda")
        if conda_path is None or not os.path.exists(conda_path):
            raise FileNotFoundError(
                "Conda not found in PATH. "
                "Ensure conda is installed and its 'Scripts' / 'bin' "
                "directory is on PATH."
            )
        return conda_path

    @property
    def conda_folder(self):
        """Root directory of the conda installation (parent of the directory containing the executable)."""
        return os.path.dirname(os.path.dirname(self.conda_path))

    def configure_conda_in_server_yaml(self) -> None:
        """
        Auto-detect the conda installation and update the placeholder
        ``~/miniconda3`` path in all server YAML files.

        Conda is located via :attr:`conda_folder`, which uses
        ``shutil.which("conda")`` and therefore works on all supported
        platforms — Linux, macOS, Windows Git Bash and Anaconda / Miniconda
        PowerShell Prompt.

        If conda cannot be found a helpful message is logged instead of
        raising an exception.
        """
        try:
            update_yaml_files(
                self.chemsmart_server, "~/miniconda3", self.conda_folder
            )
        except FileNotFoundError:
            logger.info(
                "Conda not found in PATH. "
                "Add conda to your PATH and re-run 'chemsmart config server', "
                "or update the conda path in ~/.chemsmart/server/*.yaml manually."
            )

    # ------------------------------------------------------------------ #
    # 3. POSIX shell management (bash / zsh / sh)                          #
    # ------------------------------------------------------------------ #

    @property
    def shell_config(self):
        """
        Return the shell startup file path for the active POSIX shell.

        Returns ``None`` on native Windows when no POSIX shell is active
        (i.e. the ``SHELL`` environment variable is not set).  On Windows
        Git Bash / MSYS2 the ``SHELL`` variable *is* set, so the shell rc
        file is managed exactly as on Linux / macOS.
        """
        if platform.system() == "Windows" and not os.environ.get("SHELL"):
            return None

        shell = os.environ.get("SHELL", "")
        # Use Path.stem to strip any .exe extension (e.g. /usr/bin/bash.exe
        # on some Windows Git Bash installations).
        shell_name = Path(shell).stem
        if shell_name == "bash":
            # Always target ~/.bashrc; update_shell_config creates it if absent.
            return Path.home() / ".bashrc"
        if shell_name == "zsh":
            return Path.home() / ".zshrc"
        return Path.home() / ".profile"

    @property
    def env_vars(self):
        """
        Unix-style ``export`` lines to append to the active shell rc file.

        Returns ``[]`` on native Windows (no POSIX shell active).  On
        Windows Git Bash / MSYS2 the Unix-style exports are returned so
        that the shell rc file is updated correctly.
        """
        if platform.system() == "Windows" and not os.environ.get("SHELL"):
            return []

        return [
            f'export PATH="{self.chemsmart_package_path}:$PATH"',
            f'export PATH="{self.chemsmart_package_path}/chemsmart/cli:$PATH"',
            f'export PATH="{self.chemsmart_package_path}/chemsmart/scripts:'
            f'$PATH"',
            f'export PYTHONPATH="{self.chemsmart_package_path}:$PYTHONPATH"',
        ]

    def _update_shell_config(self, shell_file: Path) -> None:
        """
        Append chemsmart ``export`` lines to *shell_file* (idempotent).

        Delegates to :func:`chemsmart.utils.io.update_shell_config`.
        """
        update_shell_config(shell_file, self.env_vars)

    # ------------------------------------------------------------------ #
    # 4. Windows PowerShell management                                      #
    # ------------------------------------------------------------------ #

    @property
    def powershell_profiles(self):
        """
        Return profile paths to update when running inside a PowerShell
        session (e.g. Anaconda / Miniconda PowerShell Prompt).

        PowerShell is detected via the ``PSModulePath`` environment variable,
        which is always set by PowerShell regardless of how it was launched.

        Both Windows PowerShell 5.x and PowerShell 7+ profile paths are
        returned so the update works across both versions.

        Returns ``[]`` on non-Windows platforms or outside a PowerShell
        session.
        """
        if platform.system() != "Windows":
            return []
        if not os.environ.get("PSModulePath"):
            return []
        home = Path.home()
        return [
            # Windows PowerShell 5.x (Anaconda/Miniconda default)
            home
            / "Documents"
            / "WindowsPowerShell"
            / "Microsoft.PowerShell_profile.ps1",
            # PowerShell 7+
            home
            / "Documents"
            / "PowerShell"
            / "Microsoft.PowerShell_profile.ps1",
        ]

    @property
    def ps_env_vars(self):
        """
        PowerShell ``$env:PATH`` / ``$env:PYTHONPATH`` assignment lines.

        These are the PowerShell equivalent of the Unix ``export`` lines
        written to ``~/.bashrc``.
        """
        pkg_path = str(self.chemsmart_package_path)
        cli_path = str(Path(self.chemsmart_package_path) / "chemsmart" / "cli")
        scripts_path = str(
            Path(self.chemsmart_package_path) / "chemsmart" / "scripts"
        )
        return [
            f'$env:PATH = "{pkg_path};$env:PATH"',
            f'$env:PATH = "{cli_path};$env:PATH"',
            f'$env:PATH = "{scripts_path};$env:PATH"',
            f'$env:PYTHONPATH = "{pkg_path};$env:PYTHONPATH"',
        ]

    def _update_powershell_profiles(self, profiles) -> None:
        """
        Append chemsmart ``$env:PATH`` lines to each PowerShell profile
        (idempotent).

        Delegates to :func:`chemsmart.utils.io.update_powershell_profiles`.
        """
        update_powershell_profiles(profiles, self.ps_env_vars)

    # ------------------------------------------------------------------ #
    # 5. Windows registry management (CMD / plain Windows fallback)         #
    # ------------------------------------------------------------------ #

    def _update_windows_env(self) -> None:
        """
        Add chemsmart paths to the Windows user PATH and PYTHONPATH via
        the registry.

        This is the fallback for plain CMD sessions that have neither a
        POSIX shell (Git Bash) nor an active PowerShell session.
        Delegates to :func:`chemsmart.utils.io.update_windows_env`.
        """
        pkg_path = str(self.chemsmart_package_path)
        paths_to_add = [
            pkg_path,
            str(Path(self.chemsmart_package_path) / "chemsmart" / "cli"),
            str(Path(self.chemsmart_package_path) / "chemsmart" / "scripts"),
        ]
        update_windows_env(paths_to_add, pkg_path)

    # ------------------------------------------------------------------ #
    # High-level orchestration                                              #
    # ------------------------------------------------------------------ #

    def setup_environment(self):
        """
        Copy bundled templates to ``~/.chemsmart`` and register the
        ``chemsmart`` command in the active shell environment.

        Dispatch logic:

        * **POSIX shell** (Linux, macOS, Git Bash on Windows): append
          ``export`` lines to the shell rc file (``~/.bashrc``,
          ``~/.zshrc``, …).
        * **Windows PowerShell** (Anaconda / Miniconda PS Prompt, detected
          via ``PSModulePath``): append ``$env:PATH`` lines to the PS
          profile.
        * **Windows CMD / other**: update PATH / PYTHONPATH in the Windows
          user registry.
        """
        # -- Copy templates -----------------------------------------------
        if not self.chemsmart_dest.exists():
            with resources.as_file(self.chemsmart_template) as src_dir:
                shutil.copytree(src_dir, self.chemsmart_dest)
            logger.info(f"Copied templates to {self.chemsmart_dest}")
        else:
            logger.info(
                f"Config directory already exists: {self.chemsmart_dest}"
            )

        # -- Register chemsmart in the active shell environment -----------
        shell_file = self.shell_config
        if shell_file is not None:
            # POSIX shell (Linux / macOS / Git Bash on Windows)
            self._update_shell_config(shell_file)
            return

        ps_profiles = self.powershell_profiles
        if ps_profiles:
            # Anaconda / Miniconda PowerShell Prompt
            self._update_powershell_profiles(ps_profiles)
        else:
            # Plain CMD or other Windows environment — use the registry
            self._update_windows_env()
            logger.info(
                "Environment variables updated in the Windows registry.\n"
                "Please restart your terminal for the changes to take "
                "effect.\n"
                "To refresh PATH in the current PowerShell session, run:\n"
                "  $env:PATH = "
                "[System.Environment]::GetEnvironmentVariable("
                "'PATH', 'Machine') + ';' + "
                "[System.Environment]::GetEnvironmentVariable("
                "'PATH', 'User')"
            )


def update_yaml_files(target_directory, value_in_file, user_value):
    """
    Update YAML files in ~/.chemsmart/server to replace value_in_file
    with the provided user_value.
    """
    target_dir = Path.home() / ".chemsmart" / target_directory
    if not target_dir.exists():
        logger.info(f"Server directory not found: {target_dir}")
        return

    for yaml_file in target_dir.glob("*.yaml"):
        logger.info(f"Processing YAML file: {yaml_file}")
        new_lines = []
        with open(yaml_file, "r") as f:
            lines = f.readlines()
            for line in lines:
                if value_in_file in line:
                    updated_line = line.replace(value_in_file, user_value)
                    new_lines.append(updated_line)
                else:
                    new_lines.append(line)

        with open(yaml_file, "w") as g:
            g.writelines(new_lines)
        logger.info(f"Updated YAML file: {yaml_file}")


def add_lines_in_yaml_files(
    target_directory, lines_in_positions, lines_to_add, prepend_string=""
):
    """
    Add lines (lines_to_add) to at specific positions after given lines
    (lines_in_positions) in all yaml files in a target directory.
    """
    if not target_directory.exists():
        logger.info(f"Target directory not found: {target_directory}")
        return

    for yaml_file in target_directory.glob("*.yaml"):
        logger.info(f"Processing YAML file: {yaml_file}")
        try:
            with open(yaml_file, "r") as f:
                yaml_content = f.readlines()  # Read file line by line

            updated_content = []
            skip_addition = (
                False  # To ensure lines are not added multiple times
            )

            for line in yaml_content:
                updated_content.append(line)

                # Check if the current line matches any line in
                # `lines_in_positions`
                if any(pos_line in line for pos_line in lines_in_positions):
                    if not skip_addition:  # Avoid duplicate additions
                        for new_line in lines_to_add:
                            updated_content.append(
                                prepend_string + new_line + "\n"
                            )
                        skip_addition = True  # Set to true after adding lines

            # Write back the updated content
            with open(yaml_file, "w") as f:
                for line in updated_content:
                    f.write(line)
                # f.writelines(updated_content)
            logger.info(f"Updated YAML file: {yaml_file}")

        except yaml.YAMLError as e:
            logger.info(f"Error reading {yaml_file}: {e}")
        except Exception as e:
            logger.error(f"Unexpected error while processing {yaml_file}: {e}")


@click.group(name="config", invoke_without_command=True)
@click.pass_context
def config(ctx):
    """Set up configuration files and environment variables."""
    cfg = Config()
    ctx.ensure_object(
        dict
    )  # Initialize the Click context object if not already initialized
    ctx.obj["cfg"] = cfg
    if ctx.invoked_subcommand is None:
        # Run the default environment setup when no subcommand is provided
        cfg.setup_environment()


@config.command()
@click.pass_context
def server(ctx):
    """
    Configure server settings in ~/.chemsmart/server/*.yaml files.

    Adds conda environment variables after the lines:
    EXTRA_COMMANDS: |
    # extra commands to activate chemsmart environment in submission script
    in the *.yaml file.

    Conda is auto-detected via ``which conda`` (works on Linux, macOS,
    Windows Git Bash and Anaconda / Miniconda PowerShell Prompt).  If
    conda is not found a helpful message is logged.

    Examples:
        chemsmart config server
    """
    cfg = ctx.obj["cfg"]
    logger.info("Configuring servers in ~/.chemsmart/server/*yaml files.")
    add_lines_in_yaml_files(
        cfg.chemsmart_server,
        [
            "#extra commands to activate chemsmart environment in "
            "submission script"
        ],
        cfg.env_vars,
        prepend_string=" " * 8,
    )

    # Auto-detect conda and update server YAML files.
    cfg.configure_conda_in_server_yaml()


@config.command()
@click.pass_context
@click.option(
    "-f",
    "--folder",
    type=str,
    required=True,
    help="Path to the Gaussian g16 folder.",
)
def gaussian(ctx, folder):
    """
    Configure paths to the g16 folder.

    Replaces '~/bin/g16' with the specified folder in YAML files.

    Examples:
        chemsmart config gaussian --folder <G16FOLDER>
    """
    cfg = ctx.obj["cfg"]
    if "~" in folder:
        g16_folder = os.path.expanduser(folder)
        assert os.path.exists(
            os.path.abspath(g16_folder)
        ), f"Gaussian folder not found: {g16_folder}"
    logger.info(f"Configuring Gaussian with folder: {folder}")
    update_yaml_files(cfg.chemsmart_server, "~/bin/g16", folder)


@config.command()
@click.pass_context
@click.option(
    "-f",
    "--folder",
    type=str,
    required=True,
    help="Path to the ORCA folder.",
)
def orca(ctx, folder):
    """
    Configure paths to the ORCA folder.

    Replaces '~/bin/orca' with the specified folder in YAML files.

    Examples:
        chemsmart config orca --folder <ORCAFOLDER>
    """
    cfg = ctx.obj["cfg"]
    if "~" in folder:
        orca_folder = os.path.expanduser(folder)
        assert os.path.exists(
            os.path.abspath(orca_folder)
        ), f"ORCA folder not found: {orca_folder}"
    logger.info(f"Configuring ORCA with folder: {folder}")
    update_yaml_files(cfg.chemsmart_server, "~/bin/orca_6_0_0", folder)


@config.command()
@click.pass_context
@click.option(
    "-f",
    "--folder",
    type=str,
    required=True,
    help="Path to the NCIPLOT folder.",
)
def nciplot(ctx, folder):
    """
    Configure paths to the NCIPLOT folder.

    Replaces '~/bin/nciplot' with the specified folder in YAML files.

    Examples:
        chemsmart config nciplot --folder <NCIPLOTFOLDER>
    """
    cfg = ctx.obj["cfg"]
    if "~" in folder:
        nciplot_folder = os.path.expanduser(folder)
        assert os.path.exists(
            os.path.abspath(nciplot_folder)
        ), f"NCIPLOT folder not found: {nciplot_folder}"
    logger.info(f"Configuring NCIPLOT with folder: {folder}")
    update_yaml_files(cfg.chemsmart_server, "~/bin/nciplot", folder)


config.add_command(server)
config.add_command(gaussian)
config.add_command(orca)
config.add_command(nciplot)

if __name__ == "__main__":
    config()
