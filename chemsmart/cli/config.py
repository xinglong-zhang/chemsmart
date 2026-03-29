import logging
import os
import platform
import shutil
from importlib import resources
from pathlib import Path

import click
import yaml

from chemsmart.utils.io import windows_update_env
from chemsmart.utils.logger import create_logger

logger = logging.getLogger(__name__)

create_logger(debug=True, stream=True)


class Config:

    @property
    def chemsmart_template(self):
        """
        Return a Traversable pointing to the bundled .chemsmart template
        directory.

        Uses importlib.resources so that the path is resolved correctly
        after installation on Windows, macOS, and Linux alike.
        """
        template_dir = (
            resources.files("chemsmart.settings") / "templates" / ".chemsmart"
        )
        return template_dir

    @property
    def chemsmart_dest(self):
        """
        Define the destination path for the .chemsmart configuration.
        """
        return Path.home() / ".chemsmart"

    @property
    def shell_config(self):
        """
        Return the shell startup file on Unix-like systems.

        Returns None on native Windows (no POSIX shell active), where shell rc
        files are not managed here.  On Windows Git Bash / MSYS2 the ``SHELL``
        environment variable is set, so the shell rc file is managed like on
        Linux/macOS.
        """
        # On native Windows (no POSIX shell), skip shell rc management.
        if platform.system() == "Windows" and not os.environ.get("SHELL"):
            return None

        shell = os.environ.get("SHELL", "")
        if shell.endswith("bash"):
            bashrc_filepath = Path.home() / ".bashrc"
            if bashrc_filepath.exists():
                return bashrc_filepath
            else:
                logger.info("Bashrc not found. Trying bash_profile.")
                return Path.home() / ".bash_profile"
        if shell.endswith("zsh"):
            return Path.home() / ".zshrc"

        return Path.home() / ".profile"

    @property
    def chemsmart_package_path(self):
        """
        Define the path for the chemsmart package.
        """
        chemsmart_path = Path(__file__).resolve().parent / ".." / ".."
        chemsmart_path = os.path.abspath(chemsmart_path)
        logger.debug(f"chemsmart package path: {chemsmart_path}")
        return chemsmart_path

    @property
    def chemsmart_server(self):
        """Define the path for the chemsmart server configuration."""
        return self.chemsmart_dest / "server"

    @property
    def chemsmart_gaussian(self):
        """Define the path for the chemsmart gaussian configuration."""
        return self.chemsmart_dest / "gaussian"

    @property
    def chemsmart_orca(self):
        """Define the path for the chemsmart orca configuration."""
        return self.chemsmart_dest / "orca"

    @property
    def conda_path(self):
        """
        Define the path for the conda environment.
        """
        conda_path = shutil.which("conda")
        if conda_path is None or not os.path.exists(conda_path):
            raise FileNotFoundError(
                "Conda not found in PATH. Please install conda!"
            )
        return conda_path

    @property
    def conda_folder(self):
        """
        Define the path to the conda folder.
        """
        # Go up 2 directories from the conda path
        return os.path.dirname(os.path.dirname(self.conda_path))

    @property
    def env_vars(self):
        """
        Define the environment variables to be added to the shell config.

        Returns an empty list on native Windows (no POSIX shell), where shell
        rc files are not managed by chemsmart.  On Windows Git Bash / MSYS2
        (``SHELL`` env var is set) the Unix-style exports are returned so that
        the shell rc file is updated correctly.
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

    @property
    def powershell_profiles(self):
        """
        Return a list of PowerShell profile paths to update when running in a
        Windows PowerShell session (e.g. Anaconda / Miniconda PowerShell
        prompt).

        PowerShell is detected via the ``PSModulePath`` environment variable,
        which is always set by PowerShell regardless of how the session was
        launched.

        Both Windows PowerShell 5.x
        (``~/Documents/WindowsPowerShell/Microsoft.PowerShell_profile.ps1``)
        and PowerShell 7+
        (``~/Documents/PowerShell/Microsoft.PowerShell_profile.ps1``) profiles
        are returned so that the update works across both versions.

        Returns an empty list on non-Windows platforms or when not running
        inside a PowerShell session.
        """
        if platform.system() != "Windows":
            return []
        if not os.environ.get("PSModulePath"):
            return []
        home = Path.home()
        return [
            # Windows PowerShell 5.x (used by Anaconda/Miniconda Prompt)
            home / "Documents" / "WindowsPowerShell"
            / "Microsoft.PowerShell_profile.ps1",
            # PowerShell 7+ (cross-platform)
            home / "Documents" / "PowerShell"
            / "Microsoft.PowerShell_profile.ps1",
        ]

    @property
    def ps_env_vars(self):
        """
        PowerShell-formatted ``$env:PATH`` / ``$env:PYTHONPATH`` lines to
        append to the PowerShell profile file.

        These are the PowerShell equivalent of the Unix ``export`` lines
        written to ``~/.bashrc``.
        """
        pkg_path = str(self.chemsmart_package_path)
        cli_path = str(
            Path(self.chemsmart_package_path) / "chemsmart" / "cli"
        )
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
        Write chemsmart PATH / PYTHONPATH entries into each PowerShell profile
        file — the PowerShell equivalent of appending ``export`` lines to
        ``~/.bashrc``.

        The update is idempotent: if the marker comment is already present the
        file is not modified again.
        """
        for ps_profile in profiles:
            ps_profile.parent.mkdir(parents=True, exist_ok=True)
            if not ps_profile.exists():
                ps_profile.touch()
            with ps_profile.open("r+", encoding="utf-8") as f:
                lines = f.readlines()
                if not any(
                    "Added by chemsmart installer" in line for line in lines
                ):
                    f.write("\n# Added by chemsmart installer\n")
                    for var in self.ps_env_vars:
                        f.write(f"{var}\n")
                    f.write("\n")
                    logger.info(f"Updated PowerShell profile: {ps_profile}")
                else:
                    logger.info(
                        f"PowerShell profile already updated: {ps_profile}"
                    )
        logger.info(
            "PowerShell profiles updated.\n"
            "To apply changes in the current PowerShell session, run:\n"
            "  . $PROFILE"
        )

    def _windows_update_env(self) -> None:
        """
        Add chemsmart paths to the Windows user PATH and PYTHONPATH via
        the registry — the Windows equivalent of appending ``export``
        lines to ``~/.bashrc``.

        Delegates to :func:`chemsmart.utils.utils.windows_update_env`.
        """
        pkg_path = str(self.chemsmart_package_path)
        paths_to_add = [
            pkg_path,
            str(Path(self.chemsmart_package_path) / "chemsmart" / "cli"),
            str(Path(self.chemsmart_package_path) / "chemsmart" / "scripts"),
        ]
        windows_update_env(paths_to_add, pkg_path)

    def setup_environment(self):
        """
        Set up configuration files and environment variables.
        """
        # Copy templates to ~/.chemsmart using importlib.resources so that
        # the source path is resolved correctly after installation on all
        # platforms.
        if not self.chemsmart_dest.exists():
            with resources.as_file(self.chemsmart_template) as src_dir:
                shutil.copytree(src_dir, self.chemsmart_dest)
            logger.info(f"Copied templates to {self.chemsmart_dest}")
        else:
            logger.info(
                f"Config directory already exists: {self.chemsmart_dest}"
            )

        shell_file = self.shell_config
        if shell_file is None:
            # Not a POSIX shell (no SHELL env var).
            # Check whether we are running inside a PowerShell session
            # (e.g. Anaconda / Miniconda PowerShell Prompt).
            ps_profiles = self.powershell_profiles
            if ps_profiles:
                self._update_powershell_profiles(ps_profiles)
            else:
                # Plain CMD or other Windows environment — update the user
                # environment via the registry.
                self._windows_update_env()
                logger.info(
                    "Environment variables updated in the Windows registry.\n"
                    "Please restart your terminal for the changes to take "
                    "effect.\n"
                    "To refresh PATH in the current PowerShell session, "
                    "run:\n"
                    "  $env:PATH = "
                    "[System.Environment]::GetEnvironmentVariable("
                    "'PATH', 'Machine') + ';' + "
                    "[System.Environment]::GetEnvironmentVariable("
                    "'PATH', 'User')"
                )
            return

        # Create shell config file if it does not exist yet
        if not shell_file.exists():
            shell_file.touch()

        # Prevent duplicate additions
        with shell_file.open("r+", encoding="utf-8") as f:
            lines = f.readlines()
            if not any(
                "Added by chemsmart installer" in line for line in lines
            ):
                f.write("\n# Added by chemsmart installer\n")
                for var in self.env_vars:
                    f.write(f"{var}\n")
                f.write("\n")
                logger.info(f"Updated shell config: {shell_file}")
            else:
                logger.info(f"Shell config already updated: {shell_file}")

        logger.info(
            f"Please restart your terminal or run " f"'source {shell_file}'."
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
@click.option(
    "--conda-path",
    "-cp",
    default=None,
    show_default=True,
    help=(
        "Path to the conda installation on the remote cluster "
        "(e.g. ~/miniconda3 or /opt/conda). "
        "Required on Windows where auto-detection is not possible. "
        "On Linux/macOS this overrides the auto-detected path."
    ),
)
def server(ctx, conda_path):
    """
    Configure server settings in ~/.chemsmart/server/*.yaml files.

    Adds conda environment variables after the lines:
    EXTRA_COMMANDS: |
    # extra commands to activate chemsmart environment in submission script
    in the *.yaml file.

    On Windows, use --conda-path to supply the Unix-style conda path for
    the remote HPC cluster (e.g. ~/miniconda3), since auto-detection would
    return a Windows-style path that is incorrect for remote clusters.

    Examples:
        chemsmart config server
        chemsmart config server --conda-path ~/miniconda3
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

    # Update the conda path in server YAML files.
    # Try auto-detection first (works on Linux/macOS and Windows Git Bash).
    # If conda is not in PATH, log a helpful message instead of failing.
    if conda_path is not None:
        # Explicit override — works on all platforms
        update_yaml_files(cfg.chemsmart_server, "~/miniconda3", conda_path)
    else:
        try:
            update_yaml_files(
                cfg.chemsmart_server, "~/miniconda3", cfg.conda_folder
            )
        except FileNotFoundError:
            logger.info(
                "Conda not found in PATH. To configure the conda path in "
                "server YAML files, run:\n"
                "  chemsmart config server --conda-path <path/to/conda>"
            )


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
