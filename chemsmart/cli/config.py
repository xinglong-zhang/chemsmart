import os
import shutil
from doctest import debug
from pathlib import Path

import click
import logging
import yaml
from chemsmart.utils.logger import create_logger

logger = logging.getLogger(__name__)

create_logger(debug=True, stream=True)


class Config:

    @property
    def chemsmart_template(self):
        """Define the source path for the .chemsmart templates."""
        return Path(__file__).resolve().parent / ".." / "settings" / "templates" / ".chemsmart"

    @property
    def chemsmart_dest(self):
        """Define the destination path for the .chemsmart configuration."""
        return Path.home() / ".chemsmart"

    @property
    def shell_config(self):
        """Define the shell configuration file path."""
        return Path.home() / (
            ".bashrc" if os.environ.get("SHELL", "").endswith("bash") else ".zshrc"
        )

    @property
    def chemsmart_path(self):
        """Define the path for the chemsmart package."""
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
    def env_vars(self):
        """Define the environment variables to be added to the shell config."""
        return [
            f'export PATH="{self.chemsmart_path}:$PATH"',
            f'export PATH="{self.chemsmart_path}/chemsmart/cli:$PATH"',
            f'export PATH="{self.chemsmart_path}/chemsmart/scripts:$PATH"',
            f'export PYTHONPATH="{self.chemsmart_path}:$PYTHONPATH"',
        ]

    def setup_environment(self):
        """Set up configuration files and environment variables."""
        # Copy templates to ~/.chemsmart
        if not self.chemsmart_dest.exists():
            shutil.copytree(self.chemsmart_template, self.chemsmart_dest)
            logger.info(f"Copied templates to {self.chemsmart_dest}")
        else:
            logger.info(f"Config directory already exists: {self.chemsmart_dest}")

        # Prevent duplicate additions
        with self.shell_config.open("r+") as f:
            lines = f.readlines()
            if not any("Added by chemsmart installer" in line for line in lines):
                f.write("\n# Added by chemsmart installer\n")
                for var in self.env_vars:
                    f.write(f"{var}\n")
                f.write("\n")
                logger.info(f"Updated shell config: {self.shell_config}")
            else:
                logger.info(f"Shell config already updated: {self.shell_config}")

        logger.info(
            "Please restart your terminal or run 'source ~/.bashrc' (or 'source ~/.zshrc')."
        )

        # Update PATH and PYTHONPATH in all yaml files after the line
        # EXTRA_COMMANDS: |
        # # extra commands to activate chemsmart environment in submission script


    # Update PATH and PYTHONPATH in all yaml files after the line
    # EXTRA_COMMANDS: |
    # # extra commands to activate chemsmart environment in submission script



def update_yaml_files(target_directory, value_in_file, user_value):
    """Update YAML files in ~/.chemsmart/server to replace
    value_in_file with the provided user_value."""
    target_dir = Path.home() / ".chemsmart" / target_directory
    if not target_dir.exists():
        logger.info(f"Server directory not found: {target_dir}")
        return

    for yaml_file in target_dir.glob("*.yaml"):
        logger.info(f"Processing YAML file: {yaml_file}")
        with yaml_file.open("r") as f:
            try:
                yaml_content = yaml.safe_load(f)  # Load YAML content
            except yaml.YAMLError as e:
                logger.info(f"Error reading {yaml_file}: {e}")
                continue

        # Check and replace any occurrences of '~/bin/g16'
        updated = False
        if yaml_content:
            yaml_str = str(yaml_content)
            if value_in_file in yaml_str:
                yaml_str = yaml_str.replace(value_in_file, user_value)
                updated = True

        # Write updated content back to the file
        if updated:
            with yaml_file.open("w") as f:
                f.write(yaml_str)
            logger.info(
                f"Updated {yaml_file} to replace {value_in_file} with '{user_value}'"
            )
        else:
            logger.info(f"No changes made to {yaml_file}")

def add_lines_in_yaml_files(target_directory, lines_in_positions, lines_to_add):
    """Add lines (lines_to_add) to at specific positions after given lines
    (lines_in_positions) in all yaml files in a target directory."""
    if not target_directory.exists():
        logger.info(f"Target directory not found: {target_directory}")
        return

    for yaml_file in target_directory.glob("*.yaml"):
        logger.info(f"Processing YAML file: {yaml_file}")
        try:
            with yaml_file.open("r") as f:
                yaml_content = f.readlines()  # Read file line by line

            updated_content = []
            skip_addition = False  # To ensure lines are not added multiple times

            for line in yaml_content:
                updated_content.append(line)

                # Check if the current line matches any line in `lines_in_positions`
                if any(pos_line in line for pos_line in lines_in_positions):
                    if not skip_addition:  # Avoid duplicate additions
                        for new_line in lines_to_add:
                            updated_content.append(new_line + "\n")
                        skip_addition = True  # Set to true after adding lines

            # Write back the updated content
            with yaml_file.open("w") as f:
                f.writelines(updated_content)
            logger.info(f"Updated YAML file: {yaml_file}")

        except yaml.YAMLError as e:
            logger.info(f"Error reading {yaml_file}: {e}")
        except Exception as e:
            logger.error(f"Unexpected error while processing {yaml_file}: {e}")


@click.group(name="config", invoke_without_command=True)
@click.pass_context
def config(ctx):
    cfg = Config()
    ctx.ensure_object(dict)  # Initialize the Click context object if not already initialized
    ctx.obj["cfg"] = cfg
    """Set up configuration files and environment variables."""
    if ctx.invoked_subcommand is None:
        # Run the default environment setup when no subcommand is provided
        cfg.setup_environment()
        ctx["env_vars"] = cfg.env_vars

@config.command()
@click.pass_context
@click.option(
    "-f",
    "--folder",
    type=click.Path(exists=True),
    required=True,
    help="Path to the Gaussian g16 folder.",
)
def server(ctx):
    """Configures server settings in ~/.chemsmart/server/*yaml files.

    Add conda env vars after the lines
    EXTRA_COMMANDS: |
    # extra commands to activate chemsmart environment in submission script
    in the *yaml file.

    Examples:
        chemsmart config server
    """
    cfg = ctx.obj["cfg"]
    logger.info(f"Configuring servers in ~/.chemsmart/server/*yaml files.")
    add_lines_in_yaml_files(cfg.chemsmart_server, ["# extra commands to activate chemsmart environment in submission script"], ctx["env_vars"])


@config.command()
@click.pass_context
@click.option(
    "-f",
    "--folder",
    type=click.Path(exists=True),
    required=True,
    help="Path to the Gaussian g16 folder.",
)
def gaussian(ctx, folder):
    """Configures paths to g16 folder.

    Replaces '~/bin/g16' with the specified folder in YAML files.

    Examples:
        chemsmart config gaussian --folder <G16FOLDER>
    """
    cfg = ctx.obj["cfg"]
    logger.info(f"Configuring Gaussian with folder: {folder}")
    update_yaml_files(cfg.chemsmart_server, "~/bin/g16", folder)


@config.command()
@click.pass_context
@click.option(
    "-f",
    "--folder",
    type=click.Path(exists=True),
    required=True,
    help="Path to the ORCA folder.",
)
def orca(ctx, folder):
    """Configures paths to g16 folder.

    Replaces '~/bin/g16' with the specified folder in YAML files.

    Examples:
        chemsmart config gaussian --folder <G16FOLDER>
    """
    cfg = ctx.obj["cfg"]
    logger.info(f"Configuring Gaussian with folder: {folder}")
    update_yaml_files(cfg.chemsmart_server, "~/bin/orca_6_0_0", folder)

config.add_command(server)
config.add_command(gaussian)
config.add_command(orca)

if __name__ == "__main__":
    config()
