import os
import shutil
from pathlib import Path

import click
import yaml


def setup_environment():
    """Set up configuration files and environment variables."""
    # Define source and destination paths
    src_templates = (
        Path(__file__).resolve().parent
        / ".."
        / "settings"
        / "templates"
        / ".chemsmart"
    )
    print(f"Source templates: {Path(__file__).resolve().parent}")

    dest_config = Path.home() / ".chemsmart"

    # Copy templates to ~/.chemsmart
    if not dest_config.exists():
        shutil.copytree(src_templates, dest_config)
        print(f"Copied templates to {dest_config}")
    else:
        print(f"Config directory already exists: {dest_config}")

    # Update PATH and PYTHONPATH in ~/.bashrc or ~/.zshrc
    shell_config = Path.home() / (
        ".bashrc" if os.environ.get("SHELL", "").endswith("bash") else ".zshrc"
    )

    chemsmart_path = Path(__file__).resolve().parent / ".." / ".."
    chemsmart_path = os.path.abspath(chemsmart_path)
    print(f"chemsmart package path: {chemsmart_path}")
    env_vars = [
        f'export PATH="{chemsmart_path}:$PATH"',
        f'export PATH="{chemsmart_path}/chemsmart/cli:$PATH"',
        f'export PATH="{chemsmart_path}/chemsmart/scripts:$PATH"',
        f'export PYTHONPATH="{chemsmart_path}:$PYTHONPATH"',
    ]

    # Prevent duplicate additions
    with shell_config.open("r+") as f:
        lines = f.readlines()
        if not any("Added by chemsmart installer" in line for line in lines):
            f.write("\n# Added by chemsmart installer\n")
            for var in env_vars:
                f.write(f"{var}\n")
            f.write("\n")
            print(f"Updated shell config: {shell_config}")
        else:
            print(f"Shell config already updated: {shell_config}")

    print(
        "Please restart your terminal or run 'source ~/.bashrc' (or 'source ~/.zshrc')."
    )


def update_yaml_files(target_directory, value_in_file, user_value):
    """Update YAML files in ~/.chemsmart/server to replace
    value_in_file with the provided user_value."""
    target_dir = Path.home() / ".chemsmart" / target_directory
    if not target_dir.exists():
        print(f"Server directory not found: {target_dir}")
        return

    for yaml_file in target_dir.glob("*.yaml"):
        print(f"Processing YAML file: {yaml_file}")
        with yaml_file.open("r") as f:
            try:
                yaml_content = yaml.safe_load(f)  # Load YAML content
            except yaml.YAMLError as e:
                print(f"Error reading {yaml_file}: {e}")
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
            print(
                f"Updated {yaml_file} to replace {value_in_file} with '{user_value}'"
            )
        else:
            print(f"No changes made to {yaml_file}")


@click.group(name="config", invoke_without_command=True)
@click.pass_context
def config(ctx):
    """Set up configuration files and environment variables."""
    if ctx.invoked_subcommand is None:
        # Run the default environment setup when no subcommand is provided
        setup_environment()


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
    print(f"Configuring Gaussian with folder: {folder}")
    update_yaml_files("server", "~/bin/g16", folder)


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
    print(f"Configuring Gaussian with folder: {folder}")
    update_yaml_files("server", "~~/bin/orca_6_0_0", folder)


config.add_command(gaussian)
config.add_command(orca)

if __name__ == "__main__":
    config()
