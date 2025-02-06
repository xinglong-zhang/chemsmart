import logging
import os
import subprocess
import tempfile
from pathlib import Path

import click
import tomlkit

from chemsmart.utils.logger import create_logger

logger = logging.getLogger(__name__)

create_logger(debug=True, stream=True)


class Updater:

    @property
    def chemsmart_package_path(self):
        """Define the path for the chemsmart package."""
        chemsmart_path = Path(__file__).resolve().parent / ".." / ".."
        chemsmart_path = os.path.abspath(chemsmart_path)
        logger.debug(f"chemsmart package path: {chemsmart_path}")
        return chemsmart_path

    @property
    def pyproject_toml_path(self):
        """Define the path for the pyproject.toml file."""
        return self.chemsmart_package_path / "pyproject.toml"

    def update_pyproject_toml(self):
        """Update the pyproject.toml file."""
        self._update_toml(self._obtain_dependencies())

    def _obtain_dependencies(self, ignore_dirs=None):
        """Obtain dependencies from the pyproject.toml file."""
        # create temporary filepath
        tmp_requirements_path = tempfile.NamedTemporaryFile().name
        if ignore_dirs is None:
            ignore_dirs = [".git", ".github", "tests"]

        ignore_arg = ",".join(ignore_dirs)
        cmd = [
            "pipreqs",
            self.chemsmart_package_path,
            "--force",
            "--ignore",
            ignore_arg,
            "--savepath",
            tmp_requirements_path,
        ]

        process = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        stdout, stderr = process.communicate()

        if process.returncode != 0:
            logger.info(f"Error running pipreqs: {stderr.decode()}")
        else:
            logger.info(stdout.decode())
        return tmp_requirements_path

    def _existing_dependencies(self):
        """Get existing dependencies from the pyproject.toml file."""
        with open(self.pyproject_toml_path, "r") as f:
            pyproject_toml = tomlkit.parse(f.read())
        return pyproject_toml["dependencies"]

    def _missing_dependencies(self, requirements_path):
        """Get missing dependencies from the pyproject.toml file."""
        with open(requirements_path, "r") as f:
            requirements = f.read().split("\n")
        existing_dependencies = self._existing_dependencies()
        missing_dependencies = []
        for req in requirements:
            if req and req not in existing_dependencies:
                missing_dependencies.append(req)
        return missing_dependencies

    def _update_toml(self, requirements_path):
        """Update the pyproject.toml file with missing dependencies
        that are not already present."""
        missing_dependencies = self._missing_dependencies(requirements_path)
        with open(self.pyproject_toml_path, "r") as f:
            pyproject_toml = tomlkit.parse(f.read())
        for dep in missing_dependencies:
            pyproject_toml["tool"]["poetry"]["dependencies"][dep] = "*"
        with open(self.pyproject_toml_path, "w") as f:
            f.write(pyproject_toml.as_string())


@click.group(name="update", invoke_without_command=True)
@click.pass_context
def update(ctx):
    updater = Updater()
    ctx.ensure_object(
        dict
    )  # Initialize the Click context object if not already initialized
    ctx.obj["updater"] = updater
    """Set up configuration files and environment variables."""


@update.command()
@click.pass_context
def deps(ctx):
    """Automatically update dependencies in the chemsmart package.

    Examples:
        chemsmart update deps
    """
    updater = ctx.obj["updater"]
    logger.info("Update dependencies in the chemsmart package.")
    updater.update_pyproject_toml()


update.add_command(deps)

if __name__ == "__main__":
    update()
