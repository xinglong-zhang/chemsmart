import logging
import subprocess
import tempfile
from pathlib import Path

import click
import tomlkit

from chemsmart.utils.logger import create_logger

# Initialize logger
logger = logging.getLogger(__name__)
create_logger(debug=True, stream=True)


class Updater:
    """Handles dependency updates in pyproject.toml."""

    def __init__(self):
        self._package_path = Path(__file__).resolve().parent.parent.parent
        self._pyproject_path = self._package_path / "pyproject.toml"

    @property
    def package_path(self) -> Path:
        logger.debug(f"Package path: {self._package_path}")
        return self._package_path

    @property
    def pyproject_path(self) -> Path:
        return self._pyproject_path

    def update_pyproject_toml(self):
        """Update dependencies in pyproject.toml."""
        requirements_path = self._generate_requirements()
        if requirements_path:
            self._update_toml(requirements_path)

    def _generate_requirements(self, ignore_dirs=None) -> Path:
        """Run pipreqs and extract dependencies."""
        ignore_dirs = ignore_dirs or [".git", ".github"]
        ignore_arg = ",".join(ignore_dirs)

        with tempfile.NamedTemporaryFile(delete=False) as tmp_file:
            tmp_path = Path(tmp_file.name)
            logger.debug(f"Temporary requirements file: {tmp_path}")

        cmd = [
            "pipreqs",
            str(self.package_path),
            "--force",
            "--ignore",
            ignore_arg,
            "--savepath",
            str(tmp_path),
        ]
        logger.info(f"Running pipreqs: {' '.join(cmd)}")

        process = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        stdout, stderr = process.communicate()

        if process.returncode != 0:
            logger.error(f"Error running pipreqs: {stderr.decode()}")
            return None
        logger.info(stdout.decode())
        return tmp_path

    def _get_existing_dependencies(self):
        """Extract dependencies from pyproject.toml."""
        if not self.pyproject_path.exists():
            logger.error("pyproject.toml not found.")
            return set()

        with self.pyproject_path.open("r") as f:
            pyproject_toml = tomlkit.parse(f.read())

        return set(pyproject_toml["project"]["dependencies"])

    def _get_missing_dependencies(self, requirements_path):
        """Identify dependencies missing from pyproject.toml, ignoring version numbers.
        Only considers dependencies found in the codebase, not the entire environment.
        """
        if not requirements_path.exists():
            logger.error("Generated requirements file not found.")
            return set()

        # Read detected dependencies from pipreqs output
        with requirements_path.open("r") as f:
            detected_deps = {line.strip() for line in f if line.strip()}

        # Extract existing dependencies from pyproject.toml
        existing_deps = self._get_existing_dependencies()

        # Function to extract only package names
        def extract_pkg_name(dep):
            # return re.split(r"[=<>~!]", dep, maxsplit=1)[0].strip()  # ignore version
            return dep

        # Convert dependencies to package names only
        detected_pkgs = {extract_pkg_name(dep) for dep in detected_deps}
        logger.debug(f"Detected packages: {detected_pkgs}")
        existing_pkgs = {extract_pkg_name(dep) for dep in existing_deps}

        # Identify missing packages (only from code, not the entire environment)
        missing_packages = detected_pkgs - existing_pkgs
        return missing_packages

    def _update_toml(self, requirements_path):
        """Update pyproject.toml with missing dependencies."""
        missing_dependencies = self._get_missing_dependencies(
            requirements_path
        )
        if not missing_dependencies:
            logger.info("No new dependencies to add.")
            return

        with self.pyproject_path.open("r") as f:
            pyproject_toml = tomlkit.parse(f.read())

        for dep in missing_dependencies:
            pyproject_toml["project"]["dependencies"].append(dep)

        with self.pyproject_path.open("w") as f:
            f.write(tomlkit.dumps(pyproject_toml))

        logger.info(
            f"Added missing dependencies: {', '.join(missing_dependencies)}"
        )


@click.group(name="update")
@click.pass_context
def update(ctx):
    """Manage updates in the chemsmart package."""
    ctx.ensure_object(dict)
    ctx.obj["updater"] = Updater()


@update.command()
@click.pass_context
def deps(ctx):
    """Automatically update dependencies in pyproject.toml."""
    logger.info("Updating dependencies...")
    ctx.obj["updater"].update_pyproject_toml()
    logger.info("Update complete.")


if __name__ == "__main__":
    update()
