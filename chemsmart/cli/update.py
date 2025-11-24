import logging
import re
import subprocess
import tempfile
from pathlib import Path

import click
import tomlkit
import yaml

from chemsmart.utils.logger import create_logger

# Initialize logger
logger = logging.getLogger(__name__)
create_logger(debug=True, stream=True)


class Updater:
    """
    Handles dependency updates in pyproject.toml.
    """

    def __init__(self):
        self._package_path = Path(__file__).resolve().parent.parent.parent
        self._pyproject_path = self._package_path / "pyproject.toml"
        self._version_file_path = self.package_path / "chemsmart" / "VERSION"
        self._docs_conf_file_path = (
            self.package_path / "docs" / "source" / "conf.py"
        )

    @property
    def package_path(self) -> Path:
        logger.debug(f"Package path: {self._package_path}")
        return self._package_path

    @property
    def pyproject_path(self) -> Path:
        return self._pyproject_path

    @property
    def version_file_path(self) -> Path:
        return self._version_file_path

    @property
    def docs_conf_file_path(self) -> Path:
        return self._docs_conf_file_path

    def update_pyproject_toml(self):
        """
        Update dependencies in pyproject.toml.
        """
        requirements_path = self._generate_requirements()
        if requirements_path:
            self._update_toml(requirements_path)

    def _generate_requirements(self, ignore_dirs=None) -> Path:
        """
        Run pipreqs and extract dependencies.
        """
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
        """
        Extract dependencies from pyproject.toml and environment.yml.
        """
        dependencies = set()
        dependencies.update(self._get_existing_dependencies_from_pyproject())
        dependencies.update(
            self._get_existing_dependencies_from_environment_yml()
        )
        return dependencies

    def _get_existing_dependencies_from_pyproject(self):
        """
        Extract dependencies from pyproject.toml.
        """
        if not self.pyproject_path.exists():
            logger.error("pyproject.toml not found.")
            return set()

        with self.pyproject_path.open("r") as f:
            pyproject_toml = tomlkit.parse(f.read())

        # Safely handle missing 'project' or 'dependencies' keys
        return set(pyproject_toml.get("project", {}).get("dependencies", []))

    def _get_existing_dependencies_from_environment_yml(self):
        """
        Extract dependencies from environment.yml.
        """
        env_yml_path = Path("environment.yml")
        if not env_yml_path.exists():
            logger.info(
                "environment.yml not found, skipping Conda dependencies."
            )
            return set()

        with env_yml_path.open("r") as f:
            env_data = yaml.safe_load(f)

        dependencies = set()
        conda_deps = env_data.get("dependencies", [])
        for dep in conda_deps:
            if isinstance(dep, str):
                dependencies.add(dep)
            elif isinstance(dep, dict) and "pip" in dep:
                dependencies.update(dep["pip"])

        return dependencies

    def _get_missing_dependencies(self, requirements_path):
        """
        Identify dependencies missing from pyproject.toml, preserving
        versions.
        """
        if not requirements_path.exists():
            logger.error("Generated requirements file not found.")
            return set()

        # Read detected dependencies from pipreqs output (with versions)
        with requirements_path.open("r") as f:
            detected_deps = {line.strip() for line in f if line.strip()}

        # Extract existing dependencies
        existing_deps = self._get_existing_dependencies()

        # Package name normalization mapping (for comparison)
        package_mapping = {
            "pymol-open-source": "pymol",  # Treat pymol-open-source as
            # pymol for matching
        }
        reverse_mapping = {
            v: k for k, v in package_mapping.items()
        }  # For display

        # Function to extract package name (for comparison)
        def extract_pkg_name(dep):
            pkg_name = re.split(r"[=<>~!]", dep, maxsplit=1)[0].lower().strip()
            return package_mapping.get(pkg_name, pkg_name)

        # Normalize package names for comparison
        detected_pkgs = {extract_pkg_name(dep) for dep in detected_deps}
        existing_pkgs = {extract_pkg_name(dep) for dep in existing_deps}

        # Create display version of detected_pkgs with reverse mapping
        detected_pkgs_display = {
            reverse_mapping.get(pkg, pkg) for pkg in detected_pkgs
        }
        logger.debug(f"Detected packages: {detected_pkgs_display}")
        logger.debug(f"Existing packages: {existing_pkgs}")

        # Identify missing dependencies with their full strings
        missing_dependencies = {
            dep
            for dep in detected_deps
            if extract_pkg_name(dep) not in existing_pkgs
        }
        return missing_dependencies

    def _update_toml(self, requirements_path):
        """
        Update pyproject.toml with missing dependencies, including versions,
        in lowercase.
        """
        missing_dependencies = self._get_missing_dependencies(
            requirements_path
        )
        if not missing_dependencies:
            logger.info("No new dependencies to add.")
            return

        with self.pyproject_path.open("r") as f:
            pyproject_toml = tomlkit.parse(f.read())

        # Function to normalize package name to lowercase while preserving
        # version
        def normalize_dep(dep):
            pkg_name, *version = re.split(r"([=<>~!])", dep, maxsplit=1)
            if version:  # If there's a version specifier
                return pkg_name.lower() + "".join(version)
            return pkg_name.lower()  # No version, just lowercase the name

        # Adjust and normalize dependencies
        adjusted_deps = {
            normalize_dep(dep).replace("==", ">=")
            for dep in missing_dependencies
        }
        for dep in adjusted_deps:
            pyproject_toml["project"]["dependencies"].append(dep)

        with self.pyproject_path.open("w") as f:
            f.write(tomlkit.dumps(pyproject_toml))

        logger.info(f"Added missing dependencies: {', '.join(adjusted_deps)}")

    def update_version_number(self, version_number):
        """Update version number in chemsmart/VERSION, pyproject.toml,
        and docs/source/conf.py"""
        from chemsmart.utils.repattern import release_pattern, version_pattern

        logger.info(f"Updating to version number: {version_number}")

        # 1) Update chemsmart/VERSION
        if not self.version_file_path.exists():
            logger.error(
                f"Version file not found at {self.version_file_path}."
            )
        else:
            logger.debug(
                f"Writing version {version_number} to {self.version_file_path}"
            )
            self.version_file_path.write_text(
                f"{version_number}\n", encoding="utf-8"
            )

        # 2) Update pyproject.toml
        # version = "0.1.9"
        if not self.pyproject_path.exists():
            logger.error(f"pyproject.toml not found at: {self.pyproject_path}")
        else:
            logger.debug(f"Updating pyproject.toml at {self.pyproject_path}")
            pyproject_text = self.pyproject_path.read_text(encoding="utf-8")

            # Replace: version = "0.1.9" (preserve spacing and comments)
            version_replacement = r"\g<1>" + version_number + r"\g<2>"

            new_pyproject_text, n_subs = re.subn(
                version_pattern,
                version_replacement,
                pyproject_text,
                count=1,
            )

            if n_subs == 0:
                logger.warning(
                    "No 'version = \"...\"' assignment found in pyproject.toml; "
                    "project version not updated."
                )
            else:
                self.pyproject_path.write_text(
                    new_pyproject_text, encoding="utf-8"
                )
                logger.info(
                    f"Updated project version to {version_number} "
                    f"in {self.pyproject_path}"
                )

        # 3) Update docs/source/conf.py: release = "x.y.z"
        if not self.docs_conf_file_path.exists():
            logger.warning(
                f"Sphinx conf.py not found at: {self.docs_conf_file_path} "
                "- skipping docs version update."
            )
        else:
            logger.debug(
                f"Updating Sphinx release in {self.docs_conf_file_path}"
            )
            text = self.docs_conf_file_path.read_text(encoding="utf-8")

            # Replace: release = "0.1.9"
            replacement = r"\g<1>" + version_number + r"\g<2>"

            new_text, n_subs = re.subn(
                release_pattern, replacement, text, count=1, flags=re.MULTILINE
            )

            if n_subs == 0:
                logger.warning(
                    "No 'release = \"...\"' assignment found in conf.py; "
                    "Sphinx version not updated."
                )
            else:
                self.docs_conf_file_path.write_text(new_text, encoding="utf-8")
                logger.info(
                    f"Updated Sphinx release to {version_number} "
                    f"in {self.docs_conf_file_path}"
                )

        logger.info(f"Version update to {version_number} completed.")


@click.group(name="update")
@click.pass_context
def update(ctx):
    """
    Manage updates in the chemsmart package.
    """
    ctx.ensure_object(dict)
    ctx.obj["updater"] = Updater()


@update.command()
@click.pass_context
def deps(ctx):
    """
    Automatically update dependencies in pyproject.toml.
    """
    logger.info("Updating dependencies...")
    ctx.obj["updater"].update_pyproject_toml()
    logger.info("Update complete.")


@update.command()
@click.option(
    "-v",
    "--version-number",
    type=str,
    required=True,
    help="Version number to be updated to.",
)
@click.pass_context
def version(ctx, version_number: str):
    """
    Automatically update chemsmart version in chemsmart/VERSION,
    pyproject.toml and docs/source/conf.py.
    """
    logger.info("Updating dependencies...")
    ctx.obj["updater"].update_version_number(version_number)
    logger.info("Update complete.")


if __name__ == "__main__":
    update()
