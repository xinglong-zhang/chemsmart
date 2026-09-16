import os

import yaml

from chemsmart.settings.user import CHEMSMARTUserSettings

user_settings = CHEMSMARTUserSettings()


class ChainProjectSettings:
    """YAML chain-project aliases for per-program project names.

    Top-level keys ``crest``, ``xtb``, ``gaussian``, and ``orca`` are aliases
    to existing per-program project names.
    """

    PROGRAMS = ("crest", "xtb", "gaussian", "orca")
    PROJECT_NAME = "chain"

    def __init__(self, aliases, project_name):
        self._aliases = dict(aliases)
        self.PROJECT_NAME = project_name

    def project_for(self, program):
        """Return the per-program project alias for ``program``.

        Args:
            program (str): One of ``crest``, ``xtb``, ``gaussian``, ``orca``.

        Returns:
            str: Project name stem, e.g. ``test`` for
            ``~/.chemsmart/gaussian/test.yaml``.

        Raises:
            ValueError: If ``program`` is unknown or has no alias in this
                chain file.
        """
        if program not in self.PROGRAMS:
            allowed = ", ".join(self.PROGRAMS)
            raise ValueError(
                f"Unknown chain program {program!r}. "
                f"Allowed programs: {allowed}."
            )
        if program not in self._aliases:
            raise ValueError(
                f"Chain project {self.PROJECT_NAME!r} has no "
                f"{program} section."
            )
        return self._aliases[program]

    @classmethod
    def from_project(cls, project):
        # Match Gaussian/xTB precedence: user projects override packaged
        # templates, and missing names fail with a discoverable config path.
        user_project_settings = cls._from_user_project_name(project)
        if user_project_settings is not None:
            return user_project_settings

        packaged_project_settings = cls._from_packaged_project_name(project)
        if packaged_project_settings is not None:
            return packaged_project_settings

        templates_path = os.path.join(os.path.dirname(__file__), "templates")
        raise FileNotFoundError(
            f"No chain project settings implemented for {project}.\n\n"
            f"Place new chain project settings .yaml file in "
            f"{user_settings.user_chain_settings_dir}.\n\n"
            f"Templates for such settings.yaml files are available at "
            f"{templates_path}\n\n "
            f"Currently available projects: "
            f"{user_settings.all_available_chain_projects}"
        )

    @classmethod
    def from_yaml(cls, filename):
        filename = os.path.abspath(filename)
        with open(filename) as handle:
            config = yaml.safe_load(handle)
        if config is None:
            config = {}
        if not isinstance(config, dict):
            raise ValueError(
                f"Chain project settings in {filename} must be a mapping."
            )
        project_name = os.path.basename(filename).removesuffix(".yaml")
        return cls._from_config(config, project_name=project_name)

    @classmethod
    def _from_config(cls, config, project_name):
        unknown_keys = sorted(set(config) - set(cls.PROGRAMS))
        if unknown_keys:
            allowed = ", ".join(cls.PROGRAMS)
            unknown = ", ".join(unknown_keys)
            raise ValueError(
                f"Unknown chain project keys: {unknown}. "
                f"Allowed keys: {allowed}."
            )
        aliases = cls._parse_aliases(config)
        return cls(aliases=aliases, project_name=project_name)

    @classmethod
    def _parse_aliases(cls, config):
        aliases = {}
        for program in cls.PROGRAMS:
            if program not in config:
                continue
            value = config[program]
            if value is None:
                continue
            if not isinstance(value, str) or not value.strip():
                raise ValueError(
                    f"Chain project {program} alias must be a "
                    f"project name string."
                )
            aliases[program] = value.strip()
        return aliases

    @classmethod
    def _from_user_project_name(cls, project_name):
        if project_name is None:
            return None
        project_name_yaml_path = os.path.join(
            CHEMSMARTUserSettings().user_chain_settings_dir,
            f"{project_name}.yaml",
        )
        try:
            return cls.from_yaml(project_name_yaml_path)
        except FileNotFoundError:
            return None

    @classmethod
    def _from_packaged_project_name(cls, project_name):
        if project_name is None:
            return None
        current_file_dir = os.path.dirname(os.path.abspath(__file__))
        project_name_yaml_path = os.path.join(
            current_file_dir,
            "templates",
            ".chemsmart",
            "chain",
            f"{project_name}.yaml",
        )
        try:
            return cls.from_yaml(project_name_yaml_path)
        except FileNotFoundError:
            return None
