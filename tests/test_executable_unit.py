"""
Direct unit tests for ``chemsmart.settings.executable``: the ``Executable``
base class and its ``GaussianExecutable``/``ORCAExecutable``/
``NCIPLOTExecutable`` subclasses.

``from_servername`` end-to-end parsing against a real server YAML is
already covered by ``test_server.py``; these tests focus on the
constructor-level properties (``scratch_dir``, ``env``,
``available_servers``) and each subclass's ``get_executable``.
"""

import os

from chemsmart.settings.executable import (
    Executable,
    GaussianExecutable,
    NCIPLOTExecutable,
    ORCAExecutable,
)


class TestExecutableScratchDir:
    def test_scratch_dir_none_without_envars(self):
        exe = Executable(executable_folder="/opt/prog")
        assert exe.scratch_dir is None

    def test_scratch_dir_none_when_not_present_in_envars(self):
        exe = Executable(envars="export FOO=bar\n")
        assert exe.scratch_dir is None

    def test_scratch_dir_extracted_from_envars(self):
        exe = Executable(envars="export SCRATCH=/tmp/scratch\n")
        assert exe.scratch_dir == "/tmp/scratch"

    def test_scratch_dir_strips_inline_comment(self):
        exe = Executable(
            envars="export SCRATCH=/tmp/scratch  # local scratch\n"
        )
        assert exe.scratch_dir == "/tmp/scratch"

    def test_scratch_dir_finds_line_among_others(self):
        exe = Executable(
            envars=(
                "export PATH=/opt/bin:$PATH\n"
                "export SCRATCH=/tmp/scratch\n"
                "export OTHER=1\n"
            )
        )
        assert exe.scratch_dir == "/tmp/scratch"


class TestExecutableEnv:
    def test_env_none_without_envars(self):
        exe = Executable()
        assert exe.env is None

    def test_env_parses_export_lines(self):
        exe = Executable(
            envars=(
                "export PATH=/opt/bin:$PATH\n" "export SCRATCH=/tmp/scratch\n"
            )
        )
        env = exe.env
        assert env == {
            "PATH": "/opt/bin:$PATH",
            "SCRATCH": "/tmp/scratch",
        }

    def test_env_ignores_non_export_lines(self):
        exe = Executable(envars=("# a comment line\n" "export FOO=bar\n" "\n"))
        assert exe.env == {"FOO": "bar"}

    def test_env_strips_inline_comment(self):
        exe = Executable(envars="export FOO=bar  # a comment\n")
        assert exe.env == {"FOO": "bar"}


class TestExecutableAvailableServers:
    def test_available_servers_delegates_to_user_settings(self):
        exe = Executable()
        # Just verify it returns whatever the (real) user settings expose,
        # without asserting specific contents (environment-dependent).
        assert isinstance(exe.available_servers, list)


class TestGaussianExecutable:
    def test_get_executable_returns_none_without_folder(self):
        exe = GaussianExecutable()
        assert exe.get_executable() is None

    def test_get_executable_joins_folder_and_g16(self):
        exe = GaussianExecutable(executable_folder="/opt/g16")
        assert exe.get_executable() == os.path.join("/opt/g16", "g16")

    def test_program_identifier(self):
        assert GaussianExecutable.PROGRAM == "GAUSSIAN"


class TestORCAExecutable:
    def test_get_executable_returns_none_without_folder(self):
        exe = ORCAExecutable()
        assert exe.get_executable() is None

    def test_get_executable_joins_folder_and_orca(self):
        exe = ORCAExecutable(executable_folder="/opt/orca")
        assert exe.get_executable() == os.path.join("/opt/orca", "orca")

    def test_program_identifier(self):
        assert ORCAExecutable.PROGRAM == "ORCA"


class TestNCIPLOTExecutable:
    def test_get_executable_returns_none_without_folder(self):
        exe = NCIPLOTExecutable()
        assert exe.get_executable() is None

    def test_get_executable_joins_folder_and_nciplot(self):
        exe = NCIPLOTExecutable(executable_folder="/opt/nciplot")
        assert exe.get_executable() == os.path.join("/opt/nciplot", "nciplot")

    def test_program_identifier(self):
        assert NCIPLOTExecutable.PROGRAM == "NCIPLOT"


class TestFromServername:
    def test_from_servername_accepts_name_without_yaml_suffix(
        self, server_yaml_file
    ):
        # server_yaml_file already ends with .yaml; strip it to exercise
        # the "append .yaml" branch of from_servername.
        assert server_yaml_file.endswith(".yaml")
        bare_name = server_yaml_file[: -len(".yaml")]
        executable = GaussianExecutable.from_servername(bare_name)
        assert executable.executable_folder == os.path.expanduser(
            "~/programs/g16"
        )

    def test_from_servername_accepts_name_with_yaml_suffix(
        self, server_yaml_file
    ):
        executable = GaussianExecutable.from_servername(server_yaml_file)
        assert executable.executable_folder == os.path.expanduser(
            "~/programs/g16"
        )
