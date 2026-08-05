"""
Direct unit tests for the ``chemsmart.cli.main`` entry point group,
beyond the ``--version`` flag already covered in ``test_version.py``.

Focuses on the ``verbose``/banner-printing branches, which are
otherwise unexercised since every other test invokes CLI groups
directly (``gaussian``, ``orca``, ...) rather than through
``entry_point``.
"""

from click.testing import CliRunner

from chemsmart.cli.main import entry_point


class TestMainEntryPoint:
    def test_group_callback_runs_and_prints_banner(self):
        # ``--help`` alone short-circuits before the group callback runs,
        # so invoke a real subcommand (itself also just showing --help)
        # to actually execute entry_point()'s body, including the
        # ASCII banner logging.
        runner = CliRunner()
        result = runner.invoke(entry_point, ["update", "--help"])
        assert result.exit_code == 0
        # ASCII-art banner lines (no literal "CHEMSMART" text appears).
        assert "____ _   _ _____ __  __ ____" in result.output

    def test_registers_expected_subcommands(self):
        runner = CliRunner()
        result = runner.invoke(entry_point, ["--help"])
        assert result.exit_code == 0
        for subcommand in ["run", "sub", "config", "update"]:
            assert subcommand in result.output
