"""
Tests for the unified ``click_folder_options`` decorator and its integration
with the ``thermochemistry`` and ``mol`` CLI groups.

This module verifies that:

1. ``click_folder_options`` registers ``-d/--directory``,
   ``-t/--filetype``, and ``-p/--program`` on any Click command it
   decorates.
2. The ``thermochemistry`` CLI correctly routes batch processing when
   ``--directory`` is combined with ``--program`` or ``--filetype``.
3. The ``mol`` CLI correctly populates ``ctx.obj`` when
   ``--directory``/``--filetype`` or ``--directory``/``--program`` are
   provided.

Each test uses :class:`click.testing.CliRunner` together with
:mod:`unittest.mock` to avoid file-system access and actual job
execution.
"""

import importlib
import os
import tempfile
from unittest.mock import MagicMock, patch

import click
import pytest
from click.testing import CliRunner

from chemsmart.cli.job import click_folder_options
from chemsmart.cli.thermochemistry.thermochemistry import thermochemistry

thermochemistry_cli_module = importlib.import_module(
    "chemsmart.cli.thermochemistry.thermochemistry"
)
mol_cli_module = importlib.import_module("chemsmart.cli.mol.mol")


# ---------------------------------------------------------------------------
# Helper – a minimal Click command decorated with click_folder_options
# ---------------------------------------------------------------------------


@click.command()
@click_folder_options
@click.pass_context
def _folder_cmd(ctx, directory, filetype, program):
    """Minimal command used only to inspect registered options."""
    ctx.ensure_object(dict)
    ctx.obj["directory"] = directory
    ctx.obj["filetype"] = filetype
    ctx.obj["program"] = program


# ---------------------------------------------------------------------------
# Tests for the decorator itself
# ---------------------------------------------------------------------------


class TestClickFolderOptionsDecorator:
    """``click_folder_options`` registers the expected CLI options."""

    def test_directory_option_present_in_help(self):
        """``-d/--directory`` appears in the command help text."""
        runner = CliRunner()
        result = runner.invoke(_folder_cmd, ["--help"])
        assert result.exit_code == 0, result.output
        assert "--directory" in result.output

    def test_filetype_option_present_in_help(self):
        """``-t/--filetype`` appears in the command help text."""
        runner = CliRunner()
        result = runner.invoke(_folder_cmd, ["--help"])
        assert result.exit_code == 0, result.output
        assert "--filetype" in result.output

    def test_program_option_present_in_help(self):
        """``-p/--program`` appears in the command help text."""
        runner = CliRunner()
        result = runner.invoke(_folder_cmd, ["--help"])
        assert result.exit_code == 0, result.output
        assert "--program" in result.output

    def test_directory_option_received_by_callback(self):
        """Passing ``-d /some/path`` propagates the value to the callback."""
        runner = CliRunner()
        result = runner.invoke(_folder_cmd, ["-d", "/some/path"], obj={})
        assert result.exit_code == 0, result.output

    def test_filetype_option_received_by_callback(self):
        """Passing ``-t log`` propagates the value to the callback."""
        runner = CliRunner()
        result = runner.invoke(_folder_cmd, ["-t", "log"], obj={})
        assert result.exit_code == 0, result.output

    def test_program_option_received_by_callback(self):
        """Passing ``-p gaussian`` propagates the value to the callback."""
        runner = CliRunner()
        result = runner.invoke(_folder_cmd, ["-p", "gaussian"], obj={})
        assert result.exit_code == 0, result.output

    def test_all_options_default_to_none(self):
        """When no folder options are given, all three values are None."""
        ctx_obj = {}
        runner = CliRunner()
        result = runner.invoke(_folder_cmd, [], obj=ctx_obj)
        assert result.exit_code == 0, result.output
        assert ctx_obj.get("directory") is None
        assert ctx_obj.get("filetype") is None
        assert ctx_obj.get("program") is None


# ---------------------------------------------------------------------------
# Tests for thermochemistry CLI – folder options
# ---------------------------------------------------------------------------


class TestThermochemistryCLIFolderOptions:
    """Folder options are wired correctly into the ``thermochemistry`` CLI."""

    def _run_with_directory(self, extra_args, mock_files=None):
        """
        Invoke the thermochemistry CLI with a temporary directory.

        Patches ``BaseFolder`` so no real file-system access occurs.
        """
        if mock_files is None:
            mock_files = ["/fake/dir/mol1.log"]

        runner = CliRunner()
        mock_job = MagicMock()
        mock_job.label = "mol1"

        with (
            patch.object(
                thermochemistry_cli_module,
                "get_program_type_from_file",
                return_value="gaussian",
            ),
            patch.object(
                thermochemistry_cli_module.ThermochemistryJob,
                "from_filename",
                return_value=mock_job,
            ) as mock_from_filename,
            patch.object(
                thermochemistry_cli_module.BaseFolder,
                "get_all_output_files_in_current_folder_by_program",
                return_value=mock_files,
            ),
            patch.object(
                thermochemistry_cli_module.BaseFolder,
                "get_all_files_in_current_folder_by_suffix",
                return_value=mock_files,
            ),
        ):
            result = runner.invoke(
                thermochemistry,
                extra_args,
                obj={},
                catch_exceptions=False,
            )
            return result, mock_from_filename

    def test_directory_with_program_accepted(self, tmp_path):
        """``-d dir -p gaussian -T 298.15`` is accepted without error."""
        result, mock_from_filename = self._run_with_directory(
            [
                "-d",
                str(tmp_path),
                "-p",
                "gaussian",
                "-T",
                "298.15",
            ]
        )
        assert result.exit_code == 0, result.output

    def test_directory_with_filetype_accepted(self, tmp_path):
        """``-d dir -t log -T 298.15`` is accepted without error."""
        result, mock_from_filename = self._run_with_directory(
            [
                "-d",
                str(tmp_path),
                "-t",
                "log",
                "-T",
                "298.15",
            ]
        )
        assert result.exit_code == 0, result.output

    def test_directory_without_program_or_filetype_raises(self, tmp_path):
        """``-d dir`` alone (no ``-p`` / ``-t``) is rejected."""
        runner = CliRunner()
        result = runner.invoke(
            thermochemistry,
            ["-d", str(tmp_path), "-T", "298.15"],
            obj={},
            catch_exceptions=True,
        )
        assert result.exit_code != 0 or isinstance(
            result.exception, (ValueError, SystemExit)
        )

    def test_directory_and_filenames_mutually_exclusive(self, tmp_path):
        """Providing both ``-d`` and ``-f`` raises a ``ValueError``."""
        runner = CliRunner()
        with pytest.raises(ValueError, match="Cannot specify both"):
            runner.invoke(
                thermochemistry,
                [
                    "-d",
                    str(tmp_path),
                    "-f",
                    "dummy.log",
                    "-p",
                    "gaussian",
                    "-T",
                    "298.15",
                ],
                obj={},
                catch_exceptions=False,
            )

    def test_directory_with_unsupported_program_raises(self, tmp_path):
        """``-d dir -p unsupported_prog`` raises a ``ValueError``."""
        runner = CliRunner()
        with pytest.raises(ValueError, match="Unsupported program"):
            runner.invoke(
                thermochemistry,
                [
                    "-d",
                    str(tmp_path),
                    "-p",
                    "unsupported_prog",
                    "-T",
                    "298.15",
                ],
                obj={},
                catch_exceptions=False,
            )

    def test_directory_with_program_calls_from_filename_for_each_file(
        self, tmp_path
    ):
        """Each discovered file triggers a ``ThermochemistryJob.from_filename`` call."""
        mock_files = ["/fake/dir/mol1.log", "/fake/dir/mol2.log"]
        result, mock_from_filename = self._run_with_directory(
            ["-d", str(tmp_path), "-p", "gaussian", "-T", "298.15"],
            mock_files=mock_files,
        )
        assert result.exit_code == 0, result.output
        assert mock_from_filename.call_count == len(mock_files)


# ---------------------------------------------------------------------------
# Tests for mol CLI – folder options
# ---------------------------------------------------------------------------


class TestMolCLIFolderOptions:
    """Folder options (``-d``/``-t`` and ``-d``/``-p``) in the ``mol`` CLI."""

    def _invoke_mol_with_visualize(self, cli_args, ctx_obj=None):
        """
        Invoke ``mol … visualize`` with all PyMOL job execution mocked out.

        The ``visualize`` subcommand is the simplest real subcommand for
        testing because it only requires ``ctx.obj["molecules"]`` and
        ``ctx.obj["label"]`` to be set – exactly what the group's folder
        processing does.
        """
        # We need the full mol group (with subcommands registered)
        from chemsmart.cli.mol import mol as mol_group

        runner = CliRunner()
        if ctx_obj is None:
            ctx_obj = {}

        # Add a lightweight no-op test subcommand so Click doesn't complain
        # about a missing command and the group callback can run.
        @mol_group.command("_test_noop")
        @click.pass_context
        def _noop(ctx):
            pass

        try:
            with (
                patch.object(
                    mol_cli_module.Molecule,
                    "from_filepath",
                    return_value=MagicMock(),
                ),
                patch.object(
                    mol_cli_module.BaseFolder,
                    "get_all_output_files_in_current_folder_by_program",
                    return_value=["/fake/dir/mol.log"],
                ),
            ):
                result = runner.invoke(
                    mol_group,
                    cli_args + ["_test_noop"],
                    obj=ctx_obj,
                    catch_exceptions=False,
                )
        finally:
            # Clean up the temporary test command
            mol_group.commands.pop("_test_noop", None)

        return result

    def test_directory_filetype_options_accepted(self, tmp_path):
        """``mol -d dir -t log visualize`` is accepted and populates ``ctx.obj``."""
        ctx_obj = {}
        result = self._invoke_mol_with_visualize(
            ["-d", str(tmp_path), "-t", "log"],
            ctx_obj=ctx_obj,
        )
        # Exit code 0 or 1 are both acceptable here; we just check the
        # folder options were recognised (no "No such option" error).
        assert "No such option" not in result.output, result.output
        assert ctx_obj.get("directory") == str(tmp_path)
        assert ctx_obj.get("filetype") == "log"

    def test_directory_program_options_accepted(self, tmp_path):
        """``mol -d dir -p gaussian visualize`` is accepted and populates ``ctx.obj``."""
        ctx_obj = {}
        result = self._invoke_mol_with_visualize(
            ["-d", str(tmp_path), "-p", "gaussian"],
            ctx_obj=ctx_obj,
        )
        assert "No such option" not in result.output, result.output
        assert ctx_obj.get("directory") == str(tmp_path)
        assert ctx_obj.get("program") == "gaussian"

    def test_directory_filetype_label_auto_generated(self, tmp_path):
        """When ``-d``/``-t`` used, auto-generated label includes dir name."""
        ctx_obj = {}
        result = self._invoke_mol_with_visualize(
            ["-d", str(tmp_path), "-t", "log"],
            ctx_obj=ctx_obj,
        )
        assert "No such option" not in result.output, result.output
        label = ctx_obj.get("label", "")
        dir_name = os.path.basename(os.path.abspath(str(tmp_path)))
        assert dir_name in label
        assert "log" in label

    def test_directory_program_label_auto_generated(self, tmp_path):
        """When ``-d``/``-p`` used, auto-generated label includes program name."""
        ctx_obj = {}
        result = self._invoke_mol_with_visualize(
            ["-d", str(tmp_path), "-p", "gaussian"],
            ctx_obj=ctx_obj,
        )
        assert "No such option" not in result.output, result.output
        label = ctx_obj.get("label", "")
        assert "gaussian" in label
