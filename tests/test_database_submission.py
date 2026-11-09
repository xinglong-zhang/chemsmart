from importlib import import_module
from unittest.mock import patch

from click.testing import CliRunner

from chemsmart.cli.main import entry_point
from chemsmart.settings.server import Server
from chemsmart.settings.submitters import SLURMSubmitter


def test_database_assembly_submission_defers_work(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    server = Server(name="dummy")
    with (
        patch.object(Server, "from_servername", return_value=server),
        patch.object(server, "submit") as submit,
        patch.object(
            import_module("chemsmart.cli.database.assemble"), "BaseFolder"
        ) as folder,
    ):
        result = CliRunner().invoke(
            entry_point,
            [
                "sub",
                "-s",
                "dummy",
                "--test",
                "database",
                "assemble",
                "-d",
                ".",
                "-p",
                "xtb",
                "-o",
                "xtb.db",
            ],
        )
    assert result.exit_code == 0, result.output + repr(result.exception)
    folder.assert_not_called()
    submit.assert_called_once()
    args = submit.call_args.kwargs["cli_args"]
    assert str(tmp_path) in args
    assert str(tmp_path / "xtb.db") in args
    assert not (tmp_path / "xtb.db").exists()


def test_database_submission_writes_scripts_without_executable_config(
    tmp_path, monkeypatch, slurm_test_server
):
    from chemsmart.cli.database.assemble import DatabaseAssemblyJob

    monkeypatch.chdir(tmp_path)
    job = DatabaseAssemblyJob("xtb.db", None)
    submitter = SLURMSubmitter(job=job, server=slurm_test_server)
    submitter.write(["database", "assemble", "-d", str(tmp_path), "-p", "xtb"])
    assert (tmp_path / submitter.run_script).exists()
    assert (tmp_path / submitter.submit_script).exists()


def test_assembly_writes_bounded_batches(tmp_path):
    from unittest.mock import MagicMock

    from chemsmart.cli.run import run

    module = import_module("chemsmart.cli.database.assemble")
    folder = MagicMock()
    folder.get_all_output_files_in_current_folder_and_subfolders_by_program.return_value = [
        str(tmp_path / f"{i}.log") for i in range(205)
    ]
    assembler = MagicMock()
    assembler.assemble_data = {"record": "example"}
    database = MagicMock()
    batches = []
    database.insert_records.side_effect = lambda records: (
        batches.append(len(records)) or len(records)
    )
    with (
        patch.object(module, "BaseFolder", return_value=folder),
        patch.object(module, "SingleFileAssembler", return_value=assembler),
        patch.object(module, "Database", return_value=database),
    ):
        result = CliRunner().invoke(
            run,
            ["database", "assemble", "-d", str(tmp_path), "-p", "gaussian"],
        )
    assert result.exit_code == 0, result.output + repr(result.exception)
    assert batches == [100, 100, 5]
