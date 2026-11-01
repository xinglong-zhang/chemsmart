"""CLI invocation and capture fixtures."""

import pytest
from click.testing import CliRunner

from chemsmart.cli.run import run
from chemsmart.cli.sub import sub
from chemsmart.settings.server import Server
from tests.helpers import write_text_file

############ CLI Fixtures ##################


@pytest.fixture()
def write_tmp_file(tmp_path):
    """Return a callable that writes a file into the test's ``tmp_path``.

    Program-agnostic: use it for any synthetic input/output file (ORCA
    ``.inp``/``.out``, Gaussian ``.com``/``.log``, ``.xyz``, ...). The
    returned path is a ``str``, which every chemsmart file reader accepts.

    Usage in tests::

        def test_something(write_tmp_file):
            path = write_tmp_file("water.inp", "! hf def2-svp\\n...")
    """

    def _write(name, content):
        return str(write_text_file(tmp_path / name, content))

    return _write


@pytest.fixture()
def invoke_gaussian_opt_cli(
    monkeypatch, single_molecule_xyz_file, gaussian_project_config_dir
):
    """Invoke ``chemsmart run``/``sub`` (``command``) with the given
    top-level ``options`` for a Gaussian ``opt`` job on
    ``single_molecule_xyz_file`` (project ``test``, charge 0, mult 1),
    using the ``gaussian_project_config_dir`` user config. Returns the
    Click result, e.g. ``invoke_gaussian_opt_cli(sub, ["--test"])``."""
    monkeypatch.setenv(
        "CHEMSMART_CONFIG_DIR", str(gaussian_project_config_dir)
    )

    def _invoke(command, options=()):
        args = [
            *options,
            "gaussian",
            "-p",
            "test",
            "-f",
            single_molecule_xyz_file,
            "-c",
            "0",
            "-m",
            "1",
            "opt",
        ]
        return CliRunner().invoke(command, args)

    return _invoke


@pytest.fixture()
def run_jobrunner_and_capture_results(monkeypatch, invoke_gaussian_opt_cli):
    """Run ``chemsmart run --fake [options] gaussian ... opt`` with
    ``JobRunner.from_job`` and ``Job.run`` stubbed out, capturing the
    ``server``/``scratch``/``fake`` arguments the CLI builds the runner
    with::

        result, captured = run_jobrunner_and_capture_results(["--scratch"])
        assert captured["scratch"] is True

    ``captured["scratch"]`` stays ``"unset"`` if ``from_job`` is never
    called.
    """

    def _run(options=()):
        captured = {"scratch": "unset"}

        def _from_job(cls, job, server, scratch=None, fake=False, **kwargs):
            captured.update(server=server, scratch=scratch, fake=fake)
            return type("R", (), {"scratch": scratch})()

        monkeypatch.setattr(
            "chemsmart.jobs.runner.JobRunner.from_job",
            classmethod(_from_job),
        )
        monkeypatch.setattr("chemsmart.jobs.job.Job.run", lambda self: None)
        result = invoke_gaussian_opt_cli(run, ["--fake", *options])
        return result, captured

    return _run


@pytest.fixture()
def submit_job_and_capture_results(monkeypatch, invoke_gaussian_opt_cli):
    """Run ``chemsmart sub --test --server dummy [options] gaussian ... opt``
    against a fake ``dummy`` :class:`Server` whose ``submit`` records the
    forwarded ``cli_args`` instead of submitting::

        result, server, captured = submit_job_and_capture_results(
            ["--queue", "gpu"]
        )
        assert server.queue_name == "gpu"
        assert "--scratch" not in captured["cli_args"]
    """

    def _submit(options=()):
        captured = {"cli_args": None}
        fake_server = Server(name="dummy")
        fake_server.submit = (
            lambda job, test=False, cli_args=None, **kw: captured.update(
                cli_args=cli_args
            )
        )
        monkeypatch.setattr(
            "chemsmart.settings.server.Server.from_servername",
            lambda _name: fake_server,
        )
        result = invoke_gaussian_opt_cli(
            sub, ["--test", "--server", "dummy", *options]
        )
        return result, fake_server, captured

    return _submit
