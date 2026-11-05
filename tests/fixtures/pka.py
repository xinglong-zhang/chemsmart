"""pKa CLI fixtures."""

from pathlib import Path

import pytest

from tests.helpers import FakeThermochemistry, write_test_backend_project


@pytest.fixture()
def install_fake_thermochemistry(monkeypatch):
    """Replace ``chemsmart.cli.pka.Thermochemistry`` with a stub that records
    the basename of every file it is constructed with; calling the factory
    installs it and returns that (optionally caller-supplied) list::

        constructed = install_fake_thermochemistry()
    """

    def _install(constructed=None):
        constructed = [] if constructed is None else constructed

        class _TrackingFakeThermochemistry(FakeThermochemistry):
            def __init__(self, filename, **kwargs):
                constructed.append(Path(filename).name)
                super().__init__(filename, **kwargs)

        monkeypatch.setattr(
            "chemsmart.cli.pka.Thermochemistry",
            _TrackingFakeThermochemistry,
        )
        return constructed

    return _install


@pytest.fixture()
def setup_sub_pka_batch_test(tmp_path, monkeypatch, captured):
    """Prepare a ``sub ... pka batch`` submission test for ``backend``: two
    acid xyz files plus a ``pka_scale.csv`` table in ``tmp_path``, a
    ``test`` project under a temporary ``CHEMSMART_CONFIG_DIR``, and a fake
    ``dummy`` server whose ``submit`` appends ``(job, test, cli_args)`` to
    ``captured["submissions"]``. Returns ``(table, captured)``::

        table, captured = setup_sub_pka_batch_test("gaussian")
    """

    def _setup(backend):
        acid1 = tmp_path / "acid1.xyz"
        acid1.write_text("2\nacid1\nC 0.0 0.0 0.0\nH 0.0 0.0 1.0\n")
        acid2 = tmp_path / "acid2.xyz"
        acid2.write_text("2\nacid2\nN 0.0 0.0 0.0\nH 0.0 0.0 1.0\n")

        table = tmp_path / "pka_scale.csv"
        table.write_text(
            "structure,filepath,proton_index,charge,multiplicity\n"
            f"acid1,{acid1},2,0,1\n"
            f"acid2,{acid2},2,1,2\n"
        )

        config_root = write_test_backend_project(tmp_path, backend)
        monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_root))

        from chemsmart.settings.server import Server

        fake_server = Server(name="dummy")
        captured["submissions"] = []
        fake_server.submit = (
            lambda job, test=False, cli_args=None, **kw: captured[
                "submissions"
            ].append((job, test, cli_args))
        )
        monkeypatch.setattr(
            "chemsmart.settings.server.Server.from_servername",
            lambda _name: fake_server,
        )
        return table, captured

    return _setup
