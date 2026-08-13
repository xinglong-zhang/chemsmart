"""Compute-free checks for the hardened PySCF execution boundary."""

import ast
import io
import json
import os
import subprocess
import sys
from types import ModuleType, SimpleNamespace
from unittest.mock import PropertyMock, patch

import numpy as np
import pytest

from chemsmart.cli.pyscf import pyscf
from chemsmart.cli.run import run
from chemsmart.cli.sub import sub
from chemsmart.io.molecules.structure import Molecule
from chemsmart.io.pyscf.output import (
    PySCFArtifactBindingError,
    read_pyscf_h5,
)
from chemsmart.jobs.pyscf import environment as environment_module
from chemsmart.jobs.pyscf.environment import (
    canonical_sha256,
    environment_blockers,
    sha256_file,
)
from chemsmart.jobs.pyscf.job import PySCFJob
from chemsmart.jobs.pyscf.runner import (
    FakePySCFJobRunner,
    PySCFJobRunner,
    PySCFPreflightError,
    PySCFResultValidationError,
    _PREVIEW_DEFERRED_RULES,
    _receipt_file_sha256,
    _run_receipt_state,
)
from chemsmart.jobs.pyscf.settings import PySCFJobSettings
from chemsmart.jobs.pyscf.writer import (
    PySCFScriptWriter,
    pyscf_td_response_materialization,
    write_pyscf_h5,
)
from chemsmart.jobs.runner import JobRunner
from chemsmart.settings.submitters import Submitter


def _source_function(source, name, namespace=None):
    """Compile one named function without executing its surrounding script."""
    tree = ast.parse(source)
    node = next(
        item
        for item in tree.body
        if isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef))
        and item.name == name
    )
    module = ast.Module(body=[node], type_ignores=[])
    ast.fix_missing_locations(module)
    scope = dict(namespace or {})
    exec(compile(module, f"<{name}-only>", "exec"), scope)
    return scope[name]


def test_run_and_sub_expose_the_same_pyscf_command_tree():
    assert run.commands["pyscf"] is pyscf
    assert sub.commands["pyscf"] is pyscf
    assert tuple(sorted(run.commands["pyscf"].commands)) == (
        "hess",
        "opt",
        "sp",
        "td",
    )
    assert tuple(sorted(sub.commands["pyscf"].commands)) == (
        "hess",
        "opt",
        "sp",
        "td",
    )


def _td_preview_job():
    return SimpleNamespace(
        settings=PySCFJobSettings(
            jobtype="td",
            functional="b3lyp",
            basis="def2-svp",
            response_method="tda",
            state_manifold="singlet",
            nstates=4,
            charge=0,
            multiplicity=1,
            engine="cpu",
        ),
        molecule=SimpleNamespace(
            chemical_symbols=["O", "H", "H"],
            positions=[
                [0.0, 0.0, 0.0],
                [0.7586, 0.0, 0.5043],
                [-0.7586, 0.0, 0.5043],
            ],
            charge=0,
            multiplicity=1,
            info={},
        ),
    )


def _prepare_preflight_runner(runner, tmp_path):
    runner.num_cores = 1
    runner.num_gpus = 0
    runner.mem_gb = 1
    runner.environment_receiptfile = str(tmp_path / "td.environment.json")
    return runner


def test_fake_td_defers_only_preview_capability_rule(tmp_path):
    runner = _prepare_preflight_runner(
        object.__new__(FakePySCFJobRunner), tmp_path
    )
    job = _td_preview_job()

    with (
        patch.object(runner, "_assign_variables"),
        patch.object(runner, "_set_receipt_paths"),
        patch.object(runner, "_current_input_artifact_binding", return_value=None),
        patch.object(runner, "_quarantine_stale_targets"),
    ):
        runner._prerun(job)

    receipt = json.loads(
        (tmp_path / "td.environment.json").read_text(encoding="utf-8")
    )
    assert receipt["deferred_rule_ids"] == [
        "pyscf.td.preview_only_capability"
    ]
    assert [item["rule_id"] for item in receipt["preflight_findings"]] == [
        "pyscf.td.preview_only_capability"
    ]
    assert "pyscf.td.preview_only_capability" in _PREVIEW_DEFERRED_RULES


def test_real_td_blocks_before_input_or_engine_launch(tmp_path):
    runner = _prepare_preflight_runner(object.__new__(PySCFJobRunner), tmp_path)
    job = _td_preview_job()

    with (
        patch.object(runner, "_assign_variables"),
        patch.object(runner, "_set_receipt_paths"),
        patch.object(runner, "_current_input_artifact_binding", return_value=None),
        patch.object(runner, "_quarantine_stale_targets"),
        patch.object(runner, "_update_os_environ", return_value={}),
        patch.object(
            PySCFJobRunner,
            "executable",
            new_callable=PropertyMock,
            return_value=SimpleNamespace(get_executable=lambda: sys.executable),
        ),
        patch(
            "chemsmart.jobs.pyscf.runner.probe_compute_environment",
            return_value={"status": "available"},
        ),
        patch(
            "chemsmart.jobs.pyscf.runner.environment_blockers",
            return_value=[],
        ),
    ):
        with pytest.raises(PySCFPreflightError) as error:
            runner._prerun(job)

    assert [item["rule_id"] for item in error.value.findings] == [
        "pyscf.td.preview_only_capability"
    ]


def test_td_preview_driver_is_inert_even_when_run_directly():
    source = PySCFScriptWriter.render(
        {
            "schema_version": "2.0",
            "label": "td-preview-only",
            "jobtype": "td",
            "preview_only": True,
        }
    )

    compile(source, "<pyscf-td-preview-driver>", "exec")
    assert 'if CONFIG.get("preview_only")' in source
    assert ".TDA()" not in source
    assert ".TDDFT()" not in source
    assert "response.kernel()" not in source


@pytest.mark.parametrize(
    ("response_method", "factory_api"),
    (
        ("tda", "pyscf.tdscf.rks.TDA"),
        ("tddft", "pyscf.tdscf.rks.TDDFT"),
    ),
)
def test_td_preview_materializes_a_typed_non_executable_response_plan(
    response_method, factory_api
):
    settings = PySCFJobSettings(
        jobtype="td",
        functional="b3lyp",
        basis="def2-svp",
        response_method=response_method,
        state_manifold="singlet",
        nstates=4,
        charge=0,
        multiplicity=1,
        engine="cpu",
    )

    materialization = pyscf_td_response_materialization(settings)

    assert materialization["response_factory_api"] == factory_api
    assert materialization["nstates"] == 4
    assert materialization["state_manifold"] == "singlet"
    assert materialization["execution_policy"] == "preview_only"
    body = dict(materialization)
    observed = body.pop("receipt_sha256")
    assert observed == canonical_sha256(body)


def test_fake_result_is_explicitly_not_evaluated_and_has_no_fake_chemistry(
    tmp_path,
):
    runner = object.__new__(FakePySCFJobRunner)
    runner.job_outputfile = str(tmp_path / "water.out")
    runner.job_resultsfile = str(tmp_path / "water.h5")
    runner._input_receipt = {
        "script_sha256": "a" * 64,
        "receipt_sha256": "b" * 64,
    }
    runner._environment_receipt = {"receipt_sha256": "c" * 64}
    runner._input_config = {
        "run_id": "preview-run",
        "run_nonce": "preview-nonce",
        "program": "pyscf",
        "jobtype": "hess",
        "engine": "cpu",
        "stages": ["scf", "hess"],
        "symbols": ["O", "H", "H"],
        "positions": [
            [0.0, 0.0, 0.1174],
            [-0.757, 0.0, -0.4696],
            [0.757, 0.0, -0.4696],
        ],
        "unit": "Angstrom",
        "charge": 0,
        "multiplicity": 1,
        "num_threads": 4,
        "settings_digest": "d" * 64,
        "requested_settings_sha256": "e" * 64,
        "input_geometry_sha256": "f" * 64,
        "input_artifact_kind": "geometry_xyz",
        "input_artifact_sha256": "1" * 64,
        "project_yaml_digest": "2" * 64,
    }

    runner._write_fake_results(SimpleNamespace())
    _spec, _provenance, status, results = read_pyscf_h5(
        runner.job_resultsfile
    )

    assert status["evaluation_state"] == "not_evaluated_preview"
    assert status["normal_termination"] is False
    assert status["engine_complete"] is False
    assert status["stages"]["scf"]["converged"] is None
    assert status["stages"]["hess"]["converged"] is None
    assert tuple(results) == ("positions",)


def test_unclassified_direct_hessian_preserves_engine_completion_only():
    state = _run_receipt_state(
        fake=False,
        findings=[],
        engine_complete=True,
        result_validation_state="unclassified",
    )

    assert state == "engine_complete"
    assert state != "validated"


def test_nonzero_pyscf_child_status_is_preserved():
    runner = object.__new__(PySCFJobRunner)
    command = ("/opt/pyscf/bin/python", "water.py")
    validation_error = PySCFResultValidationError(
        ({"rule_id": "pyscf.process.nonzero_exit"},)
    )
    with (
        patch.object(runner, "_reset_run_state"),
        patch.object(runner, "_begin_run_identity"),
        patch.object(runner, "_prerun"),
        patch.object(runner, "_write_input"),
        patch.object(runner, "_get_command", return_value=command),
        patch.object(runner, "_update_os_environ", return_value={}),
        patch.object(runner, "_create_process", return_value=object()),
        patch.object(runner, "_run", return_value=23),
        patch.object(runner, "_postrun"),
        patch.object(
            runner,
            "_postrun_cleanup",
            side_effect=validation_error,
        ),
    ):
        with pytest.raises(subprocess.CalledProcessError) as error:
            runner.run(object())

    assert error.value.returncode == 23
    assert error.value.cmd == command


def test_base_runner_preserves_nonzero_child_status_after_cleanup():
    runner = object.__new__(JobRunner)
    runner.server = SimpleNamespace(name="test-server")
    command = ("/opt/engine/bin/program", "input.dat")
    with (
        patch.object(runner, "_prerun"),
        patch.object(runner, "_write_input"),
        patch.object(runner, "_get_command", return_value=command),
        patch.object(runner, "_update_os_environ", return_value={}),
        patch.object(runner, "_create_process", return_value=object()),
        patch.object(runner, "_run", return_value=19),
        patch.object(runner, "_postrun") as postrun,
        patch.object(runner, "_postrun_cleanup") as cleanup,
    ):
        with pytest.raises(subprocess.CalledProcessError) as error:
            runner.run(object())

    assert postrun.call_count == 1
    assert cleanup.call_count == 1
    assert error.value.returncode == 19
    assert error.value.cmd == command


class _SubmitterProbe(Submitter):
    def _write_scheduler_options(self, _handle):
        return None

    def _write_change_to_job_directory(self, _handle):
        return None


def test_scheduler_replaces_shell_and_propagates_child_exit():
    submitter = object.__new__(_SubmitterProbe)
    submitter.job = SimpleNamespace(label="pyscf_sp")
    handle = io.StringIO()

    submitter._write_job_command(handle)

    assert handle.getvalue() == (
        "chmod +x ./chemsmart_run_pyscf_sp.py\n"
        "exec ./chemsmart_run_pyscf_sp.py\n"
    )


def test_pyscf_process_uses_typed_argv_with_space_containing_paths(tmp_path):
    runner = object.__new__(PySCFJobRunner)
    runner.job_inputfile = str(tmp_path / "input folder" / "job.py")
    runner.job_errfile = str(tmp_path / "job.err")
    runner.job_outputfile = str(tmp_path / "job.out")
    runner.job_resultsfile = str(tmp_path / "job.h5")
    runner.running_directory = str(tmp_path)
    runner._environment_receipt = {
        "interpreter": str(tmp_path / "python env" / "python")
    }
    argv = runner._get_command(None)

    with patch("chemsmart.jobs.pyscf.runner.subprocess.Popen") as popen:
        runner._create_process(None, argv, {"PATH": "bounded"})

    assert argv == (
        str(tmp_path / "python env" / "python"),
        str(tmp_path / "input folder" / "job.py"),
    )
    assert popen.call_args.args[0] == argv
    assert popen.call_args.kwargs.get("shell", False) is False


def test_environment_receipt_rejects_interpreter_substitution():
    receipt = {
        "status": "available",
        "interpreter": sys.executable,
        "interpreter_observed": os.devnull,
        "dependencies": {
            "pyscf": {"available": True, "version": "2.14.0"},
            "numpy": {"available": True},
            "h5py": {"available": True},
        },
    }

    findings = environment_blockers(receipt, engine="cpu")

    assert any(
        finding["rule_id"] == "pyscf.environment.interpreter_mismatch"
        for finding in findings
    )


def test_environment_receipt_rejects_unloadable_basis_and_functional():
    receipt = {
        "status": "available",
        "interpreter": sys.executable,
        "interpreter_observed": sys.executable,
        "dependencies": {
            "pyscf": {"available": True, "version": "2.14.0"},
            "numpy": {"available": True},
            "h5py": {"available": True},
        },
        "request": {"basis": "invented-basis", "functional": "invented-xc"},
        "basis_available": {"invented-basis": False},
        "functional_available": {"invented-xc": False},
    }

    rule_ids = {
        finding["rule_id"]
        for finding in environment_blockers(receipt, engine="cpu")
    }

    assert "pyscf.environment.basis_unavailable" in rule_ids
    assert "pyscf.environment.functional_unavailable" in rule_ids


def test_embedded_basis_probe_returns_known_maximum_angular_momentum():
    fake_pyscf = ModuleType("pyscf")
    fake_pyscf.gto = SimpleNamespace(
        basis=SimpleNamespace(
            load=lambda name, symbol: (
                [[0, [1.0, 1.0]], [2, [0.5, 1.0]]]
                if (name, symbol) == ("def2-svp", "O")
                else [[0, [1.0, 1.0]]]
            )
        )
    )
    basis_max_l = _source_function(
        environment_module._PROBE_SCRIPT, "basis_max_l"
    )

    with patch.dict(sys.modules, {"pyscf": fake_pyscf}):
        assert basis_max_l("def2-svp", ["H", "O"]) == 2


def test_environment_blocks_basis_that_requires_unmaterialized_ecp():
    receipt = {
        "status": "available",
        "interpreter": sys.executable,
        "interpreter_observed": sys.executable,
        "dependencies": {
            "pyscf": {"available": True, "version": "2.14.0"},
            "numpy": {"available": True},
            "h5py": {"available": True},
        },
        "request": {"basis": "def2-svp"},
        "basis_available": {"def2-svp": True},
        "basis_ecp_required_elements": ["Xe"],
    }

    findings = environment_blockers(receipt, engine="cpu")

    assert any(
        finding["rule_id"] == "pyscf.environment.ecp_unmaterialized"
        and finding["observed"]["elements_requiring_ecp"] == ["Xe"]
        for finding in findings
    )


def _complete_gpu_receipt():
    return {
        "status": "available",
        "interpreter": sys.executable,
        "interpreter_observed": sys.executable,
        "dependencies": {
            "pyscf": {"available": True, "version": "2.14.0"},
            "numpy": {"available": True, "version": np.__version__},
            "h5py": {"available": True, "version": "test"},
            "gpu4pyscf": {"available": True, "version": "1.8.0"},
            "cupy": {"available": True, "version": "13.3.0"},
        },
        "gpu4pyscf_distribution": {
            "name": "gpu4pyscf-cuda12x",
            "version": "1.8.0",
        },
        "cupy_distribution": {
            "name": "cupy-cuda12x",
            "version": "13.3.0",
        },
        "cutensor_distribution": {
            "name": "cutensor-cu12",
            "version": "2.0.2",
        },
        "device_count": 1,
        "device_name": "test GPU",
        "device_uuid": "GPU-test-uuid",
        "cuda_driver_version": 12040,
        "cuda_runtime_version": 12040,
        "cutensor_runtime": 20002,
        "cutensor_compatible": True,
    }


def test_complete_pinned_gpu_environment_has_no_environment_blocker():
    assert environment_blockers(_complete_gpu_receipt(), engine="gpu") == []


@pytest.mark.parametrize(
    ("mutation", "rule_id"),
    [
        (
            lambda receipt: receipt["dependencies"]["pyscf"].update(
                version="2.13.0"
            ),
            "pyscf.environment.version_mismatch",
        ),
        (
            lambda receipt: receipt.update(device_uuid=None),
            "pyscf.gpu.runtime_evidence_incomplete",
        ),
        (
            lambda receipt: receipt["cupy_distribution"].update(
                name="cupy-cuda11x"
            ),
            "pyscf.gpu.cuda_suffix_mismatch",
        ),
    ],
)
def test_pinned_gpu_environment_fails_closed_on_evidence_drift(
    mutation, rule_id
):
    receipt = _complete_gpu_receipt()
    mutation(receipt)

    rule_ids = {
        finding["rule_id"]
        for finding in environment_blockers(receipt, engine="gpu")
    }

    assert rule_id in rule_ids


def test_environment_probe_and_driver_sources_are_syntax_bounded():
    compile(
        environment_module._PROBE_SCRIPT,
        "<pyscf-environment-probe>",
        "exec",
    )
    source = PySCFScriptWriter.render(
        {"schema_version": "2.0", "label": "source-only"}
    )
    compile(source, "<pyscf-driver>", "exec")

    assert "gpu4pyscf-cuda" in environment_module._PROBE_SCRIPT
    assert "libxc.libxc_version()" in environment_module._PROBE_SCRIPT
    assert '"libxc_version": libxc_runtime_version' in (
        environment_module._PROBE_SCRIPT
    )
    assert "post-construction wrapper returned" in source
    assert 'status["properties"]' in source
    assert '"status": "unavailable"' in source


def test_force_free_sp_driver_does_not_launch_an_undeclared_gradient():
    source = PySCFScriptWriter.render(
        {"schema_version": "2.0", "label": "force-free-sp"}
    )

    assert "nuc_grad_method" not in source
    assert 'status["properties"]["forces"]' not in source


def test_generated_host_transfer_calls_gpu_array_get_before_numpy_conversion():
    class FakeGpuArray:
        def __init__(self):
            self.get_calls = 0

        def get(self):
            self.get_calls += 1
            return [[1.0, 2.0], [3.0, 4.0]]

    source = PySCFScriptWriter.render(
        {"schema_version": "2.0", "label": "host-transfer-only"}
    )
    to_host = _source_function(source, "_to_host_array", {"np": np})
    array = FakeGpuArray()

    observed = to_host(array)

    assert array.get_calls == 1
    np.testing.assert_array_equal(
        observed, np.asarray([[1.0, 2.0], [3.0, 4.0]])
    )
    assert source.index("def _to_host_array") < source.index("def main")
    assert "_to_host_array(mf.Hessian().kernel())" in source


def test_generated_driver_materializes_target_libxc_definition():
    class FakeLibxc:
        @staticmethod
        def parse_xc(value):
            assert value == "b3lypg"
            return (0.2, 0.2, 0.0), ((402, 1.0),)

        @staticmethod
        def libxc_version():
            return "7.0.0"

        @staticmethod
        def hybrid_coeff(value):
            assert value == "b3lypg"
            return 0.2

        @staticmethod
        def rsh_coeff(value):
            assert value == "b3lypg"
            return (0.0, 0.2, 0.0)

    source = PySCFScriptWriter.render(
        {"schema_version": "2.0", "label": "functional-materialization"}
    )
    materialize = _source_function(source, "_functional_definition")

    record = materialize(
        FakeLibxc,
        "b3lypg",
        pyscf_version="2.14.0",
        environment_receipt_sha256="e" * 64,
    )

    assert record == {
        "schema_version": "chemsmart.pyscf-functional-definition.v3",
        "field": "xc",
        "source": "pyscf.dft.libxc.parse_xc",
        "source_key": "b3lypg",
        "pyscf_version": "2.14.0",
        "libxc_version": "7.0.0",
        "environment_receipt_sha256": "e" * 64,
        "parser_hybrid_decomposition": [0.2, 0.2, 0.0],
        "exact_exchange_fraction": 0.2,
        "range_separation_coefficients": [0.0, 0.2, 0.0],
        "functional_ids": [402],
        "functional_factors": [1.0],
    }


def test_generated_driver_records_scalar_spin_diagnostic():
    class FakeMeanField:
        @staticmethod
        def spin_square():
            return np.asarray(2.3), np.asarray(np.sqrt(10.2))

    source = PySCFScriptWriter.render(
        {"schema_version": "2.0", "label": "spin-diagnostic"}
    )
    to_host = _source_function(source, "_to_host_array", {"np": np})
    diagnostic = _source_function(
        source,
        "_spin_diagnostic",
        {"_to_host_array": to_host},
    )

    observed = diagnostic(FakeMeanField())

    assert observed == pytest.approx((2.3, np.sqrt(10.2)))


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("semiempirical", "gfn2-xtb"),
        ("modred", "B 1 2 F"),
        ("forces", True),
        ("input_string", "arbitrary native input"),
    ],
)
def test_inherited_noop_settings_fail_closed(field, value):
    settings = PySCFJobSettings(
        functional="pbe",
        basis="def2-svp",
        **{field: value},
    )

    with pytest.raises(ValueError, match="not supported"):
        settings.validate()


def test_unknown_setting_is_not_silently_discarded():
    with pytest.raises(ValueError, match="Unknown PySCF setting"):
        PySCFJobSettings(
            functional="pbe",
            basis="def2-svp",
            invented_option="silent-noop",
        )


def test_absent_hdf5_is_neither_complete_nor_validated(tmp_path):
    job = object.__new__(PySCFJob)
    job.folder = str(tmp_path)
    job.label = "historical"
    job.settings = PySCFJobSettings(functional="pbe", basis="def2-svp")

    assert job.is_complete() is False
    assert job.is_validated() is False


def _write_bound_receipt(path, body):
    payload = dict(body)
    payload["receipt_sha256"] = canonical_sha256(payload)
    path.write_text(json.dumps(payload), encoding="utf-8")
    return payload


def _write_upstream_pyscf_artifact(tmp_path, label="upstream"):
    logfile = tmp_path / f"{label}.out"
    resultsfile = tmp_path / f"{label}.h5"
    logfile.write_text("PySCF version 2.14.0\n", encoding="utf-8")
    positions = np.asarray(
        [[0.0, 0.0, 0.1], [0.0, 0.7, -0.4], [0.0, -0.7, -0.4]],
        dtype=float,
    )
    write_pyscf_h5(
        resultsfile,
        spec={
            "program": "pyscf",
            "jobtype": "opt",
            "stages": ["scf", "opt"],
            "symbols": ["O", "H", "H"],
            "charge": 0,
            "multiplicity": 1,
        },
        provenance={"engine": "cpu"},
        status={"normal_termination": True, "failure": None},
        results={"positions": positions, "energies": [-76.0]},
    )
    return logfile, resultsfile, positions


def test_legacy_normal_hdf5_is_complete_without_modern_receipts(tmp_path):
    _write_upstream_pyscf_artifact(tmp_path, label="legacy")
    job = object.__new__(PySCFJob)
    job.folder = str(tmp_path)
    job.label = "legacy"
    job.settings = PySCFJobSettings(
        functional="pbe",
        basis="def2-svp",
        jobtype="opt",
        charge=0,
        multiplicity=1,
    )

    assert job.is_complete() is True
    assert job.is_validated() is False


def test_out_entrypoint_binds_hdf5_and_substitution_fails_preflight(tmp_path):
    logfile, resultsfile, _ = _write_upstream_pyscf_artifact(tmp_path)
    molecule = Molecule.from_filepath(logfile)
    job = SimpleNamespace(
        molecule=molecule,
        inputfile=str(tmp_path / "downstream.py"),
        outputfile=str(tmp_path / "downstream.out"),
        resultsfile=str(tmp_path / "downstream.h5"),
        errfile=str(tmp_path / "downstream.err"),
    )
    runner = object.__new__(PySCFJobRunner)
    runner.job_inputfile = job.inputfile
    runner.job_outputfile = job.outputfile
    runner.job_resultsfile = job.resultsfile
    runner.job_errfile = job.errfile

    binding = runner._current_input_artifact_binding(job)

    assert binding["path"] == str(resultsfile)
    assert binding["sha256"] == sha256_file(resultsfile)
    resultsfile.write_bytes(resultsfile.read_bytes() + b"substituted")
    with pytest.raises(
        PySCFArtifactBindingError,
        match="pyscf.input_artifact.hash_mismatch",
    ):
        runner._current_input_artifact_binding(job)


def test_source_target_collision_is_blocked_without_quarantining_source(
    tmp_path,
):
    _, resultsfile, _ = _write_upstream_pyscf_artifact(tmp_path)
    molecule = Molecule.from_filepath(resultsfile)
    job = SimpleNamespace(
        molecule=molecule,
        folder=str(tmp_path),
        label="upstream",
        inputfile=str(tmp_path / "upstream.py"),
        outputfile=str(tmp_path / "upstream.out"),
        resultsfile=str(resultsfile),
        errfile=str(tmp_path / "upstream.err"),
    )
    runner = object.__new__(PySCFJobRunner)
    runner._run_id = "collision-run"
    runner.job_inputfile = job.inputfile
    runner.job_outputfile = job.outputfile
    runner.job_resultsfile = job.resultsfile
    runner.job_errfile = job.errfile
    runner._set_receipt_paths(job)

    with pytest.raises(
        PySCFArtifactBindingError,
        match="pyscf.input_artifact.target_collision",
    ):
        runner._current_input_artifact_binding(job)

    original_sha256 = sha256_file(resultsfile)
    runner._quarantine_stale_targets(
        job,
        protected_paths=runner._declared_input_artifact_paths(job),
    )
    assert resultsfile.is_file()
    assert sha256_file(resultsfile) == original_sha256


def test_receipt_digest_rejects_file_substitution(tmp_path):
    path = tmp_path / "bound.input.json"
    receipt = _write_bound_receipt(path, {"state": "previewed"})
    assert _receipt_file_sha256(path) == receipt["receipt_sha256"]

    path.write_text(
        json.dumps(dict(receipt, state="substituted")), encoding="utf-8"
    )
    assert _receipt_file_sha256(path) is None


def test_json_receipt_writer_refuses_symlink_destination(tmp_path):
    victim = tmp_path / "victim.json"
    victim.write_text("unchanged\n", encoding="utf-8")
    destination = tmp_path / "receipt.json"
    destination.symlink_to(victim)

    with pytest.raises(ValueError, match="symlink receipt destination"):
        environment_module.write_json_receipt(destination, {"state": "red"})

    assert victim.read_text(encoding="utf-8") == "unchanged\n"
    assert not tuple(tmp_path.glob(".receipt.json.*.tmp"))


def test_validated_state_is_invalidated_by_current_geometry_change(tmp_path):
    job = object.__new__(PySCFJob)
    job.folder = str(tmp_path)
    job.label = "bound"
    job.settings = PySCFJobSettings(
        functional="pbe",
        basis="def2-svp",
        charge=0,
        multiplicity=1,
    )
    job.molecule = SimpleNamespace(
        chemical_symbols=["H", "H"],
        positions=[[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]],
        charge=0,
        multiplicity=1,
    )
    (tmp_path / "bound.py").write_text("# immutable script\n")
    (tmp_path / "bound.h5").write_bytes(b"immutable result")
    run_id = "run-bound"
    run_nonce = "nonce-bound"
    environment = _write_bound_receipt(
        tmp_path / "bound.environment.json",
        {
            "run_id": run_id,
            "run_nonce": run_nonce,
            "status": "available",
        },
    )
    input_receipt = _write_bound_receipt(
        tmp_path / "bound.input.json",
        {
            "run_id": run_id,
            "run_nonce": run_nonce,
            "state": "previewed",
        },
    )
    _write_bound_receipt(
        tmp_path / "bound.receipt.json",
        {
            "run_id": run_id,
            "run_nonce": run_nonce,
            "state": "validated",
            "fake": False,
            "child_returncode": 0,
            "findings": [],
            "script_sha256": sha256_file(tmp_path / "bound.py"),
            "result_sha256": sha256_file(tmp_path / "bound.h5"),
            "input_receipt_sha256": input_receipt["receipt_sha256"],
            "environment_receipt_sha256": environment["receipt_sha256"],
            "input_geometry_sha256": job._current_geometry_sha256(),
            "requested_settings_sha256": "requested",
            "project_yaml_sha256": None,
        },
    )

    with patch(
        "chemsmart.jobs.pyscf.job.verify_provenance", return_value=[]
    ):
        assert job.is_validated() is True
        job.molecule.positions[1][2] = 0.75
        assert job.is_validated() is False


def test_validated_state_is_invalidated_by_source_hdf5_substitution(tmp_path):
    _, source_results, _ = _write_upstream_pyscf_artifact(
        tmp_path, label="optimized"
    )
    job = object.__new__(PySCFJob)
    job.folder = str(tmp_path)
    job.label = "downstream_hess"
    job.settings = PySCFJobSettings(
        functional="pbe",
        basis="def2-svp",
        jobtype="hess",
        freq=True,
        charge=0,
        multiplicity=1,
    )
    job.molecule = Molecule.from_filepath(source_results)
    (tmp_path / "downstream_hess.py").write_text(
        "# immutable script\n", encoding="utf-8"
    )
    (tmp_path / "downstream_hess.h5").write_bytes(b"immutable result")
    run_id = "run-source-bound"
    run_nonce = "nonce-source-bound"
    environment = _write_bound_receipt(
        tmp_path / "downstream_hess.environment.json",
        {"run_id": run_id, "run_nonce": run_nonce, "status": "available"},
    )
    input_receipt = _write_bound_receipt(
        tmp_path / "downstream_hess.input.json",
        {"run_id": run_id, "run_nonce": run_nonce, "state": "previewed"},
    )
    binding = job.molecule.info["chemsmart_source_artifact"]
    _write_bound_receipt(
        tmp_path / "downstream_hess.receipt.json",
        {
            "run_id": run_id,
            "run_nonce": run_nonce,
            "state": "validated",
            "fake": False,
            "child_returncode": 0,
            "findings": [],
            "script_sha256": sha256_file(tmp_path / "downstream_hess.py"),
            "result_sha256": sha256_file(tmp_path / "downstream_hess.h5"),
            "input_receipt_sha256": input_receipt["receipt_sha256"],
            "environment_receipt_sha256": environment["receipt_sha256"],
            "input_geometry_sha256": job._current_geometry_sha256(),
            "input_artifact_kind": binding["kind"],
            "input_artifact_path": binding["path"],
            "input_artifact_sha256": binding["sha256"],
            "requested_settings_sha256": "requested",
            "project_yaml_sha256": None,
        },
    )

    with patch(
        "chemsmart.jobs.pyscf.job.verify_provenance", return_value=[]
    ):
        assert job.is_validated() is True
        source_results.write_bytes(
            source_results.read_bytes() + b"substituted"
        )
        assert job.is_validated() is False


def test_validated_state_rejects_receipt_from_a_different_run_nonce(tmp_path):
    job = object.__new__(PySCFJob)
    job.folder = str(tmp_path)
    job.label = "bound"
    job.settings = PySCFJobSettings(
        functional="pbe",
        basis="def2-svp",
        jobtype="sp",
        charge=0,
        multiplicity=1,
    )
    job.molecule = SimpleNamespace(
        chemical_symbols=["H", "H"],
        positions=[[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]],
        charge=0,
        multiplicity=1,
    )
    (tmp_path / "bound.py").write_text("# immutable script\n")
    (tmp_path / "bound.h5").write_bytes(b"immutable result")
    environment = _write_bound_receipt(
        tmp_path / "bound.environment.json",
        {"run_id": "run-a", "run_nonce": "nonce-b"},
    )
    input_receipt = _write_bound_receipt(
        tmp_path / "bound.input.json",
        {"run_id": "run-a", "run_nonce": "nonce-a"},
    )
    _write_bound_receipt(
        tmp_path / "bound.receipt.json",
        {
            "run_id": "run-a",
            "run_nonce": "nonce-a",
            "state": "validated",
            "fake": False,
            "child_returncode": 0,
            "findings": [],
            "script_sha256": sha256_file(tmp_path / "bound.py"),
            "result_sha256": sha256_file(tmp_path / "bound.h5"),
            "input_receipt_sha256": input_receipt["receipt_sha256"],
            "environment_receipt_sha256": environment["receipt_sha256"],
            "input_geometry_sha256": job._current_geometry_sha256(),
        },
    )

    assert job.is_validated() is False


def test_stale_quarantine_moves_only_exact_run_targets(tmp_path):
    job = SimpleNamespace(
        folder=str(tmp_path),
        label="water",
        inputfile=str(tmp_path / "water.py"),
        outputfile=str(tmp_path / "water.out"),
        resultsfile=str(tmp_path / "water.h5"),
        errfile=str(tmp_path / "water.err"),
    )
    runner = object.__new__(PySCFJobRunner)
    runner._run_id = "run-identity"
    runner.job_inputfile = job.inputfile
    runner.job_outputfile = job.outputfile
    runner.job_resultsfile = job.resultsfile
    runner.job_errfile = job.errfile
    runner._set_receipt_paths(job)
    exact = tmp_path / "water.h5"
    exact.write_bytes(b"stale exact result")
    near_collision = tmp_path / "water.h5.backup"
    near_collision.write_bytes(b"must remain")

    runner._quarantine_stale_targets(job)

    assert not exact.exists()
    assert near_collision.read_bytes() == b"must remain"
    assert len(runner._quarantined_targets) == 1
    quarantined = tmp_path / runner._quarantined_targets[0]["quarantined_to"]
    assert quarantined.read_bytes() == b"stale exact result"
    assert runner._quarantined_targets[0]["source_sha256"] == sha256_file(
        quarantined
    )


@pytest.mark.parametrize(
    ("num_cores", "num_gpus", "mem_gb", "engine", "message"),
    [
        (0, 0, 1, "cpu", "num_cores"),
        (1, -1, 1, "cpu", "num_gpus"),
        (1, 0, 0, "cpu", "mem_gb"),
        (1, 0, 1, "gpu", "positive resolved num_gpus"),
    ],
)
def test_runner_resources_fail_closed(
    num_cores, num_gpus, mem_gb, engine, message
):
    runner = object.__new__(PySCFJobRunner)
    runner.num_cores = num_cores
    runner.num_gpus = num_gpus
    runner.mem_gb = mem_gb
    job = SimpleNamespace(settings=SimpleNamespace(engine=engine))

    with pytest.raises(ValueError, match=message):
        runner._validate_resources(job)


def test_fake_gpu_resources_allow_preview_without_device_allocation():
    runner = object.__new__(FakePySCFJobRunner)
    runner.num_cores = 1
    runner.num_gpus = 0
    runner.mem_gb = 1
    job = SimpleNamespace(settings=SimpleNamespace(engine="gpu"))

    runner._validate_resources(job)


def test_optional_property_omission_is_reported_but_not_required():
    job = SimpleNamespace(required_properties=(), kwargs={})
    status = {
        "properties": {
            "dipole_moment": {
                "status": "unavailable",
                "failure": {"type": "NotImplementedError"},
            }
        }
    }

    findings = PySCFJobRunner._property_findings(job, status)

    assert findings == [
        {
            "rule_id": "pyscf.property.unavailable",
            "field": "properties.dipole_moment",
            "expected": "optional",
            "observed": status["properties"]["dipole_moment"],
            "required": False,
            "evidence_ref": "h5:/status/properties/dipole_moment",
        }
    ]


def test_spin_diagnostic_unavailability_uses_existing_required_property_gate():
    status = {
        "properties": {
            "spin_square": {
                "status": "unavailable",
                "failure": {"type": "NotImplementedError"},
            }
        }
    }
    optional = SimpleNamespace(required_properties=(), kwargs={})
    required = SimpleNamespace(
        required_properties=("spin_square",), kwargs={}
    )

    optional_findings = PySCFJobRunner._property_findings(optional, status)
    required_findings = PySCFJobRunner._property_findings(required, status)

    assert optional_findings[0]["required"] is False
    assert required_findings[0]["required"] is True


def test_postrun_environment_reloads_receipt_and_rehashes_interpreter(tmp_path):
    interpreter = tmp_path / "python"
    interpreter.write_bytes(b"pinned interpreter bytes")
    receipt_path = tmp_path / "job.environment.json"
    _write_bound_receipt(
        receipt_path,
        {
            "interpreter": str(interpreter),
            "interpreter_observed": str(interpreter),
            "interpreter_sha256": sha256_file(interpreter),
            "libxc_version": "7.0.0",
            "dependencies": {
                "pyscf": {"version": "2.14.0"},
                "numpy": {"version": "1.26.4"},
                "h5py": {"version": "3.16.0"},
            },
        },
    )
    provenance = {
        "interpreter": str(interpreter),
        "pyscf_version": "2.14.0",
        "libxc_version": "7.0.0",
        "numpy_version": "1.26.4",
        "h5py_version": "3.16.0",
        "engine": "cpu",
    }
    runner = object.__new__(PySCFJobRunner)
    runner.environment_receiptfile = str(receipt_path)

    assert runner._runtime_environment_findings(provenance) == []

    interpreter.write_bytes(b"substituted interpreter bytes")
    findings = runner._runtime_environment_findings(provenance)
    assert any(
        finding["rule_id"] == "pyscf.environment.runtime_mismatch"
        and finding["field"] == "interpreter_sha256"
        for finding in findings
    )


def test_postrun_environment_rejects_libxc_runtime_mismatch(tmp_path):
    interpreter = tmp_path / "python"
    interpreter.write_bytes(b"pinned interpreter bytes")
    receipt_path = tmp_path / "job.environment.json"
    _write_bound_receipt(
        receipt_path,
        {
            "interpreter": str(interpreter),
            "interpreter_observed": str(interpreter),
            "interpreter_sha256": sha256_file(interpreter),
            "libxc_version": "7.0.0",
            "dependencies": {
                "pyscf": {"version": "2.14.0"},
                "numpy": {"version": "1.26.4"},
                "h5py": {"version": "3.16.0"},
            },
        },
    )
    runner = object.__new__(PySCFJobRunner)
    runner.environment_receiptfile = str(receipt_path)
    provenance = {
        "interpreter": str(interpreter),
        "pyscf_version": "2.14.0",
        "libxc_version": "6.2.2",
        "numpy_version": "1.26.4",
        "h5py_version": "3.16.0",
        "engine": "cpu",
    }

    findings = runner._runtime_environment_findings(provenance)

    assert any(
        finding["field"] == "provenance.libxc_version"
        for finding in findings
    )
