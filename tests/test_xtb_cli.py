"""Focused contracts for the xTB CLI execution-plane slice."""

import json
import math
import shutil
import subprocess
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest
import yaml
from click.testing import CliRunner

from chemsmart.cli.run import run
from chemsmart.cli.sub import sub
from chemsmart.cli.xtb.xtb import xtb
from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.xtb.hess import XTBHessJob
from chemsmart.jobs.xtb.opt import XTBOptJob
from chemsmart.jobs.xtb.runner import FakeXTBJobRunner, XTBJobRunner
from chemsmart.jobs.xtb.settings import XTBJobSettings
from chemsmart.jobs.xtb.singlepoint import XTBSinglePointJob
from chemsmart.jobs.xtb.validation import (
    _openblas_openmp_warning_count,
    audit_xtb_result_receipt,
    bind_xtb_execution_input,
    canonical_sha256,
    finalize_receipt,
    probe_xtb_environment,
    sha256_file,
    validate_xtb_result,
)
from chemsmart.settings.capabilities import program_capability
from chemsmart.settings.executable import XTBExecutable
from chemsmart.settings.server import Server
from chemsmart.settings.submitters import RunScript, SLURMSubmitter
from chemsmart.settings.xtb import XTBProjectSettings


def test_openblas_openmp_warning_amplification_is_counted_from_stderr(
    tmp_path,
):
    stderr = tmp_path / "xtb.err"
    warning = (
        "OpenBLAS Warning: Detect OpenMP Loop and this application may hang."
    )
    stderr.write_text(f"{warning}\nordinary diagnostic\n{warning}\n")

    assert _openblas_openmp_warning_count(stderr) == 2


def test_run_and_sub_expose_the_same_xtb_command_tree():
    assert run.commands["xtb"] is xtb
    assert sub.commands["xtb"] is xtb
    assert tuple(sorted(run.commands["xtb"].commands)) == (
        "hess",
        "opt",
        "sp",
    )
    assert tuple(sorted(sub.commands["xtb"].commands)) == (
        "hess",
        "opt",
        "sp",
    )


def test_xtb_capability_is_optional_project_cpu_only():
    capability = program_capability("xtb")
    assert capability is not None
    assert capability.requires_project_configuration is False
    assert capability.supports_project_configuration is True
    assert capability.jobtypes == ("hess", "opt", "sp")
    assert capability.engines == ("cpu",)


def test_registry_exposes_loader_bounded_parameter_domains():
    pyscf = program_capability("pyscf")
    xtb_capability = program_capability("xtb")

    assert dict(pyscf.project_parameter_domains)["ab_initio"] == ("hf",)
    assert dict(pyscf.project_parameter_domains)["opt_solver"] == (
        "ase",
        "berny",
        "geometric",
    )
    assert dict(xtb_capability.project_parameter_domains)["gfn_version"] == (
        "gfn0",
        "gfn1",
        "gfn2",
        "gfnff",
    )


def test_xtb_agent_audit_rejects_status_only_receipt(tmp_path):
    receipt_path = tmp_path / "forged.xtb-result-receipt.json"
    finalize_receipt(
        receipt_path,
        {
            "engine_state": "completed",
            "validation_state": "validated",
            "ready": True,
        },
    )
    _, findings = audit_xtb_result_receipt(
        receipt_path,
        expected_jobtype="sp",
        expected_charge=0,
        expected_multiplicity=1,
    )
    assert "xtb.result.schema_version_mismatch" in findings
    assert "xtb.result.requested_settings_missing" in findings
    assert "xtb.result.environment_receipt_unavailable" in findings


def _artifact_manifest_record(path):
    return {
        "path": path.name,
        "size": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def _write_auditable_scratch_xtb_receipt(tmp_path, *, keep_original=False):
    receipt_dir = tmp_path / "durable"
    receipt_dir.mkdir()
    scratch_dir = tmp_path / "scratch"
    scratch_dir.mkdir()
    original_input = scratch_dir / "case.xyz"
    original_input.write_text("2\nhydrogen\nH 0 0 0\nH 0 0 0.7\n")
    durable_input = receipt_dir / original_input.name
    durable_input.write_bytes(original_input.read_bytes())
    output = receipt_dir / "case.out"
    output.write_text("archived xTB output\n")
    executable = tmp_path / "xtb-6.7.1"
    executable.write_bytes(b"synthetic pinned executable identity\n")
    executable.chmod(0o700)
    environment_path = receipt_dir / "case.xtb-environment-receipt.json"
    environment = finalize_receipt(
        environment_path,
        {
            "schema_version": "chemsmart.xtb-environment.v1",
            "required_version": "6.7.1",
            "status": "available",
            "preflight_state": "ready",
            "execution_ready": True,
            "observed_version": "6.7.1",
            "executable": str(executable.resolve()),
            "executable_sha256": sha256_file(executable),
            "findings": [],
        },
    )
    settings = XTBJobSettings(jobtype="sp")
    requested = {
        name: getattr(settings, name) for name in sorted(settings.FIELDS)
    }
    molecule = {
        "symbols": ["H", "H"],
        "positions_angstrom": [[0.0, 0.0, 0.0], [0.0, 0.0, 0.7]],
        "charge": 0,
        "multiplicity": 1,
    }
    execution_input = bind_xtb_execution_input(original_input)
    provenance = {
        "source_artifact": None,
        "project_artifact": None,
        "execution_input_artifact": execution_input,
    }
    command = [str(executable.resolve()), str(original_input.resolve())]
    receipt_path = receipt_dir / "case.xtb-result-receipt.json"
    receipt = finalize_receipt(
        receipt_path,
        {
            "schema_version": "chemsmart.xtb-result-validation.v1",
            "state": "validated",
            "engine_state": "completed",
            "validation_state": "validated",
            "engine_completion": {
                "state": "completed",
                "returncode": 0,
                "normal_termination": True,
            },
            "program": "xtb",
            "jobtype": "sp",
            "ready": True,
            "returncode": 0,
            "canonical_argv": command,
            "command_sha256": canonical_sha256(command),
            "requested_settings": requested,
            "requested_settings_sha256": canonical_sha256(requested),
            "requested_molecule": molecule,
            "requested_molecule_sha256": canonical_sha256(molecule),
            "command_settings": requested,
            "command_settings_sha256": canonical_sha256(requested),
            "applied_settings": requested,
            "applied_settings_sha256": canonical_sha256(requested),
            "results": {"final_energy_hartree": -1.0},
            **provenance,
            "provenance_binding_sha256": canonical_sha256(provenance),
            "environment_receipt_sha256": environment["receipt_sha256"],
            "artifacts": {
                durable_input.name: _artifact_manifest_record(durable_input),
                output.name: _artifact_manifest_record(output),
            },
            "findings": [],
        },
    )
    if not keep_original:
        original_input.unlink()
    return receipt_path, receipt, original_input, durable_input


def _audit_synthetic_xtb_receipt(receipt_path):
    return audit_xtb_result_receipt(
        receipt_path,
        expected_jobtype="sp",
        expected_charge=0,
        expected_multiplicity=1,
        expected_settings={
            name: getattr(XTBJobSettings(jobtype="sp"), name)
            for name in sorted(XTBJobSettings.FIELDS)
        },
    )


def test_xtb_agent_audit_accepts_exact_durable_input_after_scratch_cleanup(
    tmp_path,
):
    receipt_path, receipt, original_input, _ = (
        _write_auditable_scratch_xtb_receipt(tmp_path)
    )

    observation, findings = _audit_synthetic_xtb_receipt(receipt_path)

    assert findings == (), observation
    assert receipt["canonical_argv"][1] == str(original_input.resolve())
    assert receipt["execution_input_artifact"]["declared_path"] == str(
        original_input.resolve()
    )


def test_xtb_agent_audit_rejects_missing_durable_input_fallback(tmp_path):
    receipt_path, _, _, durable_input = _write_auditable_scratch_xtb_receipt(
        tmp_path
    )
    durable_input.unlink()

    _, findings = _audit_synthetic_xtb_receipt(receipt_path)

    assert "xtb.provenance.artifact_mismatch" in findings
    assert f"xtb.result.artifact.{durable_input.name}_mismatch" in findings


def test_xtb_agent_audit_rejects_mutated_durable_input_fallback(tmp_path):
    receipt_path, _, _, durable_input = _write_auditable_scratch_xtb_receipt(
        tmp_path
    )
    durable_input.write_text("mutated durable input\n")

    _, findings = _audit_synthetic_xtb_receipt(receipt_path)

    assert "xtb.provenance.artifact_mismatch" in findings
    assert f"xtb.result.artifact.{durable_input.name}_mismatch" in findings


def test_xtb_agent_audit_never_falls_back_for_existing_mutated_original(
    tmp_path,
):
    receipt_path, _, original_input, durable_input = (
        _write_auditable_scratch_xtb_receipt(tmp_path, keep_original=True)
    )
    original_input.write_text("mutated original input\n")

    _, findings = _audit_synthetic_xtb_receipt(receipt_path)

    assert "xtb.provenance.artifact_mismatch" in findings
    assert f"xtb.result.artifact.{durable_input.name}_mismatch" not in findings


def test_xtb_executable_prefers_configured_folder(tmp_path):
    binary = tmp_path / "xtb"
    binary.write_text("placeholder")
    binary.chmod(0o700)
    executable = XTBExecutable(executable_folder=str(tmp_path), local_run=True)
    assert executable.get_executable() == str(binary)


def _invoke_xtb_and_capture_settings(job_path, args):
    mock_job = MagicMock()
    with patch(job_path, return_value=mock_job) as mock_job_cls:
        result = CliRunner().invoke(
            xtb,
            args,
            obj={"jobrunner": MagicMock()},
            catch_exceptions=False,
        )
        settings = mock_job_cls.call_args.kwargs["settings"]
    return result, settings


class TestXTBSettingsContract:
    def test_exact_job_and_method_space(self):
        assert XTBJobSettings.JOBTYPES == ("sp", "opt", "hess")
        assert XTBJobSettings.GFN_VERSIONS == (
            "gfn0",
            "gfn1",
            "gfn2",
            "gfnff",
        )

    def test_unknown_key_is_rejected(self):
        with pytest.raises(ValueError, match="Unknown xTB setting"):
            XTBJobSettings.from_dict({"gfn_version": "gfn2", "typo": 1})

    def test_unknown_merge_key_is_rejected(self):
        with pytest.raises(ValueError, match="Unknown xTB merge"):
            XTBJobSettings.default().merge({"charg": -1}, keywords=("charg",))

    @pytest.mark.parametrize("charge", [True, 0.0, "0"])
    def test_non_integer_charge_is_rejected(self, charge):
        with pytest.raises((TypeError, ValueError)):
            XTBJobSettings(charge=charge)

    @pytest.mark.parametrize("multiplicity", [True, 0, -1, 1.0, "1"])
    def test_invalid_multiplicity_is_rejected(self, multiplicity):
        with pytest.raises((TypeError, ValueError)):
            XTBJobSettings(multiplicity=multiplicity)

    def test_half_specified_solvent_is_rejected(self):
        with pytest.raises(ValueError, match="requires solvent_model"):
            XTBJobSettings(solvent_model="alpb")
        with pytest.raises(ValueError, match="requires solvent_model"):
            XTBJobSettings(solvent_id="water")

    def test_numeric_cosmo_identifier_is_preserved_as_argv_text(self):
        settings = XTBJobSettings(
            jobtype="sp", solvent_model="cosmo", solvent_id=78.4
        )
        assert settings.solvent_id == "78.4"
        assert XTBJobRunner._settings_args(settings)[-2:] == [
            "--cosmo",
            "78.4",
        ]

    def test_unknown_textual_solvent_is_retained_for_preflight(self):
        settings = XTBJobSettings(
            jobtype="sp",
            solvent_model="cpcmx",
            solvent_id="custom-solvent",
        )
        capability = settings.solvent_capability()
        assert capability.status == "unknown"
        assert settings.requires_solvent_clarification is True
        assert XTBJobRunner._settings_args(settings)[-2:] == [
            "--cpcmx",
            "custom-solvent",
        ]

    def test_confirmed_unsupported_gbsa_method_is_rejected(self):
        with pytest.raises(ValueError, match="Confirmed unsupported"):
            XTBJobSettings(
                jobtype="sp",
                gfn_version="gfn0",
                solvent_model="gbsa",
                solvent_id="water",
            )

    @pytest.mark.parametrize("jobtype", ["opt", "hess"])
    def test_grad_is_rejected_outside_single_point(self, jobtype):
        with pytest.raises(ValueError, match="supported only for the sp"):
            XTBJobSettings(jobtype=jobtype, grad=True)

    def test_title_is_not_an_execution_setting(self):
        with pytest.raises(ValueError, match="Unknown xTB setting"):
            XTBJobSettings.from_dict({"jobtype": "sp", "title": "ignored"})

    @pytest.mark.parametrize(
        "solvent_id",
        ["", "   ", "--md", "-input", "water\n--md", "water\tinput", "a\x00b"],
    )
    def test_solvent_id_rejects_empty_option_like_or_control_data(
        self, solvent_id
    ):
        with pytest.raises(ValueError):
            XTBJobSettings(
                jobtype="sp",
                solvent_model="alpb",
                solvent_id=solvent_id,
            )

    @pytest.mark.parametrize(
        ("jobtype", "expected"),
        [
            ("sp", ["--gfn", "2", "--chrg", "0", "--uhf", "0"]),
            (
                "opt",
                [
                    "--gfn",
                    "2",
                    "--opt",
                    "vtight",
                    "--chrg",
                    "0",
                    "--uhf",
                    "0",
                ],
            ),
            (
                "hess",
                ["--gfn", "2", "--hess", "--chrg", "0", "--uhf", "0"],
            ),
        ],
    )
    def test_exact_sp_opt_hess_argv(self, jobtype, expected):
        assert (
            XTBJobRunner._settings_args(XTBJobSettings(jobtype=jobtype))
            == expected
        )

    def test_gfnff_has_dedicated_flag(self):
        settings = XTBJobSettings(jobtype="sp", gfn_version="gfnff")
        assert XTBJobRunner._settings_args(settings)[:1] == ["--gfnff"]

    def test_electron_parity_is_checked_before_job_execution(
        self, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        water = Molecule(
            symbols=["O", "H", "H"],
            positions=[[0.0, 0.0, 0.0], [0.8, 0.0, 0.5], [-0.8, 0.0, 0.5]],
        )
        with pytest.raises(ValueError, match="parity"):
            XTBSinglePointJob(
                molecule=water,
                settings=XTBJobSettings(jobtype="sp", multiplicity=2),
                label="invalid_water_doublet",
            )

    def test_job_copy_uses_resolved_state_without_mutating_source(
        self, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        source = Molecule(
            symbols=["O", "H", "H"],
            positions=[[0.0, 0.0, 0.0], [0.8, 0.0, 0.5], [-0.8, 0.0, 0.5]],
            charge=-1,
            multiplicity=2,
        )
        job = XTBSinglePointJob(
            molecule=source,
            settings=XTBJobSettings(jobtype="sp", charge=0, multiplicity=1),
            label="resolved_water",
        )
        assert (job.molecule.charge, job.molecule.multiplicity) == (0, 1)
        assert (source.charge, source.multiplicity) == (-1, 2)


class TestXTBProjectSchema:
    def _write(self, tmp_path, payload):
        path = tmp_path / "project.yaml"
        path.write_text(yaml.safe_dump(payload))
        return path

    def test_project_is_optional_and_defaults_to_gfn2(self):
        settings = XTBProjectSettings.from_project(None).opt_settings()
        assert settings.gfn_version == "gfn2"
        assert settings.jobtype == "opt"
        assert settings.optimization_level == "vtight"

    def test_distributed_template_does_not_force_neutral_singlet(self):
        template = (
            Path(__file__).resolve().parents[1]
            / "chemsmart"
            / "settings"
            / "templates"
            / "template_xtb_simple.yaml"
        )
        payload = yaml.safe_load(template.read_text())
        for section in ("sp", "opt", "hess"):
            assert "charge" not in payload[section]
            assert "multiplicity" not in payload[section]

    def test_explicit_project_path_loads(self, tmp_path):
        path = self._write(
            tmp_path,
            {"sp": {"gfn_version": "gfn1", "charge": -1, "multiplicity": 2}},
        )
        project = XTBProjectSettings.from_project(path)
        settings = project.sp_settings()
        assert settings.gfn_version == "gfn1"
        assert settings.charge == -1
        assert settings.multiplicity == 2
        assert project.source_file == str(path.resolve())

    def test_unknown_top_level_section_is_rejected(self, tmp_path):
        path = self._write(tmp_path, {"md": {"gfn_version": "gfn2"}})
        with pytest.raises(ValueError, match="Unknown xTB project section"):
            XTBProjectSettings.from_project(path)

    def test_section_jobtype_mismatch_is_rejected(self, tmp_path):
        path = self._write(tmp_path, {"sp": {"jobtype": "hess"}})
        with pytest.raises(ValueError, match="incompatible jobtype"):
            XTBProjectSettings.from_project(path)

    def test_unknown_section_setting_is_rejected(self, tmp_path):
        path = self._write(tmp_path, {"sp": {"gf_version": "gfn2"}})
        with pytest.raises(ValueError, match="Unknown xTB setting"):
            XTBProjectSettings.from_project(path)

    @pytest.mark.parametrize(
        "payload",
        [
            {"sp": {"charge": -1}},
            {"sp": {"multiplicity": 2}},
        ],
    )
    def test_project_electronic_state_must_be_declared_together(
        self, tmp_path, payload
    ):
        path = self._write(tmp_path, payload)
        with pytest.raises(
            ValueError, match="charge and multiplicity together"
        ):
            XTBProjectSettings.from_project(path)

    @pytest.mark.parametrize("section", ["sp", "hess"])
    def test_optimization_level_is_scoped_to_opt(self, tmp_path, section):
        path = self._write(
            tmp_path, {section: {"optimization_level": "tight"}}
        )
        with pytest.raises(ValueError, match="valid only in the opt"):
            XTBProjectSettings.from_project(path)


class TestXTBCLI:
    def _state_project(self, tmp_path):
        path = tmp_path / "charged-doublet.yaml"
        path.write_text(
            yaml.safe_dump(
                {
                    "sp": {
                        "gfn_version": "gfn2",
                        "charge": -1,
                        "multiplicity": 1,
                    }
                }
            )
        )
        return path

    def test_yaml_charge_and_multiplicity_survive_omitted_cli_state(
        self, tmp_path, single_molecule_xyz_file
    ):
        project = self._state_project(tmp_path)
        result, settings = _invoke_xtb_and_capture_settings(
            "chemsmart.jobs.xtb.singlepoint.XTBSinglePointJob",
            ["-p", str(project), "-f", single_molecule_xyz_file, "sp"],
        )
        assert result.exit_code == 0, result.output
        assert (settings.charge, settings.multiplicity) == (-1, 1)

    def test_explicit_cli_state_overrides_yaml_state(
        self, tmp_path, single_molecule_xyz_file
    ):
        project = self._state_project(tmp_path)
        result, settings = _invoke_xtb_and_capture_settings(
            "chemsmart.jobs.xtb.singlepoint.XTBSinglePointJob",
            [
                "-p",
                str(project),
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "2",
                "sp",
            ],
        )
        assert result.exit_code == 0, result.output
        assert (settings.charge, settings.multiplicity) == (0, 2)

    def test_source_state_survives_default_project(
        self, tmp_path, single_molecule_xyz_file
    ):
        source = Molecule(
            symbols=["O", "H", "H"],
            positions=[
                [0.0, 0.0, 0.0],
                [0.8, 0.0, 0.5],
                [-0.8, 0.0, 0.5],
            ],
            charge=-1,
            multiplicity=2,
        )
        with patch.object(
            Molecule,
            "from_filepath",
            return_value=[source],
        ):
            result, settings = _invoke_xtb_and_capture_settings(
                "chemsmart.jobs.xtb.singlepoint.XTBSinglePointJob",
                ["-f", single_molecule_xyz_file, "sp"],
            )
        assert result.exit_code == 0, result.output
        assert (settings.charge, settings.multiplicity) == (-1, 2)

    def test_sp_preserves_scientific_overrides(self, single_molecule_xyz_file):
        result, settings = _invoke_xtb_and_capture_settings(
            "chemsmart.jobs.xtb.singlepoint.XTBSinglePointJob",
            [
                "-p",
                "test",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "-1",
                "-m",
                "1",
                "-g",
                "gfn1",
                "-sm",
                "alpb",
                "-si",
                "water",
                "--grad",
                "sp",
            ],
        )
        assert result.exit_code == 0, result.output
        assert settings.jobtype == "sp"
        assert settings.charge == -1
        assert settings.multiplicity == 1
        assert settings.gfn_version == "gfn1"
        assert settings.solvent_model == "alpb"
        assert settings.solvent_id == "water"
        assert settings.grad is True

    def test_opt_override_is_scoped_to_opt(self, single_molecule_xyz_file):
        result, settings = _invoke_xtb_and_capture_settings(
            "chemsmart.jobs.xtb.opt.XTBOptJob",
            [
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "2",
                "opt",
                "--optimization-level",
                "loose",
            ],
        )
        assert result.exit_code == 0, result.output
        assert settings.jobtype == "opt"
        assert settings.optimization_level == "loose"

    def test_partial_solvent_override_fails_before_job_creation(
        self, single_molecule_xyz_file
    ):
        result = CliRunner().invoke(
            xtb,
            [
                "-f",
                single_molecule_xyz_file,
                "-sm",
                "alpb",
                "sp",
            ],
            obj={"jobrunner": MagicMock()},
        )
        assert result.exit_code != 0
        assert "must be supplied together" in result.output

    def test_leaf_help_does_not_require_coordinates(self):
        result = CliRunner().invoke(xtb, ["opt", "--help"])
        assert result.exit_code == 0, result.output


class TestXTBRunnerSafety:
    def test_stale_output_is_rejected(self, tmp_path):
        (tmp_path / "previous.out").write_text("old")
        with pytest.raises(RuntimeError, match="stale output artifacts"):
            XTBJobRunner._reject_stale_artifacts(tmp_path)

    def test_fresh_run_workspace_preserves_old_outputs(
        self, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        runner = XTBJobRunner(server=Server("cpu"), scratch=False)
        molecule = Molecule(
            symbols=["H", "H"],
            positions=[[0.0, 0.0, 0.0], [0.0, 0.0, 0.7]],
        )
        job = XTBSinglePointJob(
            molecule=molecule,
            settings=XTBJobSettings(jobtype="sp"),
            label="h2_sp",
            jobrunner=runner,
        )
        old_output = Path(job.execution_root) / "previous.out"
        old_output.write_text("preserved")
        workspace = runner._prepare_job_workspace(job)
        XTBJobRunner._reject_stale_artifacts(workspace)
        assert old_output.read_text() == "preserved"
        assert Path(job.folder) == workspace
        assert workspace.parent.name == ".chemsmart-runs"

    def test_unrelated_xyz_is_not_treated_as_output(self, tmp_path):
        (tmp_path / "input.xyz").write_text("0\n\n")
        XTBJobRunner._reject_stale_artifacts(tmp_path)

    @pytest.mark.parametrize(
        "basename",
        [
            "xtbopt.xyz",
            "xtbopt.sdf",
            "xtbopt.pdb",
            "xtbopt.coord",
            "xtbopt.poscar",
            "xtbopt.gen",
            "xtbopt.EIn",
        ],
    )
    def test_every_parser_selected_optimized_geometry_is_stale(
        self, tmp_path, basename
    ):
        (tmp_path / basename).write_text("old")
        with pytest.raises(RuntimeError, match=basename):
            XTBJobRunner._reject_stale_artifacts(tmp_path)

    @pytest.mark.parametrize(
        "basename",
        [
            ".xtboptok",
            "charges",
            "energy",
            "job.engrad",
            "g98.out",
            "gradient",
            "hessian",
            "vibspectrum",
            "wbo",
            "xtbopt.log",
            "xtbrestart",
            "xtbtopo.mol",
        ],
    )
    def test_every_parser_selected_auxiliary_artifact_is_stale(
        self, tmp_path, basename
    ):
        (tmp_path / basename).write_text("old")
        with pytest.raises(RuntimeError):
            XTBJobRunner._reject_stale_artifacts(tmp_path)

    def test_explicit_gpu_request_is_rejected(self):
        with pytest.raises(ValueError, match="CPU-only"):
            XTBJobRunner(
                server=Server("cpu"),
                scratch=False,
                num_gpus=1,
            )

    def test_gpu_host_inventory_does_not_become_an_xtb_gpu_request(self):
        runner = XTBJobRunner(
            server=Server("gpu", NUM_GPUS=1),
            scratch=False,
        )
        assert runner.detected_server_num_gpus == 1
        assert runner.num_gpus == 0

    def test_job_boundary_rejects_gpu_runner(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        water = Molecule(
            symbols=["O", "H", "H"],
            positions=[[0.0, 0.0, 0.0], [0.8, 0.0, 0.5], [-0.8, 0.0, 0.5]],
        )
        with pytest.raises(ValueError, match="CPU-only"):
            XTBSinglePointJob(
                molecule=water,
                settings=XTBJobSettings(jobtype="sp"),
                label="gpu_water",
                jobrunner=SimpleNamespace(num_gpus=1),
            )

    def test_create_process_passes_argv_without_shell(self, tmp_path):
        runner = object.__new__(XTBJobRunner)
        runner.job_outputfile = str(tmp_path / "job.out")
        runner.job_errfile = str(tmp_path / "job.err")
        runner.running_directory = str(tmp_path)
        with patch("subprocess.Popen") as popen:
            runner._create_process(
                MagicMock(),
                ["xtb", "input.xyz", "--gfn", "2"],
                {"PATH": "/bin"},
            )
        assert popen.call_args.args[0] == [
            "xtb",
            "input.xyz",
            "--gfn",
            "2",
        ]
        assert "shell" not in popen.call_args.kwargs

    def test_execution_environment_binds_omp_to_requested_cores(self):
        runner = object.__new__(XTBJobRunner)
        runner.num_cores = 4
        with patch.object(
            XTBJobRunner,
            "_update_os_environ",
            return_value={"PATH": "/bin", "OMP_NUM_THREADS": "99"},
        ):
            env = runner._execution_environment(MagicMock())
        assert env["OMP_NUM_THREADS"] == "4"
        assert env["MKL_NUM_THREADS"] == "4"
        assert env["OPENBLAS_NUM_THREADS"] == "4"

    def test_nonzero_child_status_is_preserved(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        molecule = Molecule(
            symbols=["H", "H"],
            positions=[[0.0, 0.0, 0.0], [0.0, 0.0, 0.7]],
        )
        job = XTBSinglePointJob(
            molecule=molecule,
            settings=XTBJobSettings(jobtype="sp"),
            label="failed_h2",
        )
        runner = object.__new__(XTBJobRunner)
        runner.num_cores = 2
        runner.num_gpus = 0
        runner.mem_gb = 1
        runner.scratch = False
        runner.delete_scratch = False
        execution_input = tmp_path / "execution-input.xyz"
        execution_input.write_text("2\nh2\nH 0 0 0\nH 0 0 0.7\n")
        runner.job_xyzfile = str(execution_input)
        command = ["/opt/xtb", job.inputfile, "--gfn", "2"]
        green_environment = {
            "execution_ready": True,
            "receipt_sha256": "environment-digest",
            "findings": [],
        }
        with (
            patch.object(runner, "_prerun"),
            patch.object(runner, "_write_input"),
            patch.object(runner, "_get_command", return_value=command),
            patch.object(
                runner,
                "_execution_environment",
                return_value={"OMP_NUM_THREADS": "2"},
            ),
            patch.object(runner, "_create_process", return_value=MagicMock()),
            patch.object(runner, "_run", return_value=23),
            patch.object(runner, "_postrun"),
            patch(
                "chemsmart.jobs.xtb.runner.probe_xtb_environment",
                return_value=green_environment,
            ),
        ):
            with pytest.raises(subprocess.CalledProcessError) as error:
                runner.run(job)
        assert error.value.returncode == 23
        receipt = json.loads(Path(job.result_receiptfile).read_text())
        assert receipt["returncode"] == 23
        assert receipt["engine_state"] == "failed"
        assert receipt["validation_state"] == "failed"
        assert receipt["ready"] is False


class TestXTBSubmissionReconstruction:
    def _job(self, tmp_path, monkeypatch):
        source = tmp_path / "source.xyz"
        source.write_text("2\nsource\nH 0 0 0\nH 0 0 0.7\n")
        monkeypatch.chdir(tmp_path)
        molecule = Molecule.from_filepath(source)
        return XTBSinglePointJob(
            molecule=molecule,
            settings=XTBJobSettings(jobtype="sp"),
            label="h2_sp",
            source_filename=source,
            project_reference=None,
            source_index=1,
        )

    def test_reconstructed_submission_binds_source_index_and_label(
        self, tmp_path, monkeypatch
    ):
        job = self._job(tmp_path, monkeypatch)
        args = job.rewrite_submission_cli_args(
            [
                "--num-cores",
                "4",
                "xtb",
                "-f",
                "source.xyz",
                "--append-label",
                "trial",
                "sp",
            ]
        )
        assert args[args.index("--filename") + 1] == str(
            (tmp_path / "source.xyz").resolve()
        )
        assert args[args.index("--index") + 1] == "1"
        assert args[args.index("--label") + 1] == "h2_sp"
        assert args[args.index("--num-cores") + 1] == "4"
        assert "--append-label" not in args
        assert job.submission_execution_cwd == str(tmp_path.resolve())

    def test_runscript_is_written_in_job_folder_and_restores_parent_cwd(
        self, tmp_path, monkeypatch
    ):
        job = self._job(tmp_path, monkeypatch)
        server = Server(
            "slurm",
            SCHEDULER="SLURM",
            NUM_CORES=4,
            NUM_GPUS=0,
        )
        submitter = SLURMSubmitter(job=job, server=server)
        submitter._write_runscript(
            ["xtb", "--filename", job.source_filename, "sp"]
        )
        script = Path(job.folder) / submitter.run_script
        contents = script.read_text()
        assert script.is_file()
        assert f"os.chdir({str(tmp_path.resolve())!r})" in contents
        assert "except subprocess.CalledProcessError" in contents

    def test_runscript_maps_child_failure_to_exact_process_exit(
        self, tmp_path
    ):
        script = tmp_path / "run.py"
        RunScript(script, ["xtb", "sp"], execution_cwd=tmp_path).write()
        contents = script.read_text()
        assert "sys.exit(error.returncode)" in contents


class TestXTBPreviewAndValidationReceipts:
    def _molecule(self):
        return Molecule(
            symbols=["O", "H", "H"],
            positions=[
                [0.0, 0.0, 0.0],
                [0.8, 0.0, 0.5],
                [-0.8, 0.0, 0.5],
            ],
        )

    def test_fake_preview_cannot_alias_or_complete_real_job(
        self, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        runner = FakeXTBJobRunner(
            server=Server("preview", NUM_CORES=4, NUM_GPUS=0),
            num_cores=4,
            num_gpus=0,
        )
        job = XTBOptJob(
            molecule=self._molecule(),
            settings=XTBJobSettings(jobtype="opt"),
            label="water_opt",
            jobrunner=runner,
        )
        assert runner.run(job) == 0
        receipt = json.loads(Path(job.preview_receiptfile).read_text())
        assert receipt["state"] == "previewed"
        assert receipt["executed"] is False
        assert receipt["execution_ready"] is False
        assert not Path(job.outputfile).exists()
        assert not (Path(job.folder) / "xtbopt.xyz").exists()
        assert job.is_complete() is False

    def test_environment_receipt_requires_exact_671_and_core_parity(
        self, tmp_path
    ):
        binary = tmp_path / "xtb"
        binary.write_text("binary bytes")
        binary.chmod(0o700)
        completed = subprocess.CompletedProcess(
            [str(binary), "--version"],
            0,
            stdout="xtb version 6.7.1\n",
            stderr="",
        )
        receipt_path = tmp_path / "environment.json"
        with patch(
            "chemsmart.jobs.xtb.validation.subprocess.run",
            return_value=completed,
        ):
            receipt = probe_xtb_environment(
                executable=binary,
                num_cores=4,
                num_gpus=0,
                mem_gb=2,
                env={"OMP_NUM_THREADS": "4"},
                receipt_path=receipt_path,
            )
        assert receipt["execution_ready"] is True
        assert receipt["observed_version"] == "6.7.1"
        assert len(receipt["executable_sha256"]) == 64

    def test_environment_receipt_rejects_other_xtb_version(self, tmp_path):
        binary = tmp_path / "xtb"
        binary.write_text("binary bytes")
        binary.chmod(0o700)
        completed = subprocess.CompletedProcess(
            [str(binary), "--version"],
            0,
            stdout="xtb version 6.7.0\n",
            stderr="",
        )
        with patch(
            "chemsmart.jobs.xtb.validation.subprocess.run",
            return_value=completed,
        ):
            receipt = probe_xtb_environment(
                executable=binary,
                num_cores=2,
                num_gpus=0,
                mem_gb=1,
                env={"OMP_NUM_THREADS": "2"},
                receipt_path=tmp_path / "wrong-version.json",
            )
        assert receipt["execution_ready"] is False
        assert {item["rule_id"] for item in receipt["findings"]} == {
            "xtb.environment.version_mismatch"
        }

    def test_unknown_solvent_is_localized_as_scientific_preflight(
        self, tmp_path
    ):
        binary = tmp_path / "xtb"
        binary.write_text("binary bytes")
        binary.chmod(0o700)
        completed = subprocess.CompletedProcess(
            [str(binary), "--version"],
            0,
            stdout="xtb version 6.7.1\n",
            stderr="",
        )
        settings = XTBJobSettings(
            jobtype="sp",
            solvent_model="cpcmx",
            solvent_id="custom-solvent",
        )
        with patch(
            "chemsmart.jobs.xtb.validation.subprocess.run",
            return_value=completed,
        ):
            receipt = probe_xtb_environment(
                executable=binary,
                num_cores=2,
                num_gpus=0,
                mem_gb=1,
                env={"OMP_NUM_THREADS": "2"},
                receipt_path=tmp_path / "unknown-solvent.json",
                settings=settings,
            )
        assert receipt["status"] == "available"
        assert receipt["execution_ready"] is False
        assert receipt["preflight_state"] == "needs_clarification"
        assert {
            item["rule_id"] for item in receipt["scientific_findings"]
        } == {"xtb.solvent.capability_unknown"}

    def test_environment_receipt_rejects_omp_core_mismatch(self, tmp_path):
        binary = tmp_path / "xtb"
        binary.write_text("binary bytes")
        binary.chmod(0o700)
        receipt = probe_xtb_environment(
            executable=binary,
            num_cores=4,
            num_gpus=0,
            mem_gb=1,
            env={"OMP_NUM_THREADS": "1"},
            receipt_path=tmp_path / "wrong-cores.json",
        )
        assert receipt["execution_ready"] is False
        assert "xtb.resources.omp_mismatch" in {
            item["rule_id"] for item in receipt["findings"]
        }

    @pytest.mark.parametrize(
        (
            "fixture_name",
            "job_cls",
            "jobtype",
            "settings_kwargs",
            "source_name",
        ),
        [
            (
                "p_benzyne_sp_alpb_toluene",
                XTBSinglePointJob,
                "sp",
                {
                    "multiplicity": 3,
                    "solvent_model": "alpb",
                    "solvent_id": "toluene",
                },
                "p_benzyne.xyz",
            ),
            (
                "p_benzyne_opt_alpb_toluene",
                XTBOptJob,
                "opt",
                {
                    "multiplicity": 3,
                    "optimization_level": "loose",
                    "solvent_model": "alpb",
                    "solvent_id": "toluene",
                },
                "p_benzyne.xyz",
            ),
            (
                "acetaldehyde_hess",
                XTBHessJob,
                "hess",
                {},
                "acetaldehyde.xyz",
            ),
        ],
    )
    def test_archived_sp_opt_hess_require_green_bound_receipt(
        self,
        tmp_path,
        monkeypatch,
        fixture_name,
        job_cls,
        jobtype,
        settings_kwargs,
        source_name,
    ):
        fixture = (
            Path(__file__).parent
            / "data"
            / "XTBTests"
            / "outputs"
            / fixture_name
        )
        source = fixture / source_name
        molecule = Molecule.from_filepath(source)
        if "multiplicity" in settings_kwargs:
            molecule.multiplicity = settings_kwargs["multiplicity"]
        monkeypatch.chdir(tmp_path)
        settings = XTBJobSettings(jobtype=jobtype, **settings_kwargs)
        job = job_cls(
            molecule=molecule,
            settings=settings,
            label=fixture_name,
            source_filename=source,
        )
        shutil.copytree(fixture, job.folder, dirs_exist_ok=True)
        rendered_input = Path(job.inputfile)
        rendered_input.write_bytes(source.read_bytes())
        copied_source = Path(job.folder) / source_name
        if copied_source != rendered_input:
            copied_source.unlink()

        if jobtype == "opt":
            output_text = Path(job.outputfile).read_text()
            Path(job.outputfile).write_text(
                output_text.replace(" --grad --json", " --json")
            )

        binary = tmp_path / "xtb-6.7.1"
        binary.write_bytes(b"synthetic pinned executable identity\n")
        binary.chmod(0o700)
        environment = finalize_receipt(
            job.environment_receiptfile,
            {
                "schema_version": "chemsmart.xtb-environment.v1",
                "required_version": "6.7.1",
                "status": "available",
                "preflight_state": "ready",
                "execution_ready": True,
                "observed_version": "6.7.1",
                "executable": str(binary.resolve()),
                "executable_sha256": sha256_file(binary),
                "findings": [],
            },
        )
        command = [str(binary), job.inputfile]
        command.extend(XTBJobRunner._settings_args(settings))
        provenance = dict(job.declared_provenance_binding)
        provenance["execution_input_artifact"] = bind_xtb_execution_input(
            job.inputfile
        )
        receipt = validate_xtb_result(
            job=job,
            command=command,
            environment_receipt=environment,
            provenance_binding=provenance,
            returncode=0,
            receipt_path=job.result_receiptfile,
        )
        assert receipt["ready"] is True, receipt["findings"]
        assert receipt["engine_state"] == "completed"
        assert receipt["validation_state"] == "validated"
        observation, audit_findings = audit_xtb_result_receipt(
            job.result_receiptfile,
            expected_jobtype=jobtype,
            expected_charge=settings.charge,
            expected_multiplicity=settings.multiplicity,
            expected_settings={
                name: getattr(settings, name)
                for name in sorted(settings.FIELDS)
            },
            expected_source_sha256=sha256_file(source),
        )
        assert audit_findings == (), observation
        assert math.isfinite(observation["final_energy_hartree"])
        assert observation["final_energy_evidence"] in {
            "main_output",
            "receipt",
        }
        assert job.is_complete() is True

        if jobtype == "sp":
            Path(job.errfile).write_text(
                "OpenBLAS Warning: Detect OpenMP Loop and this application "
                "may hang.\n"
            )
            warned = validate_xtb_result(
                job=job,
                command=command,
                environment_receipt=environment,
                provenance_binding=provenance,
                returncode=0,
                receipt_path=job.result_receiptfile,
            )
            finding = next(
                item
                for item in warned["warnings"]
                if item["rule_id"] == "xtb.environment.openblas_openmp_warning"
            )
            assert warned["ready"] is True
            assert warned["validation_state"] == "validated"
            assert warned["findings"] == []
            assert finding["observed"] == 1

        with open(job.outputfile, "a") as handle:
            handle.write("mutated\n")
        assert job.is_complete() is False


class TestXTBExamples:
    def test_example_xyz_files_exist(self):
        examples = Path(__file__).resolve().parents[1] / "examples" / "xtb"
        for filename in ("water.xyz", "methane.xyz"):
            assert (examples / filename).stat().st_size > 0
