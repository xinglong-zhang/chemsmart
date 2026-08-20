"""Focused tests for the versioned PySCF provenance artifact.

These tests construct temporary HDF5 files only. No PySCF calculation or
optional chemistry solver is executed.
"""

import json

import h5py
import numpy as np
import pytest

from chemsmart.io.pyscf.output import PySCFOutput, read_pyscf_h5
from chemsmart.jobs.pyscf.writer import (
    H5_NULL_ATTRIBUTE,
    RESULT_UNITS,
    RESULTS_SCHEMA_VERSION,
    PySCFScriptWriter,
    write_pyscf_h5,
)


def _v2_payload():
    """Return a representative B2 payload with exact numeric arrays."""
    spec = {
        "label": "water_hess",
        "program": "pyscf",
        "jobtype": "hess",
        "engine": "cpu",
        "xc": "b3lyp",
        "method": "b3lyp",
        "basis": "def2-svp",
        "charge": 0,
        "spin": 0,
        "multiplicity": 1,
        "solvent_model": "smd",
        "solvent_id": "water",
        "density_fit": True,
        "aux_basis": "def2-universal-jkfit",
        "defgrid": "defgrid3",
        "scf_tol": 1.0e-9,
        "opt_solver": "geometric",
        "opt_maxsteps": 100,
        "stages": ["scf", "hess"],
        "symbols": ["O", "H", "H"],
        "settings_digest": "settings-sha256",
    }
    provenance = {
        "pyscf_version": "2.14.0",
        "gpu4pyscf_version": None,
        "libxc_version": "7.0.0",
        "chemsmart_version": "3.1.4",
        "numpy_version": "1.26.4",
        "h5py_version": "3.16.0",
        "python_version": "3.10.20",
        "interpreter": "/opt/compute/bin/python",
        "host": "compute-node",
        "engine": "cpu",
        "num_threads": 2,
        "cuda_visible_devices": "",
        "started_at": "2026-08-02T12:00:00+0000",
        "ended_at": "2026-08-02T12:00:02+0000",
        "wall_seconds": 2.0,
        "core_seconds": 4.0,
        "settings_digest": "settings-sha256",
        "project_yaml_digest": "project-sha256",
    }
    status = {
        "engine_complete": True,
        "normal_termination": True,
        "stages": {
            "scf": {"converged": True, "iterations": 7},
            "hess": {"converged": True},
        },
        "failure": None,
    }
    results = {
        "energies": np.asarray([-76.3], dtype=np.float64),
        "positions": np.asarray(
            [[0.0, 0.0, 0.0], [0.7, 0.5, 0.0], [-0.7, 0.5, 0.0]],
            dtype=np.float64,
        ),
        "forces": np.asarray(
            [[0.0, 0.0, 0.0], [1.0e-5, 0.0, 0.0], [-1.0e-5, 0.0, 0.0]],
            dtype=np.float64,
        ),
        "mo_energy": np.asarray([-0.5, 0.1], dtype=np.float64),
        "mo_occ": np.asarray([2.0, 0.0], dtype=np.float64),
        "vibrational_frequencies": np.asarray(
            [1600.0, 3900.0, 4010.0], dtype=np.float64
        ),
        "normal_modes": np.arange(27, dtype=np.float64).reshape(3, 3, 3),
        "mulliken_charges": np.asarray([-0.3, 0.15, 0.15], dtype=np.float64),
        "dipole_moment": np.asarray([0.0, 1.8, 0.0], dtype=np.float64),
        "point_group": "C2v",
    }
    return spec, provenance, status, results


def test_v2_groups_round_trip_typed_status_nulls_and_exact_arrays(tmp_path):
    path = tmp_path / "label.h5"
    spec, provenance, status, results = _v2_payload()

    write_pyscf_h5(
        path,
        spec=spec,
        provenance=provenance,
        status=status,
        results=results,
    )

    with h5py.File(path, "r") as handle:
        assert handle.attrs["schema_version"] == RESULTS_SCHEMA_VERSION
        for name in ("spec", "provenance", "status", "results"):
            assert isinstance(handle[name], h5py.Group)
        null = handle["provenance/gpu4pyscf_version"]
        assert null.shape == (0,)
        assert bool(null.attrs[H5_NULL_ATTRIBUTE]) is True
        for name, node in handle["results"].items():
            if isinstance(node, h5py.Dataset) and node.dtype.kind in {
                "b",
                "i",
                "u",
                "f",
                "c",
            }:
                assert node.attrs["unit"] == RESULT_UNITS[name]

    observed_spec, observed_provenance, observed_status, observed_results = (
        read_pyscf_h5(path)
    )
    assert observed_spec == spec
    assert observed_provenance == provenance
    assert observed_status == status
    for key, expected in results.items():
        if isinstance(expected, str):
            assert observed_results[key] == expected
        else:
            np.testing.assert_array_equal(observed_results[key], expected)
            assert observed_results[key].dtype == expected.dtype

    output = PySCFOutput(tmp_path / "label.out")
    assert output.normal_termination is True
    assert output.engine_complete is True
    assert output.failure is None
    assert output.engine == "cpu"
    assert output.settings_digest == "settings-sha256"
    assert output.project_yaml_digest == "project-sha256"
    assert output.point_group == "C2v"
    np.testing.assert_array_equal(output.forces, results["forces"])


def test_spin_diagnostic_datasets_have_explicit_dimensionless_units(tmp_path):
    path = tmp_path / "spin-diagnostic.h5"
    spec, provenance, status, results = _v2_payload()
    results["spin_square"] = np.asarray(2.0)
    results["spin_square_effective_multiplicity"] = np.asarray(3.0)

    write_pyscf_h5(
        path,
        spec=spec,
        provenance=provenance,
        status=status,
        results=results,
    )

    with h5py.File(path, "r") as handle:
        assert handle["results/spin_square"].attrs["unit"] == "dimensionless"
        assert (
            handle["results/spin_square_effective_multiplicity"].attrs["unit"]
            == "dimensionless"
        )


def test_generated_standalone_encoder_writes_same_v2_contract(tmp_path):
    """Exercise only the generated script's encoder, never its chemistry."""
    source = PySCFScriptWriter.render(
        {"schema_version": RESULTS_SCHEMA_VERSION, "label": "encoder-only"}
    )
    namespace = {"__name__": "chemsmart_generated_encoder_test"}
    exec(compile(source, "<generated-pyscf-driver>", "exec"), namespace)

    path = tmp_path / "generated.h5"
    spec, provenance, status, results = _v2_payload()
    namespace["_write_h5"](path, spec, provenance, status, results)

    observed = read_pyscf_h5(path)
    assert observed[:3] == (spec, provenance, status)
    for key, expected in results.items():
        if isinstance(expected, str):
            assert observed[3][key] == expected
        else:
            np.testing.assert_array_equal(observed[3][key], expected)

    output = PySCFOutput(path)
    assert output.result_units["results/energies"] == "Eh"
    assert output.result_units["results/positions"] == "Angstrom"
    assert output.result_units["results/vibrational_frequencies"] == "cm^-1"


def test_generated_driver_preserves_imaginary_modes_as_negative_values():
    source = PySCFScriptWriter.render(
        {"schema_version": RESULTS_SCHEMA_VERSION, "label": "source-only"}
    )

    assert "imaginary_freq=False" in source


def test_solvent_resolution_is_deferred_to_target_interpreter():
    source = PySCFScriptWriter.render(
        {"schema_version": RESULTS_SCHEMA_VERSION, "label": "source-only"}
    )

    assert "from pyscf.solvent.smd import solvent_db" in source
    assert 'config["solvent_eps"] = eps' in source
    assert '"source": "pyscf.solvent.smd.solvent_db"' in source
    assert '"unit": "dimensionless_relative_permittivity"' in source
    assert '"pyscf_version": pyscf.__version__' in source
    assert '"environment_receipt_sha256": os.environ.get(' in source
    assert "CHEMSMART_PYSCF_ENVIRONMENT_RECEIPT_SHA256" in source
    assert "mf.with_solvent.lebedev_order" in source


def test_optimizer_probe_accepts_expected_sentinel_rejection_shapes():
    from chemsmart.jobs.pyscf.environment import _PROBE_SCRIPT

    # PySCF adapters reject the non-computing sentinel at different API
    # boundaries: geomeTRIC uses NotImplementedError, Berny commonly reaches
    # a missing Mole attribute, and ASE raises RuntimeError.  All demonstrate
    # that the optional dependency and real adapter entry point were reached.
    assert "AttributeError," in _PROBE_SCRIPT
    assert "NotImplementedError," in _PROBE_SCRIPT
    assert "RuntimeError," in _PROBE_SCRIPT
    assert "TypeError," in _PROBE_SCRIPT
    assert '"call_probe": "entrypoint_reached"' in _PROBE_SCRIPT


def test_gpu_solvent_quadrature_is_explicit_and_receipted():
    from chemsmart.jobs.pyscf.settings import PySCFJobSettings

    settings = PySCFJobSettings(
        jobtype="sp",
        functional="pbe0",
        basis="def2-svp",
        engine="gpu",
        charge=0,
        multiplicity=1,
        solvent_model="smd",
        solvent_id="water",
    )
    molecule = type(
        "Molecule",
        (),
        {
            "chemical_symbols": ["H", "H"],
            "positions": [[0.0, 0.0, 0.0], [0.0, 0.0, 0.7]],
            "charge": 0,
            "multiplicity": 1,
        },
    )()
    job = type(
        "Job",
        (),
        {
            "settings": settings,
            "molecule": molecule,
            "jobrunner": type(
                "Runner",
                (),
                {
                    "num_cores": 1,
                    "mem_gb": 1,
                    "_run_id": "run-gpu-solvent",
                    "_run_nonce": "nonce-gpu-solvent",
                },
            )(),
            "label": "h2_gpu_smd",
            "TYPE": "pyscf_sp",
            "stages": ("scf",),
        },
    )()

    config = PySCFScriptWriter(job).build_config()

    assert config["solvent_lebedev_order"] == 29


def test_non_opt_spec_does_not_claim_optimizer_settings():
    settings = type(
        "Settings",
        (),
        {
            "jobtype": "sp",
            "title": "Hydrogen receipt",
            "functional": "pbe0",
            "ab_initio": None,
            "basis": "def2-svp",
            "aux_basis": None,
            "defgrid": None,
            "dispersion": None,
            "density_fit": True,
            "scf_tol": None,
            "scf_maxiter": None,
            "solvent_model": None,
            "solvent_id": None,
            "opt_solver": "geometric",
            "opt_maxsteps": 100,
            "engine": "cpu",
            "charge": 0,
            "multiplicity": 1,
            "xc": "pbe0",
            "method_name": "pbe0",
            "project_yaml_digest": None,
            "validate": lambda self: self,
        },
    )()
    molecule = type(
        "Molecule",
        (),
        {
            "chemical_symbols": ["H", "H"],
            "positions": [[0.0, 0.0, 0.0], [0.0, 0.0, 0.7]],
            "charge": 0,
            "multiplicity": 1,
        },
    )()
    job = type(
        "Job",
        (),
        {
            "settings": settings,
            "molecule": molecule,
            "jobrunner": type(
                "Runner",
                (),
                {
                    "num_cores": 1,
                    "mem_gb": None,
                    "_run_id": "run-config",
                    "_run_nonce": "nonce-config",
                },
            )(),
            "label": "h2_sp",
            "TYPE": "pyscf_sp",
            "stages": ("scf",),
        },
    )()

    config = PySCFScriptWriter(job).build_config()

    assert config["opt_solver"] is None
    assert config["opt_maxsteps"] is None
    assert config["title"] == "Hydrogen receipt"
    assert config["run_id"] == "run-config"
    assert config["run_nonce"] == "nonce-config"


@pytest.mark.parametrize(
    ("field", "expected", "injected"),
    [
        ("xc", "b3lyp", "pbe"),
        ("spin", 0, 2),
        ("solvent_id", "water", "ethanol"),
    ],
)
def test_v2_round_trip_does_not_mask_injected_spec_mismatch(
    tmp_path, field, expected, injected
):
    """An applied-setting mismatch must survive exact HDF5 round-trip."""
    path = tmp_path / f"wrong-{field}.h5"
    spec, provenance, status, results = _v2_payload()
    assert spec[field] == expected
    spec[field] = injected

    write_pyscf_h5(
        path,
        spec=spec,
        provenance=provenance,
        status=status,
        results=results,
    )

    observed_spec, _, _, _ = read_pyscf_h5(path)
    assert observed_spec[field] == injected
    assert observed_spec[field] != expected


def test_v1_json_payload_remains_readable(tmp_path):
    path = tmp_path / "legacy.h5"
    spec = {
        "label": "legacy",
        "basis": "def2-svp",
        "point_group": "C1",
        "stages": ["scf"],
    }
    provenance = {
        "engine": "cpu",
        "settings_digest": "legacy-digest",
    }
    status = {
        "normal_termination": False,
        "stages": {"scf": {"converged": False, "cycles": 50}},
        "failure": {"type": "RuntimeError", "message": "not converged"},
    }
    energy = np.asarray([-1.25], dtype=np.float64)

    with h5py.File(path, "w") as handle:
        handle.attrs["schema_version"] = "1.0"
        handle.create_dataset("spec", data=json.dumps(spec))
        handle.create_dataset("provenance", data=json.dumps(provenance))
        handle.create_dataset("status", data=json.dumps(status))
        handle.create_group("results").create_dataset("energies", data=energy)

    observed_spec, observed_provenance, observed_status, observed_results = (
        read_pyscf_h5(path)
    )
    assert observed_spec == spec
    assert observed_provenance == provenance
    assert observed_status == status
    np.testing.assert_array_equal(observed_results["energies"], energy)


def test_schema_version_layout_mismatch_is_rejected(tmp_path):
    path = tmp_path / "wrong-layout.h5"
    with h5py.File(path, "w") as handle:
        handle.attrs["schema_version"] = RESULTS_SCHEMA_VERSION
        handle.create_dataset("spec", data="{}")
        handle.create_group("provenance")
        handle.create_group("status")
        handle.create_group("results")

    with pytest.raises(ValueError, match=r"2\.0 requires /spec to be a group"):
        read_pyscf_h5(path)


def test_settings_digest_ignores_provenance_only_fields_but_not_engine():
    config = {
        "schema_version": "2.0",
        "chemsmart_version": "3.1.4",
        "label": "label-a",
        "program": "pyscf",
        "jobtype": "sp",
        "engine": "cpu",
        "stages": ["scf"],
        "symbols": ["H", "H"],
        "positions": [[0.0, 0.0, 0.0], [0.0, 0.0, 0.7]],
        "basis": "def2-svp",
        "xc": "b3lyp",
        "num_threads": 2,
        "max_memory_mb": 1024.0,
        "project_yaml_digest": "project-a",
    }
    baseline = PySCFScriptWriter.settings_digest(config)

    provenance_change = dict(
        config,
        run_id="another-run",
        run_nonce="another-nonce",
        schema_version="3.0",
        chemsmart_version="3.2.0",
        label="label-b",
        num_threads=8,
        max_memory_mb=4096.0,
        project_yaml_digest="project-b",
        input_artifact_kind="pyscf_hdf5",
        input_artifact_sha256="source-artifact-b",
    )
    assert PySCFScriptWriter.settings_digest(provenance_change) == baseline

    engine_change = dict(config, engine="gpu")
    assert PySCFScriptWriter.settings_digest(engine_change) != baseline
