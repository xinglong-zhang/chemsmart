"""Generate the standalone PySCF driver script for a job.

The emitted ``label.py`` is a **fixed skeleton plus one configuration dict**.
The skeleton runs the ground-state SCF, an optimisation (on the ground
state, on an excited root, or on a correlated surface), a TDA/TDDFT
response stage, a correlated-energy stage and an analytic Hessian, as the
resolved stages ask. The driver logic is invariant, so
nothing is gained by templating source text -- and free-form generation adds a
quoting, escaping and injection surface for no benefit. Only ``CONFIG``
varies, and it holds scalars, strings and a geometry array.

The script imports ``pyscf``, ``numpy``, ``h5py`` and the standard library
**only**. It must never import ``chemsmart``: that would drag ChemSmart's
``numpy~=1.26.4`` pin and the rdkit/pymatgen/ase tree into every compute
process, and would make the per-program ``CONDA_ENV`` -- the mechanism that
lets a GPU job run in its own environment -- unusable.
"""

import hashlib
import json
import logging
import math
import os
from numbers import Integral, Real

from ase.data import atomic_numbers as ASE_ATOMIC_NUMBERS

from chemsmart import __version__ as chemsmart_version
from chemsmart.io.pyscf.output import pyscf_source_artifact_binding
from chemsmart.jobs.pyscf.settings import (
    PYSCF_DEFGRIDS,
    PYSCF_SOLVENT_MODELS,
    PYSCF_UNRESTRICTED_MANIFOLD,
)

logger = logging.getLogger(__name__)

#: Bumped when the layout of ``label.h5`` changes in a non-additive way.
#: Version 1.0 stored ``spec``, ``provenance`` and ``status`` as JSON string
#: datasets. Version 2.0 makes all four sections real HDF5 groups so an agent
#: can inspect individual fields without first decoding an opaque blob.
RESULTS_SCHEMA_VERSION = "2.0"
LEGACY_RESULTS_SCHEMA_VERSION = "1.0"

#: Additive contract marker carried inside ``spec``.  The HDF5 container
#: remains schema 2.0 so historical artifacts stay readable; this marker
#: identifies records that satisfy the stricter state, status, and runtime
#: reference checks required for new execution/data-edge admission.
RESULT_CONTRACT_VERSION = "chemsmart.pyscf-result-contract.v5"
#: Contract versions this ChemSmart still reads as executed evidence.  v3 and
#: v4 share one applied-spec vocabulary; v4 adds datasets (forces at the
#: Hessian geometry, spin populations) and status facts (the mass convention
#: behind the frequencies, the symmetry tolerance behind the point group, the
#: optimiser's criteria, the continuation SCF) without touching any v3 field,
#: so a v3 artifact stays analysis-ready and stays a legal geometry source.
#: v5 makes the response stage executable, adds the excited-state and
#: correlated-method datasets and status blocks, and widens the applied-spec
#: vocabulary by four fields -- under its own version only, because the
#: applied-settings digest of every archived artifact is reconstructed from
#: the vocabulary of *that* artifact's contract.
PREVIOUS_RESULT_CONTRACT_VERSIONS = (
    "chemsmart.pyscf-result-contract.v3",
    "chemsmart.pyscf-result-contract.v4",
)
SUPPORTED_RESULT_CONTRACT_VERSIONS = PREVIOUS_RESULT_CONTRACT_VERSIONS + (
    RESULT_CONTRACT_VERSION,
)
TD_RESPONSE_MATERIALIZATION_SCHEMA_VERSION = (
    "chemsmart.pyscf-td-response-materialization.v1"
)

#: Marker used for an explicit JSON-style null. HDF5 has no native null scalar,
#: so a marked empty uint8 dataset distinguishes "known to be unavailable" from
#: a missing field.
H5_NULL_ATTRIBUTE = "chemsmart_is_null"

#: Fields copied from the resolved driver configuration into ``spec/``.
#: Fake and real artifacts use the same applied-settings vocabulary.
LEGACY_APPLIED_SPEC_FIELDS = (
    "run_id",
    "run_nonce",
    "label",
    "title",
    "program",
    "jobtype",
    "engine",
    "stages",
    "symbols",
    "positions",
    "unit",
    "xc",
    "ab_initio",
    "method",
    "basis",
    "charge",
    "spin",
    "multiplicity",
    "dispersion",
    "density_fit",
    "aux_basis",
    "defgrid",
    "atom_grid",
    "scf_tol",
    "scf_maxiter",
    "solvent_model",
    "solvent_call",
    "solvent_method",
    "solvent_id",
    "solvent_eps",
    "solvent_lebedev_order",
    "opt_solver",
    "opt_maxsteps",
    "response_method",
    "state_manifold",
    "nstates",
    "preview_only",
    "materializations",
    "num_threads",
    "max_memory_mb",
    "input_geometry_sha256",
    "input_artifact_kind",
    "input_artifact_sha256",
    "requested_settings_sha256",
    "applied_settings_sha256",
    "settings_digest",
)

#: The v3/v4 vocabulary, frozen: the digest stored in every archived v3/v4
#: artifact was computed over exactly these fields.
APPLIED_SPEC_FIELDS_V4 = LEGACY_APPLIED_SPEC_FIELDS + (
    "result_contract_version",
    "reference_family",
)

#: The current (v5) vocabulary: the excited-root, frozen-core and iteration
#: controls the driver now applies.
APPLIED_SPEC_FIELDS = APPLIED_SPEC_FIELDS_V4 + (
    "excited_state_root",
    "frozen_core",
    "td_max_cycle",
    "cc_max_cycle",
)

#: Digest vocabulary per contract version.  Extending the current tuple in
#: place would silently change the reconstruction for every archived
#: artifact, so each version names its own.
APPLIED_SPEC_FIELDS_BY_CONTRACT = {
    "chemsmart.pyscf-result-contract.v3": APPLIED_SPEC_FIELDS_V4,
    "chemsmart.pyscf-result-contract.v4": APPLIED_SPEC_FIELDS_V4,
    RESULT_CONTRACT_VERSION: APPLIED_SPEC_FIELDS,
}

# Units are part of the machine contract rather than prose in a reader.  Every
# numeric dataset the fixed driver can emit has one explicit unit, including
# quantities that are dimensionless.
RESULT_UNITS = {
    "atomic_numbers": "dimensionless",
    "ccsd_correlation_energy": "Eh",
    "correlation_energy": "Eh",
    "dipole_moment": "Debye",
    "energies": "Eh",
    "excitation_energies": "Eh",
    "excited_state_converged": "dimensionless",
    "excited_state_multiplicities": "dimensionless",
    "force_constants": "Dyne/Angstrom",
    "forces": "Eh/Bohr",
    "hessian": "Eh/Bohr^2",
    "mo_energy": "Eh",
    "mo_occ": "electron",
    "mulliken_charges": "elementary_charge",
    "mulliken_spin_populations": "electron",
    "normal_modes": "atomic_mass_unit^-1/2",
    "oscillator_strengths": "dimensionless",
    "positions": "Angstrom",
    "reduced_masses": "atomic_mass_unit",
    "reference_energy": "Eh",
    "scf_energy": "Eh",
    "spin_square": "dimensionless",
    "spin_square_effective_multiplicity": "dimensionless",
    "total_energy": "Eh",
    "transition_dipole_moments": "Debye",
    "triples_correction": "Eh",
    "vibrational_frequencies": "cm^-1",
}


def applied_pyscf_spec(config):
    """Return the resolved fields written under ``spec/``."""
    return {key: config.get(key) for key in APPLIED_SPEC_FIELDS}


def applied_pyscf_spec_fields(spec):
    """Return the digest vocabulary of *this* artifact's contract version.

    A v4 artifact's ``applied_settings_sha256`` was computed over the v4
    fields; reconstructing it from the v5 vocabulary would mark every
    archived artifact tampered.  A contract-less artifact keeps the legacy
    vocabulary.
    """

    version = spec.get("result_contract_version")
    if version in APPLIED_SPEC_FIELDS_BY_CONTRACT:
        return APPLIED_SPEC_FIELDS_BY_CONTRACT[version]
    return LEGACY_APPLIED_SPEC_FIELDS


def pyscf_td_response_materialization(settings, *, reference_family=None):
    """Map typed TDA/TDDFT intent to the PySCF response construction.

    The manifest names the exact classes the driver instantiates, so the
    review packet, the artifact and the validator agree on the response
    object by construction.  It covers the ``td`` jobtype and an ``opt``
    on an excited root; ``reference_family`` (``rks``/``uks``) selects the
    restricted or unrestricted factory.
    """

    jobtype = str(getattr(settings, "jobtype", "")).strip().lower()
    excited_root = getattr(settings, "excited_state_root", None)
    if jobtype != "td" and not (jobtype == "opt" and excited_root):
        return None
    response_method = str(settings.response_method).strip().lower()
    manifold = str(settings.state_manifold).strip().lower()
    family = str(reference_family or "rks").strip().lower()
    module = "uks" if family == "uks" else "rks"
    factory_api = {
        "tda": f"pyscf.tdscf.{module}.TDA",
        "tddft": f"pyscf.tdscf.{module}.TDDFT",
    }[response_method]
    operations = ["ground_state_scf", "response_construct"]
    if manifold != PYSCF_UNRESTRICTED_MANIFOLD:
        operations.append(f"set_{manifold}_channel")
    operations.append("set_nstates")
    if getattr(settings, "td_max_cycle", None) is not None:
        operations.append("set_max_cycle")
    if excited_root:
        operations.extend(
            (
                "gradient_scanner_on_root",
                "geometry_optimisation",
                "final_scf",
            )
        )
    operations.append("vertical_excitation_kernel")
    body = {
        "schema_version": TD_RESPONSE_MATERIALIZATION_SCHEMA_VERSION,
        "ground_state_reference_family": module,
        "ground_state_stage": "scf",
        "response_method": response_method,
        "response_factory_api": factory_api,
        "state_manifold": manifold,
        "nstates": int(settings.nstates),
        "excited_state_root": (
            int(excited_root) if excited_root is not None else None
        ),
        "max_cycle": (
            int(settings.td_max_cycle)
            if getattr(settings, "td_max_cycle", None) is not None
            else None
        ),
        "operation_order": tuple(operations),
        "execution_policy": "executable",
    }
    payload = json.dumps(body, sort_keys=True, separators=(",", ":"))
    return {
        **body,
        "receipt_sha256": hashlib.sha256(payload.encode("utf-8")).hexdigest(),
    }


def pyscf_reference_family(*, symbols, charge, multiplicity, xc):
    """Derive the PySCF reference family from host-owned molecular state.

    ``pyscf.scf.HF`` selects its concrete class from electron count and spin;
    in particular, the one-electron open-shell path is ``HF1e``/ROHF rather
    than UHF.  Recording that decision before script generation lets result
    validation detect a runtime reference substitution independently.
    """

    if isinstance(charge, bool) or not isinstance(charge, Integral):
        raise ValueError("charge must be an integer, not a boolean")
    if isinstance(multiplicity, bool) or not isinstance(
        multiplicity, Integral
    ):
        raise ValueError("multiplicity must be an integer, not a boolean")
    try:
        electron_count = sum(
            int(ASE_ATOMIC_NUMBERS[str(symbol)]) for symbol in symbols
        ) - int(charge)
    except (KeyError, TypeError) as exc:
        raise ValueError("symbols must resolve to atomic numbers") from exc
    spin = int(multiplicity) - 1
    if electron_count < 0 or spin < 0 or spin > electron_count:
        raise ValueError("electron count and multiplicity are inconsistent")
    if (electron_count - spin) % 2:
        raise ValueError("electron count and multiplicity have invalid parity")
    if xc is not None:
        return "rks" if spin == 0 else "uks"
    if spin == 0:
        return "rhf"
    if electron_count == 1:
        return "rohf"
    return "uhf"


def write_pyscf_h5(
    filename,
    *,
    spec,
    provenance,
    status,
    results,
    schema_version=RESULTS_SCHEMA_VERSION,
):
    """Write a versioned, machine-operable PySCF results artifact.

    This public helper is also the contract for synthetic/fake runners. Numeric
    result arrays are passed directly to h5py, preserving their dtype and shape;
    mappings are represented as nested groups and scalar datasets. Unsupported
    Python objects raise rather than being silently converted to strings.
    """
    import h5py

    with h5py.File(filename, "w") as handle:
        handle.attrs["schema_version"] = schema_version
        spec_group = handle.create_group("spec")
        _write_mapping(spec_group, spec)
        _attach_supplied_geometry_unit(spec_group, spec)
        _write_mapping(handle.create_group("provenance"), provenance)
        _write_mapping(handle.create_group("status"), status)
        results_group = handle.create_group("results")
        _write_mapping(results_group, results)
        _attach_result_units(results_group)


def _attach_supplied_geometry_unit(spec_group, spec):
    """State the unit of the supplied geometry on the dataset itself.

    ``spec/positions`` is the structure the calculation was handed and
    ``spec/unit`` its unit; a reader serving the supplied structure as a
    typed quantity must find the unit where every other numeric dataset
    carries it, not in a sibling field.
    """

    import h5py

    node = spec_group.get("positions")
    if not isinstance(node, h5py.Dataset):
        return
    if bool(node.attrs.get(H5_NULL_ATTRIBUTE, False)):
        return
    if node.dtype.kind not in {"f", "i", "u"}:
        return
    unit = spec.get("unit")
    if isinstance(unit, str) and unit:
        node.attrs["unit"] = unit


def _attach_result_units(group):
    """Attach mandatory units to every numeric result dataset."""

    import h5py

    for key, node in group.items():
        if isinstance(node, h5py.Group):
            _attach_result_units(node)
            continue
        if bool(node.attrs.get(H5_NULL_ATTRIBUTE, False)):
            continue
        if node.dtype.kind not in {"b", "i", "u", "f", "c"}:
            continue
        unit = RESULT_UNITS.get(key)
        if unit is None:
            raise ValueError(
                f"Numeric PySCF result {node.name} has no declared unit."
            )
        node.attrs["unit"] = unit


def _write_mapping(group, mapping):
    """Write *mapping* below an open h5py group."""
    for key in sorted(mapping):
        _write_value(group, str(key), mapping[key])


def _write_value(group, key, value):
    """Write one supported scalar, array, null or nested mapping."""
    import h5py
    import numpy as np

    if isinstance(value, dict):
        _write_mapping(group.create_group(key), value)
        return

    if value is None:
        dataset = group.create_dataset(
            key, data=np.empty((0,), dtype=np.uint8)
        )
        dataset.attrs[H5_NULL_ATTRIBUTE] = True
        return

    string_dtype = h5py.string_dtype(encoding="utf-8")
    if isinstance(value, str):
        group.create_dataset(key, data=value, dtype=string_dtype)
        return

    array = np.asarray(value)
    if array.dtype.kind in {"U", "S"}:
        group.create_dataset(
            key, data=array.astype(object), dtype=string_dtype
        )
        return
    if array.dtype.kind == "O":
        flat = array.reshape(-1).tolist()
        if all(isinstance(item, str) for item in flat):
            group.create_dataset(key, data=array, dtype=string_dtype)
            return
        raise TypeError(
            f"Cannot store object-valued field {group.name}/{key} in the "
            "PySCF HDF5 schema"
        )
    group.create_dataset(key, data=value)


class PySCFScriptWriter:
    """Write the generated driver script for a :class:`PySCFJob`."""

    def __init__(self, job):
        self.job = job

    # ------------------------------------------------------------------
    # configuration
    # ------------------------------------------------------------------

    def _solvent_config(self, settings):
        """Resolve the solvent model without importing controller PySCF.

        Returns a ``(call, method, solvent_id, eps)`` tuple.

        The dielectric is deliberately left unresolved here.  The standalone
        child may use another ``CONDA_ENV``/interpreter, so only that target
        environment's solvent database is authoritative.  The generated
        driver resolves PCM dielectric data before the first SCF cycle.
        """
        if settings.solvent_model is None:
            return None, None, None, None

        model = str(settings.solvent_model).lower()
        call, method = PYSCF_SOLVENT_MODELS[model]
        solvent_id = str(settings.solvent_id).strip().lower()
        return call, method, solvent_id, None

    def build_config(self):
        """Return the configuration dict embedded in the generated script."""
        job = self.job
        settings = job.settings
        settings.validate()

        molecule = job.molecule
        input_artifact = pyscf_source_artifact_binding(molecule)
        symbols = list(molecule.chemical_symbols)
        positions = [
            [float(x), float(y), float(z)] for x, y, z in molecule.positions
        ]

        charge = settings.charge
        if charge is None:
            charge = molecule.charge
        multiplicity = settings.multiplicity
        if multiplicity is None:
            multiplicity = molecule.multiplicity

        atom_grid = None
        if settings.xc is not None and settings.defgrid is not None:
            atom_grid = list(PYSCF_DEFGRIDS[str(settings.defgrid).lower()])

        call, method, solvent_id, eps = self._solvent_config(settings)

        jobrunner = job.jobrunner
        num_threads = getattr(jobrunner, "num_cores", None)
        if (
            isinstance(num_threads, bool)
            or not isinstance(num_threads, Integral)
            or int(num_threads) <= 0
        ):
            raise ValueError(
                f"num_cores must be a positive integer, got {num_threads!r}"
            )
        mem_gb = getattr(jobrunner, "mem_gb", None)
        # PySCF's max_memory is per-process in MB and bounds integral
        # buffering, not total allocation. Leave PySCF's own default in
        # place when the runner has no explicit budget.
        if mem_gb is not None and (
            isinstance(mem_gb, bool)
            or not isinstance(mem_gb, Real)
            or not math.isfinite(float(mem_gb))
            or float(mem_gb) <= 0
        ):
            raise ValueError(
                f"mem_gb must be a positive finite scalar, got {mem_gb!r}"
            )
        max_memory_mb = float(mem_gb) * 1024 if mem_gb is not None else None

        config = {
            "schema_version": RESULTS_SCHEMA_VERSION,
            "result_contract_version": RESULT_CONTRACT_VERSION,
            "run_id": getattr(jobrunner, "_run_id", None),
            "run_nonce": getattr(jobrunner, "_run_nonce", None),
            # Captured by the ChemSmart process because the standalone driver
            # deliberately cannot import chemsmart from the compute env.
            "chemsmart_version": chemsmart_version,
            "label": job.label,
            "title": settings.title,
            "program": "pyscf",
            "jobtype": settings.jobtype or job.TYPE,
            "stages": list(job.stages),
            "symbols": symbols,
            "positions": positions,
            "unit": "Angstrom",
            "basis": settings.basis,
            "charge": int(charge),
            "multiplicity": int(multiplicity),
            # PySCF's spin is 2S = Nalpha - Nbeta, NOT the multiplicity.
            "spin": int(multiplicity) - 1,
            "xc": settings.xc,
            "reference_family": pyscf_reference_family(
                symbols=symbols,
                charge=charge,
                multiplicity=multiplicity,
                xc=settings.xc,
            ),
            "ab_initio": settings.ab_initio,
            "method": settings.method_name,
            "dispersion": settings.dispersion,
            "density_fit": bool(settings.density_fit),
            "aux_basis": settings.aux_basis if settings.density_fit else None,
            "defgrid": settings.defgrid if settings.xc is not None else None,
            "atom_grid": atom_grid,
            "scf_tol": settings.scf_tol,
            "scf_maxiter": settings.scf_maxiter,
            "solvent_model": settings.solvent_model,
            "solvent_call": call,
            "solvent_method": method,
            "solvent_id": solvent_id,
            "solvent_eps": eps,
            "solvent_lebedev_order": (
                29 if settings.engine == "gpu" and call is not None else None
            ),
            "opt_solver": (
                settings.opt_solver if "opt" in job.stages else None
            ),
            "opt_maxsteps": (
                settings.opt_maxsteps if "opt" in job.stages else None
            ),
            "response_method": (
                str(getattr(settings, "response_method", None)).strip().lower()
                if getattr(settings, "response_method", None) is not None
                else None
            ),
            "state_manifold": (
                str(getattr(settings, "state_manifold", None)).strip().lower()
                if getattr(settings, "state_manifold", None) is not None
                else None
            ),
            "nstates": (
                int(getattr(settings, "nstates", None))
                if getattr(settings, "nstates", None) is not None
                else None
            ),
            "excited_state_root": (
                int(settings.excited_state_root)
                if getattr(settings, "excited_state_root", None) is not None
                else None
            ),
            "td_max_cycle": (
                int(settings.td_max_cycle)
                if getattr(settings, "td_max_cycle", None) is not None
                else None
            ),
            "frozen_core": (
                None
                if getattr(settings, "frozen_core", None) is None
                else (
                    str(settings.frozen_core).strip().lower()
                    if isinstance(settings.frozen_core, str)
                    else int(settings.frozen_core)
                )
            ),
            "cc_max_cycle": (
                int(settings.cc_max_cycle)
                if getattr(settings, "cc_max_cycle", None) is not None
                else None
            ),
            # Kept in the vocabulary for every artifact on disk; no stage
            # is a preview any more, so it is False for every new run.
            "preview_only": False,
            "materializations": {},
            "engine": settings.engine,
            "num_threads": int(num_threads),
            "max_memory_mb": max_memory_mb,
            # The project-settings builder attaches this only when a concrete
            # source YAML exists. A null records that it was unavailable; it
            # is never inferred from a project name or regenerated content.
            "project_yaml_digest": getattr(
                settings, "project_yaml_digest", None
            ),
            "input_artifact_kind": (
                input_artifact["kind"] if input_artifact else None
            ),
            "input_artifact_sha256": (
                input_artifact["sha256"] if input_artifact else None
            ),
        }
        td_materialization = pyscf_td_response_materialization(
            settings, reference_family=config["reference_family"]
        )
        if td_materialization is not None:
            config["materializations"]["td_response_plan"] = td_materialization
        geometry_payload = {
            "symbols": config["symbols"],
            "positions": config["positions"],
            "unit": config["unit"],
            "charge": config["charge"],
            "multiplicity": config["multiplicity"],
        }
        config["input_geometry_sha256"] = self._json_digest(geometry_payload)
        requested = self.settings_digest(config)
        config["requested_settings_sha256"] = requested
        config["applied_settings_sha256"] = None
        # Historical field retained as the controller-side requested identity.
        config["settings_digest"] = requested
        return config

    @staticmethod
    def _json_digest(value):
        body = json.dumps(value, sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(body.encode("utf-8")).hexdigest()

    @staticmethod
    def settings_digest(config):
        """Return a stable digest of the scientifically meaningful settings.

        Resource knobs, artifact labels, provenance-only software/source
        identifiers and the HDF5 layout version are excluded: they do not
        change the requested calculation. ``engine`` remains included because
        changing CPU/GPU execution is a scientific change and must invalidate
        approval.
        """
        ignored = {
            "schema_version",
            "result_contract_version",
            "run_id",
            "run_nonce",
            "chemsmart_version",
            "label",
            "title",
            "num_threads",
            "max_memory_mb",
            "project_yaml_digest",
            "input_geometry_sha256",
            "input_artifact_kind",
            "input_artifact_sha256",
            "requested_settings_sha256",
            "applied_settings_sha256",
            "settings_digest",
        }
        payload = {k: v for k, v in config.items() if k not in ignored}
        body = json.dumps(payload, sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(body.encode("utf-8")).hexdigest()

    # ------------------------------------------------------------------
    # emission
    # ------------------------------------------------------------------

    def write(self, target_directory=None, config=None):
        """Write ``label.py`` and return its path."""
        directory = target_directory or self.job.folder
        os.makedirs(directory, exist_ok=True)
        path = os.path.join(directory, f"{self.job.label}.py")
        if config is None:
            config = self.build_config()
        with open(path, "w") as handle:
            handle.write(self.render(config))
        logger.debug(f"Wrote PySCF driver script: {path}")
        return path

    @staticmethod
    def render(config):
        """Return the full script text for ``config``.

        The configuration is embedded as a JSON document parsed at run time,
        not as a Python literal. It stays readable and diffable, and it is
        unambiguously *data*: nothing in CONFIG can execute, whatever a label
        or solvent name happens to contain.
        """
        payload = json.dumps(config, indent=4)
        # Well-formed JSON cannot contain three consecutive double quotes --
        # a closing quote is always followed by , : } ] or whitespace -- so
        # the triple-quoted literal below is safe. Assert it anyway.
        assert '"""' not in payload, "JSON payload would break the literal"
        return _SKELETON.replace("__CHEMSMART_CONFIG__", payload)


# The skeleton is invariant. Only the CONFIG literal is substituted, and it
# is emitted as JSON so the file stays readable and diffable.
_SKELETON = '''#!/usr/bin/env python
# ChemSmart PySCF driver
"""PySCF driver generated by ChemSmart.

DO NOT EDIT: this file is regenerated on every run, so edits are lost.
Executable artifacts may be rerun with the bound environment.

Imports pyscf, numpy, h5py and the standard library only -- never chemsmart.
"""

import json
import os
import platform
import socket
import sys
import time
import traceback

import h5py
import numpy as np

CONFIG = json.loads("""__CHEMSMART_CONFIG__""")


def _apply_threads(config):
    """Bound PySCF's thread pool.

    The BLAS/OpenMP pools are bound by environment variables the jobrunner
    sets before this interpreter starts; this only covers PySCF's own pool.
    """
    from pyscf import lib

    lib.num_threads(config["num_threads"])


def _build_mole(config, log_path):
    import pyscf

    atoms = [
        (sym, tuple(pos))
        for sym, pos in zip(config["symbols"], config["positions"])
    ]
    kwargs = dict(
        atom=atoms,
        basis=config["basis"],
        charge=config["charge"],
        # PySCF's spin is 2S = Nalpha - Nbeta, not the multiplicity.
        spin=config["spin"],
        unit=config["unit"],
        output=log_path,
        verbose=4,
    )
    if config["max_memory_mb"]:
        kwargs["max_memory"] = config["max_memory_mb"]
    return pyscf.M(**kwargs)


def _functional_definition(
    libxc,
    xc,
    *,
    pyscf_version,
    environment_receipt_sha256,
):
    """Materialize the target LibXC interpretation of one DFT literal."""
    parser_hybrid, components = libxc.parse_xc(xc)
    ordered_components = sorted(components)
    return {
        "schema_version": "chemsmart.pyscf-functional-definition.v3",
        "field": "xc",
        "source": "pyscf.dft.libxc.parse_xc",
        "source_key": xc,
        "pyscf_version": str(pyscf_version),
        "libxc_version": str(libxc.libxc_version()),
        "environment_receipt_sha256": environment_receipt_sha256,
        # ``parse_xc`` returns parser decomposition metadata.  For compound
        # LibXC aliases such as B3LYPG/PBE0 this tuple may be all zero even
        # though the actual exact-exchange fraction is non-zero, so it must
        # not be described as the physical hybrid coefficient.
        "parser_hybrid_decomposition": [
            float(value) for value in parser_hybrid
        ],
        "exact_exchange_fraction": float(libxc.hybrid_coeff(xc)),
        "range_separation_coefficients": [
            float(value) for value in libxc.rsh_coeff(xc)
        ],
        # Parallel primitive vectors remain lossless while being directly
        # serializable as numeric HDF5 datasets.  A list of dictionaries would
        # become an object array and fail before any real DFT result is saved.
        "functional_ids": [
            int(functional_id) for functional_id, _ in ordered_components
        ],
        "functional_factors": [
            float(factor) for _, factor in ordered_components
        ],
    }


def _build_method(config, mol):
    """Construct the mean-field object.

    Order matters and follows gpu4pyscf's own method_from_config: grids and
    density fitting are applied on the CPU object, `.to_gpu()` comes next,
    and the solvent is attached last. Attaching a solvent before `.to_gpu()`
    leaves a CPU-resident solvent object on a GPU method.
    """
    import pyscf
    from pyscf import dft, scf

    def _is_gpu4pyscf_object(value):
        return any(
            cls.__module__.startswith("gpu4pyscf")
            for cls in type(value).__mro__
        )

    def _solvent_dielectric(solvent_id):
        """Resolve dielectric data in the target PySCF environment."""
        from pyscf.solvent.smd import solvent_db

        key = str(solvent_id).strip().lower()
        if key not in solvent_db:
            raise ValueError(
                "Solvent %r is not in target PySCF solvent_db" % solvent_id
            )
        # PySCF SMD records unpack as
        # n, _, alpha, beta, gamma, eps, phi, psi.
        return float(solvent_db[key][5])

    if config["xc"] is None:
        mf = scf.HF(mol)
    else:
        mf = dft.KS(mol, xc=config["xc"])
        config.setdefault("materializations", {})[
            "functional_definition"
        ] = _functional_definition(
            dft.libxc,
            config["xc"],
            pyscf_version=pyscf.__version__,
            environment_receipt_sha256=os.environ.get(
                "CHEMSMART_PYSCF_ENVIRONMENT_RECEIPT_SHA256"
            ),
        )
        if config["atom_grid"]:
            mf.grids.atom_grid = tuple(config["atom_grid"])

    if config["dispersion"]:
        mf.disp = config["dispersion"]

    if config["scf_tol"] is not None:
        mf.conv_tol = float(config["scf_tol"])
    if config["scf_maxiter"] is not None:
        mf.max_cycle = int(config["scf_maxiter"])

    if config["density_fit"]:
        mf = mf.density_fit(auxbasis=config["aux_basis"])

    if config["engine"] == "gpu":
        mf = mf.to_gpu()
        if not _is_gpu4pyscf_object(mf):
            raise RuntimeError(
                "GPU engine requested but PySCF did not return a "
                "GPU4PySCF mean-field object"
            )

    call = config["solvent_call"]
    if call == "PCM":
        mf = mf.PCM()
        mf.with_solvent.method = config["solvent_method"]
        # PySCF's PCM defaults eps to water regardless of the requested
        # solvent, so it must be set explicitly or every PCM run is aqueous.
        eps = config["solvent_eps"]
        if eps is None:
            eps = _solvent_dielectric(config["solvent_id"])
        mf.with_solvent.eps = eps
        config["solvent_eps"] = eps
    elif call == "SMD":
        mf = mf.SMD()
        mf.with_solvent.solvent = config["solvent_id"]
        # Resolve now as a target-environment validation and applied-spec
        # record.  SMD still performs its own descriptor lookup at kernel
        # time; setting the name, rather than overriding eps, preserves that
        # implementation path.
        config["solvent_eps"] = _solvent_dielectric(config["solvent_id"])
    if call is not None:
        config.setdefault("materializations", {})[
            "solvent_dielectric"
        ] = {
            "schema_version": (
                "chemsmart.pyscf-solvent-materialization.v1"
            ),
            "field": "solvent_eps",
            "source": "pyscf.solvent.smd.solvent_db",
            "source_key": config["solvent_id"],
            "value": config["solvent_eps"],
            "unit": "dimensionless_relative_permittivity",
            "pyscf_version": pyscf.__version__,
            "environment_receipt_sha256": os.environ.get(
                "CHEMSMART_PYSCF_ENVIRONMENT_RECEIPT_SHA256"
            ),
        }
    if call is not None and config["solvent_lebedev_order"] is not None:
        # GPU4PySCF 1.8.0's solvent path is pinned to this quadrature order;
        # record and apply it explicitly instead of relying on a runtime
        # default that may drift between environments.
        mf.with_solvent.lebedev_order = int(
            config["solvent_lebedev_order"]
        )
    if config["engine"] == "gpu" and not _is_gpu4pyscf_object(mf):
        raise RuntimeError(
            "GPU engine requested but a post-construction wrapper returned "
            "a non-GPU4PySCF mean-field object"
        )
    return mf


def _run_opt(config, method):
    """Optimise the geometry and return (converged, optimised Mole).

    ``method`` is the mean-field object for a ground-state surface, or a
    nuclear-gradient scanner for an excited root or a correlated method.
    assert_convergence is disabled so that a failed optimisation still
    produces a results file recording the failure, rather than raising with
    no artifact. Downstream, normal_termination stays False.
    """
    solver = config["opt_solver"]
    maxsteps = config["opt_maxsteps"]
    if solver == "geometric":
        from pyscf.geomopt import geometric_solver

        return geometric_solver.kernel(
            method, assert_convergence=False, maxsteps=maxsteps
        )
    if solver == "berny":
        from pyscf.geomopt import berny_solver

        return berny_solver.kernel(
            method, assert_convergence=False, maxsteps=maxsteps
        )
    if solver == "ase":
        from pyscf.geomopt import ase_solver

        return ase_solver.kernel(method, target="atoms", max_steps=maxsteps)
    raise ValueError("Unknown opt_solver: %s" % solver)


class FollowedRootFiltered(RuntimeError):
    """The followed excited root fell below PySCF's positive-eigenvalue
    filter and vanished from the spectrum; the driver never switches roots.
    """


def _class_name(value):
    return type(value).__module__ + "." + type(value).__qualname__


def _build_response(config, mf):
    """Construct the TDA/TDDFT response object on a converged mean field.

    ``pyscf.tdscf.TDA/TDDFT`` dispatch on the mean-field class (RKS or
    UKS; a PCM-wrapped reference attaches the non-equilibrium response).
    The manifold flag is set for a restricted reference only: an
    unrestricted reference has one spin-conserving manifold and PySCF
    leaves ``singlet`` as None there.
    """
    from pyscf import tdscf

    method = str(config["response_method"]).strip().lower()
    manifold = str(config["state_manifold"]).strip().lower()
    factory = {"tda": tdscf.TDA, "tddft": tdscf.TDDFT}[method]
    td = factory(mf)
    if manifold in ("singlet", "triplet"):
        td.singlet = manifold == "singlet"
    td.nstates = int(config["nstates"])
    if config.get("td_max_cycle") is not None:
        td.max_cycle = int(config["td_max_cycle"])
    return td


def _root_convergence(td):
    """Per-root convergence flags as a 1-d boolean array."""
    flags = np.asarray(td.converged, dtype=bool).reshape(-1)
    obtained = int(np.asarray(td.e).reshape(-1).size)
    if flags.size == 1 and obtained > 1:
        flags = np.repeat(flags, obtained)
    return flags[:obtained]


def _response_solvent_record(mf, td):
    """What the solvent model applied to the response, if any."""
    ground = getattr(mf, "with_solvent", None)
    response = getattr(td, "with_solvent", None)
    if ground is None and response is None:
        return None
    return {
        "equilibrium_solvation": (
            None
            if response is None
            else bool(getattr(response, "equilibrium_solvation", False))
        ),
        "static_eps_applied": (
            None if ground is None else float(getattr(ground, "eps", 0.0))
        ),
        # PySCF applies a fixed optical dielectric to the fast response
        # for every solvent (solvent/_attach_solvent.py); the number it
        # used is recorded here so a non-aqueous run shows the divergence.
        "response_eps_applied": (
            None if response is None else float(getattr(response, "eps", 0.0))
        ),
    }


def _run_td(config, mf, results, status, runtime):
    """Vertical excitations at the mean field's current geometry.

    Roots are ascending within the requested manifold at this geometry.
    PySCF drops eigenvalues below ``positive_eig_threshold`` before
    reporting, so ``nstates_obtained`` may be smaller than the request;
    the count of filtered roots is a stage fact, never a silent shift.
    """
    from pyscf.data import nist

    td = _build_response(config, mf)
    td.kernel()
    runtime["response_class"] = _class_name(td)
    excitations = np.asarray(td.e, dtype=float).reshape(-1)
    obtained = int(excitations.size)
    requested = int(config["nstates"])
    converged = _root_convergence(td)
    manifold = str(config["state_manifold"]).strip().lower()
    results["excitation_energies"] = excitations
    results["excited_state_converged"] = converged
    if manifold in ("singlet", "triplet"):
        results["excited_state_multiplicities"] = np.full(
            obtained, 1 if manifold == "singlet" else 3, dtype=int
        )
    stage = {
        "converged": bool(obtained > 0 and converged.all()),
        "all_converged": bool(obtained > 0 and converged.all()),
        "unconverged_roots": [
            int(index + 1) for index, flag in enumerate(converged) if not flag
        ],
        "nstates_requested": requested,
        "nstates_obtained": obtained,
        "roots_filtered": int(max(requested - obtained, 0)),
        "positive_eig_threshold_applied": float(
            getattr(td, "positive_eig_threshold", float("nan"))
        ),
        "max_cycle_applied": int(td.max_cycle),
        "response_method_applied": str(config["response_method"]),
        "state_manifold_applied": manifold,
        "singlet_flag_applied": (
            None if td.singlet is None else bool(td.singlet)
        ),
        "transition_dipole_au_to_debye": float(nist.AU2DEBYE),
        "solvent": _response_solvent_record(mf, td),
    }
    try:
        results["oscillator_strengths"] = np.asarray(
            td.oscillator_strength(), dtype=float
        ).reshape(-1)[:obtained]
        status["properties"]["oscillator_strengths"] = {"status": "ok"}
    except Exception as exc:
        status["properties"]["oscillator_strengths"] = {
            "status": "unavailable",
            "failure": {"type": type(exc).__name__, "message": str(exc)},
        }
    try:
        dipoles = np.asarray(td.transition_dipole(), dtype=float)
        results["transition_dipole_moments"] = (
            dipoles.reshape(obtained, -1)[:, :3] * float(nist.AU2DEBYE)
        )
        status["properties"]["transition_dipole_moments"] = {"status": "ok"}
    except Exception as exc:
        status["properties"]["transition_dipole_moments"] = {
            "status": "unavailable",
            "failure": {"type": type(exc).__name__, "message": str(exc)},
        }
    status["stages"]["td"] = stage
    return td, stage


def _frozen_core_count(method):
    frozen = getattr(method, "frozen", None)
    if frozen is None:
        return 0
    if isinstance(frozen, (int, np.integer)):
        return int(frozen)
    return int(len(frozen))


def _build_correlated(config, mf):
    """Construct the MP2 or coupled-cluster object on the HF reference."""
    from pyscf import cc, mp

    method = str(config["ab_initio"]).strip().lower()
    if method == "mp2":
        obj = mp.MP2(mf)
    elif method in ("ccsd", "ccsd(t)"):
        obj = cc.CCSD(mf)
        if config.get("cc_max_cycle") is not None:
            obj.max_cycle = int(config["cc_max_cycle"])
    else:
        raise ValueError("Unknown correlated method: %s" % method)
    frozen = config.get("frozen_core")
    if frozen is None:
        pass
    elif isinstance(frozen, str) and frozen.strip().lower() == "auto":
        setter = getattr(obj, "set_frozen", None)
        if not callable(setter):
            raise RuntimeError(
                "this PySCF cannot resolve frozen_core 'auto' for %s"
                % _class_name(obj)
            )
        setter(method="auto")
    else:
        obj.frozen = int(frozen)
    return obj


def _run_corr(config, mf, results, status, runtime):
    """The correlated components at the mean field's current geometry."""
    method = str(config["ab_initio"]).strip().lower()
    obj = _build_correlated(config, mf)
    runtime["correlated_class"] = _class_name(obj)
    obj.kernel()
    reference = float(mf.e_tot)
    correlation = float(obj.e_corr)
    stage = {
        "converged": bool(getattr(obj, "converged", True)),
        "method": method,
        "frozen_core_requested": config.get("frozen_core"),
        "frozen_core_applied": _frozen_core_count(obj),
        "max_cycle_applied": (
            int(obj.max_cycle) if method != "mp2" else None
        ),
        "amplitudes_converged": (
            bool(getattr(obj, "converged", True)) if method != "mp2" else None
        ),
        "triples_computed": False,
    }
    results["reference_energy"] = reference
    if method in ("ccsd", "ccsd(t)"):
        results["ccsd_correlation_energy"] = correlation
    if method == "ccsd(t)":
        triples = float(obj.ccsd_t())
        results["triples_correction"] = triples
        correlation += triples
        stage["triples_computed"] = True
    results["correlation_energy"] = correlation
    status["stages"]["corr"] = stage
    return reference + correlation


def _distribution_versions():
    from importlib import metadata

    prefixes = ("gpu4pyscf", "cupy", "cutensor")
    found = {}
    for distribution in metadata.distributions():
        name = str(distribution.metadata.get("Name") or "").lower()
        if name.startswith(prefixes):
            found[name] = distribution.version
    return dict(sorted(found.items()))


def _runtime_provenance(config, mf):
    evidence = {
        "mean_field_class": (
            type(mf).__module__ + "." + type(mf).__qualname__
        ),
        "mean_field_mro": [
            cls.__module__ + "." + cls.__qualname__
            for cls in type(mf).__mro__
        ],
        "packages": _distribution_versions(),
    }
    gpu_distribution = None
    for name in sorted(evidence["packages"]):
        if name == "gpu4pyscf" or name.startswith("gpu4pyscf-cuda"):
            gpu_distribution = {
                "name": name,
                "version": evidence["packages"][name],
            }
            break
    evidence["gpu4pyscf_distribution"] = gpu_distribution
    if config["engine"] != "gpu":
        return evidence
    try:
        import cupy

        evidence["cupy_version"] = cupy.__version__
        evidence["device_count"] = int(
            cupy.cuda.runtime.getDeviceCount()
        )
        evidence["cuda_driver_version"] = int(
            cupy.cuda.runtime.driverGetVersion()
        )
        evidence["cuda_runtime_version"] = int(
            cupy.cuda.runtime.runtimeGetVersion()
        )
        if evidence["device_count"]:
            properties = cupy.cuda.runtime.getDeviceProperties(0)
            name = properties.get("name")
            if isinstance(name, bytes):
                name = name.decode("utf-8", errors="replace")
            evidence["device_name"] = name
            uuid = properties.get("uuid")
            if isinstance(uuid, bytes):
                uuid = uuid.hex()
            elif isinstance(uuid, (list, tuple)):
                uuid = "".join("%02x" % int(item) for item in uuid)
            evidence["device_uuid"] = uuid
        try:
            from cupy_backends.cuda.libs import cutensor

            evidence["cutensor_runtime"] = int(cutensor.get_version())
        except Exception as exc:
            evidence["cutensor_error_type"] = type(exc).__name__
    except Exception as exc:
        evidence["gpu_runtime_error_type"] = type(exc).__name__
    return evidence


def _provenance(config, started_at, ended_at, wall_seconds, runtime):
    import pyscf

    gpu_distribution = runtime.get("gpu4pyscf_distribution")
    gpu_version = (
        gpu_distribution.get("version")
        if isinstance(gpu_distribution, dict)
        else None
    )

    try:
        from pyscf.dft import libxc

        libxc_version = str(libxc.libxc_version())
    except Exception:
        libxc_version = None

    return {
        "run_id": config.get("run_id"),
        "run_nonce": config.get("run_nonce"),
        "pyscf_version": pyscf.__version__,
        "gpu4pyscf_version": gpu_version,
        "libxc_version": libxc_version,
        "chemsmart_version": config.get("chemsmart_version"),
        "numpy_version": np.__version__,
        "h5py_version": h5py.__version__,
        "python_version": platform.python_version(),
        "interpreter": sys.executable,
        "host": socket.gethostname(),
        "platform": platform.platform(),
        "engine": config["engine"],
        "num_threads": config["num_threads"],
        "cuda_visible_devices": os.environ.get("CUDA_VISIBLE_DEVICES"),
        "omp_num_threads": os.environ.get("OMP_NUM_THREADS"),
        "mkl_num_threads": os.environ.get("MKL_NUM_THREADS"),
        "openblas_num_threads": os.environ.get("OPENBLAS_NUM_THREADS"),
        "settings_digest": config["settings_digest"],
        "requested_settings_sha256": config["requested_settings_sha256"],
        "applied_settings_sha256": config.get("applied_settings_sha256"),
        "input_geometry_sha256": config["input_geometry_sha256"],
        "input_artifact_kind": config.get("input_artifact_kind"),
        "input_artifact_sha256": config.get("input_artifact_sha256"),
        "project_yaml_digest": config.get("project_yaml_digest"),
        "script_sha256": os.environ.get("CHEMSMART_PYSCF_SCRIPT_SHA256"),
        "input_receipt_sha256": os.environ.get(
            "CHEMSMART_PYSCF_INPUT_RECEIPT_SHA256"
        ),
        "environment_receipt_sha256": os.environ.get(
            "CHEMSMART_PYSCF_ENVIRONMENT_RECEIPT_SHA256"
        ),
        "started_at": started_at,
        "ended_at": ended_at,
        "wall_seconds": wall_seconds,
        "core_seconds": wall_seconds * config["num_threads"],
        "core_seconds_kind": "wall_times_configured_threads",
        "runtime": runtime,
    }


def _settings_digest(config):
    ignored = {
        "schema_version", "result_contract_version",
        "run_id", "run_nonce", "chemsmart_version",
        "label", "title", "num_threads", "max_memory_mb",
        "project_yaml_digest", "input_geometry_sha256",
        "input_artifact_kind", "input_artifact_sha256",
        "requested_settings_sha256", "applied_settings_sha256",
        "settings_digest",
    }
    payload = {key: value for key, value in config.items() if key not in ignored}
    body = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    import hashlib

    return hashlib.sha256(body.encode("utf-8")).hexdigest()


_H5_NULL_ATTRIBUTE = "chemsmart_is_null"


def _write_mapping(group, mapping):
    """Write a mapping as nested groups and typed datasets."""
    for key in sorted(mapping):
        _write_value(group, str(key), mapping[key])


def _write_value(group, key, value):
    """Write one scalar, array, explicit null or nested mapping."""
    if isinstance(value, dict):
        _write_mapping(group.create_group(key), value)
        return

    if value is None:
        dataset = group.create_dataset(
            key, data=np.empty((0,), dtype=np.uint8)
        )
        dataset.attrs[_H5_NULL_ATTRIBUTE] = True
        return

    string_dtype = h5py.string_dtype(encoding="utf-8")
    if isinstance(value, str):
        group.create_dataset(key, data=value, dtype=string_dtype)
        return

    array = np.asarray(value)
    if array.dtype.kind in {"U", "S"}:
        group.create_dataset(
            key, data=array.astype(object), dtype=string_dtype
        )
        return
    if array.dtype.kind == "O":
        flat = array.reshape(-1).tolist()
        if all(isinstance(item, str) for item in flat):
            group.create_dataset(key, data=array, dtype=string_dtype)
            return
        raise TypeError(
            "Cannot store object-valued field %s/%s" % (group.name, key)
        )
    group.create_dataset(key, data=value)


def _write_h5(path, spec, provenance, status, results):
    """Write the versioned machine contract."""
    with h5py.File(path, "w") as handle:
        handle.attrs["schema_version"] = CONFIG["schema_version"]
        spec_group = handle.create_group("spec")
        _write_mapping(spec_group, spec)
        positions_node = spec_group.get("positions")
        if isinstance(positions_node, h5py.Dataset) and isinstance(
            spec.get("unit"), str
        ):
            positions_node.attrs["unit"] = spec["unit"]
        _write_mapping(handle.create_group("provenance"), provenance)
        _write_mapping(handle.create_group("status"), status)
        results_group = handle.create_group("results")
        _write_mapping(results_group, results)
        _attach_result_units(results_group)


def _attach_result_units(group):
    """Attach the controller-declared unit to each numeric result dataset."""
    units = __CHEMSMART_RESULT_UNITS__
    for key, node in group.items():
        if isinstance(node, h5py.Group):
            _attach_result_units(node)
            continue
        if bool(node.attrs.get(_H5_NULL_ATTRIBUTE, False)):
            continue
        if node.dtype.kind not in {"b", "i", "u", "f", "c"}:
            continue
        unit = units.get(key)
        if unit is None:
            raise ValueError(
                "Numeric PySCF result %s has no declared unit" % node.name
            )
        node.attrs["unit"] = unit


def _to_host_array(array):
    """Materialize a CPU NumPy array from NumPy or a GPU array wrapper."""
    getter = getattr(array, "get", None)
    return np.asarray(getter() if callable(getter) else array)


def _symmetrize_cartesian_hessian(hessian):
    """Return a symmetric Cartesian Hessian and the raw antisymmetry.

    PySCF stores analytic Hessians as ``(atom, atom, xyz, xyz)``.  Finite
    integration grids can leave the two mixed-derivative evaluations different
    at roundoff-to-quadrature scale.  ``numpy.linalg.eigh`` assumes a symmetric
    matrix and otherwise consumes only one triangle, so make that assumption
    explicit and preserve the size of the correction in the stage status.
    """

    values = np.asarray(hessian, dtype=float)
    if values.ndim == 4:
        natm = values.shape[0]
        expected = (natm, natm, 3, 3)
        if values.shape != expected:
            raise ValueError(
                "PySCF Cartesian Hessian has unexpected shape %r" %
                (values.shape,)
            )
        matrix = values.transpose(0, 2, 1, 3).reshape(3 * natm, 3 * natm)
        maximum_antisymmetry = float(np.max(np.abs(matrix - matrix.T)))
        symmetric = 0.5 * (matrix + matrix.T)
        restored = symmetric.reshape(natm, 3, natm, 3).transpose(0, 2, 1, 3)
        return restored, maximum_antisymmetry
    if values.ndim == 2 and values.shape[0] == values.shape[1]:
        maximum_antisymmetry = float(np.max(np.abs(values - values.T)))
        return 0.5 * (values + values.T), maximum_antisymmetry
    raise ValueError(
        "PySCF Cartesian Hessian must be square or atom-blocked"
    )


def _optimizer_criteria(config):
    """Record the convergence standard the selected optimiser applies.

    ``optimizer_converged`` means a different physical test for each
    ``opt_solver``; only the solver's name reached the artifact before.  The
    values are read from the installed optimiser rather than typed here, so
    the record follows the environment that ran.
    """
    solver = config.get("opt_solver")
    record = {"solver": solver, "source": None}
    try:
        if solver == "geometric":
            from geometric.params import OptParams

            params = OptParams()
            record["source"] = "geometric.params.OptParams (pyscf passes the same defaults)"
            for key in ("energy", "grms", "gmax", "drms", "dmax"):
                record["convergence_" + key] = float(
                    getattr(params, "Convergence_" + key)
                )
        elif solver == "berny":
            record["source"] = "pyberny defaults (not recorded)"
        elif solver == "ase":
            record["source"] = "ase BFGS fmax default 0.05 eV/Angstrom"
    except Exception as exc:
        record["error_type"] = type(exc).__name__
    return record


def _spin_populations(mf, config):
    """Per-atom Mulliken spin populations, open-shell references only."""
    if int(config.get("spin") or 0) == 0:
        return None
    evaluator = getattr(mf, "mulliken_spin_pop", None)
    if not callable(evaluator):
        raise AttributeError("mean-field object has no mulliken_spin_pop")
    _per_ao, per_atom = evaluator(verbose=0)
    return _to_host_array(per_atom).astype(float)


def _spin_diagnostic(mf):
    """Return PySCF's final-state ``(<S^2>, effective multiplicity)``."""
    evaluator = getattr(mf, "spin_square", None)
    if not callable(evaluator):
        raise AttributeError("mean-field object has no spin_square evaluator")
    observed = evaluator()
    if not isinstance(observed, (tuple, list)) or len(observed) != 2:
        raise ValueError("spin_square must return two scalar values")
    values = []
    for value in observed:
        materialized = _to_host_array(value)
        if materialized.size != 1:
            raise ValueError("spin_square returned a non-scalar value")
        values.append(float(materialized.reshape(-1)[0]))
    return tuple(values)


def main():
    label = CONFIG["label"]
    log_path = label + ".out"
    results_path = label + ".h5"

    started_at = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    t0 = time.time()

    spec_fields = __CHEMSMART_APPLIED_SPEC_FIELDS__
    spec = {key: CONFIG.get(key) for key in spec_fields}
    status = {
        "stages": {},
        "engine_complete": False,
        "normal_termination": False,
        "failure": None,
        "properties": {},
    }
    results = {}
    runtime = {}
    current_stage = "initialization"

    try:
        _apply_threads(CONFIG)
        mol = _build_mole(CONFIG, log_path)
        mf = _build_method(CONFIG, mol)
        runtime = _runtime_provenance(CONFIG, mf)
        # ``_build_method`` resolves environment-owned values such as the PCM
        # dielectric.  Echo the value actually applied, not the controller's
        # unresolved placeholder.
        spec["solvent_eps"] = CONFIG.get("solvent_eps")
        spec["solvent_lebedev_order"] = CONFIG.get(
            "solvent_lebedev_order"
        )
        spec["materializations"] = CONFIG.get("materializations", {})
        CONFIG["applied_settings_sha256"] = _settings_digest(spec)
        spec["applied_settings_sha256"] = CONFIG[
            "applied_settings_sha256"
        ]

        energies = []
        total_energy = None
        excited_root = CONFIG.get("excited_state_root")
        correlated = CONFIG.get("ab_initio") in ("mp2", "ccsd", "ccsd(t)")
        for stage in CONFIG["stages"]:
            current_stage = stage
            if stage == "scf":
                energy = mf.kernel()
                energies.append(float(energy))
                stage_status = {"converged": bool(mf.converged)}
                cycles = getattr(mf, "cycles", None)
                if cycles is not None:
                    stage_status["iterations"] = int(cycles)
                    # Historical alias retained for existing consumers.
                    stage_status["cycles"] = int(cycles)
                status["stages"]["scf"] = stage_status
            elif stage == "opt":
                surface = mf
                surface_record = None
                if excited_root is not None:
                    # The surface is root k of the response at every step.
                    # The scanner's own ``converged`` property indexes the
                    # per-root flags with the 1-based root (PySCF 2.14,
                    # grad/tdrhf.py) and raises for root == nstates, so
                    # the driver reads the response object directly.
                    td_start = _build_response(CONFIG, mf)
                    td_start.kernel()
                    runtime["response_class"] = _class_name(td_start)
                    start_roots = int(np.asarray(td_start.e).reshape(-1).size)
                    root = int(excited_root)
                    if start_roots < root:
                        raise FollowedRootFiltered(
                            "excited_state_root %d is not in the %d roots "
                            "PySCF kept at the supplied geometry (%d "
                            "requested; roots below positive_eig_threshold "
                            "%g Eh are dropped)"
                            % (
                                root,
                                start_roots,
                                int(CONFIG["nstates"]),
                                float(td_start.positive_eig_threshold),
                            )
                        )
                    surface = td_start.nuc_grad_method().as_scanner(state=root)
                    surface_record = {
                        "root": root,
                        "response_method": str(CONFIG["response_method"]),
                        "state_manifold": str(CONFIG["state_manifold"]),
                        "nstates": int(CONFIG["nstates"]),
                        "start_root_converged": bool(
                            _root_convergence(td_start)[root - 1]
                        ),
                        "start_roots_obtained": start_roots,
                        "followed_root_total_energy_start": float(
                            td_start.e_tot[root - 1]
                        ),
                        "gradient_scanner_class": _class_name(surface),
                    }
                elif correlated:
                    # The surface is the correlated method's energy; the
                    # scanner rebuilds SCF, amplitudes and (for CC) the
                    # lambda equations at every geometry.
                    corr_start = _build_correlated(CONFIG, mf)
                    runtime["correlated_class"] = _class_name(corr_start)
                    surface = corr_start.nuc_grad_method().as_scanner()
                    surface_record = {
                        "method": str(CONFIG["ab_initio"]),
                        "frozen_core_applied": _frozen_core_count(corr_start),
                        "gradient_scanner_class": _class_name(surface),
                    }
                try:
                    optimizer_converged, mol_eq = _run_opt(CONFIG, surface)
                except IndexError as exc:
                    if excited_root is None:
                        raise
                    raise FollowedRootFiltered(
                        "excited_state_root %d vanished from the spectrum "
                        "during the optimisation (PySCF drops roots below "
                        "positive_eig_threshold): %s" % (int(excited_root), exc)
                    ) from exc
                # Re-converge on the optimised geometry so that every
                # reported property belongs to the same structure.
                mf.reset(mol_eq)
                energy = mf.kernel()
                energies.append(float(energy))
                mol = mol_eq
                final_scf_converged = bool(mf.converged)
                stage_status = {
                    "converged": bool(
                        optimizer_converged and final_scf_converged
                    ),
                    "optimizer_converged": bool(optimizer_converged),
                    "final_scf_converged": final_scf_converged,
                }
                cycles = getattr(mf, "cycles", None)
                if cycles is not None:
                    stage_status["final_scf_iterations"] = int(cycles)
                # The final SCF starts from the optimiser's last density
                # (mf.reset keeps the orbitals), so the electronic state is
                # the one that path reached; say so on the record.
                stage_status["final_scf_from_optimizer_density"] = True
                stage_status["convergence_criteria"] = _optimizer_criteria(
                    CONFIG
                )
                if surface_record is not None:
                    # The last gradient the scanner evaluated is the
                    # number behind ``optimizer_converged``.
                    last_gradient = getattr(surface, "de", None)
                    if last_gradient is not None:
                        surface_record["final_gradient_max_eh_per_bohr"] = (
                            float(np.max(np.abs(_to_host_array(last_gradient))))
                        )
                        surface_record["gradient_source"] = (
                            "last scanner evaluation (scanner.de)"
                        )
                    base = getattr(surface, "base", None)
                    if correlated and base is not None:
                        surface_record["amplitudes_converged"] = (
                            None
                            if not hasattr(base, "converged")
                            else bool(base.converged)
                        )
                        surface_record["lambda_converged"] = (
                            None
                            if not hasattr(base, "converged_lambda")
                            else bool(base.converged_lambda)
                        )
                    stage_status[
                        "excited_state" if excited_root is not None
                        else "correlated"
                    ] = surface_record
                status["stages"]["opt"] = stage_status
            elif stage == "td":
                td, td_stage = _run_td(CONFIG, mf, results, status, runtime)
                if excited_root is not None:
                    root = int(excited_root)
                    obtained = int(td_stage["nstates_obtained"])
                    if obtained < root:
                        raise FollowedRootFiltered(
                            "excited_state_root %d is not in the %d roots "
                            "PySCF kept at the reached geometry"
                            % (root, obtained)
                        )
                    hartree_to_ev = 27.211386245988
                    excitations = np.asarray(td.e, dtype=float).reshape(-1)
                    total_energy = float(td.e_tot[root - 1])
                    neighbours = []
                    if root >= 2:
                        neighbours.append(
                            float(excitations[root - 1] - excitations[root - 2])
                        )
                    if root < obtained:
                        neighbours.append(
                            float(excitations[root] - excitations[root - 1])
                        )
                    record = status["stages"]["opt"].setdefault(
                        "excited_state", {}
                    )
                    record["followed_root_total_energy_end"] = total_energy
                    record["followed_root_converged_end"] = bool(
                        _root_convergence(td)[root - 1]
                    )
                    record["root_gap_to_ground_end_ev"] = float(
                        excitations[root - 1] * hartree_to_ev
                    )
                    record["root_gap_to_neighbour_end_ev"] = (
                        float(min(neighbours) * hartree_to_ev)
                        if neighbours
                        else None
                    )
                    record["roots_filtered_end"] = int(
                        td_stage["roots_filtered"]
                    )
            elif stage == "corr":
                total_energy = _run_corr(CONFIG, mf, results, status, runtime)
            elif stage == "hess":
                # GPU4PySCF returns a CuPy-like Hessian. PySCF's CPU thermo
                # helper and HDF5 writer must never receive that device array.
                raw_hessian = _to_host_array(mf.Hessian().kernel()).astype(float)
                hessian, raw_hessian_antisymmetry = (
                    _symmetrize_cartesian_hessian(raw_hessian)
                )
                from pyscf.hessian import thermo

                # Preserve imaginary modes as negative real wavenumbers.
                # PySCF's default returns complex values, which a float HDF5
                # dataset would otherwise truncate to zero imaginary parts.
                analysis = thermo.harmonic_analysis(
                    mol, hessian, imaginary_freq=False
                )
                results["hessian"] = hessian
                results["vibrational_frequencies"] = np.asarray(
                    analysis["freq_wavenumber"], dtype=float
                )
                results["normal_modes"] = np.asarray(
                    analysis["norm_mode"], dtype=float
                )
                results["reduced_masses"] = np.asarray(
                    analysis["reduced_mass"], dtype=float
                )
                results["force_constants"] = np.asarray(
                    analysis["force_const_dyne"], dtype=float
                )
                hess_status = {
                    "converged": True,
                    "cartesian_symmetrization_applied": True,
                    "raw_max_abs_antisymmetry_eh_per_bohr2": (
                        raw_hessian_antisymmetry
                    ),
                    # harmonic_analysis takes mol.atom_mass_list(
                    # isotope_avg=True) when no mass is passed; the host's
                    # thermochemistry applies its own table to rotation and
                    # translation, so the table behind the frequencies is
                    # stated where the frequencies are.
                    "mass_convention": "isotope_averaged",
                    "mass_source": (
                        "pyscf.gto.Mole.atom_mass_list(isotope_avg=True)"
                    ),
                }
                # The projected spectrum can be all-real at a geometry that
                # is not stationary, so the gradient at the Hessian geometry
                # is the only fact that says how stationary it was.  It is
                # a stage fact, never a refusal: a Hessian off a stationary
                # point is a legitimate request.
                try:
                    gradient = _to_host_array(
                        mf.nuc_grad_method().kernel()
                    ).astype(float)
                    results["forces"] = -gradient
                    hess_status["gradient_computed"] = True
                    hess_status["max_abs_gradient_eh_per_bohr"] = float(
                        np.max(np.abs(gradient))
                    )
                except Exception as exc:
                    hess_status["gradient_computed"] = False
                    hess_status["gradient_failure"] = {
                        "type": type(exc).__name__,
                        "message": str(exc),
                    }
                status["stages"]["hess"] = hess_status
            else:
                raise ValueError("Unknown stage: %s" % stage)

        results["energies"] = np.asarray(energies, dtype=float)
        # The SCF energy at the final geometry, stated once, and the
        # total energy of the surface the job computed on: the SCF for
        # HF/DFT, the correlated total, the followed root's total for an
        # excited-state optimisation, the reference for a td spectrum.
        results["scf_energy"] = float(energies[-1])
        results["total_energy"] = float(
            energies[-1] if total_energy is None else total_energy
        )
        results["positions"] = np.asarray(
            mol.atom_coords(unit="Angstrom"), dtype=float
        )
        results["atomic_numbers"] = np.asarray(
            mol.atom_charges(), dtype=int
        )
        results["mo_energy"] = _to_host_array(mf.mo_energy).astype(float)
        results["mo_occ"] = _to_host_array(mf.mo_occ).astype(float)

        # Optional properties are explicit status records. Missing values must
        # never be indistinguishable from a property that was not attempted.
        # ``forces`` is deliberately absent: the public PySCF settings reject
        # force requests, so an ordinary SP must not pay for an undeclared
        # nuclear-gradient calculation after its requested SCF has completed.
        try:
            spin_square, effective_multiplicity = _spin_diagnostic(mf)
            results["spin_square"] = spin_square
            results[
                "spin_square_effective_multiplicity"
            ] = effective_multiplicity
            status["properties"]["spin_square"] = {"status": "ok"}
        except Exception as exc:
            status["properties"]["spin_square"] = {
                "status": "unavailable",
                "failure": {
                    "type": type(exc).__name__,
                    "message": str(exc),
                },
            }
        try:
            spin_populations = _spin_populations(mf, CONFIG)
            if spin_populations is None:
                status["properties"]["mulliken_spin_populations"] = {
                    "status": "not_applicable",
                    "reason": "closed-shell reference carries no spin",
                }
            else:
                results["mulliken_spin_populations"] = spin_populations
                status["properties"]["mulliken_spin_populations"] = {
                    "status": "ok"
                }
        except Exception as exc:
            status["properties"]["mulliken_spin_populations"] = {
                "status": "unavailable",
                "failure": {
                    "type": type(exc).__name__,
                    "message": str(exc),
                },
            }
        try:
            results["mulliken_charges"] = _to_host_array(
                mf.mulliken_pop(verbose=0)[1]
            ).astype(float)
            status["properties"]["mulliken_charges"] = {"status": "ok"}
        except Exception as exc:
            status["properties"]["mulliken_charges"] = {
                "status": "unavailable",
                "failure": {
                    "type": type(exc).__name__,
                    "message": str(exc),
                },
            }
        try:
            results["dipole_moment"] = _to_host_array(
                mf.dip_moment(unit="Debye", verbose=0)
            ).astype(float)
            status["properties"]["dipole_moment"] = {"status": "ok"}
        except Exception as exc:
            status["properties"]["dipole_moment"] = {
                "status": "unavailable",
                "failure": {
                    "type": type(exc).__name__,
                    "message": str(exc),
                },
            }
        try:
            from pyscf import symm
            from pyscf.symm import geom as symm_geom

            results["point_group"] = symm.detect_symm(mol._atom)[0]
            # The group is a threshold decision (symm_geom_tol, scaled by
            # 1/sqrt(1+natm) inside detect_symm), far tighter than any
            # optimiser's displacement criterion; the tolerance rides beside
            # the label so a symmetry number derived from it is a stated
            # convention and not a printed fact.
            status["properties"]["point_group"] = {
                "status": "ok",
                "source": "pyscf.symm.detect_symm",
                "tolerance_bohr": float(symm_geom.TOLERANCE),
            }
        except Exception as exc:
            results["point_group"] = None
            status["properties"]["point_group"] = {
                "status": "unavailable",
                "failure": {
                    "type": type(exc).__name__,
                    "message": str(exc),
                },
            }

        spec["num_basis_functions"] = int(mol.nao)
        spec["num_shells"] = int(mol.nbas)
        spec["num_electrons"] = int(mol.nelectron)
        spec["nelec"] = [int(n) for n in mol.nelec]

        status["engine_complete"] = (
            len(status["stages"]) == len(CONFIG["stages"])
        )
        status["normal_termination"] = (
            status["engine_complete"]
            and all(
                entry.get("converged", False)
                for entry in status["stages"].values()
            )
        )
    except Exception as exc:  # noqa: BLE001 - recorded, then re-raised
        status["failure"] = {
            "type": type(exc).__name__,
            "message": str(exc),
            "stage": current_stage,
            "traceback": traceback.format_exc(),
        }

    wall_seconds = time.time() - t0
    ended_at = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    provenance = _provenance(
        CONFIG, started_at, ended_at, wall_seconds, runtime
    )

    _write_h5(results_path, spec, provenance, status, results)

    if status["failure"] is not None:
        sys.stderr.write(status["failure"]["traceback"])
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
'''

# Keep the standalone skeleton invariant while deriving its applied-spec field
# tuple from the same declaration used by the fake runner.
_SKELETON = _SKELETON.replace(
    "__CHEMSMART_APPLIED_SPEC_FIELDS__", repr(APPLIED_SPEC_FIELDS)
)
_SKELETON = _SKELETON.replace("__CHEMSMART_RESULT_UNITS__", repr(RESULT_UNITS))
