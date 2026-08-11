"""Generate the standalone PySCF driver script for a job.

The emitted ``label.py`` is a **fixed skeleton plus one configuration dict**.
Executable v1 workflows cover only sp/opt/hess; td is a declarative preview
that the emitted driver refuses to run. The driver logic is invariant, so
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
RESULT_CONTRACT_VERSION = "chemsmart.pyscf-result-contract.v3"
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

APPLIED_SPEC_FIELDS = LEGACY_APPLIED_SPEC_FIELDS + (
    "result_contract_version",
    "reference_family",
)

# Units are part of the machine contract rather than prose in a reader.  Every
# numeric dataset the fixed driver can emit has one explicit unit, including
# quantities that are dimensionless.
RESULT_UNITS = {
    "atomic_numbers": "dimensionless",
    "dipole_moment": "Debye",
    "energies": "Eh",
    "excitation_energies": "Eh",
    "force_constants": "Dyne/Angstrom",
    "forces": "Eh/Bohr",
    "hessian": "Eh/Bohr^2",
    "mo_energy": "Eh",
    "mo_occ": "electron",
    "mulliken_charges": "elementary_charge",
    "normal_modes": "atomic_mass_unit^-1/2",
    "oscillator_strengths": "dimensionless",
    "positions": "Angstrom",
    "reduced_masses": "atomic_mass_unit",
    "spin_square": "dimensionless",
    "spin_square_effective_multiplicity": "dimensionless",
    "vibrational_frequencies": "cm^-1",
}


def applied_pyscf_spec(config):
    """Return the resolved fields written under ``spec/``."""
    return {key: config.get(key) for key in APPLIED_SPEC_FIELDS}


def applied_pyscf_spec_fields(spec):
    """Return the digest vocabulary for a current or historical artifact."""

    if spec.get("result_contract_version") == RESULT_CONTRACT_VERSION:
        return APPLIED_SPEC_FIELDS
    return LEGACY_APPLIED_SPEC_FIELDS


def pyscf_td_response_materialization(settings):
    """Map typed TD/TDA intent to a non-executable PySCF response plan.

    The preview artifact remains inert. This manifest makes the intended
    target-library construction explicit and testable without embedding a
    dormant engine call that could later become an execution escape hatch.
    """

    if str(getattr(settings, "jobtype", "")).strip().lower() != "td":
        return None
    response_method = str(settings.response_method).strip().lower()
    factory_api = {
        "tda": "pyscf.tdscf.rks.TDA",
        "tddft": "pyscf.tdscf.rks.TDDFT",
    }[response_method]
    body = {
        "schema_version": TD_RESPONSE_MATERIALIZATION_SCHEMA_VERSION,
        "ground_state_reference_family": "rks",
        "ground_state_stage": "scf",
        "response_method": response_method,
        "response_factory_api": factory_api,
        "state_manifold": "singlet",
        "nstates": int(settings.nstates),
        "operation_order": (
            "ground_state_scf",
            "response_construct",
            "set_singlet_channel",
            "set_nstates",
            "vertical_excitation_kernel",
        ),
        "execution_policy": "preview_only",
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
        _write_mapping(handle.create_group("spec"), spec)
        _write_mapping(handle.create_group("provenance"), provenance)
        _write_mapping(handle.create_group("status"), status)
        results_group = handle.create_group("results")
        _write_mapping(results_group, results)
        _attach_result_units(results_group)


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
                29
                if settings.engine == "gpu" and call is not None
                else None
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
            # TD/TDA is deliberately a typed preview, not a latent engine
            # escape hatch. The generated artifact carries the requested
            # semantics but refuses direct execution even outside ChemSmart.
            "preview_only": settings.jobtype == "td",
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
        td_materialization = pyscf_td_response_materialization(settings)
        if td_materialization is not None:
            config["materializations"]["td_response_plan"] = (
                td_materialization
            )
        geometry_payload = {
            "symbols": config["symbols"],
            "positions": config["positions"],
            "unit": config["unit"],
            "charge": config["charge"],
            "multiplicity": config["multiplicity"],
        }
        config["input_geometry_sha256"] = self._json_digest(
            geometry_payload
        )
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
Executable SP/OPT/HESS artifacts may be rerun with the bound environment.
Preview-only artifacts refuse direct execution.

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


def _run_opt(config, mf):
    """Optimise the geometry and return (converged, optimised Mole).

    assert_convergence is disabled so that a failed optimisation still
    produces a results file recording the failure, rather than raising with
    no artifact. Downstream, normal_termination stays False.
    """
    solver = config["opt_solver"]
    maxsteps = config["opt_maxsteps"]
    if solver == "geometric":
        from pyscf.geomopt import geometric_solver

        return geometric_solver.kernel(
            mf, assert_convergence=False, maxsteps=maxsteps
        )
    if solver == "berny":
        from pyscf.geomopt import berny_solver

        return berny_solver.kernel(
            mf, assert_convergence=False, maxsteps=maxsteps
        )
    if solver == "ase":
        from pyscf.geomopt import ase_solver

        return ase_solver.kernel(mf, target="atoms", max_steps=maxsteps)
    raise ValueError("Unknown opt_solver: %s" % solver)


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
        _write_mapping(handle.create_group("spec"), spec)
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

    if CONFIG.get("preview_only"):
        raise RuntimeError(
            "This ChemSmart PySCF artifact is preview-only and cannot launch "
            "a chemistry engine."
        )

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
                optimizer_converged, mol_eq = _run_opt(CONFIG, mf)
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
                status["stages"]["opt"] = stage_status
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
                status["stages"]["hess"] = {
                    "converged": True,
                    "cartesian_symmetrization_applied": True,
                    "raw_max_abs_antisymmetry_eh_per_bohr2": (
                        raw_hessian_antisymmetry
                    ),
                }
            elif stage == "td":
                raise RuntimeError(
                    "PySCF TD/TDA is a ChemSmart preview-only capability."
                )
            else:
                raise ValueError("Unknown stage: %s" % stage)

        results["energies"] = np.asarray(energies, dtype=float)
        results["positions"] = np.asarray(
            mol.atom_coords(unit="Angstrom"), dtype=float
        )
        results["atomic_numbers"] = np.asarray(
            mol.atom_charges(), dtype=int
        )
        results["mo_energy"] = _to_host_array(mf.mo_energy).astype(float)
        results["mo_occ"] = _to_host_array(mf.mo_occ).astype(float)

        # PySCF gradients are dE/dR in Hartree/Bohr. Molecule and database
        # contracts expose forces, so negate the exact final-geometry gradient.
        # Optional properties are explicit status records. Missing values must
        # never be indistinguishable from a property that was not attempted.
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
            results["forces"] = -_to_host_array(
                mf.nuc_grad_method().kernel()
            ).astype(float)
            status["properties"]["forces"] = {"status": "ok"}
        except Exception as exc:
            status["properties"]["forces"] = {
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

            results["point_group"] = symm.detect_symm(mol._atom)[0]
            status["properties"]["point_group"] = {"status": "ok"}
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
_SKELETON = _SKELETON.replace(
    "__CHEMSMART_RESULT_UNITS__", repr(RESULT_UNITS)
)
