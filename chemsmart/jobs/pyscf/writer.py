"""Generate the standalone PySCF driver script for a job.

The emitted ``label.py`` is the source of ``chemsmart.io.pyscf.driver`` with
one JSON configuration substituted in. Because v1 covers only sp/opt/hess the
driver logic is invariant, so nothing is gained by templating source text --
and free-form generation would add a quoting, escaping and injection surface
for no benefit. Only the configuration varies, and it holds scalars, strings
and a geometry array.

The driver imports ``pyscf``, ``numpy``, ``h5py`` and the standard library
**only**. It must never import ``chemsmart``: that would drag CHEMSMART's
``numpy~=1.26.4`` pin and the rdkit/pymatgen/ase tree into every compute
process, and would make the per-program ``CONDA_ENV`` -- the mechanism that
lets a GPU job run in its own environment -- unusable.
"""

import hashlib
import json
import logging
import os
from pathlib import Path

from chemsmart import __version__ as chemsmart_version
from chemsmart.io.pyscf import driver
from chemsmart.io.pyscf.driver import (  # noqa: F401
    APPLIED_SPEC_FIELDS,
    H5_NULL_ATTRIBUTE,
    RESULTS_SCHEMA_VERSION,
    write_pyscf_h5,
)
from chemsmart.jobs.pyscf.settings import (
    PYSCF_DEFGRIDS,
    PYSCF_SOLVENT_MODELS,
)

logger = logging.getLogger(__name__)


def applied_pyscf_spec(config):
    """Return the resolved fields written under ``spec/``."""
    return {key: config.get(key) for key in APPLIED_SPEC_FIELDS}


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
        # A transition-state search is an optimisation towards a saddle
        # point, and an IRC a sequence of constrained optimisations.
        optimises = bool({"opt", "ts", "irc"} & set(job.stages))

        # A job can be previewed before a runner is attached.
        jobrunner = job.jobrunner
        num_threads = (jobrunner.num_cores if jobrunner else None) or 1
        mem_gb = jobrunner.mem_gb if jobrunner else None
        # PySCF's max_memory is per-process in MB and bounds integral
        # buffering, not total allocation. Leave PySCF's own default in
        # place when the runner has no explicit budget.
        max_memory_mb = float(mem_gb) * 1024 if mem_gb else None

        config = {
            "schema_version": RESULTS_SCHEMA_VERSION,
            # Captured by the CHEMSMART process because the standalone driver
            # deliberately cannot import chemsmart from the compute env.
            "chemsmart_version": chemsmart_version,
            "label": job.label,
            "program": "pyscf",
            # job.JOBTYPE ("ts"), not job.TYPE ("pyscf_ts"): thermochemistry
            # recognises a transition state by the former.
            "jobtype": settings.jobtype or job.JOBTYPE,
            "stages": list(job.stages),
            "symbols": symbols,
            "positions": positions,
            "unit": "Angstrom",
            "basis": settings.basis,
            "charge": int(charge),
            "multiplicity": int(multiplicity),
            "spin": settings.multiplicity_to_spin(multiplicity),
            "xc": settings.xc,
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
            "opt_solver": settings.opt_solver if optimises else None,
            "opt_maxsteps": settings.opt_maxsteps if optimises else None,
            # A single point needs no gradient unless asked for; an optimised
            # geometry's residual forces are part of its result.
            "forces": bool(settings.forces or optimises),
            "point_group_tol": float(settings.point_group_tol),
            "correlation": settings.resolved_correlation,
            "active_space": settings.resolved_active_space,
            "scf_stability": settings.scf_stability,
            "irc_direction": (
                settings.irc_direction if "irc" in job.stages else None
            ),
            "engine": settings.resolved_engine,
            "num_threads": int(num_threads),
            "max_memory_mb": max_memory_mb,
            # The project-settings builder attaches this only when a concrete
            # source YAML exists. A null records that it was unavailable; it
            # is never inferred from a project name or regenerated content.
            "project_yaml_digest": settings.project_yaml_digest,
        }
        config["settings_digest"] = self.settings_digest(config)
        return config

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
            "chemsmart_version",
            "label",
            "num_threads",
            "max_memory_mb",
            "project_yaml_digest",
            "settings_digest",
        }
        payload = {k: v for k, v in config.items() if k not in ignored}
        body = json.dumps(payload, sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(body.encode("utf-8")).hexdigest()

    # ------------------------------------------------------------------
    # emission
    # ------------------------------------------------------------------

    def write(self, target_directory=None):
        """Write ``label.py`` and return its path."""
        directory = target_directory or self.job.folder
        os.makedirs(directory, exist_ok=True)
        path = os.path.join(directory, f"{self.job.label}.py")
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
        source = Path(driver.__file__).read_text(encoding="utf-8")
        return source.replace("__CHEMSMART_CONFIG__", payload, 1)
