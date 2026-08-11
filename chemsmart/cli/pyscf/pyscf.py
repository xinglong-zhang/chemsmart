"""PySCF CLI group.

Group-level options are execution controls and per-invocation overrides only.
Method, basis, response settings and phase are governed by project YAML, as
for Gaussian and ORCA.

Option names deliberately reuse the existing ChemSmart vocabulary
(``--functional``, ``--basis``, ``--aux-basis``, ``--defgrid``, ``--scf-tol``)
rather than inventing PySCF-specific spellings. One normalised vocabulary
across programs is what lets a downstream tool compare a PySCF invocation
with an ORCA one without a per-program translation table.
"""

import functools
import logging
import os

import click

from chemsmart.cli.job import (
    click_file_label_and_index_options,
    click_filename_options,
    click_pubchem_options,
)
from chemsmart.io.molecules.structure import Molecule
from chemsmart.utils.cli import MyGroup
from chemsmart.utils.io import clean_label
from chemsmart.utils.utils import return_objects_and_indices_from_string_index

logger = logging.getLogger(__name__)


def click_pyscf_options(f):
    """Project selection for PySCF jobs."""

    @click.option(
        "--project", "-p", type=str, default=None, help="Project settings."
    )
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper


def click_pyscf_settings_options(f):
    """Per-invocation overrides of project-YAML settings.

    Every option here must also be appended to the ``keywords`` tuple in the
    group body, or it is silently dropped by ``settings.merge()``.
    """

    @click.option(
        "-t", "--title", type=str, default=None, help="PySCF job title."
    )
    @click.option(
        "-c", "--charge", type=int, default=None, help="Molecular charge."
    )
    @click.option(
        "-m",
        "--multiplicity",
        type=int,
        default=None,
        help="Spin multiplicity (2S+1). PySCF's own `spin` is 2S; ChemSmart "
        "converts, so always give the multiplicity here.",
    )
    @click.option(
        "-x",
        "--functional",
        type=str,
        default=None,
        help="DFT functional, resolved through libxc (e.g. b3lyp, m062x).",
    )
    @click.option(
        "-A",
        "--ab-initio",
        type=str,
        default=None,
        help="Ab initio method. Use 'hf' for a Hartree-Fock reference.",
    )
    @click.option(
        "-b",
        "--basis",
        type=str,
        default=None,
        help="Basis set, hyphenated PySCF spelling (e.g. def2-svp).",
    )
    @click.option(
        "--aux-basis",
        "-ab",
        type=str,
        default=None,
        help="Auxiliary basis for density fitting. Omit to let PySCF pick.",
    )
    @click.option(
        "-d",
        "--dispersion",
        type=str,
        default=None,
        help="Dispersion correction (e.g. d3bj). Requires pyscf-dispersion.",
    )
    @click.option(
        "--defgrid",
        type=click.Choice(
            ["defgrid1", "defgrid2", "defgrid3"], case_sensitive=False
        ),
        default=None,
        help="DFT integration grid density.",
    )
    @click.option(
        "--scf-tol",
        type=float,
        default=None,
        help="SCF convergence tolerance (mf.conv_tol).",
    )
    @click.option(
        "--scf-maxiter",
        type=int,
        default=None,
        help="Maximum SCF cycles (mf.max_cycle).",
    )
    @click.option(
        "--density-fit/--no-density-fit",
        default=None,
        help="Use resolution-of-the-identity density fitting.",
    )
    @click.option(
        "--opt-solver",
        type=click.Choice(["geometric", "berny", "ase"], case_sensitive=False),
        default=None,
        help="Geometry-optimisation backend. Recorded in the results file, "
        "since the backends do not share convergence criteria.",
    )
    @click.option(
        "--opt-maxsteps",
        type=int,
        default=None,
        help="Maximum geometry-optimisation steps.",
    )
    @click.option(
        "--gpu/--no-gpu",
        default=None,
        help="Override the project engine. Real --gpu execution requires a "
        "positive resolved NUM_GPUS; --fake may preview the GPU script without "
        "a device. --no-gpu selects CPU. Omission preserves project YAML, and "
        "an unavailable GPU request never falls back to CPU.",
    )
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper


def click_pyscf_solvent_options(f):
    """Implicit-solvation options for PySCF jobs."""

    @click.option(
        "-sm",
        "--solvent-model",
        type=click.Choice(
            ["pcm", "iefpcm", "cpcm", "cosmo", "ssvpe", "smd"],
            case_sensitive=False,
        ),
        default=None,
        help="Implicit solvation model.",
    )
    @click.option(
        "-si",
        "--solvent-id",
        type=str,
        default=None,
        help="Solvent name; must exist in pyscf.solvent.smd.solvent_db. "
        "Its dielectric is resolved explicitly, because PySCF's PCM "
        "otherwise defaults to water regardless of the requested solvent.",
    )
    @click.option(
        "--remove-solvent/--no-remove-solvent",
        default=False,
        help="Drop the solvent model from the project settings.",
    )
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper


@click.group(cls=MyGroup)
@click_pyscf_options
@click_filename_options
@click_file_label_and_index_options
@click_pyscf_settings_options
@click_pyscf_solvent_options
@click_pubchem_options
@click.pass_context
def pyscf(
    ctx,
    project,
    filename,
    label,
    append_label,
    index,
    title,
    charge,
    multiplicity,
    functional,
    ab_initio,
    basis,
    aux_basis,
    dispersion,
    defgrid,
    scf_tol,
    scf_maxiter,
    density_fit,
    opt_solver,
    opt_maxsteps,
    gpu,
    solvent_model,
    solvent_id,
    remove_solvent,
    pubchem,
):
    """CLI subcommand for running PySCF jobs using the chemsmart framework.

    Resolves the project settings, the molecules and the per-invocation
    overrides, then publishes them on ``ctx.obj`` for the sp/opt/hess/td leaves.
    """
    from chemsmart.jobs.pyscf.settings import PySCFJobSettings
    from chemsmart.settings.pyscf import PySCFProjectSettings

    project_settings = PySCFProjectSettings.from_project(project)
    logger.debug(f"Loaded PySCF project settings: {project_settings}")
    project_explicit_fields = project_settings.explicit_fields(
        ctx.invoked_subcommand
    )

    job_settings = PySCFJobSettings.default()

    # `keywords` is a merge whitelist: an attribute set below only survives
    # `project_settings.<jobtype>_settings().merge(job_settings, keywords)`
    # if its name is in this tuple. Forgetting an entry is a silent no-op.
    keywords = ()
    explicit_state_fields = set(
        project_explicit_fields & {"charge", "multiplicity"}
    )

    if (charge is None) != (multiplicity is None):
        raise ValueError(
            "PySCF CLI charge and multiplicity overrides must be supplied "
            "together so the electronic state cannot be partially replaced."
        )

    if title is not None:
        job_settings.title = title
        keywords += ("title",)
    if charge is not None:
        job_settings.charge = charge
        keywords += ("charge",)
        explicit_state_fields.add("charge")
    if multiplicity is not None:
        job_settings.multiplicity = multiplicity
        keywords += ("multiplicity",)
        explicit_state_fields.add("multiplicity")
    if functional is not None and ab_initio is not None:
        raise ValueError(
            "Specify either --functional or --ab-initio, not both."
        )
    if ab_initio is not None and defgrid is not None:
        raise ValueError("--defgrid is DFT-only and cannot be used with HF.")
    if density_fit is False and aux_basis is not None:
        raise ValueError(
            "--aux-basis cannot be used when --no-density-fit is set."
        )
    if functional is not None:
        job_settings.functional = functional
        # A functional override must also replace a project-owned HF method.
        job_settings.ab_initio = None
        keywords += ("functional", "ab_initio")
    if ab_initio is not None:
        job_settings.ab_initio = ab_initio
        # An HF override must also replace a project-owned DFT functional.
        job_settings.functional = None
        # A DFT grid inherited from that project would otherwise survive the
        # whitelist merge and then be recorded despite being inapplicable.
        job_settings.defgrid = None
        keywords += ("ab_initio", "functional", "defgrid")
    if basis is not None:
        job_settings.basis = basis
        keywords += ("basis",)
    if aux_basis is not None:
        job_settings.aux_basis = aux_basis
        keywords += ("aux_basis",)
    if dispersion is not None:
        job_settings.dispersion = dispersion
        keywords += ("dispersion",)
    if defgrid is not None:
        job_settings.defgrid = defgrid
        keywords += ("defgrid",)
    if scf_tol is not None:
        job_settings.scf_tol = scf_tol
        keywords += ("scf_tol",)
    if scf_maxiter is not None:
        job_settings.scf_maxiter = scf_maxiter
        keywords += ("scf_maxiter",)
    if density_fit is not None:
        job_settings.density_fit = density_fit
        keywords += ("density_fit",)
        if density_fit is False:
            # Clear an inherited auxiliary basis together with the operation
            # that would consume it.
            job_settings.aux_basis = None
            keywords += ("aux_basis",)
    if opt_solver is not None:
        job_settings.opt_solver = opt_solver
        keywords += ("opt_solver",)
    if opt_maxsteps is not None:
        job_settings.opt_maxsteps = opt_maxsteps
        keywords += ("opt_maxsteps",)

    # An explicit flag overrides the project engine.  Omission deliberately
    # leaves ``engine`` out of the merge whitelist so project YAML remains
    # authoritative.  Allocation conformance is resolved later in the leaf.
    if gpu is not None:
        job_settings.engine = "gpu" if gpu else "cpu"
        keywords += ("engine",)

    if remove_solvent:
        job_settings.solvent_model = None
        job_settings.solvent_id = None
        keywords += ("solvent_model", "solvent_id")
    else:
        if solvent_model is not None:
            job_settings.solvent_model = solvent_model
            keywords += ("solvent_model",)
        if solvent_id is not None:
            job_settings.solvent_id = solvent_id
            keywords += ("solvent_id",)

    if filename is None and pubchem is None:
        raise ValueError(
            "[filename] or [pubchem] has not been specified!\n"
            "Please specify one of them!"
        )
    if filename and pubchem:
        raise ValueError(
            "Both [filename] and [pubchem] have been specified!\n"
            "Please specify only one of them."
        )

    molecules = None
    if filename:
        molecules = Molecule.from_filepath(
            filepath=filename, index=":", return_list=True
        )
        assert (
            molecules is not None
        ), f"Could not obtain molecule from {filename}!"
    if pubchem:
        molecules = Molecule.from_pubchem(identifier=pubchem, return_list=True)
        assert (
            molecules is not None
        ), f"Could not obtain molecule from PubChem {pubchem}!"

    if label is not None and append_label is not None:
        raise ValueError(
            "Only give PySCF input filename or name to be appended, "
            "but not both!"
        )
    if append_label is not None:
        label = os.path.splitext(os.path.basename(filename))[0]
        label = f"{label}_{append_label}"
    if label is None and append_label is None:
        if filename:
            label = os.path.splitext(os.path.basename(filename))[0]
        else:
            label = "output"
        label = f"{label}_{ctx.invoked_subcommand}"
    label = clean_label(label)

    molecule_indices = None
    if index is not None:
        molecules, molecule_indices = (
            return_objects_and_indices_from_string_index(
                list_of_objects=molecules, index=index
            )
        )

    if not isinstance(molecules, list):
        molecules = [molecules]
        if molecule_indices is not None and not isinstance(
            molecule_indices, list
        ):
            molecule_indices = [molecule_indices]

    ctx.obj["project_settings"] = project_settings
    ctx.obj["job_settings"] = job_settings
    ctx.obj["keywords"] = keywords
    ctx.obj["molecules"] = molecules
    ctx.obj["molecule_indices"] = molecule_indices
    ctx.obj["label"] = label
    ctx.obj["filename"] = filename
    ctx.obj["explicit_state_fields"] = frozenset(explicit_state_fields)


@pyscf.result_callback()
@click.pass_context
def pyscf_process_pipeline(ctx, *args, **kwargs):
    """Hand the constructed job(s) back up to ``run``/``sub``."""
    kwargs.update({"subcommand": ctx.invoked_subcommand})
    ctx.obj[ctx.info_name] = kwargs
    logger.debug(
        f"Pipeline completed for subcommand: {ctx.invoked_subcommand}"
    )
    return args[0]
