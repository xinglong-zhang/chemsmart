"""Shared reaction-job submit options.

Used by ``chemsmart sub/run chain … reaction``. Output analysis is not
part of this command.
"""

import logging

import click

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.reaction import PATH_SEARCH_SINGLE_STRUCTURE
from chemsmart.utils.datasets import ReactionTableEntry
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


def click_reaction_shared_options(f):
    """Submit options for ``chain … reaction``."""
    f = click.option(
        "-ts",
        "--ts-guess",
        type=click.Path(exists=True, dir_okay=False),
        default=None,
        help="QST3/NEB intermediate when -f is the reactant.",
    )(f)
    f = click.option(
        "-p",
        "--product",
        "products",
        multiple=True,
        type=click.Path(exists=True, dir_okay=False),
        help=(
            "Product geometry file. Repeatable for extra fragments. "
            "A single product without extra reactants selects path search."
        ),
    )(f)
    f = click.option(
        "-r",
        "--reactant",
        "reactants",
        multiple=True,
        type=click.Path(exists=True, dir_okay=False),
        help=(
            "Reactant geometry file. Repeatable for extra fragments. "
            "With --product, parent -f is the TS guess."
        ),
    )(f)
    return f


def merge_reaction_options(
    ctx,
    reactants=(),
    products=(),
    ts_guess=None,
):
    """Merge submit-level options with values stored on the reaction group."""
    shared = ctx.obj.get("reaction_shared") or {}
    reactants = reactants or shared.get("reactants") or ()
    products = products or shared.get("products") or ()
    ts_guess = ts_guess or shared.get("ts_guess")
    return reactants, products, ts_guess


def store_reaction_shared(ctx, reactants, products, ts_guess):
    """Record group-level reaction options for submit."""
    ctx.ensure_object(dict)
    ctx.obj["reaction_shared"] = dict(
        reactants=reactants or (),
        products=products or (),
        ts_guess=ts_guess,
    )


def load_reaction_molecule(filepath):
    """Load a single molecule from *filepath*."""
    molecule = Molecule.from_filepath(str(filepath))
    if molecule is None:
        raise click.UsageError(f"Could not read molecule from {filepath}.")
    return molecule


def apply_settings_charge_multiplicity(molecule, settings):
    """Copy charge/multiplicity from *settings* when the molecule lacks them."""
    molecule = molecule.copy()
    if molecule.charge is None and settings.charge is not None:
        molecule.charge = settings.charge
    if molecule.multiplicity is None and settings.multiplicity is not None:
        molecule.multiplicity = settings.multiplicity
    return molecule


def molecule_from_reaction_entry(entry):
    """Load a table row and apply its charge/multiplicity."""
    molecule = load_reaction_molecule(entry.filepath).copy()
    molecule.charge = int(entry.charge)
    molecule.multiplicity = int(entry.multiplicity)
    return molecule


def resolve_reaction_structures(
    parent_molecule,
    reactant_molecules=(),
    product_molecules=(),
    ts_guess_molecule=None,
):
    """Map parent/-f and role files onto reaction-job constructor kwargs.

    Dispatch:
    - ``-f`` only, or ``-f`` plus reactants without products: case 1, ``-f`` is TS
    - ``-f`` plus products without reactants: case 2, ``-f`` is reactant
    - reactants and products: case 2; ``-f`` is the QST3/NEB intermediate
    - extra fragments with a TS parent skip path search
    """
    reactant_molecules = tuple(reactant_molecules or ())
    product_molecules = tuple(product_molecules or ())
    if ts_guess_molecule is not None and not product_molecules:
        raise ValueError(
            "--ts-guess requires a product structure for path search."
        )
    kwargs = dict(
        molecule=parent_molecule,
        reactants=reactant_molecules,
        products=product_molecules,
        ts_guess=ts_guess_molecule if product_molecules else None,
    )
    if not product_molecules:
        return kwargs
    if not reactant_molecules and len(product_molecules) > 1:
        raise ValueError(PATH_SEARCH_SINGLE_STRUCTURE)
    if reactant_molecules and (
        len(reactant_molecules) > 1 or len(product_molecules) > 1
    ):
        kwargs["no_path_search"] = True
    return kwargs


def resolve_reaction_group_from_entries(entries):
    """Map one reaction_id group of table rows onto job constructor kwargs.

    A group with both reactant and product roles is case 2 when each role
    has a single geometry. Extra fragments with a ``ts`` row skip Guess.
    Multiple reactant or product rows without a ``ts`` row are an error.
    """
    ts_entries = [entry for entry in entries if entry.role == "ts"]
    reactant_entries = [entry for entry in entries if entry.role == "reactant"]
    product_entries = [entry for entry in entries if entry.role == "product"]
    if len(ts_entries) > 1:
        raise ValueError(
            f"Reaction {entries[0].reaction_id!r} has {len(ts_entries)} "
            "ts rows; at most one is allowed."
        )
    if not ts_entries and not reactant_entries:
        raise ValueError(
            f"Reaction {entries[0].reaction_id!r} needs a ts or reactant row."
        )

    reactant_molecules = tuple(
        molecule_from_reaction_entry(entry) for entry in reactant_entries
    )
    product_molecules = tuple(
        molecule_from_reaction_entry(entry) for entry in product_entries
    )
    ts_molecule = (
        molecule_from_reaction_entry(ts_entries[0]) if ts_entries else None
    )
    multi_fragment = len(reactant_molecules) > 1 or len(product_molecules) > 1
    uses_path_search = (
        bool(reactant_molecules)
        and bool(product_molecules)
        and not multi_fragment
    )

    if multi_fragment and ts_molecule is None:
        raise ValueError(PATH_SEARCH_SINGLE_STRUCTURE)

    if not uses_path_search:
        parent = (
            ts_molecule if ts_molecule is not None else reactant_molecules[0]
        )
        extra_reactants = (
            reactant_molecules
            if ts_molecule is not None
            else reactant_molecules[1:]
        )
        job_kwargs = dict(
            molecule=parent,
            reactants=extra_reactants,
            products=product_molecules,
            ts_guess=None,
        )
        if product_molecules:
            job_kwargs["no_path_search"] = True
        return job_kwargs, {
            "filepath": (
                ts_entries[0].filepath
                if ts_entries
                else reactant_entries[0].filepath
            ),
            "charge": parent.charge,
            "multiplicity": parent.multiplicity,
        }

    if ts_molecule is not None:
        return dict(
            molecule=ts_molecule,
            reactants=reactant_molecules,
            products=product_molecules,
            ts_guess=None,
        ), {
            "filepath": ts_entries[0].filepath,
            "charge": ts_molecule.charge,
            "multiplicity": ts_molecule.multiplicity,
        }

    return dict(
        molecule=reactant_molecules[0],
        reactants=reactant_molecules[1:],
        products=product_molecules,
        ts_guess=None,
    ), {
        "filepath": reactant_entries[0].filepath,
        "charge": reactant_molecules[0].charge,
        "multiplicity": reactant_molecules[0].multiplicity,
    }


def reaction_child_settings(ctx, include_neb=False):
    """Copy project opt/ts/sp (and optional neb) settings, merged with CLI."""
    project_settings = ctx.obj["project_settings"]
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]
    opt_settings = project_settings.opt_settings().merge(
        job_settings, keywords=keywords
    )
    ts_settings = project_settings.ts_settings().merge(
        job_settings, keywords=keywords
    )
    sp_settings = project_settings.sp_settings().merge(
        job_settings, keywords=keywords
    )
    settings = dict(
        opt_settings=opt_settings,
        ts_settings=ts_settings,
        sp_settings=sp_settings,
    )
    if include_neb:
        settings["neb_settings"] = project_settings.neb_settings().merge(
            job_settings, keywords=keywords
        )
    return settings


def build_reaction_job(
    ctx,
    job_cls,
    skip_completed,
    include_neb=False,
    reactants=(),
    products=(),
    ts_guess=None,
    **kwargs,
):
    """Create one reaction job from parent ``-f`` and role files."""
    filename = ctx.obj.get("filename")
    if ReactionTableEntry.is_submission_table(filename):
        raise click.UsageError(
            "Table input detected for reaction job submission. "
            "Use the 'batch' subcommand to process the table."
        )

    molecules = ctx.obj.get("molecules")
    if not molecules:
        raise click.UsageError(
            "Reaction submit requires a structure file via parent -f/--filename."
        )

    child_settings = reaction_child_settings(ctx, include_neb=include_neb)
    parent_settings = child_settings["ts_settings"]
    check_charge_and_multiplicity(parent_settings)
    parent_molecule = apply_settings_charge_multiplicity(
        molecules[-1], parent_settings
    )

    try:
        reactant_mols = tuple(
            load_reaction_molecule(path) for path in reactants
        )
        product_mols = tuple(load_reaction_molecule(path) for path in products)
        ts_guess_mol = (
            load_reaction_molecule(ts_guess) if ts_guess is not None else None
        )
        job_kwargs = resolve_reaction_structures(
            parent_molecule,
            reactant_molecules=reactant_mols,
            product_molecules=product_mols,
            ts_guess_molecule=ts_guess_mol,
        )
        label = ctx.obj["label"]
        logger.info("Creating %s job %s", job_cls.__name__, label)
        return job_cls(
            settings=parent_settings,
            label=label,
            jobrunner=ctx.obj["jobrunner"],
            skip_completed=skip_completed,
            **child_settings,
            **job_kwargs,
            **kwargs,
        )
    except ValueError as exc:
        raise click.UsageError(str(exc)) from exc


def build_reaction_batch_jobs(
    ctx,
    job_cls,
    skip_completed,
    include_neb=False,
    **kwargs,
):
    """Create one reaction job per ``reaction_id`` group in the parent table."""
    input_table_path = ctx.obj.get("filename")
    if not input_table_path:
        raise click.UsageError(
            "Batch mode requires the parent -f/--filename to specify the "
            "table file path."
        )
    try:
        entries = ReactionTableEntry.parse_reaction_table(input_table_path)
    except (FileNotFoundError, ValueError) as exc:
        raise click.UsageError(str(exc)) from exc

    child_settings = reaction_child_settings(ctx, include_neb=include_neb)
    parent_settings = child_settings["ts_settings"]
    jobrunner = ctx.obj["jobrunner"]
    jobs = []
    for reaction_id, group in ReactionTableEntry.group_by_reaction_id(
        entries
    ).items():
        try:
            job_kwargs, batch_meta = resolve_reaction_group_from_entries(group)
            job = job_cls(
                settings=parent_settings,
                label=str(reaction_id),
                jobrunner=jobrunner,
                skip_completed=skip_completed,
                **child_settings,
                **job_kwargs,
                **kwargs,
            )
        except ValueError as exc:
            raise click.UsageError(str(exc)) from exc
        parent_path = str(batch_meta["filepath"])
        job._batch_entry = {
            "kind": "reaction",
            "filepath": parent_path,
            "charge": batch_meta["charge"],
            "multiplicity": batch_meta["multiplicity"],
            "label": str(reaction_id),
            "reactants": [
                str(entry.filepath)
                for entry in group
                if entry.role == "reactant"
                and str(entry.filepath) != parent_path
            ],
            "products": [
                str(entry.filepath)
                for entry in group
                if entry.role == "product"
            ],
            "ts_guess": None,
        }
        jobs.append(job)
        logger.info("Created reaction job %s from table", reaction_id)

    logger.info("Created %s reaction jobs from table", len(jobs))
    return jobs


def replace_reaction_batch_tokens(cli_args, batch_entry):
    """Rewrite table-mode args to a single-reaction ``submit`` invocation."""
    args = list(cli_args)
    for idx in range(len(args) - 1):
        if args[idx] in {"-f", "--filename"}:
            args[idx + 1] = str(batch_entry["filepath"])
            break
    if "batch" in args:
        args[args.index("batch")] = "submit"

    def _drop_option_pair(tokens, option_names):
        idx = 0
        while idx < len(tokens):
            if tokens[idx] in option_names:
                del tokens[idx : idx + 2]
                continue
            idx += 1

    def _set_option(tokens, long_opt, short_opt, value, insert_before=None):
        if long_opt in tokens:
            pos = tokens.index(long_opt)
            if pos + 1 < len(tokens):
                tokens[pos + 1] = value
            return
        if short_opt is not None and short_opt in tokens:
            pos = tokens.index(short_opt)
            if pos + 1 < len(tokens):
                tokens[pos + 1] = value
            return
        insert_idx = len(tokens)
        if insert_before in tokens:
            insert_idx = tokens.index(insert_before)
        tokens[insert_idx:insert_idx] = [long_opt, value]

    _drop_option_pair(args, {"--reactant", "-r"})
    _drop_option_pair(args, {"--product"})
    _drop_option_pair(args, {"--ts-guess"})

    if batch_entry.get("charge") is not None:
        _set_option(
            args,
            "--charge",
            "-c",
            str(batch_entry["charge"]),
            insert_before="reaction",
        )
    if batch_entry.get("multiplicity") is not None:
        _set_option(
            args,
            "--multiplicity",
            "-m",
            str(batch_entry["multiplicity"]),
            insert_before="reaction",
        )
    if batch_entry.get("label") is not None:
        _set_option(
            args,
            "--label",
            "-l",
            str(batch_entry["label"]),
            insert_before="reaction",
        )

    extra = []
    for path in batch_entry.get("reactants") or ():
        extra.extend(["--reactant", str(path)])
    for path in batch_entry.get("products") or ():
        extra.extend(["--product", str(path)])
    ts_guess = batch_entry.get("ts_guess")
    if ts_guess:
        extra.extend(["--ts-guess", str(ts_guess)])
    insert_idx = args.index("submit") if "submit" in args else len(args)
    args[insert_idx:insert_idx] = extra
    return args
