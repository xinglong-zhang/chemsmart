import functools
import logging
import os

import click

from chemsmart.io.molecules.structure import Molecule
from chemsmart.utils.cli import MyGroup
from chemsmart.utils.utils import get_list_from_string_range

logger = logging.getLogger(__name__)


def click_pymol_options(f):
    """Common click options for PyMOL CLI."""

    @click.option(
        "-f",
        "--filename",
        type=str,
        default=None,
        help="filename from which new Gaussian input is prepared.",
    )
    @click.option(
        "-l",
        "--label",
        type=str,
        default=None,
        help="write user input filename for the job (without extension)",
    )
    @click.option(
        "-a",
        "--append-label",
        type=str,
        default=None,
        help="name to be appended to file for the job",
    )
    @click.option(
        "-i",
        "--index",
        type=str,
        default=None,
        help="Index of molecules to use; 1-based indices. "
        "Default to the last molecule structure. 1-based index.",
    )
    @click.option(
        "-P",
        "--pubchem",
        type=str,
        default=None,
        help="Queries structure from PubChem using name, smiles, cid and conformer information.",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


def click_pymol_visualization_options(f):
    """Common click options for PyMOL visualization."""

    @click.option(
        "-s",
        "--style-file",
        type=str,
        default=None,
        help="PyMOL script or style file. If not specified, defaults to zhang_group_pymol_style.py.",
    )
    @click.option(
        "-q",
        "--quiet",
        is_flag=True,
        default=True,
        help="Run PyMOL in quiet mode. Default to True.",
    )
    @click.option(
        "-c",
        "--command-line-only",
        is_flag=True,
        default=True,
        help="Run PyMOL in command line only. Default to True.",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


@click.group(cls=MyGroup)
@click_pymol_options
@click.pass_context
def mol(
    ctx,
    filename,
    label,
    append_label,
    index,
    pubchem,
):
    # obtain molecule structure
    if filename is None and pubchem is None:
        raise ValueError(
            "[filename] or [pubchem] has not been specified!\nPlease specify one of them!"
        )
    if filename and pubchem:
        raise ValueError(
            "Both [filename] and [pubchem] have been specified!\nPlease specify only one of them."
        )

    if filename:
        molecules = Molecule.from_filepath(
            filepath=filename, index=":", return_list=True
        )
        assert (
            molecules is not None
        ), f"Could not obtain molecule from {filename}!"
        logger.debug(f"Obtained molecule {molecules} from {filename}")

    if pubchem:
        molecules = Molecule.from_pubchem(identifier=pubchem, return_list=True)
        assert (
            molecules is not None
        ), f"Could not obtain molecule from PubChem {pubchem}!"
        logger.debug(f"Obtained molecule {molecules} from PubChem {pubchem}")

    # update labels
    if label is not None and append_label is not None:
        raise ValueError(
            "Only give Gaussian input filename or name to be be appended, but not both!"
        )
    if append_label is not None:
        label = os.path.splitext(os.path.basename(filename))[0]
        label = f"{label}_{append_label}"
    if label is None and append_label is None:
        label = os.path.splitext(os.path.basename(filename))[0]
        label = f"{label}_{ctx.invoked_subcommand}"

    logger.debug(f"Obtained molecules: {molecules} before applying indices")

    # if user has specified an index to use to access particular structure
    # then return that structure as a list
    if index is not None:
        logger.debug(f"Using molecule with index: {index}")
        try:
            # try to get molecule using python style string indexing,
            # but in 1-based
            from chemsmart.utils.utils import string2index_1based

            index = string2index_1based(index)
            molecules = molecules[index]
            if not isinstance(molecules, list):
                molecules = [molecules]
        except ValueError:
            # except user defined indices such as s='[1-3,28-31,34-41]'
            # or s='1-3,28-31,34-41' which cannot be parsed by string2index_1based
            index = get_list_from_string_range(index)
            molecules = [molecules[i - 1] for i in index]

    logger.debug(f"Obtained molecules: {molecules}")

    # store objects
    ctx.obj["molecules"] = (
        molecules  # molecules as a list, as some jobs requires all structures to be used
    )
    ctx.obj["label"] = label
    ctx.obj["filename"] = filename


@mol.result_callback()
@click.pass_context
def mol_process_pipeline(ctx, *args, **kwargs):
    kwargs.update({"subcommand": ctx.invoked_subcommand})
    ctx.obj[ctx.info_name] = kwargs
    return args[0]
