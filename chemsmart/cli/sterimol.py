"""
Sterimol parameter CLI.

Registered as ``chemsmart run sterimol`` and ``chemsmart sub sterimol``.
Selects a substituent from a structure file and returns a SterimolJob, so a
wrong selection fails here, before any job is submitted.
"""

import logging
import os

import click

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.sterimol import SterimolJob
from chemsmart.utils.cli import MyCommand
from chemsmart.utils.io import clean_label
from chemsmart.utils.utils import get_list_from_string_range

logger = logging.getLogger(__name__)


def _read_molecule(filename, index):
    """Return the single molecule that ``index`` selects from ``filename``."""
    try:
        # Without return_list, a slice such as ":" silently gives its first
        # structure instead of reporting that it selects several.
        molecules = Molecule.from_filepath(
            filename, index=index, return_list=True
        )
    except IndexError as err:
        raise ValueError(
            f"--index {index} does not select a structure in {filename}."
        ) from err
    if molecules is None:
        raise ValueError(f"No structure could be read from {filename}.")
    if len(molecules) != 1:
        raise ValueError(
            f"--index {index} selects {len(molecules)} structures. "
            "Select a single structure."
        )
    return molecules[0]


def _parse_explicit_atoms(molecule, atoms, atom1, atom2):
    """Return the sorted atoms of an explicit ``--atoms`` selection."""
    selected = sorted(set(get_list_from_string_range(atoms)))
    invalid = [
        idx
        for idx in (atom1, atom2, *selected)
        if not 1 <= idx <= molecule.num_atoms
    ]
    if invalid:
        raise ValueError(
            f"Atom indices {invalid} are outside the valid 1-based range "
            f"1..{molecule.num_atoms}."
        )
    if atom1 in selected:
        raise ValueError(
            f"--atoms cannot contain the first axis atom ({atom1}), because "
            "it is not part of the substituent."
        )
    if atom2 not in selected:
        raise ValueError(
            f"--atoms must contain the second axis atom ({atom2}), because "
            "it is the first atom of the substituent."
        )
    return selected


def _select_auto_fragment(molecule, atom1, atom2):
    """Return the atoms that cutting the ``atom1``-``atom2`` bond separates."""
    try:
        return molecule.get_substituent_atoms(atom1, atom2)
    except ValueError as err:
        raise ValueError(f"{err} Select the atoms with --atoms.") from err


def _default_label(filename, index, atom1, atom2):
    """Return a label that differs between axes and selected structures."""
    stem = os.path.splitext(os.path.basename(filename))[0]
    label = f"{stem}_{atom1}_{atom2}"
    if index != "-1":
        label = f"{label}_i{index}"
    return clean_label(label)


@click.command(name="sterimol", cls=MyCommand)
@click.option(
    "-f",
    "--filename",
    "--file",
    "filename",
    type=click.Path(exists=True, dir_okay=False),
    required=True,
    help="Geometry or calculation output file.",
)
@click.option(
    "--axis",
    nargs=2,
    type=click.IntRange(min=1),
    required=True,
    metavar="ATOM1 ATOM2",
    help="1-based parent atom and first substituent atom. "
    "The axis points from ATOM1 to ATOM2.",
)
@click.option(
    "--auto-fragment",
    is_flag=True,
    help="Select the substituent by cutting the ATOM1-ATOM2 bond.",
)
@click.option(
    "--atoms",
    type=str,
    default=None,
    help="Explicit 1-based substituent atoms, including ATOM2, "
    "e.g. '6,7,12-18'. Needed when the bond is in a ring or is not detected.",
)
@click.option(
    "-i",
    "--index",
    type=str,
    default="-1",
    show_default=True,
    help="1-based index of the structure to use; -1 is the last one.",
)
@click.option(
    "-l",
    "--label",
    type=str,
    default=None,
    help="Job label. Defaults to the file name, the axis, and any --index.",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(dir_okay=False),
    default=None,
    help="JSON report to write. Defaults to {label}_sterimol.json.",
)
def sterimol(filename, axis, auto_fragment, atoms, index, label, output):
    """
    Calculate Sterimol parameters L, B1 and B5 of a substituent.

    \b
    Run on this machine:
        chemsmart run sterimol -f molecule.xyz --axis 5 6 --auto-fragment
        chemsmart run sterimol -f molecule.xyz --axis 5 6 --atoms "6,7,12-18"
        chemsmart run sterimol -f traj.xyz -i 3 --axis 5 6 --auto-fragment \\
            -o sterimol.json

    \b
    Submit to a compute node:
        chemsmart sub -s SERVER sterimol -f molecule.xyz --axis 5 6 \\
            --auto-fragment
    """
    if auto_fragment == (atoms is not None):
        raise click.UsageError(
            "Use exactly one of --auto-fragment or --atoms."
        )

    atom1, atom2 = axis
    try:
        molecule = _read_molecule(filename, index)
        if auto_fragment:
            selection = "auto-fragment"
            selected = _select_auto_fragment(molecule, atom1, atom2)
        else:
            selection = "explicit"
            selected = _parse_explicit_atoms(molecule, atoms, atom1, atom2)
    except (ValueError, IndexError, TypeError) as err:
        raise click.ClickException(str(err)) from err

    return SterimolJob(
        molecule=molecule,
        atom1=atom1,
        atom2=atom2,
        substituent_atoms=selected,
        filename=filename,
        index=index,
        selection=selection,
        output=output,
        label=label or _default_label(filename, index, atom1, atom2),
    )
