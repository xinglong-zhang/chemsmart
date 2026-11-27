"""Condensed Fukui index analysis from charge-state outputs.

Post-processes finished neutral / radical-cation / radical-anion calculations.
Job submission lives under ``chemsmart sub … gaussian … fukui`` or
``chemsmart sub … orca … fukui``; analysis is exposed as
``chemsmart run fukui``.
"""

from __future__ import annotations

import logging
from pathlib import Path

from chemsmart.io.gaussian.output import Gaussian16WBIOutput
from chemsmart.io.orca.output import ORCAOutput
from chemsmart.utils.io import get_program_type_from_file

logger = logging.getLogger(__name__)

FUKUI_MODES = ("mulliken", "nbo", "hirshfeld", "cm5")
ORCA_FUKUI_MODES = ("mulliken", "hirshfeld")

_FUKUI_OUTPUT_EXTENSIONS = (".log", ".out", ".LOG", ".OUT")


def radical_ion_charge_and_multiplicity(charge, multiplicity, delta_charge):
    """Return charge and multiplicity after changing molecular charge.

    ``delta_charge`` is the change in molecular charge (``+1`` for the
    radical cation, ``-1`` for the radical anion). Closed-shell singlets
    become doublets; other multiplicities decrease by one.
    """
    ion_charge = charge + delta_charge
    if multiplicity == 1:
        ion_multiplicity = 2
    else:
        ion_multiplicity = multiplicity - 1
    return ion_charge, ion_multiplicity


def discover_fukui_companion_outputs(neutral_filename):
    """Discover radical-cation / radical-anion outputs from a neutral file.

    Expects job labels ``<base>_n``, ``<base>_rc``, and ``<base>_ra`` as written
    by ``GaussianFukuiJob`` and ``ORCAFukuiJob``.
    """
    path = Path(neutral_filename)
    stem = path.stem
    if stem.endswith("_n"):
        base = stem[: -len("_n")]
    else:
        base = stem

    directory = path.parent
    cation = _find_companion_output(directory, f"{base}_rc")
    anion = _find_companion_output(directory, f"{base}_ra")
    return {"radical_cation": cation, "radical_anion": anion}


def _find_companion_output(directory, stem):
    for ext in _FUKUI_OUTPUT_EXTENSIONS:
        candidate = directory / f"{stem}{ext}"
        if candidate.is_file():
            return str(candidate)
    return None


def _load_output(filename):
    program = get_program_type_from_file(filename)
    if program == "gaussian":
        return Gaussian16WBIOutput(filename), program
    if program == "orca":
        return ORCAOutput(filename), program
    raise TypeError(f"File {filename} is of unknown filetype.")


def _charges_for_mode(output, mode, program):
    if mode == "mulliken":
        return output.mulliken_atomic_charges
    if mode == "nbo":
        if program != "gaussian":
            raise ValueError(
                "NBO charges are only available for Gaussian outputs."
            )
        return output.natural_charges
    if mode == "hirshfeld":
        return output.hirshfeld_charges
    if mode == "cm5":
        if program != "gaussian":
            raise ValueError(
                "CM5 charges are only available for Gaussian outputs."
            )
        return output.hirshfeld_cm5_charges
    raise ValueError(
        f"Unknown mode {mode}. Supported modes are: "
        f"{', '.join(FUKUI_MODES)}."
    )


def _emit(message, lines):
    if lines is None:
        logger.info(message)
    else:
        lines.append(message)


def _log_charges(title, charges, lines):
    _emit(f"\n{title}", lines)
    for key, value in charges.items():
        _emit(f"{key:<6}  :  {value:>8.3f}", lines)
    _emit("\n", lines)


def _write_fukui_output(output, lines):
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def analyze_fukui(
    neutral_filename,
    radical_cation_filename=None,
    radical_anion_filename=None,
    mode="mulliken",
    output=None,
):
    """Compute and log condensed Fukui indices from existing outputs.

    Parameters
    ----------
    neutral_filename
        Output file for the neutral system.
    radical_cation_filename, radical_anion_filename
        Optional ion outputs. At least one must be provided.
    mode
        Charge partitioning: ``mulliken``, ``nbo``, ``hirshfeld``, or ``cm5``.
    output
        Optional path to write results. When set, results are written to this
        file and not logged.
    """
    mode = mode.lower()
    if mode not in FUKUI_MODES:
        raise ValueError(
            f"Unknown mode {mode}. Supported modes are: "
            f"{', '.join(FUKUI_MODES)}."
        )
    if radical_cation_filename is None and radical_anion_filename is None:
        raise ValueError(
            "At least one of radical cation or radical anion files must "
            "be provided."
        )

    lines = [] if output is not None else None
    neutral_output, program = _load_output(neutral_filename)

    radical_cation_output = None
    if radical_cation_filename is not None:
        radical_cation_output, _ = _load_output(radical_cation_filename)

    radical_anion_output = None
    if radical_anion_filename is not None:
        radical_anion_output, _ = _load_output(radical_anion_filename)

    if None not in (
        neutral_output,
        radical_cation_output,
        radical_anion_output,
    ):
        ionization_energy = (
            radical_cation_output.energies[-1] - neutral_output.energies[-1]
        )
        affinity_energy = (
            neutral_output.energies[-1] - radical_anion_output.energies[-1]
        )
        chemical_potential = -0.5 * (ionization_energy + affinity_energy)
        chemical_hardness = ionization_energy - affinity_energy
        if abs(chemical_hardness) < 1e-12:
            warning = (
                "Chemical hardness is effectively zero; global "
                "electrophilicity index cannot be computed to avoid "
                "division by zero."
            )
            if lines is None:
                logger.warning(warning)
            else:
                lines.append(warning)
            global_electrophilicity_index = None
        else:
            global_electrophilicity_index = chemical_potential**2 / (
                2 * chemical_hardness
            )
        _emit(f"Ionization energy = {ionization_energy}", lines)
        _emit(f"Electron affinity energy = {affinity_energy}", lines)
        _emit(f"Chemical potential = {chemical_potential}", lines)
        _emit(f"Chemical hardness = {chemical_hardness}", lines)
        _emit(
            f"Global electrophilicity index = {global_electrophilicity_index}",
            lines,
        )

    mode_labels = {
        "mulliken": "Mulliken",
        "nbo": "NBO",
        "hirshfeld": "Hirshfeld",
        "cm5": "CM5",
    }
    _emit(
        f"\nUsing {mode_labels[mode]} Charges for computing Fukui "
        "Reactivity Indices.",
        lines,
    )

    charge_for_neutral = _charges_for_mode(neutral_output, mode, program)
    charge_for_radical_cation = None
    if radical_cation_output is not None:
        charge_for_radical_cation = _charges_for_mode(
            radical_cation_output, mode, program
        )
    charge_for_radical_anion = None
    if radical_anion_output is not None:
        charge_for_radical_anion = _charges_for_mode(
            radical_anion_output, mode, program
        )

    _log_charges("Neutral System Charges:", charge_for_neutral, lines)
    if charge_for_radical_cation is not None:
        _log_charges(
            "Radical Cationic System Charges:",
            charge_for_radical_cation,
            lines,
        )
    if charge_for_radical_anion is not None:
        _log_charges(
            "Radical Anionic System Charges:",
            charge_for_radical_anion,
            lines,
        )

    _emit(
        "\nAtom        Fukui Minus (f-)   Fukui Plus(f+)    "
        "Fukui Zero(f0)    Fukui Dual Descriptor(f(2))",
        lines,
    )
    for key, value in charge_for_neutral.items():
        if charge_for_radical_cation is not None:
            fukui_minus = charge_for_radical_cation[key] - value
        else:
            fukui_minus = 0.0
        if charge_for_radical_anion is not None:
            fukui_plus = value - charge_for_radical_anion[key]
        else:
            fukui_plus = 0.0
        if (
            charge_for_radical_cation is not None
            and charge_for_radical_anion is not None
        ):
            fukui_zero = 0.5 * (
                charge_for_radical_cation[key] - charge_for_radical_anion[key]
            )
            fukui_dual = fukui_plus - fukui_minus
        else:
            fukui_zero = 0.0
            fukui_dual = 0.0
        _emit(
            f"{key:<4}        {fukui_minus:>12.3f}     {fukui_plus:>12.3f}     "
            f"{fukui_zero:>12.3f}     {fukui_dual:>12.3f}",
            lines,
        )

    if output is not None:
        _write_fukui_output(output, lines)
