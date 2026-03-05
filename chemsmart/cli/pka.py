import logging

import click

from chemsmart.utils.io import get_program_type_from_file

logger = logging.getLogger(__name__)


def _detect_program_from_outputs(filepaths):
    """Detect and validate a unique QC program across all provided outputs.

    Delegates per-file detection to ``get_program_type_from_file`` so
    output-format parsing remains centralized in ``chemsmart.utils.io``.
    """
    detected = {}
    programs = set()
    for fp in filepaths:
        p = get_program_type_from_file(fp)
        detected[fp] = p
        if p != "unknown":
            programs.add(p)

    if not programs:
        raise click.UsageError(
            "Could not detect output-file program type from supplied files. "
            "Supported for pKa analysis: Gaussian and ORCA output files."
        )

    if len(programs) > 1:
        pairs = ", ".join([f"{k}: {v}" for k, v in detected.items()])
        raise click.UsageError(
            "Supplied files contain mixed program types. "
            "Please use outputs from a single QC program (Gaussian or ORCA).\n"
            f"Detected: {pairs}"
        )

    program = next(iter(programs))
    if program not in {"gaussian", "orca"}:
        raise click.UsageError(
            f"Detected unsupported program '{program}'. "
            "Only Gaussian and ORCA are supported for pKa output analysis."
        )
    return program


@click.command(name="pka")
@click.option(
    "-rp",
    "--reference-pka",
    type=float,
    required=True,
    help="Experimental pKa of the reference acid HB used in the proton-exchange scheme.",
)
@click.option(
    "-ha",
    "--ha-output",
    type=click.Path(exists=True),
    required=True,
    help="Path to HA (protonated acid) gas-phase optimization/frequency output file.",
)
@click.option(
    "-a",
    "--a-output",
    type=click.Path(exists=True),
    required=True,
    help="Path to A- (conjugate base) gas-phase optimization/frequency output file.",
)
@click.option(
    "-hb",
    "--hb-output",
    type=click.Path(exists=True),
    required=True,
    help="Path to HB (reference acid) gas-phase optimization/frequency output file.",
)
@click.option(
    "-b",
    "--b-output",
    type=click.Path(exists=True),
    required=True,
    help="Path to B- (reference conjugate base) gas-phase optimization/frequency output file.",
)
@click.option(
    "-has",
    "--ha-solv-output",
    type=click.Path(exists=True),
    required=True,
    help="Path to HA solvent single-point output file.",
)
@click.option(
    "-as",
    "--a-solv-output",
    type=click.Path(exists=True),
    required=True,
    help="Path to A- solvent single-point output file.",
)
@click.option(
    "-hbs",
    "--hb-solv-output",
    type=click.Path(exists=True),
    required=True,
    help="Path to HB solvent single-point output file.",
)
@click.option(
    "-bs",
    "--b-solv-output",
    type=click.Path(exists=True),
    required=True,
    help="Path to B- solvent single-point output file.",
)
@click.option(
    "-T",
    "--temperature",
    type=float,
    default=298.15,
    show_default=True,
    help="Temperature in Kelvin used in thermochemistry corrections and pKa calculation.",
)
@click.option(
    "-conc",
    "--concentration",
    type=float,
    default=1.0,
    show_default=True,
    help="Concentration in mol/L used in quasi-harmonic thermochemistry.",
)
@click.option(
    "-csg",
    "--cutoff-entropy-grimme",
    type=float,
    default=100.0,
    show_default=True,
    help="Cutoff frequency (cm^-1) for entropy using Grimme's quasi-RRHO method.",
)
@click.option(
    "-ch",
    "--cutoff-enthalpy",
    type=float,
    default=100.0,
    show_default=True,
    help="Cutoff frequency (cm^-1) for enthalpy using Head-Gordon's quasi-RRHO method.",
)
@click.pass_context
def pka(
    ctx,
    reference_pka,
    ha_output,
    a_output,
    hb_output,
    b_output,
    ha_solv_output,
    a_solv_output,
    hb_solv_output,
    b_solv_output,
    temperature,
    concentration,
    cutoff_entropy_grimme,
    cutoff_enthalpy,
):
    """Analyze pKa from existing outputs without selecting Gaussian/ORCA subcommands.

    This mode is strictly post-processing and does not submit jobs.
    The QC backend is auto-detected from output-file signatures.
    """
    output_files = [
        ha_output,
        a_output,
        hb_output,
        b_output,
        ha_solv_output,
        a_solv_output,
        hb_solv_output,
        b_solv_output,
    ]

    program = _detect_program_from_outputs(output_files)
    logger.info(f"Detected output program: {program}")

    kwargs = dict(
        ha_gas_file=ha_output,
        a_gas_file=a_output,
        hb_gas_file=hb_output,
        b_gas_file=b_output,
        ha_solv_file=ha_solv_output,
        a_solv_file=a_solv_output,
        hb_solv_file=hb_solv_output,
        b_solv_file=b_solv_output,
        pka_reference=reference_pka,
        temperature=temperature,
        concentration=concentration,
        cutoff_entropy_grimme=cutoff_entropy_grimme,
        cutoff_enthalpy=cutoff_enthalpy,
    )

    if program == "gaussian":
        from chemsmart.io.gaussian.output import Gaussian16pKaOutput

        Gaussian16pKaOutput.print_pka_summary(**kwargs)
    else:
        from chemsmart.io.orca.output import ORCApKaOutput

        ORCApKaOutput.print_pka_summary(**kwargs)

    # Return None so run.process_pipeline skips jobrunner execution.
    return None
