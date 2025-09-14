#!/usr/bin/env python
import logging
import os

import click

from chemsmart.io.gaussian.output import Gaussian16Output
from chemsmart.utils.logger import create_logger

logger = logging.getLogger(__name__)
os.environ["OMP_NUM_THREADS"] = "1"


@click.command()
@click.option(
    "-d",
    "--directory",
    default=None,
    help="Directory in which to convert files.",
)
@click.option(
    "-t",
    "--type",
    default=None,
    help="Type of file to be converted, if direcotry is specified.",
)
@click.option(
    "-f", "--filename", default=None, help="Input filename to be converted."
)
@click.option(
    "-v",
    "--vib-mode",
    default=1,
    type=int,
    help="Vibrational mode to visualize (1-indexed).",
)
@click.option("--amp", default=0.6, type=float, help="Peak displacement in Å.")
@click.option(
    "--nframes", default=40, type=int, help="Number of animation frames."
)
@click.option(
    "--interval-ms",
    default=40,
    type=int,
    help="Frame interval (ms) in HTML animation.",
)
@click.option(
    "--html",
    "out_html",
    type=click.Path(dir_okay=False),
    help="Write interactive HTML to this path.",
)
@click.option(
    "--xyz",
    "out_xyz",
    type=click.Path(dir_okay=False),
    help="Write multi-model XYZ trajectory to this path.",
)
@click.option(
    "--open",
    "open_after",
    is_flag=True,
    help="Open the HTML in your default browser.",
)
def entry_point(
    directory,
    type,
    filename,
    vib_mode,
    amp,
    nframes,
    interval_ms,
    out_html,
    out_xyz,
    open_after,
):
    """Script for visualizing vibrational modes in a QM calculation outputfile.
    Can write either an interactive HTML file (using 3Dmol.js) or an XYZ trajectory file.
    Interactive HTML (recommended from terminal):
        file_visualizer.py -f CN_attack_quintet_from_triplet_ts.log -v 1 --html vib_mode01.html --open
    This writes vib_mode01.html and opens it in your browser (loads 3Dmol.js from CDN).

    XYZ trajectory (for Avogadro/VMD/ChimeraX):
        file_visualizer.py -f CN_attack_quintet_from_triplet_ts.log -v 1 --xyz mode_001.xyz
    Then open the file and play the trajectory. structures in different formats.
    """
    create_logger()
    g16 = Gaussian16Output(filename) if filename else None
    if g16 is None:
        raise click.ClickException("Please provide -f/--filename")

    # Prefer HTML if requested; otherwise XYZ; otherwise emit a hint
    if out_html:
        path = g16.save_vibration_html(
            vib_mode,
            amp=amp,
            nframes=nframes,
            interval_ms=interval_ms,
            outfile=out_html,
        )
        click.echo(f"Wrote interactive viewer: {path}")
        if open_after:
            import os
            import webbrowser

            webbrowser.open(f"file://{os.path.abspath(path)}")
    elif out_xyz:
        path = g16.write_mode_xyz(
            vib_mode, amp=amp, nframes=nframes, outfile=out_xyz
        )
        click.echo(f"Wrote XYZ trajectory: {path}")
    else:
        click.echo(
            "Tip: use --html path.html (interactive) or --xyz path.xyz (for Avogadro/VMD)."
        )


if __name__ == "__main__":
    entry_point()
