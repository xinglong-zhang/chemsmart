#!/usr/bin/env python
"""
Render derived ChemSmart PyMOL styles for documentation.

For each ``scientific_styles.py`` style, runs::

    chemsmart run mol -f INPUT.xyz visualize -s STYLE

loads the resulting ``.pse`` session in PyMOL, ray-traces, and exports PNG
images to ``docs/source/_static/pymol_styles/``. Figures are referenced from
``docs/source/pymol-visualization.rst``.

Run with the ChemSmart Python environment (not inside PyMOL)::

    python chemsmart/scripts/render_structure_styles.py INPUT.xyz

Regenerate all documentation figures::

    python chemsmart/scripts/render_structure_styles.py \\
        tests/data/StructuresTests/xyz/1-mer.xyz \\
        --output-dir docs/source/_static/pymol_styles \\
        --comic-coordinates '[[1,2],[1,5],[1,36],[1,3],[1,15],[1,8]]'

Pass ``--coordinates`` for every rendered style, or ``--comic-coordinates`` for
comic highlighting only.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path


def _ensure_chemsmart_python() -> None:
    """Require the ChemSmart interpreter, not PyMOL's bundled Python."""
    if "pymol" in sys.modules:
        raise SystemExit(
            "Run this script with the ChemSmart Python environment, not PyMOL:\n"
            "  python chemsmart/scripts/render_structure_styles.py INPUT.xyz"
        )
    try:
        import ase  # noqa: F401
    except ImportError as exc:
        raise SystemExit(
            "ChemSmart dependencies are missing in this Python environment.\n"
            "Activate the chemsmart conda/venv environment, then rerun with:\n"
            "  python chemsmart/scripts/render_structure_styles.py INPUT.xyz\n"
            "Do not invoke this script with `pymol -cq -r ...`."
        ) from exc


def _get_derived_styles() -> frozenset[str]:
    from chemsmart.jobs.mol.runner import PYMOL_SCIENTIFIC_STYLE_COMMANDS

    return frozenset(PYMOL_SCIENTIFIC_STYLE_COMMANDS)


def _default_derived_styles() -> list[str]:
    """Return derived style keys in registry order."""
    from chemsmart.jobs.mol.runner import PYMOL_SCIENTIFIC_STYLE_COMMANDS

    return list(PYMOL_SCIENTIFIC_STYLE_COMMANDS)


def _clean_label(label: str) -> str:
    from chemsmart.utils.io import clean_label

    return clean_label(label)


def _style_cli_name(style_key: str) -> str:
    """Convert internal style key to ``visualize -s`` CLI value."""
    return style_key.replace("_", "-")


def _resolve_chemsmart_command(explicit: str | None) -> list[str]:
    if explicit:
        return [explicit]
    return [sys.executable, "-m", "chemsmart.cli.main"]


def _resolve_pymol_command(explicit: str | None) -> list[str]:
    if explicit:
        return [explicit]
    found = shutil.which("pymol")
    if found:
        return [found]
    raise SystemExit(
        "PyMOL executable not found on PATH. Install PyMOL or pass --pymol."
    )


def _expected_pse_path(work_dir: Path, label: str, style_key: str) -> Path:
    return work_dir / f"{label}_{style_key}_visualization.pse"


def run_chemsmart_visualize(
    chemsmart_cmd: list[str],
    xyz_path: Path,
    style_key: str,
    work_dir: Path,
    coordinates: str | None = None,
) -> None:
    """Run ``chemsmart run mol -f … visualize -s …`` for one derived style."""
    cli_style = _style_cli_name(style_key)
    command = [
        *chemsmart_cmd,
        "run",
        "--no-scratch",
        "mol",
        "-f",
        str(xyz_path),
        "visualize",
        "-s",
        cli_style,
        "--no-trace",
    ]
    if coordinates is not None:
        command.extend(["-c", coordinates])

    print(f"Running: {' '.join(command)}")
    subprocess.run(command, cwd=work_dir, check=True)


def export_pse_png(
    pymol_cmd: list[str],
    pse_path: Path,
    png_path: Path,
    width: int = 1200,
    height: int = 1200,
    dpi: int = 300,
) -> None:
    """Load a ChemSmart ``.pse`` session in PyMOL, ray-trace, and save a PNG."""
    pse_path = pse_path.resolve()
    png_path = png_path.resolve()
    png_path.parent.mkdir(parents=True, exist_ok=True)

    command = [
        *pymol_cmd,
        "-cq",
        str(pse_path),
        "-d",
        f"ray {int(width)}, {int(height)}",
        "-d",
        f"png {png_path}, dpi={int(dpi)}",
        "-d",
        "quit",
    ]
    subprocess.run(command, check=True)
    if not png_path.is_file():
        raise FileNotFoundError(
            f"PyMOL did not write expected PNG: {png_path}"
        )


def render_all_styles(
    xyz_path: Path,
    output_dir: Path,
    work_dir: Path,
    label: str,
    chemsmart_cmd: list[str],
    pymol_cmd: list[str],
    styles: list[str],
    derived_styles: frozenset[str],
    coordinates: str | None = None,
    comic_coordinates: str | None = None,
    width: int = 1200,
    height: int = 1200,
    dpi: int = 300,
) -> dict[str, Path]:
    """Run ChemSmart visualize for each style and export PNG figures."""
    output_dir.mkdir(parents=True, exist_ok=True)
    work_dir.mkdir(parents=True, exist_ok=True)

    exported: dict[str, Path] = {}
    for style_key in styles:
        if style_key not in derived_styles:
            raise ValueError(f"Unsupported derived style: {style_key}")

        style_coordinates = coordinates
        if style_coordinates is None and style_key == "comic":
            style_coordinates = comic_coordinates

        run_chemsmart_visualize(
            chemsmart_cmd=chemsmart_cmd,
            xyz_path=xyz_path,
            style_key=style_key,
            work_dir=work_dir,
            coordinates=style_coordinates,
        )

        pse_path = _expected_pse_path(work_dir, label, style_key)
        if not pse_path.is_file():
            raise FileNotFoundError(
                f"Expected ChemSmart output not found: {pse_path}"
            )

        png_path = output_dir / f"style_{style_key}.png"
        export_pse_png(
            pymol_cmd=pymol_cmd,
            pse_path=pse_path,
            png_path=png_path,
            width=width,
            height=height,
            dpi=dpi,
        )
        exported[style_key] = png_path
        print(f"Wrote {png_path}")

    return exported


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Run ChemSmart derived visualize styles and export PNG figures "
            "for docs/source/pymol-visualization.rst."
        )
    )
    parser.add_argument(
        "xyz_file",
        type=Path,
        help="Input XYZ structure file passed to ``mol -f``.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("docs/source/_static/pymol_styles"),
        help="Directory for exported PNG images.",
    )
    parser.add_argument(
        "--work-dir",
        type=Path,
        default=None,
        help="ChemSmart working directory for ``.pse`` outputs (default: OUTPUT_DIR/chemsmart).",
    )
    parser.add_argument(
        "--label",
        default=None,
        help="ChemSmart job label stem (defaults to cleaned XYZ basename).",
    )
    parser.add_argument(
        "--chemsmart",
        default=None,
        help="Path to chemsmart executable (defaults to current Python module).",
    )
    parser.add_argument(
        "--pymol",
        default=None,
        help="Path to pymol executable (defaults to PATH lookup).",
    )
    parser.add_argument(
        "--styles",
        nargs="+",
        default=None,
        help="Derived style keys to render (default: all scientific styles).",
    )
    parser.add_argument(
        "--coordinates",
        default=None,
        help=(
            "Coordinate spec passed to ``visualize -c`` for every rendered "
            "style."
        ),
    )
    parser.add_argument(
        "--comic-coordinates",
        default=None,
        help=(
            "Coordinate spec passed to ``visualize -c`` for the comic style "
            "only."
        ),
    )
    parser.add_argument(
        "--width",
        type=int,
        default=1200,
        help="Ray-traced image width in pixels.",
    )
    parser.add_argument(
        "--height",
        type=int,
        default=1200,
        help="Ray-traced image height in pixels.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="PNG export DPI.",
    )
    return parser.parse_args(argv)


def main(argv=None):
    _ensure_chemsmart_python()
    args = parse_args(argv)
    xyz_path = args.xyz_file.resolve()
    if not xyz_path.is_file():
        raise FileNotFoundError(f"XYZ file not found: {xyz_path}")

    derived_styles = _get_derived_styles()
    styles = (
        args.styles if args.styles is not None else _default_derived_styles()
    )
    label = _clean_label(args.label or xyz_path.stem)
    output_dir = args.output_dir.resolve()
    work_dir = (args.work_dir or output_dir / "chemsmart").resolve()
    chemsmart_cmd = _resolve_chemsmart_command(args.chemsmart)
    pymol_cmd = _resolve_pymol_command(args.pymol)

    render_all_styles(
        xyz_path=xyz_path,
        output_dir=output_dir,
        work_dir=work_dir,
        label=label,
        chemsmart_cmd=chemsmart_cmd,
        pymol_cmd=pymol_cmd,
        styles=styles,
        derived_styles=derived_styles,
        coordinates=args.coordinates,
        comic_coordinates=args.comic_coordinates,
        width=args.width,
        height=args.height,
        dpi=args.dpi,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
