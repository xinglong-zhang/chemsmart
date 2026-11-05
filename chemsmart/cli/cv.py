"""CV reports: labelled geometry changes, optional CSV export, no QC jobs."""

import csv
import json
from pathlib import Path

import click
import numpy as np

from chemsmart.analysis.collective_variables import CollectiveVariables


def parse_geometry_pair(ctx, param, value):
    """Parse a bracketed pair while preserving comparison order."""
    if value is None:
        return None
    try:
        pair = json.loads(value)
    except (ValueError, TypeError):
        pair = None
    if (
        not isinstance(pair, list)
        or len(pair) != 2
        or any(type(index) is not int for index in pair)
    ):
        raise click.BadParameter(
            'Use two integer geometry numbers, e.g. "[50,10]".',
            ctx=ctx,
            param=param,
        )
    # Keep the order entered: [50,10] asks for geometry 10 minus 50.
    return tuple(pair)


class CVReport:
    """Build reports from stored profiles without recalculating geometries."""

    def __init__(self, analysis, geometries=None):
        self.analysis = analysis
        count = len(analysis.molecules)
        # A user-selected pair overrides the search across the trajectory.
        self.explicit = geometries is not None
        self.geometries = geometries or (1, count)
        if any(i < 1 or i > count for i in self.geometries):
            raise ValueError(f"Geometry numbers must be in 1..{count}.")

    @staticmethod
    def difference(first, second, circular):
        """Return second minus first, wrapped for dihedrals."""
        delta = second - first
        # For a torsion, 179 to -179 degrees is a change of +2 degrees.
        return (delta + 180.0) % 360.0 - 180.0 if circular else delta

    def rows(self, profile, category, trajectory=False):
        """Return one labelled comparison per coordinate, including NaNs."""
        data = profile[category]
        values = data["values_matrix"]
        circular = category == "dihedrals"
        symbols = self.analysis.molecules[0].symbols
        records = []
        # Each column follows one atom pair, angle or torsion across frames.
        for column, atoms in enumerate(data["atom_indices"]):
            # Reader positions start at 1; matrix rows start at 0.
            a, b = (i - 1 for i in self.geometries)
            series = values[:, column]
            if trajectory and not self.explicit:
                # Compare each finite value with all later finite values.
                valid = np.flatnonzero(np.isfinite(series))
                best = -1.0
                if len(valid) == 1:
                    # One valid frame gives a self-comparison with zero change.
                    a = b = int(valid[0])
                for offset, first in enumerate(valid[:-1]):
                    later = valid[offset + 1 :]
                    differences = self.difference(
                        series[first], series[later], circular
                    )
                    winner = int(np.argmax(np.abs(differences)))
                    magnitude = abs(differences[winner])
                    if magnitude > best:
                        best = magnitude
                        a, b = int(first), int(later[winner])
            start, end = series[a], series[b]
            change = self.difference(start, end, circular)
            # Use the same columns for values and connectivity. An unbonded
            # pair can still have a measurable distance in either frame.
            records.append(
                {
                    "atoms": "-".join(f"{symbols[i - 1]}{i}" for i in atoms),
                    "geometry_first": a + 1,
                    "geometry_second": b + 1,
                    "first_value": float(start),
                    "second_value": float(end),
                    "signed_change": float(change),
                    "absolute_change": float(abs(change)),
                    "unit": (
                        "angstrom" if category == "distances" else "degree"
                    ),
                    "connected_first": bool(
                        data["connected_matrix"][a, column]
                    ),
                    "connected_second": bool(
                        data["connected_matrix"][b, column]
                    ),
                }
            )
        return records

    @staticmethod
    def ranked(rows):
        """Rank finite changes by magnitude; preserve atom order for ties."""
        # Keep the signed change in each record, but rank by its magnitude.
        return sorted(
            (r for r in rows if np.isfinite(r["absolute_change"])),
            key=lambda r: -r["absolute_change"],
        )

    def reports(self, profile, options):
        """Yield requested reports; largest-per-atom uses endpoint changes."""
        # Each category offers the same four reports under its own flags.
        for category, largest, endpoint, ranked, endpoint_ranked in (
            ("distances", "dl", "epd", "rd", "eprd"),
            ("angles", "al", "epa", "ra", "epra"),
            ("dihedrals", "dihl", "epdi", "rdi", "eprdi"),
        ):
            if not any(
                options[k]
                for k in (largest, endpoint, ranked, endpoint_ranked)
            ):
                continue
            rows = self.rows(profile, category)
            if options[endpoint]:
                yield endpoint, rows
            if options[endpoint_ranked]:
                yield endpoint_ranked, self.ranked(rows)[:5]
            if options[ranked]:
                yield ranked, self.ranked(self.rows(profile, category, True))[
                    :5
                ]
            if options[largest]:
                result = []
                # A manual group restricts the atoms receiving maxima.
                atoms = profile.get(
                    "selection_indices",
                    range(1, len(self.analysis.molecules[0].symbols) + 1),
                )
                for atom in atoms:
                    label = (
                        f"{self.analysis.molecules[0].symbols[atom - 1]}{atom}"
                    )
                    # Include unbonded pairs when ranking distance changes.
                    # This lets a forming bond appear among the candidates.
                    candidates = self.ranked(
                        [r for r in rows if label in r["atoms"].split("-")]
                    )
                    if candidates:
                        # Keep every tie for this atom's largest change.
                        maximum = candidates[0]["absolute_change"]
                        result.extend(
                            dict(r, focus_atom=label)
                            for r in candidates
                            if r["absolute_change"] == maximum
                        )
                    else:
                        a, b = self.geometries
                        result.append(
                            dict(
                                focus_atom=label,
                                atoms="No valid candidate",
                                geometry_first=a,
                                geometry_second=b,
                                first_value=np.nan,
                                second_value=np.nan,
                                signed_change=np.nan,
                                absolute_change=np.nan,
                                unit=(
                                    "angstrom"
                                    if category == "distances"
                                    else "degree"
                                ),
                                connected_first=False,
                                connected_second=False,
                            )
                        )
                yield largest, result


# The single-letter file option is -f. Multi-character options use --.
@click.command()
@click.option(
    "-f",
    "--file",
    "filepath",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
)
@click.option(
    "--sg",
    "--geometries",
    "geometries",
    type=str,
    callback=parse_geometry_pair,
    metavar='"[FIRST,SECOND]"',
    default=None,
    help=(
        "Bracketed pair of 1-based reader positions; "
        "overrides every comparison."
    ),
)
@click.option(
    "--ms",
    "--manual-select",
    "manual_select",
    multiple=True,
    help=(
        "Ordered comma-separated labels, e.g. C13,C19,C28. "
        "Repeat for independent groups."
    ),
)
@click.option(
    "--central-bond-only",
    is_flag=True,
    help="For two-atom selections, require the central dihedral bond.",
)
@click.option(
    "--out",
    type=click.Path(file_okay=False),
    help="New directory for CSV reports (never overwrite a directory).",
)
@click.option(
    "--dl",
    is_flag=True,
    help="Largest absolute endpoint distance change per atom.",
)
@click.option(
    "--al",
    is_flag=True,
    help="Largest absolute endpoint angle change per atom.",
)
@click.option(
    "--dihl",
    is_flag=True,
    help="Largest absolute endpoint dihedral change per atom.",
)
@click.option("--epd", is_flag=True, help="All distance endpoint changes.")
@click.option("--epa", is_flag=True, help="All angle endpoint changes.")
@click.option("--epdi", is_flag=True, help="All dihedral endpoint changes.")
@click.option(
    "--rd",
    is_flag=True,
    help="Top five distance changes across geometry pairs.",
)
@click.option(
    "--ra", is_flag=True, help="Top five angle changes across geometry pairs."
)
@click.option(
    "--rdi",
    is_flag=True,
    help="Top five dihedral changes across geometry pairs.",
)
@click.option(
    "--eprd", is_flag=True, help="Top five endpoint distance changes."
)
@click.option("--epra", is_flag=True, help="Top five endpoint angle changes.")
@click.option(
    "--eprdi", is_flag=True, help="Top five endpoint dihedral changes."
)
def cv(filepath, geometries, manual_select, central_bond_only, out, **options):
    """Screen geometric changes. Default: all endpoint-change reports.

    Rankings use absolute changes; signed changes are second minus first.
    Distances include non-bonded pairs. Angular candidates require a complete
    connected path in at least one loaded geometry. Geometry numbers refer to
    reader output, not necessarily accepted IRC points. 
    """
    if central_bond_only and not manual_select:
        raise click.UsageError("--central-bond-only requires --ms.")
    if out and Path(out).exists():
        raise click.ClickException(
            "Output directory already exists; choose a new one."
        )
    # With no report flags, list the endpoint changes in all three categories.
    if not any(options.values()):
        options.update(epd=True, epa=True, epdi=True)
    try:
        # Build the profiles once so all reports use the same measurements.
        analysis = CollectiveVariables.from_filepath(filepath)
        report = CVReport(analysis, geometries)
        # Keep each repeated --ms group separate, preserving its atom order.
        profiles = (
            analysis.select(
                [
                    [label.strip() for label in group.split(",")]
                    for group in manual_select
                ],
                central_bond_only=central_bond_only,
            )
            if manual_select
            else [analysis.profile]
        )
    except (ValueError, TypeError, IndexError, OSError) as exc:
        raise click.ClickException(str(exc)) from exc
    if out:
        Path(out).mkdir(parents=True, exist_ok=False)
    click.echo(f"Loaded {len(analysis.molecules)} geometries from {filepath}")
    click.echo(
        f"Endpoint comparison: {report.geometries[0]} -> "
        f"{report.geometries[1]}"
    )
    click.echo(
        "Geometry numbers are 1-based reader positions, "
        "not accepted IRC point numbers."
    )
    click.echo(
        "Rankings use absolute changes. Distances include all pairs, "
        "including non-bonded pairs."
    )
    for group_number, profile in enumerate(profiles, 1):
        label = ",".join(profile.get("selection_labels", ())) or "all atoms"
        for name, rows in report.reports(profile, options):
            click.echo(f"\n[{name}] {label}")
            if not rows:
                click.echo(
                    "No eligible finite changes."
                    if name not in ("epd", "epa", "epdi")
                    else "No eligible candidates."
                )
            for row in rows:
                unit = "Å" if row["unit"] == "angstrom" else "°"
                focus = f"{row['focus_atom']}: " if "focus_atom" in row else ""
                click.echo(
                    f"{focus}{row['atoms']} "
                    f"{row['absolute_change']:.5f} {unit} "
                    f"(signed {row['signed_change']:+.5f}; "
                    f"geometries {row['geometry_first']} -> "
                    f"{row['geometry_second']})"
                )
            # Export the same rows shown in the terminal, with the two measured
            # values and connectivity flags. Write a header for empty reports.
            if out:
                fields = (
                    list(rows[0])
                    if rows
                    else [
                        "atoms",
                        "geometry_first",
                        "geometry_second",
                        "first_value",
                        "second_value",
                        "signed_change",
                        "absolute_change",
                        "unit",
                        "connected_first",
                        "connected_second",
                    ]
                )
                with (Path(out) / f"group_{group_number}_{name}.csv").open(
                    "w", newline="", encoding="utf-8"
                ) as stream:
                    writer = csv.DictWriter(stream, fieldnames=fields)
                    writer.writeheader()
                    writer.writerows(rows)
