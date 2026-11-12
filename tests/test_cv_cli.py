"""Test CV reports through Click and the multi-frame XYZ reader."""

import csv

import pytest
from click.testing import CliRunner

from chemsmart.cli.cv import cv

# The central distance is 1.45, 1.60 and 1.50 Å. The torsion crosses ±180°.
TRAJECTORY_XYZ = """4
geometry 1
C -0.75 1.30 0
C 0 0 0
C 1.45 0 0
C 2.20 -1.299802003703 0.022688128368
4
geometry 2
C -0.75 1.30 0
C 0 0 0
C 1.60 0 0
C 2.35 -1.125833024920 0.65
4
geometry 3
C -0.75 1.30 0
C 0 0 0
C 1.50 0 0
C 2.25 -1.299802003703 -0.022688128368
"""


def test_default_endpoint_reports(write_tmp_file):
    trajectory = write_tmp_file("path.xyz", TRAJECTORY_XYZ)
    result = CliRunner().invoke(cv, ["-f", trajectory])
    assert result.exit_code == 0, result.output
    headers = [
        line for line in result.output.splitlines() if line.startswith("[")
    ]
    assert headers == [
        "[epd] all atoms",
        "[epa] all atoms",
        "[epdi] all atoms",
    ]
    bond = "C2-C3 0.05000 Å (signed +0.05000; geometries 1 -> 3)"
    assert bond in result.output
    # The endpoints +179° and -179° are two degrees apart.
    assert (
        "C1-C2-C3-C4 2.00000 ° (signed +2.00000; geometries 1 -> 3)"
        in result.output
    )


def test_all_reports_and_csv_values(write_tmp_file, tmp_path):
    trajectory = write_tmp_file("path.xyz", TRAJECTORY_XYZ)
    names = "dl al dihl epd epa epdi rd ra rdi eprd epra eprdi".split()
    out = tmp_path / "reports"
    flags = [f"--{name}" for name in names]
    result = CliRunner().invoke(
        cv, ["--file", trajectory, *flags, "--out", str(out)]
    )
    assert result.exit_code == 0, result.output
    assert {path.name for path in out.iterdir()} == {
        f"group_1_{name}.csv" for name in names
    }
    reports = {}
    for name in names:
        path = out / f"group_1_{name}.csv"
        with path.open(newline="", encoding="utf-8") as stream:
            reports[name] = list(csv.DictReader(stream))

    # All six distances are eligible, including the unbonded terminal pair.
    endpoints = {row["atoms"]: row for row in reports["epd"]}
    assert len(endpoints) == 6
    assert endpoints["C1-C4"]["connected_first"] == "False"
    bond = endpoints["C2-C3"]
    assert bond["unit"] == "angstrom"
    assert (
        float(bond["first_value"]),
        float(bond["second_value"]),
    ) == pytest.approx((1.45, 1.50))
    assert {row["atoms"] for row in reports["epa"]} == {"C1-C2-C3", "C2-C3-C4"}
    assert [row["atoms"] for row in reports["epdi"]] == ["C1-C2-C3-C4"]

    # Check endpoint ranks against all endpoint values in the same category.
    for ranked, endpoint in (
        ("eprd", "epd"),
        ("epra", "epa"),
        ("eprdi", "epdi"),
    ):
        expected = sorted(
            reports[endpoint], key=lambda row: -float(row["absolute_change"])
        )
        assert reports[ranked] == expected[:5]

    # The middle frame produces a larger change than the endpoint comparison.
    assert len(reports["rd"]) == 5
    for name in ("rd", "ra", "rdi"):
        sizes = [float(row["absolute_change"]) for row in reports[name]]
        assert sizes and sizes == sorted(sizes, reverse=True)
    largest = reports["rd"][0]
    assert (
        largest["atoms"],
        largest["geometry_first"],
        largest["geometry_second"],
    ) == ("C2-C3", "1", "2")
    assert float(largest["signed_change"]) == pytest.approx(0.15)
    assert float(reports["rdi"][0]["absolute_change"]) == pytest.approx(31)

    # Each atom gets every coordinate tied for its largest finite change.
    for maxima, endpoint in (("dl", "epd"), ("al", "epa"), ("dihl", "epdi")):
        for atom in ("C1", "C2", "C3", "C4"):
            candidates = [
                row
                for row in reports[endpoint]
                if atom in row["atoms"].split("-")
            ]
            maximum = max(float(row["absolute_change"]) for row in candidates)
            expected = {
                row["atoms"]
                for row in candidates
                if float(row["absolute_change"]) == maximum
            }
            actual = {
                row["atoms"]
                for row in reports[maxima]
                if row["focus_atom"] == atom
            }
            assert actual == expected


def test_selected_frames_override_rankings_for_each_group(
    write_tmp_file, tmp_path
):
    trajectory = write_tmp_file("path.xyz", TRAJECTORY_XYZ)
    out = tmp_path / "chosen"
    result = CliRunner().invoke(
        cv,
        [
            "-f",
            trajectory,
            "--sg",
            "[3,1]",
            "--ms",
            "C3,C2",
            "--manual-select",
            "C1,C2",
            "--central-bond-only",
            "--rdi",
            "--out",
            str(out),
        ],
    )
    assert result.exit_code == 0, result.output
    assert "[rdi] C3,C2" in result.output
    assert "[rdi] C1,C2" in result.output
    path = out / "group_1_rdi.csv"
    with path.open(newline="", encoding="utf-8") as stream:
        (row,) = csv.DictReader(stream)
    assert (row["atoms"], row["geometry_first"], row["geometry_second"]) == (
        "C1-C2-C3-C4",
        "3",
        "1",
    )
    assert float(row["signed_change"]) == pytest.approx(-2)
    assert float(row["absolute_change"]) == pytest.approx(2)
    # The terminal pair has no proper dihedral with this central bond.
    empty_report = out / "group_2_rdi.csv"
    assert empty_report.read_text(encoding="utf-8").splitlines() == [
        ",".join(row)
    ]


@pytest.mark.parametrize(
    ("options", "message"),
    [
        (["-dl"], "No such option"),
        (["--sg", "1,3"], "two integer geometry numbers"),
        (["--geometries", "[0,3]"], "Geometry numbers must be in 1..3"),
        (["--ms", "H1"], "does not match C1"),
        (["--central-bond-only"], "requires --ms"),
    ],
)
def test_invalid_options(write_tmp_file, options, message):
    trajectory = write_tmp_file("path.xyz", TRAJECTORY_XYZ)
    result = CliRunner().invoke(cv, ["-f", trajectory, *options])
    assert result.exit_code != 0
    assert message in result.output


def test_existing_report_directory_is_preserved(write_tmp_file, tmp_path):
    trajectory = write_tmp_file("path.xyz", TRAJECTORY_XYZ)
    out = tmp_path / "reports"
    out.mkdir()
    previous = out / "group_1_epd.csv"
    previous.write_text("previous results\n", encoding="utf-8")
    result = CliRunner().invoke(cv, ["-f", trajectory, "--out", str(out)])
    assert result.exit_code != 0
    assert "Output directory already exists" in result.output
    assert list(out.iterdir()) == [previous]
    assert previous.read_text(encoding="utf-8") == "previous results\n"