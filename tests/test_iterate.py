"""Focused integration tests for the Iterate v2 generation workflow.

The suite deliberately tests public behavior with small molecules.  Geometry
tests run the complete YAML -> CLI -> runner -> analyzer -> XYZ path and
compare against manually generated golden files.  Combination tests stop
before optimization so traversal rules can be checked exactly and quickly.
"""

from __future__ import annotations

from collections import Counter
from itertools import product
from pathlib import Path

import numpy as np
import pytest
import yaml
from click.testing import CliRunner

from chemsmart.cli.iterate.yaml_cmd import yaml_cmd
from chemsmart.jobs.iterate.job import IterateJob
from chemsmart.jobs.iterate.report import ERROR_CODE_INPUT, ERROR_CODE_TIMEOUT
from chemsmart.jobs.iterate.settings import (
    IterateJobSettings,
    resolve_algorithm_config,
)
from chemsmart.utils.iterate import (
    generate_yaml_template,
    validate_yaml_config,
)

pytestmark = pytest.mark.usefixtures("chemsmart_templates_config")

TESTS_DIR = Path(__file__).resolve().parent
DATA_DIR = TESTS_DIR / "data" / "IterateTests"
CONFIG_DIR = DATA_DIR / "configs"
INPUT_DIR = DATA_DIR / "input"
EXPECTED_DIR = DATA_DIR / "expected_output"
TEMPLATE_GOLDEN = CONFIG_DIR / "iterate_template.yaml"

GENERATION_LABELS = {
    "benzene_1Me",
    "benzene_1OH",
    "benzene_3Me",
    "benzene_3OH",
    "benzene_1Me_3Me",
    "benzene_1Me_3OH",
    "benzene_1OH_3Me",
    "benzene_1OH_3OH",
}


def _read_xyz(path: Path) -> dict[str, tuple[list[str], np.ndarray]]:
    """Parse a merged or single-structure XYZ file by comment label."""
    lines = path.read_text().splitlines()
    structures = {}
    cursor = 0

    while cursor < len(lines):
        if not lines[cursor].strip():
            cursor += 1
            continue

        atom_count = int(lines[cursor].strip())
        cursor += 1
        label = lines[cursor].strip()
        cursor += 1

        symbols = []
        coordinates = []
        for line in lines[cursor : cursor + atom_count]:
            fields = line.split()
            assert (
                len(fields) >= 4
            ), f"Malformed XYZ atom line in {path}: {line}"
            symbols.append(fields[0])
            coordinates.append([float(value) for value in fields[1:4]])
        cursor += atom_count

        assert label, f"Missing XYZ comment label in {path}"
        assert label not in structures, f"Duplicate XYZ label: {label}"
        structures[label] = (symbols, np.asarray(coordinates, dtype=float))

    return structures


def _assert_structure_maps_equal(
    actual: dict[str, tuple[list[str], np.ndarray]],
    expected: dict[str, tuple[list[str], np.ndarray]],
) -> None:
    """Compare labels, atom order, and coordinates for two structure maps."""
    assert set(actual) == set(expected)
    for label in sorted(expected):
        actual_symbols, actual_positions = actual[label]
        expected_symbols, expected_positions = expected[label]
        assert actual_symbols == expected_symbols, label
        np.testing.assert_allclose(
            actual_positions,
            expected_positions,
            rtol=0,
            atol=5e-5,
            err_msg=label,
        )


def _load_validated_config(name: str) -> tuple[Path, dict]:
    config_path = CONFIG_DIR / name
    raw_config = yaml.safe_load(config_path.read_text())
    return config_path, validate_yaml_config(raw_config, str(config_path))


def _build_job(
    config_name: str,
    jobrunner,
    tmp_path: Path,
    *,
    combination_mode: str = "independent",
    timeout: float = 120,
) -> IterateJob:
    """Build a real IterateJob from one of the small v2 YAML fixtures."""
    config_path, config = _load_validated_config(config_name)
    algorithm_config = resolve_algorithm_config(
        yaml_algorithm=config.get("algorithm")
    )
    settings = IterateJobSettings(
        config_file=str(config_path),
        algorithm_config=algorithm_config,
        combination_mode=combination_mode,
    )
    settings.skeleton_list = config["skeletons"]
    settings.substituent_list = config["substituents"]
    return IterateJob(
        settings=settings,
        jobrunner=jobrunner,
        nprocs=1,
        timeout=timeout,
        outputfile=str(tmp_path / config_path.stem),
    )


@pytest.mark.parametrize(
    ("config_name", "expected_name"),
    [
        pytest.param(
            "lagrange_generation.yaml",
            "lagrange_generation.xyz",
            id="joint-lagrange",
        ),
        pytest.param(
            "etkdg_generation.yaml",
            "etkdg_generation.xyz",
            id="etkdg",
        ),
    ],
)
def test_iterate_cli_generation_matches_golden(
    config_name: str,
    expected_name: str,
    tmp_path: Path,
):
    """Both algorithms reproduce their manually generated eight structures."""
    output_base = tmp_path / Path(config_name).stem
    result = CliRunner().invoke(
        yaml_cmd,
        [
            "-f",
            str(CONFIG_DIR / config_name),
            "-cm",
            "global",
            "-np",
            "1",
            "-o",
            str(output_base),
        ],
        obj={},
    )

    assert result.exit_code == 0, result.output
    assert "Total combinations:        8" in result.output
    assert "Successful combinations:   8" in result.output
    assert "Failed combinations:       0" in result.output

    output_path = output_base.with_suffix(".xyz")
    actual = _read_xyz(output_path)
    expected = _read_xyz(EXPECTED_DIR / expected_name)
    assert set(actual) == GENERATION_LABELS
    _assert_structure_maps_equal(actual, expected)

    report_path = tmp_path / f"{Path(config_name).stem}_iterate.out"
    report_text = report_path.read_text()
    assert "CHEMSMART ITERATE JOB REPORT" in report_text
    assert "Total combinations:      8" in report_text
    assert "Generated successfully:      8" in report_text
    assert "Normal termination of ChemSmart Iterate" in report_text


def test_iterate_timeout_is_reported(iterate_jobrunner, tmp_path: Path):
    """A worker exceeding its limit becomes a timed-out combination."""
    job = _build_job(
        "single_generation.yaml",
        iterate_jobrunner,
        tmp_path,
        timeout=1e-9,
    )

    summary = job.run()

    assert summary.total == 1
    assert summary.succeeded == 0
    assert summary.failed == 0
    assert summary.timed_out == 1
    assert summary.structures_written == 0
    assert summary.output_paths == []
    assert summary.error_codes == [ERROR_CODE_TIMEOUT]
    assert summary.exit_code == 1
    assert summary.summary_path is not None
    assert Path(summary.summary_path).exists()


def test_iterate_template_matches_golden(tmp_path: Path):
    """The generated YAML template remains identical to its golden copy."""
    generated_path = Path(generate_yaml_template(str(tmp_path / "template")))
    generated_text = generated_path.read_text()

    assert generated_text.strip() == TEMPLATE_GOLDEN.read_text().strip()
    parsed = yaml.safe_load(generated_text)
    assert parsed["algorithm"]["name"] == "lagrange_multipliers"
    assert len(parsed["skeletons"]) == 1
    assert len(parsed["substituents"]) == 1


def test_iterate_cli_rejects_link_outside_skeleton_indices(tmp_path: Path):
    """A declared link atom must belong to the retained skeleton atoms."""
    config_path = tmp_path / "invalid_link.yaml"
    config_path.write_text("""\
skeletons:
  - file_path: skeleton.xyz
    label: skeleton
    link_index: 7
    skeleton_indices: 1-6
substituents:
  - file_path: substituent.xyz
    label: Me
    link_index: 1
    groups: [1]
""")

    result = CliRunner().invoke(yaml_cmd, ["-f", str(config_path)], obj={})

    assert result.exit_code == 2
    assert "Invalid value" in result.output
    assert (
        "link_index [7] is not included in 'skeleton_indices'" in result.output
    )


@pytest.mark.parametrize(
    ("content", "message"),
    [
        pytest.param(
            """\
skeletons:
  - label: skeleton
    link_index: 1
substituents: []
""",
            "Missing required field 'file_path'",
            id="missing-file-path",
        ),
        pytest.param(
            """\
skeletons:
  - file_path: skeleton.xyz
    label: unsafe/label
    link_index: 1
substituents: []
""",
            "Contains invalid characters",
            id="unsafe-label",
        ),
        pytest.param(
            """\
skeletons:
  - file_path: skeleton.xyz
    label: skeleton
    link_index: 1
substituents:
  - file_path: substituent.xyz
    label: Me
    link_index: 1,2
    groups: [1]
""",
            "Multiple values found in 'link_index'",
            id="multiple-substituent-links",
        ),
        pytest.param(
            """\
skeletons:
  - file_path: skeleton.xyz
    label: duplicate
    link_index: 1
  - file_path: skeleton2.xyz
    label: duplicate
    link_index: 1
substituents:
  - file_path: substituent.xyz
    label: Me
    link_index: 1
    groups: [1,2]
""",
            "duplicate label 'duplicate'",
            id="duplicate-label",
        ),
    ],
)
def test_iterate_cli_rejects_invalid_config(
    content: str,
    message: str,
    tmp_path: Path,
):
    """Representative malformed configurations fail before molecule loading."""
    config_path = tmp_path / "invalid.yaml"
    config_path.write_text(content)

    result = CliRunner().invoke(yaml_cmd, ["-f", str(config_path)], obj={})

    assert result.exit_code == 2
    assert "Invalid value" in result.output
    assert message in result.output


@pytest.mark.parametrize("field", ["link_index", "skeleton_indices"])
def test_iterate_runner_rejects_out_of_bounds_indices(
    field: str,
    iterate_jobrunner,
):
    """Runtime validation rejects indices beyond the loaded molecule size."""
    config = {
        "file_path": str(INPUT_DIR / "benzene.xyz"),
        "file_path_raw": "benzene.xyz",
        "label": "benzene",
        "link_index": [1],
        "skeleton_indices": [1, 2, 3, 4, 5, 6],
    }
    config[field] = [99]
    errors = []

    molecule, label = iterate_jobrunner._load_molecule(
        config, "skeleton", 0, errors
    )

    assert molecule is None
    assert label == "benzene"
    assert len(errors) == 1
    assert errors[0]["error_type"] == "IndexError"
    assert field in errors[0]["error_message"]


def test_iterate_cli_rejects_missing_config_file(tmp_path: Path):
    """A missing YAML file is a CLI usage/configuration error."""
    missing_path = tmp_path / "missing.yaml"

    result = CliRunner().invoke(yaml_cmd, ["-f", str(missing_path)], obj={})

    assert result.exit_code == 2
    assert "does not exist" in result.output


def test_iterate_cli_missing_molecule_file_is_input_error(tmp_path: Path):
    """Absent molecule files fail as an input error (exit 1, ITR-INPUT-001).

    The config is structurally valid, so validation passes, but every declared
    molecule file is missing.  No structure can be produced, so the run is an
    error termination that honours the input-error contract: a non-zero exit
    and the ``ITR-INPUT-001`` code in both the terminal and the run report.
    """
    config_path = tmp_path / "missing_molecule.yaml"
    config_path.write_text("""\
skeletons:
  - file_path: does_not_exist.xyz
    label: skeleton
    link_index: 1
substituents:
  - file_path: also_missing.xyz
    label: Me
    link_index: 1
    groups: [1]
""")
    output_base = tmp_path / "missing_out"

    result = CliRunner().invoke(
        yaml_cmd,
        ["-f", str(config_path), "-np", "1", "-o", str(output_base)],
        obj={},
    )

    assert result.exit_code == 1, result.output
    assert ERROR_CODE_INPUT in result.output
    # Nothing was generated, so the merged output file is never created.
    assert not output_base.with_suffix(".xyz").exists()
    # The best-effort report records the input error and error termination.
    report_path = tmp_path / "missing_molecule_iterate.out"
    assert report_path.exists()
    report_text = report_path.read_text()
    assert ERROR_CODE_INPUT in report_text
    assert "Error termination of ChemSmart Iterate" in report_text


def test_iterate_cli_partial_failure_retains_successes(tmp_path: Path):
    """A partial input failure exits 1 yet keeps the successful structure.

    One substituent loads and is placed while a second substituent's file is
    missing.  The run is an error termination (exit 1, ITR-INPUT-001), but the
    structure that did generate is still written to the merged output file.
    """
    benzene = INPUT_DIR / "benzene.xyz"
    methane = INPUT_DIR / "methane.xyz"
    missing = tmp_path / "missing.xyz"
    config_path = tmp_path / "partial_failure.yaml"
    config_path.write_text(
        "skeletons:\n"
        f'  - file_path: "{benzene}"\n'
        "    label: benzene\n"
        '    skeleton_indices: "1-6"\n'
        "    slots:\n"
        "      - group: 1\n"
        "        link_indices: 1\n"
        "substituents:\n"
        f'  - file_path: "{methane}"\n'
        "    label: Me\n"
        "    link_index: 1\n"
        "    groups: [1]\n"
        f'  - file_path: "{missing}"\n'
        "    label: Bad\n"
        "    link_index: 1\n"
        "    groups: [1]\n"
    )
    output_base = tmp_path / "partial_out"

    result = CliRunner().invoke(
        yaml_cmd,
        ["-f", str(config_path), "-np", "1", "-o", str(output_base)],
        obj={},
    )

    assert result.exit_code == 1, result.output
    assert ERROR_CODE_INPUT in result.output
    assert "Successful combinations:" in result.output
    # The single loadable substituent still produced and retained its output.
    output_path = output_base.with_suffix(".xyz")
    assert output_path.exists()
    assert set(_read_xyz(output_path)) == {"benzene_1Me"}


def test_iterate_three_site_combination_traversal(
    iterate_jobrunner,
    tmp_path: Path,
):
    """Two substituents traverse three sites as 6/12/8 substitutions."""
    job = _build_job(
        "combination_traversal.yaml",
        iterate_jobrunner,
        tmp_path,
        combination_mode="global",
    )

    _, combinations, input_errors, attachment_sites = (
        iterate_jobrunner._generate_combinations(job)
    )
    labels = [combination.label for combination in combinations]

    expected_labels = set()
    for choices in product((None, "Me", "OH"), repeat=3):
        assignments = [
            f"{site}{choice}"
            for site, choice in zip((1, 3, 5), choices)
            if choice is not None
        ]
        if assignments:
            expected_labels.add("_".join(["benzene", *assignments]))

    assert input_errors == []
    assert attachment_sites == 3
    assert len(combinations) == 26
    assert len(labels) == len(set(labels))
    assert set(labels) == expected_labels
    assert Counter(len(combo.assignments) for combo in combinations) == {
        1: 6,
        2: 12,
        3: 8,
    }


def test_iterate_independent_and_global_combination_modes(
    iterate_jobrunner,
    tmp_path: Path,
):
    """Global mode adds the cross-group structure to independent results."""
    labels_by_mode = {}
    for mode in ("independent", "global"):
        job = _build_job(
            "combination_modes.yaml",
            iterate_jobrunner,
            tmp_path,
            combination_mode=mode,
        )
        _, combinations, input_errors, _ = (
            iterate_jobrunner._generate_combinations(job)
        )
        assert input_errors == []
        labels_by_mode[mode] = {combo.label for combo in combinations}

    assert labels_by_mode["independent"] == {
        "benzene_1Me",
        "benzene_3OH",
    }
    assert labels_by_mode["global"] == {
        "benzene_1Me",
        "benzene_3OH",
        "benzene_1Me_3OH",
    }


def test_iterate_cli_separate_outputs(tmp_path: Path):
    """Separate-output mode writes one correctly named XYZ per combination."""
    output_directory = tmp_path / "separate"
    result = CliRunner().invoke(
        yaml_cmd,
        [
            "-f",
            str(CONFIG_DIR / "etkdg_generation.yaml"),
            "-cm",
            "global",
            "-np",
            "1",
            "--separate-outputs",
            "-d",
            str(output_directory),
        ],
        obj={},
    )

    assert result.exit_code == 0, result.output
    output_files = sorted(output_directory.glob("*.xyz"))
    assert len(output_files) == 8
    assert {path.stem for path in output_files} == GENERATION_LABELS

    actual = {}
    for output_file in output_files:
        structures = _read_xyz(output_file)
        assert set(structures) == {output_file.stem}
        actual.update(structures)

    expected = _read_xyz(EXPECTED_DIR / "etkdg_generation.xyz")
    _assert_structure_maps_equal(actual, expected)
    assert (output_directory / "etkdg_generation_iterate.out").exists()
