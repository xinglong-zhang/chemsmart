"""Focused integration tests for the Iterate v2 generation workflow.

The suite deliberately tests public behavior with small molecules.  Geometry
tests run the complete YAML -> CLI -> runner -> analyzer -> XYZ path and
compare against manually generated golden files.  Combination tests stop
before optimization so traversal rules can be checked exactly and quickly.
"""

from __future__ import annotations

import ast
import importlib
import shlex
import sys
from collections import Counter
from datetime import datetime
from itertools import product
from pathlib import Path

import click
import numpy as np
import pytest
import yaml
from click.testing import CliRunner

from chemsmart.cli.main import entry_point
from chemsmart.cli.run import run
from chemsmart.jobs.iterate.job import IterateExecutionError, IterateJob
from chemsmart.jobs.iterate.report import (
    ERROR_CODE_INPUT,
    ERROR_CODE_INTERNAL,
    ERROR_CODE_TIMEOUT,
    IterateReport,
)
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
run_cli_module = importlib.import_module("chemsmart.cli.run")

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


def _build_job_from_config_path(
    config_path: Path,
    jobrunner,
    tmp_path: Path,
    *,
    combination_mode: str = "independent",
    timeout: float = 120,
    max_substituted_sites: int | None = None,
) -> IterateJob:
    raw_config = yaml.safe_load(config_path.read_text())
    config = validate_yaml_config(raw_config, str(config_path))
    algorithm_config = resolve_algorithm_config(
        yaml_algorithm=config.get("algorithm")
    )
    settings = IterateJobSettings(
        config_file=str(config_path),
        algorithm_config=algorithm_config,
        combination_mode=combination_mode,
        max_substituted_sites=max_substituted_sites,
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


def _build_job(
    config_name: str,
    jobrunner,
    tmp_path: Path,
    *,
    combination_mode: str = "independent",
    timeout: float = 120,
    max_substituted_sites: int | None = None,
) -> IterateJob:
    """Build a real IterateJob from one of the small v2 YAML fixtures."""
    return _build_job_from_config_path(
        CONFIG_DIR / config_name,
        jobrunner,
        tmp_path,
        combination_mode=combination_mode,
        timeout=timeout,
        max_substituted_sites=max_substituted_sites,
    )


def _write_iterate_config(tmp_path: Path, name: str, content: str) -> Path:
    config_path = tmp_path / name
    config_path.write_text(content)
    return config_path


def _assignment_signature(combo) -> tuple[tuple[int, str], ...]:
    return tuple(
        sorted(
            (
                assignment.skeleton_link_index,
                assignment.substituent_label,
            )
            for assignment in combo.assignments
        )
    )


def _capture_iterate_cli_jobs(monkeypatch) -> list[IterateJob]:
    captured_jobs = []

    def fake_run_single_job(job, jobrunner):
        captured_jobs.append(job)

    monkeypatch.setattr(
        run_cli_module,
        "_run_single_job",
        fake_run_single_job,
    )
    return captured_jobs


def _minimal_iterate_report(
    max_substituted_sites: int | None,
    *,
    error: bool = False,
) -> IterateReport:
    return IterateReport(
        run_id="test",
        chemsmart_version="test",
        rdkit_version="test",
        started_at=datetime(2026, 1, 1, 0, 0, 0),
        finished_at=datetime(2026, 1, 1, 0, 0, 1),
        duration_seconds=1.0,
        working_directory=".",
        command_line="chemsmart run iterate yaml -f config.yaml",
        config_file="config.yaml",
        algorithm_name="etkdg",
        algorithm_options={},
        combination_mode="independent",
        max_substituted_sites=max_substituted_sites,
        extra_error_codes=[ERROR_CODE_INTERNAL] if error else [],
        exit_code_override=1 if error else None,
    )


def test_iterate_job_names_follow_config_stem(tmp_path: Path):
    """Job, report, and default merged output share one config-based stem."""
    settings = IterateJobSettings(config_file=str(tmp_path / "my config.yaml"))
    job = IterateJob(settings=settings)
    job.folder = str(tmp_path)

    assert job.label == "my_config_iterate"
    assert job.outputfile == "my_config_iterate.xyz"
    assert job.reportfile == str(tmp_path / "my_config_iterate.out")

    custom_output = tmp_path / "custom"
    custom_job = IterateJob(settings=settings, outputfile=str(custom_output))
    assert custom_job.label == "my_config_iterate"
    assert custom_job.outputfile == str(custom_output.with_suffix(".xyz"))
    assert custom_job.reportfile == str(tmp_path / "my_config_iterate.out")


def test_iterate_job_skips_complete_output_unless_rerun_requested(
    tmp_path: Path,
):
    """A normal report and its declared output activate skip-completed."""
    settings = IterateJobSettings(config_file=str(tmp_path / "complete.yaml"))
    job = IterateJob(settings=settings)
    job.folder = str(tmp_path)
    (tmp_path / job.outputfile).write_text("1\ncomplete\nH 0 0 0\n")
    Path(job.reportfile).write_text(
        "Structures written: 1\n"
        "Normal termination of CHEMSMART Iterate at 2026-01-01 00:00:00.\n"
    )

    class RecordingRunner:
        def __init__(self):
            self.calls = 0

        def run(self, iterate_job, progress_callback=None):
            self.calls += 1

            class Summary:
                exit_code = 0

            return Summary()

    runner = RecordingRunner()
    job.jobrunner = runner

    assert job.is_complete()
    assert job.skip_completed
    assert job.run() is None
    assert runner.calls == 0

    job.skip_completed = False
    summary = job.run()
    assert summary.exit_code == 0
    assert runner.calls == 1


@pytest.mark.parametrize(
    ("config_name", "expected_name", "algorithm_args"),
    [
        pytest.param(
            "jlgo_generation.yaml",
            "jlgo_generation.xyz",
            ["jlgo"],
            id="jlgo",
        ),
        pytest.param(
            "etkdg_generation.yaml",
            "etkdg_generation.xyz",
            ["etkdg"],
            id="etkdg",
        ),
    ],
)
def test_iterate_cli_generation_matches_golden(
    config_name: str,
    expected_name: str,
    algorithm_args: list[str],
    tmp_path: Path,
):
    """Both algorithms reproduce their manually generated eight structures."""
    output_base = tmp_path / Path(config_name).stem
    result = CliRunner().invoke(
        run,
        [
            "--num-cores",
            "1",
            "iterate",
            "yaml",
            "-f",
            str(CONFIG_DIR / config_name),
            "-cm",
            "global",
            "-o",
            str(output_base),
            *algorithm_args,
        ],
        obj={},
    )

    assert result.exit_code == 0, result.output

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
    assert "Number of processes:     1" in report_text
    assert "Normal termination of CHEMSMART Iterate" in report_text
    assert "Configuration SHA256" not in report_text


def test_iterate_timeout_is_reported(iterate_jobrunner, tmp_path: Path):
    """A worker exceeding its limit becomes a timed-out combination."""
    job = _build_job(
        "single_generation.yaml",
        iterate_jobrunner,
        tmp_path,
        timeout=1e-9,
    )

    with pytest.raises(IterateExecutionError) as exc_info:
        job.run()
    summary = exc_info.value.summary

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
    assert parsed["algorithm"]["name"] == "etkdg"
    assert len(parsed["skeletons"]) == 1
    assert len(parsed["substituents"]) == 1


def test_iterate_template_generation_is_rejected_under_sub(tmp_path: Path):
    """The real CLI hierarchy must not run template utilities under ``sub``."""
    output_path = tmp_path / "should_not_be_written.yaml"

    result = CliRunner().invoke(
        entry_point,
        [
            "sub",
            "iterate",
            "yaml",
            "--generate-template",
            str(output_path),
        ],
        obj={},
    )

    assert result.exit_code == 2
    assert "local utility" in result.output
    assert not output_path.exists()


@pytest.mark.parametrize(
    ("server_name", "resource_directive"),
    [
        (
            "PBS",
            "#PBS -l select=1:ncpus=4:mpiprocs=1:mem=375G",
        ),
        (
            "SLURM",
            "#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=375G",
        ),
    ],
)
def test_iterate_sub_test_writes_rerunnable_python_native_scripts(
    tmp_path: Path,
    monkeypatch,
    server_name: str,
    resource_directive: str,
):
    """``sub --test`` reconstructs Iterate and writes native CPU scripts."""
    monkeypatch.chdir(tmp_path)
    config_path = CONFIG_DIR / "etkdg_generation.yaml"

    result = CliRunner().invoke(
        entry_point,
        [
            "sub",
            "--server",
            server_name,
            "--num-cores",
            "4",
            "--test",
            "iterate",
            "yaml",
            "-f",
            str(config_path),
            "-ms",
            "2",
            "etkdg",
            "--num-conformers",
            "3",
        ],
        obj={},
    )

    assert result.exit_code == 0, result.output

    label = "etkdg_generation_iterate"
    run_script = tmp_path / f"chemsmart_run_{label}.py"
    submit_script = tmp_path / f"chemsmart_sub_{label}.sh"
    assert run_script.exists()
    assert submit_script.exists()

    run_tree = ast.parse(run_script.read_text())
    cli_assignment = next(
        node
        for node in ast.walk(run_tree)
        if isinstance(node, ast.Assign)
        and any(
            isinstance(target, ast.Name) and target.id == "cli_args"
            for target in node.targets
        )
    )
    reconstructed_args = ast.literal_eval(cli_assignment.value)
    assert reconstructed_args[reconstructed_args.index("--server") + 1] == (
        server_name
    )
    assert (
        reconstructed_args[reconstructed_args.index("--num-cores") + 1] == "4"
    )
    assert "--skip-completed" in reconstructed_args
    max_sites_index = reconstructed_args.index("--max-substituted-sites")
    assert reconstructed_args[max_sites_index + 1] == "2"
    assert (
        reconstructed_args[reconstructed_args.index("--num-conformers") + 1]
        == "3"
    )
    assert "sub" not in reconstructed_args
    assert (
        reconstructed_args.index("iterate")
        < reconstructed_args.index("yaml")
        < max_sites_index
        < reconstructed_args.index("etkdg")
    )
    assert (
        "sys.argv = ['chemsmart', 'run', *cli_args]" in run_script.read_text()
    )

    submit_contents = submit_script.read_text()
    assert resource_directive in submit_contents
    assert "conda activate" not in submit_contents
    assert "export OMP_NUM_THREADS=1" in submit_contents
    assert "export MKL_NUM_THREADS=1" in submit_contents
    assert (
        f"{shlex.quote(sys.executable)} ./{run_script.name} &"
        in submit_contents
    )


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

    result = CliRunner().invoke(
        run, ["iterate", "yaml", "-f", str(config_path)], obj={}
    )

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

    result = CliRunner().invoke(
        run, ["iterate", "yaml", "-f", str(config_path)], obj={}
    )

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

    result = CliRunner().invoke(
        run, ["iterate", "yaml", "-f", str(missing_path)], obj={}
    )

    assert result.exit_code == 2
    assert "does not exist" in result.output


@pytest.mark.parametrize(
    ("option_args", "expected"),
    [
        pytest.param(["-ms", "2"], 2, id="short"),
        pytest.param(
            ["--max-substituted-sites", "2"],
            2,
            id="long",
        ),
        pytest.param([], None, id="unset"),
        pytest.param(["-ms", "0"], None, id="zero"),
        pytest.param(["-ms", "2", "jlgo"], 2, id="before-algorithm"),
    ],
)
def test_iterate_cli_max_substituted_sites_option(
    option_args: list[str],
    expected: int | None,
    monkeypatch,
):
    """The YAML input command forwards the reusable Iterate limit option."""
    captured_jobs = _capture_iterate_cli_jobs(monkeypatch)

    result = CliRunner().invoke(
        run,
        [
            "iterate",
            "yaml",
            "-f",
            str(CONFIG_DIR / "combination_traversal.yaml"),
            *option_args,
        ],
        obj={},
    )

    assert result.exit_code == 0, result.output
    assert len(captured_jobs) == 1
    assert captured_jobs[0].settings.max_substituted_sites == expected


@pytest.mark.parametrize(
    "option_args",
    [
        pytest.param(["-ms", "-1"], id="negative"),
        pytest.param(["-ms", "1.5"], id="float"),
        pytest.param(["-ms", "abc"], id="nonnumeric"),
        pytest.param(["-ms"], id="missing-short-value"),
    ],
)
def test_iterate_cli_rejects_invalid_max_substituted_sites(
    option_args: list[str],
    monkeypatch,
):
    """Invalid limits fail as Click usage errors before job execution."""
    captured_jobs = _capture_iterate_cli_jobs(monkeypatch)

    result = CliRunner().invoke(
        run,
        [
            "iterate",
            "yaml",
            "-f",
            str(CONFIG_DIR / "combination_traversal.yaml"),
            *option_args,
        ],
        obj={},
    )

    assert result.exit_code == 2
    assert "Error:" in result.output
    assert captured_jobs == []


def test_iterate_cli_missing_molecule_file_is_input_error(tmp_path: Path):
    """Absent molecule files are recorded as an input error.

    The config is structurally valid, so validation passes, but every declared
    molecule file is missing.  No structure can be produced, so the run is an
    error termination recorded with ``ITR-INPUT-001`` in the run report.
    Iterate propagates its internal completion status to the CLI process so a
    scheduler can distinguish this failure from a successful calculation.
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
        run,
        [
            "--num-cores",
            "1",
            "iterate",
            "yaml",
            "-f",
            str(config_path),
            "-o",
            str(output_base),
        ],
        obj={},
    )

    assert result.exit_code == 1
    assert isinstance(result.exception, IterateExecutionError)
    # Nothing was generated, so the merged output file is never created.
    assert not output_base.with_suffix(".xyz").exists()
    # The best-effort report records the input error and error termination.
    report_path = tmp_path / "missing_molecule_iterate.out"
    assert report_path.exists()
    report_text = report_path.read_text()
    assert ERROR_CODE_INPUT in report_text
    assert "Error termination of CHEMSMART Iterate" in report_text


def test_iterate_cli_partial_failure_retains_successes(tmp_path: Path):
    """A partial input failure keeps the successful structure and report.

    One substituent loads and is placed while a second substituent's file is
    missing. The report records an error termination with ``ITR-INPUT-001``,
    while the structure that did generate is retained in the merged output.
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
        run,
        [
            "--num-cores",
            "1",
            "iterate",
            "yaml",
            "-f",
            str(config_path),
            "-o",
            str(output_base),
        ],
        obj={},
    )

    assert result.exit_code == 1
    assert isinstance(result.exception, IterateExecutionError)
    # The single loadable substituent still produced and retained its output.
    output_path = output_base.with_suffix(".xyz")
    assert output_path.exists()
    assert set(_read_xyz(output_path)) == {"benzene_1Me"}
    report_text = (tmp_path / "partial_failure_iterate.out").read_text()
    assert ERROR_CODE_INPUT in report_text
    assert "Error termination of CHEMSMART Iterate" in report_text


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


@pytest.mark.parametrize(
    ("raw_value", "expected"),
    [
        pytest.param(None, None, id="none"),
        pytest.param(0, None, id="zero"),
        pytest.param(2, 2, id="positive"),
    ],
)
def test_iterate_settings_normalizes_max_substituted_sites(
    raw_value,
    expected,
):
    settings = IterateJobSettings(max_substituted_sites=raw_value)

    assert settings.max_substituted_sites == expected


@pytest.mark.parametrize(
    "raw_value",
    [
        pytest.param(-1, id="negative"),
        pytest.param(1.0, id="float"),
        pytest.param("1", id="string"),
        pytest.param(True, id="true-bool"),
        pytest.param(False, id="false-bool"),
    ],
)
def test_iterate_settings_rejects_invalid_max_substituted_sites(raw_value):
    with pytest.raises(ValueError, match="max_substituted_sites"):
        IterateJobSettings(max_substituted_sites=raw_value)


def test_iterate_settings_copy_preserves_max_substituted_sites():
    limited = IterateJobSettings(max_substituted_sites=2)
    unlimited = IterateJobSettings(max_substituted_sites=0)

    assert limited.copy().max_substituted_sites == 2
    assert unlimited.copy().max_substituted_sites is None


def test_iterate_independent_max_substituted_sites_limits_current_group(
    iterate_jobrunner,
    tmp_path: Path,
):
    """A three-position slot is pruned to singles and doubles at K=2."""
    unlimited_job = _build_job(
        "combination_traversal.yaml",
        iterate_jobrunner,
        tmp_path,
        combination_mode="independent",
    )
    limited_job = _build_job(
        "combination_traversal.yaml",
        iterate_jobrunner,
        tmp_path,
        combination_mode="independent",
        max_substituted_sites=2,
    )

    _, unlimited, unlimited_errors, _ = (
        iterate_jobrunner._generate_combinations(unlimited_job)
    )
    _, limited, limited_errors, _ = iterate_jobrunner._generate_combinations(
        limited_job
    )

    assert unlimited_errors == []
    assert limited_errors == []
    assert len(unlimited) == 26
    assert len(limited) == 18
    assert Counter(len(combo.assignments) for combo in limited) == {
        1: 6,
        2: 12,
    }
    assert {_assignment_signature(combo) for combo in limited} == {
        _assignment_signature(combo)
        for combo in unlimited
        if len(combo.assignments) <= 2
    }


def test_iterate_global_max_substituted_sites_spans_groups(
    iterate_jobrunner,
    tmp_path: Path,
):
    benzene = INPUT_DIR / "benzene.xyz"
    methane = INPUT_DIR / "methane.xyz"
    config_path = _write_iterate_config(
        tmp_path,
        "three_group_global.yaml",
        f"""\
skeletons:
  - file_path: "{benzene}"
    label: tri
    skeleton_indices: "1-6"
    slots:
      - group: 1
        link_indices: 1
      - group: 2
        link_indices: 3
      - group: 3
        link_indices: 5
substituents:
  - file_path: "{methane}"
    label: Me1
    link_index: 1
    groups: [1]
  - file_path: "{methane}"
    label: Me2
    link_index: 1
    groups: [2]
  - file_path: "{methane}"
    label: Me3
    link_index: 1
    groups: [3]
""",
    )
    job = _build_job_from_config_path(
        config_path,
        iterate_jobrunner,
        tmp_path,
        combination_mode="global",
        max_substituted_sites=2,
    )

    _, combinations, input_errors, _ = (
        iterate_jobrunner._generate_combinations(job)
    )

    assert input_errors == []
    assert len(combinations) == 6
    assert all(len(combo.assignments) <= 2 for combo in combinations)
    assert any(
        {assignment.skeleton_link_index for assignment in combo.assignments}
        == {1, 3}
        for combo in combinations
    )


def test_iterate_max_substituted_sites_is_per_skeleton(
    iterate_jobrunner,
    tmp_path: Path,
):
    benzene = INPUT_DIR / "benzene.xyz"
    methane = INPUT_DIR / "methane.xyz"
    water = INPUT_DIR / "water.xyz"
    config_path = _write_iterate_config(
        tmp_path,
        "multi_skeleton_limit.yaml",
        f"""\
skeletons:
  - file_path: "{benzene}"
    label: two_sites
    skeleton_indices: "1-6"
    slots:
      - group: 1
        link_indices: "1,3"
  - file_path: "{benzene}"
    label: three_sites
    skeleton_indices: "1-6"
    slots:
      - group: 2
        link_indices: "1,3,5"
substituents:
  - file_path: "{methane}"
    label: Me
    link_index: 1
    groups: [1]
  - file_path: "{water}"
    label: OH
    link_index: 1
    groups: [2]
""",
    )
    unlimited_job = _build_job_from_config_path(
        config_path,
        iterate_jobrunner,
        tmp_path,
        combination_mode="global",
    )
    limited_job = _build_job_from_config_path(
        config_path,
        iterate_jobrunner,
        tmp_path,
        combination_mode="global",
        max_substituted_sites=2,
    )

    _, unlimited, unlimited_errors, _ = (
        iterate_jobrunner._generate_combinations(unlimited_job)
    )
    _, limited, limited_errors, _ = iterate_jobrunner._generate_combinations(
        limited_job
    )

    assert unlimited_errors == []
    assert limited_errors == []
    assert Counter(combo.skeleton_label for combo in unlimited) == {
        "two_sites": 3,
        "three_sites": 7,
    }
    assert Counter(combo.skeleton_label for combo in limited) == {
        "two_sites": 3,
        "three_sites": 6,
    }
    assert all(len(combo.assignments) <= 2 for combo in limited)


def test_iterate_multi_link_index_forms_one_implicit_group(
    iterate_jobrunner,
    tmp_path: Path,
):
    benzene = INPUT_DIR / "benzene.xyz"
    methane = INPUT_DIR / "methane.xyz"
    water = INPUT_DIR / "water.xyz"
    config_path = _write_iterate_config(
        tmp_path,
        "link_index_implicit_group.yaml",
        f"""\
skeletons:
  - file_path: "{benzene}"
    label: shorthand
    skeleton_indices: "1-6"
    link_index: "1,3,5"
substituents:
  - file_path: "{methane}"
    label: Me
    link_index: 1
    groups: [1]
  - file_path: "{water}"
    label: OH
    link_index: 1
    groups: [1]
""",
    )
    combinations_by_case = {}
    for case, mode, limit in (
        ("unlimited", "independent", None),
        ("independent-limited", "independent", 2),
        ("global-limited", "global", 2),
    ):
        job = _build_job_from_config_path(
            config_path,
            iterate_jobrunner,
            tmp_path,
            combination_mode=mode,
            max_substituted_sites=limit,
        )
        _, combinations, input_errors, _ = (
            iterate_jobrunner._generate_combinations(job)
        )
        assert input_errors == []
        combinations_by_case[case] = combinations

    unlimited = combinations_by_case["unlimited"]
    assert len(unlimited) == 26
    assert Counter(len(combo.assignments) for combo in unlimited) == {
        1: 6,
        2: 12,
        3: 8,
    }
    labels = {combo.label for combo in unlimited}
    assert {
        "shorthand_1Me",
        "shorthand_3OH",
        "shorthand_1Me_3OH",
        "shorthand_1OH_5Me",
        "shorthand_1Me_3OH_5Me",
    }.issubset(labels)

    limited_by_mode = {
        case: combinations_by_case[case]
        for case in ("independent-limited", "global-limited")
    }
    for limited in limited_by_mode.values():
        assert len(limited) == 18
        assert Counter(len(combo.assignments) for combo in limited) == {
            1: 6,
            2: 12,
        }
    assert {
        _assignment_signature(combo)
        for combo in limited_by_mode["independent-limited"]
    } == {
        _assignment_signature(combo)
        for combo in limited_by_mode["global-limited"]
    }


def test_iterate_multi_link_index_uses_one_implicit_group_number(
    iterate_jobrunner,
    tmp_path: Path,
):
    benzene = INPUT_DIR / "benzene.xyz"
    methane = INPUT_DIR / "methane.xyz"
    water = INPUT_DIR / "water.xyz"
    config_path = _write_iterate_config(
        tmp_path,
        "link_index_group_numbering.yaml",
        f"""\
skeletons:
  - file_path: "{benzene}"
    label: first
    skeleton_indices: "1-6"
    link_index: "1,3,5"
  - file_path: "{benzene}"
    label: second
    skeleton_indices: "1-6"
    link_index: 1
substituents:
  - file_path: "{methane}"
    label: Me
    link_index: 1
    groups: [1]
  - file_path: "{water}"
    label: OH
    link_index: 1
    groups: [2]
""",
    )
    job = _build_job_from_config_path(
        config_path,
        iterate_jobrunner,
        tmp_path,
        combination_mode="global",
    )

    _, combinations, input_errors, _ = (
        iterate_jobrunner._generate_combinations(job)
    )

    assert input_errors == []
    assert {combo.label for combo in combinations} == {
        "first_1Me",
        "first_3Me",
        "first_5Me",
        "first_1Me_3Me",
        "first_1Me_5Me",
        "first_3Me_5Me",
        "first_1Me_3Me_5Me",
        "second_1OH",
    }


@pytest.mark.parametrize(
    ("content", "message"),
    [
        pytest.param(
            """\
max_substituted_sites: 2
skeletons: []
substituents: []
""",
            "Unknown top-level key",
            id="top-level",
        ),
        pytest.param(
            """\
algorithm:
  name: etkdg
  max_substituted_sites: 2
skeletons: []
substituents: []
""",
            "Unknown key(s) in 'algorithm' block",
            id="algorithm",
        ),
        pytest.param(
            """\
skeletons:
  - file_path: skeleton.xyz
    label: skel
    link_index: 1
    max_substituted_sites: 2
substituents: []
""",
            "Unknown key(s) in skeleton entry",
            id="skeleton",
        ),
        pytest.param(
            """\
skeletons:
  - file_path: skeleton.xyz
    label: skel
    link_index: 1
substituents:
  - file_path: sub.xyz
    label: sub
    link_index: 1
    groups: [1]
    max_substituted_sites: 2
""",
            "Unknown key(s) in substituent entry",
            id="substituent",
        ),
    ],
)
def test_iterate_yaml_rejects_max_substituted_sites_key(
    content: str,
    message: str,
    tmp_path: Path,
):
    config_path = tmp_path / "bad.yaml"

    with pytest.raises(click.BadParameter) as exc_info:
        validate_yaml_config(yaml.safe_load(content), str(config_path))

    assert message in str(exc_info.value)
    assert "max_substituted_sites" in str(exc_info.value)


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


@pytest.mark.parametrize(
    ("raw_value", "expected_text"),
    [
        pytest.param(None, "Maximum substituted sites: unlimited", id="none"),
        pytest.param(2, "Maximum substituted sites: 2", id="positive"),
    ],
)
def test_iterate_report_renders_max_substituted_sites(
    raw_value,
    expected_text: str,
):
    text = _minimal_iterate_report(raw_value).render()

    assert expected_text in text
    assert "Normal termination of CHEMSMART Iterate" in text


def test_iterate_error_report_renders_max_substituted_sites():
    text = _minimal_iterate_report(2, error=True).render()

    assert "Maximum substituted sites: 2" in text
    assert ERROR_CODE_INTERNAL in text
    assert "Error termination of CHEMSMART Iterate" in text


@pytest.mark.parametrize(
    ("report_value", "current_value", "expected"),
    [
        pytest.param("unlimited", None, True, id="unlimited-unset"),
        pytest.param("unlimited", 2, False, id="unlimited-to-two"),
        pytest.param("2", None, False, id="two-to-unlimited"),
        pytest.param("2", 3, False, id="two-to-three"),
        pytest.param("2", 2, True, id="two-to-two"),
        pytest.param(None, None, True, id="legacy-unset"),
        pytest.param(None, 2, False, id="legacy-to-two"),
    ],
)
def test_iterate_skip_completed_compares_max_substituted_sites(
    report_value: str | None,
    current_value: int | None,
    expected: bool,
    tmp_path: Path,
):
    settings = IterateJobSettings(
        config_file=str(tmp_path / "complete.yaml"),
        max_substituted_sites=current_value,
    )
    job = IterateJob(settings=settings)
    job.folder = str(tmp_path)

    report_lines = ["CHEMSMART ITERATE JOB REPORT\n"]
    if report_value is not None:
        report_lines.append(f" Maximum substituted sites: {report_value}\n")
    report_lines += [
        " Structures written: 0\n",
        " Normal termination of CHEMSMART Iterate at 2026-01-01 00:00:00.\n",
    ]
    Path(job.reportfile).write_text("".join(report_lines))

    assert job.is_complete() is expected


def test_iterate_cli_separate_outputs(tmp_path: Path):
    """Separate-output mode writes one correctly named XYZ per combination."""
    output_directory = tmp_path / "separate"
    result = CliRunner().invoke(
        run,
        [
            "--num-cores",
            "1",
            "iterate",
            "yaml",
            "-f",
            str(CONFIG_DIR / "etkdg_generation.yaml"),
            "-cm",
            "global",
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
