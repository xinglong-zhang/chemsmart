import os

import numpy as np
import pytest
import yaml

from chemsmart.jobs.iterate.job import IterateJob
from chemsmart.jobs.iterate.settings import IterateJobSettings
from chemsmart.utils.iterate import validate_yaml_config

pytestmark = pytest.mark.usefixtures("chemsmart_templates_config")


def test_iterate_integration_workflow(
    tmpdir,
    iterate_integration_config_file,
    iterate_input_directory,
    iterate_expected_output_file,
    iterate_jobrunner,
):
    """
    Test the full Iterate workflow (Integration Test):
    1. Load config from YAML (integration_iterate.yaml)
    2. Run IterateJob
    3. Compare output with expected XYZ file
    """
    # Change CWD to input directory so relative paths in YAML work
    original_cwd = os.getcwd()
    os.chdir(iterate_input_directory)

    try:
        # 1. Load and validate configuration
        with open(iterate_integration_config_file, "r") as f:
            raw_config = yaml.safe_load(f)

        config = validate_yaml_config(
            raw_config, iterate_integration_config_file
        )

        # 2. Setup Job Settings
        job_settings = IterateJobSettings(
            config_file=iterate_integration_config_file,
        )
        job_settings.skeleton_list = config["skeletons"]
        job_settings.substituent_list = config["substituents"]

        # 3. Setup Job
        # Use a temporary file for output
        output_file = tmpdir / "test_output"

        jobrunner = iterate_jobrunner
        job = IterateJob(
            settings=job_settings,
            jobrunner=jobrunner,
            nprocs=4,
            outputfile=str(output_file),
        )

        # 4. Run Job
        summary = job.run()

        # 5. Verify Output
        assert summary.succeeded > 0, "No structures were generated"
        assert (
            len(summary.output_paths) == 1
        ), "Merged mode should write exactly one file"
        generated_output_path = summary.output_paths[0]
        assert os.path.exists(
            generated_output_path
        ), "Output file was not generated"

        # The run report is written next to the merged output.
        assert summary.summary_path is not None
        assert os.path.exists(summary.summary_path)
        with open(summary.summary_path) as report_handle:
            report_text = report_handle.read()
        assert "CHEMSMART ITERATE JOB REPORT" in report_text
        assert "FINAL SUMMARY" in report_text
        assert "Normal termination of ChemSmart Iterate" in report_text

        # Compare generated output with expected output
        # Semantic comparison (atoms and coordinates)
        # is preferred over byte-comparison
        # to robustly handle floating point formatting differences in XYZ files
        from chemsmart.io.xyz.xyzfile import XYZFile

        generated_xyz = XYZFile(generated_output_path)
        generated_structures = generated_xyz.get_molecules(
            index=":", return_list=True
        )
        generated_structures_comments = generated_xyz.get_comments(
            index=":", return_list=True
        )
        # attach comment to structure as molecule.info for easier comparison
        for mol, comment in zip(
            generated_structures, generated_structures_comments
        ):
            mol.info["comment"] = comment

        expected_xyz = XYZFile(iterate_expected_output_file)
        expected_structures = expected_xyz.get_molecules(
            index=":", return_list=True
        )
        expected_structures_comments = expected_xyz.get_comments(
            index=":", return_list=True
        )
        for mol, comment in zip(
            expected_structures, expected_structures_comments
        ):
            mol.info["comment"] = comment

        # avoid silent truncation by zip()
        assert len(generated_structures) == len(
            generated_structures_comments
        ), (
            f"Generated molecules/comments length mismatch: "
            f"{len(generated_structures)} != {len(generated_structures_comments)}"
        )
        assert len(expected_structures) == len(expected_structures_comments), (
            f"Expected molecules/comments length mismatch: "
            f"{len(expected_structures)} != {len(expected_structures_comments)}"
        )

        # attach comment to structure as molecule.info for easier comparison
        for mol, comment in zip(
            generated_structures, generated_structures_comments
        ):
            mol.info["comment"] = (comment or "").strip()

        for mol, comment in zip(
            expected_structures, expected_structures_comments
        ):
            mol.info["comment"] = (comment or "").strip()

        # Sort structures by comment to ensure order-independent comparison
        generated_structures.sort(key=lambda m: (m.info.get("comment") or ""))
        expected_structures.sort(key=lambda m: (m.info.get("comment") or ""))

        assert len(generated_structures) == len(expected_structures), (
            f"Number of generated structures ({len(generated_structures)}) "
            f"does not match expected ({len(expected_structures)})"
        )

        for gen, exp in zip(generated_structures, expected_structures):
            # Compare labels/comments
            assert gen.info.get("comment") == exp.info.get(
                "comment"
            ), f"Comment mismatch: {gen.info.get('comment')} != {exp.info.get('comment')}"

            # Compare atom symbols
            assert list(gen.symbols) == list(
                exp.symbols
            ), f"Atom symbols mismatch for {gen.info.get('comment')}"

            # Compare coordinates
            np.allclose(
                np.asarray(gen.positions, dtype=float),
                np.asarray(exp.positions, dtype=float),
                atol=5e-5,
            )

    finally:
        # Restore CWD
        os.chdir(original_cwd)


def test_run_combinations_progress_callback(
    iterate_integration_config_file,
    iterate_input_directory,
    iterate_jobrunner,
    tmpdir,
):
    """run_combinations reports display-only progress once per combination.

    The callback must fire an initial ``(0, total)`` tick followed by exactly
    one increasing tick per finalized combination, ending at ``(total, total)``
    with each combination counted exactly once.
    """
    original_cwd = os.getcwd()
    os.chdir(iterate_input_directory)
    try:
        with open(iterate_integration_config_file, "r") as f:
            raw_config = yaml.safe_load(f)
        config = validate_yaml_config(
            raw_config, iterate_integration_config_file
        )
        job_settings = IterateJobSettings(
            config_file=iterate_integration_config_file,
        )
        job_settings.skeleton_list = config["skeletons"]
        job_settings.substituent_list = config["substituents"]

        job = IterateJob(
            settings=job_settings,
            jobrunner=iterate_jobrunner,
            nprocs=1,
            outputfile=str(tmpdir / "progress_output"),
        )

        pool, combinations, _input_errors, _sites = (
            iterate_jobrunner._generate_combinations(job)
        )
        assert combinations, "Expected at least one combination to run"
        total = len(combinations)

        calls = []

        def _record(completed, reported_total):
            calls.append((completed, reported_total))

        results = iterate_jobrunner.run_combinations(
            pool,
            combinations,
            nprocs=1,
            progress_callback=_record,
        )

        # One result per combination, in original order.
        assert len(results) == total

        # Every tick reports the same, correct total.
        assert all(reported == total for _c, reported in calls)

        completed_values = [completed for completed, _t in calls]
        # Exactly one initial (0) tick plus one finalize tick per combination,
        # strictly increasing by one: [0, 1, 2, ..., total]. This proves each
        # combination is counted exactly once and the bar reaches 100%.
        assert completed_values == list(range(0, total + 1))
    finally:
        os.chdir(original_cwd)


def test_iterate_timeout(
    iterate_timeout_config_file,
    iterate_input_directory,
    iterate_jobrunner,
    tmpdir,
    caplog,
):
    """
    Test that the timeout mechanism works correctly.
    We set a very short timeout (e.g. 0.001s) which
    should cause the worker to fail due to timeout.
    """
    import logging

    # Change CWD to input directory for relative file paths
    original_cwd = os.getcwd()
    os.chdir(iterate_input_directory)

    # Capture logs. The per-combination timeout notice is emitted at DEBUG so
    # it stays out of normal terminal output; capture at DEBUG to observe it.
    caplog.set_level(logging.DEBUG)

    try:
        # 1. Load Config
        with open(iterate_timeout_config_file, "r") as f:
            raw_config = yaml.safe_load(f)
        config = validate_yaml_config(raw_config, iterate_timeout_config_file)

        # 2. Setup Job with very short timeout
        job_settings = IterateJobSettings(
            config_file=iterate_timeout_config_file,
        )
        job_settings.skeleton_list = config["skeletons"]
        job_settings.substituent_list = config["substituents"]

        output_file = tmpdir / "timeout_output"

        jobrunner = iterate_jobrunner
        job = IterateJob(
            settings=job_settings,
            jobrunner=jobrunner,
            nprocs=1,  # Use 1 proc to ensure we hit it
            timeout=0.00000001,  # Ultra short timeout
            outputfile=str(output_file),
        )

        # 4. Run Job
        # It should not raise an exception, but handle the timeout gracefully
        job.run()

        # 5. Verify results
        # Check logs for timeout warning
        # Expected log from runner.py: "Timeout
        # ({timeout}s) for combination: {label}"

        # We need to construct the label to search for
        # Carbene1_34OTf
        label = "Carbene1_34OTf"

        found_timeout_log = False
        for record in caplog.records:
            if "Timeout" in record.message and label in record.message:
                found_timeout_log = True
                break

        assert (
            found_timeout_log
        ), f"Timeout log not found for the combination {label}. Logs: {[r.message for r in caplog.records]}"

    finally:
        os.chdir(original_cwd)


def test_iterate_template_generation(tmpdir, iterate_template_file):
    """
    Test that the iterate configuration template is
    generated correctly and matches the golden copy.
    """
    from chemsmart.utils.iterate import generate_yaml_template

    # 1. Generate template
    generated_path = tmpdir / "test_template.yaml"
    generate_yaml_template(str(generated_path))

    # 2. Assert file exists
    assert generated_path.exists()

    # 3. Compare content with expected template
    with open(generated_path, "r") as f:
        generated_content = f.read()

    with open(iterate_template_file, "r") as f:
        expected_content = f.read()

    # Normalize newlines and strip whitespace for robust comparison
    assert (
        generated_content.strip() == expected_content.strip()
    ), "Generated template does not match expected template content."

    # 4. Verify it is valid YAML
    parsed = yaml.safe_load(generated_content)
    assert "skeletons" in parsed
    assert "substituents" in parsed
    assert len(parsed["skeletons"]) == 1  # Based on current template examples
    assert len(parsed["substituents"]) == 1
    # Algorithm layer is documented in the template with a default block.
    assert "algorithm" in parsed
    assert parsed["algorithm"]["name"] == "lagrange_multipliers"


def test_iterate_validation_fails_on_invalid_link_index(
    iterate_invalid_skeleton_link_index_config_file,
):
    """
    Test that the validation logic correctly
    identifies when a link_index is not
    included in skeleton_indices.
    """
    from click.testing import CliRunner

    from chemsmart.cli.iterate.yaml_cmd import yaml_cmd

    runner = CliRunner()
    # Pass obj={} to initialize context object, required by MyGroup middleware
    result = runner.invoke(
        yaml_cmd,
        ["-f", iterate_invalid_skeleton_link_index_config_file],
        obj={},
    )

    assert result.exit_code != 0
    # Click 8.x format for BadParameter: "Error: Invalid value for ..."
    # We check for key parts of the error message
    assert "Invalid value" in result.output
    assert (
        "The link_index [6] is not included in 'skeleton_indices'"
        in result.output
    )


def test_iterate_validation_failures_comprehensive(tmpdir):
    """
    Test various validation failure scenarios
    using dynamically generated configs.
    Covers missing required fields and forbidden keys.
    """
    from click.testing import CliRunner

    from chemsmart.cli.iterate.yaml_cmd import yaml_cmd

    test_cases = [
        # Case 1: Skeleton missing file_path
        (
            """\
skeletons:
  - label: "s1"
    link_index: "1"
substituents:
  - file_path: "sub.xyz"
    label: "sub1"
    link_index: "1"
    groups: [1]
""",
            ["Missing required field 'file_path'", "Skeleton entry 1"],
            "Skeleton missing file_path",
        ),
        # Case 2: Substituent missing file_path
        (
            """\
skeletons:
  - file_path: "skel.xyz"
    label: "s1"
    link_index: "1"
substituents:
  - label: "sub1"
    link_index: "1"
    groups: [1]
""",
            ["Missing required field 'file_path'", "Substituent entry 1"],
            "Substituent missing file_path",
        ),
        # Case 3: Skeleton has forbidden key 'smiles'
        (
            """\
skeletons:
  - file_path: "skel.xyz"
    label: "s1"
    link_index: "1"
    smiles: "C"
substituents:
  - file_path: "sub.xyz"
    label: "sub1"
    link_index: "1"
    groups: [1]
""",
            ["Unknown key(s) in skeleton entry 1", "{'smiles'}"],
            "Skeleton forbidden key smiles",
        ),
        # Case 4: Substituent has forbidden key 'pubchem'
        (
            """\
skeletons:
  - file_path: "skel.xyz"
    label: "s1"
    link_index: "1"
substituents:
  - file_path: "sub.xyz"
    label: "sub1"
    link_index: "1"
    pubchem: "100"
    groups: [1]
""",
            ["Unknown key(s) in substituent entry 1", "{'pubchem'}"],
            "Substituent forbidden key pubchem",
        ),
        # Case 5: Random garbage key
        (
            """\
skeletons:
  - file_path: "skel.xyz"
    label: "s1"
    link_index: "1"
    random_key: "garbage"
""",
            ["Unknown key(s) in skeleton entry 1", "{'random_key'}"],
            "Skeleton random garbage key",
        ),
        # Case 6: Zero link_index (S2 Check)
        (
            """\
skeletons:
  - file_path: "skel.xyz"
    label: "s1"
    link_index: "0"
""",
            ["Found invalid index <= 0", "link_index"],
            "Skeleton zero link_index",
        ),
        # Case 7: Negative skeleton_indices (S2 Check)
        (
            """\
skeletons:
  - file_path: "skel.xyz"
    label: "s1"
    link_index: "1"
    skeleton_indices: [1, -5, 3]
""",
            ["Found invalid index <= 0", "skeleton_indices"],
            "Skeleton negative skeleton_indices",
        ),
        # Case 8: Empty skeleton_indices (S2 Check - Updated)
        (
            """\
skeletons:
  - file_path: "skel.xyz"
    label: "s1"
    link_index: "1"
    skeleton_indices: []
""",
            ["Found empty list in 'skeleton_indices'", "Skeleton entry 1"],
            "Skeleton empty skeleton_indices",
        ),
        # Case 9: Empty link_index (S2 Check - Updated)
        (
            """\
skeletons:
  - file_path: "skel.xyz"
    label: "s1"
    link_index: []
""",
            ["Found empty list in 'link_index'", "Skeleton entry 1"],
            "Skeleton empty link_index",
        ),
        # Case 10: Substituent multiple link_indices (S3 Check)
        (
            """\
skeletons:
  - file_path: "skel.xyz"
    label: "s1"
    link_index: "1"
substituents:
  - file_path: "sub.xyz"
    label: "sub1"
    link_index: "1, 2"
    groups: [1]
""",
            ["Multiple values found in 'link_index'", "exactly one link atom"],
            "Substituent multiple link_index",
        ),
        # Case 11: Skeleton label with invalid
        # characters (S4 Check: Safe Label)
        (
            """\
skeletons:
  - file_path: "skel.xyz"
    label: "s1/unsafe"
    link_index: "1"
""",
            ["Contains invalid characters", "Allowed characters"],
            "Skeleton label unsafe characters",
        ),
        # Case 12: Substituent label with invalid
        # characters (S4 Check: Safe Label)
        (
            """\
skeletons:
  - file_path: "skel.xyz"
    label: "s1"
    link_index: "1"
substituents:
  - file_path: "sub.xyz"
    label: "sub .. 1"
    link_index: "1"
    groups: [1]
""",
            ["Contains invalid characters", "Allowed characters"],
            "Substituent label unsafe characters",
        ),
        # Case 13: Label with space
        (
            """\
skeletons:
  - file_path: "skel.xyz"
    label: "s1 "
    link_index: "1"
""",
            ["Contains invalid characters"],
            "Label with space",
        ),
    ]

    runner = CliRunner()

    for idx, (config_content, expected_fragments, case_name) in enumerate(
        test_cases
    ):
        config_file = tmpdir / f"test_config_{idx}.yaml"
        with open(config_file, "w") as f:
            f.write(config_content)

        result = runner.invoke(yaml_cmd, ["-f", str(config_file)], obj={})

        assert (
            result.exit_code != 0
        ), f"Case '{case_name}' failed to raise error"
        assert (
            "Invalid value" in result.output
        ), f"Case '{case_name}' missing generic invalid value message"
        for fragment in expected_fragments:
            assert fragment in result.output, (
                f"Case '{case_name}' failed. "
                f"Expected fragment '{fragment}' not found in output:\n{result.output}"
            )


def test_iterate_runner_bounds_validation(tmpdir, fake_iterate_jobrunner):
    """
    Test that IterateJobRunner correctly
    validates indices against molecule size.
    (S2) Check: indices > num_atoms checking LOGS, not exceptions.
    """
    from unittest.mock import MagicMock, patch

    from chemsmart.io.molecules.structure import Molecule

    # Setup wrapper for the test
    runner = fake_iterate_jobrunner

    # Mock Molecule.from_filepath to return a predictable molecule
    # Create a dummy molecule with 5 atoms
    mock_mol = MagicMock(spec=Molecule)
    mock_mol.num_atoms = 5
    # Configure mock to return itself when from_filepath is called

    with (
        patch(
            "chemsmart.io.molecules.structure.Molecule.from_filepath",
            return_value=mock_mol,
        ),
        patch("chemsmart.jobs.iterate.runner.logger") as mock_logger,
    ):

        # Case 1: Skeleton link_index out of bounds
        config_1 = {
            "file_path": "dummy.xyz",
            "label": "skel1",
            "link_index": [6],  # > 5
        }
        res, label = runner._load_molecule(config_1, "skeleton", 0)

        assert res is None
        assert label == "skel1"
        # Check logs for "out of bounds"
        # We need to check if ANY of the call args contain our message
        found_error = False
        for call_args in mock_logger.error.call_args_list:
            msg = call_args[0][0]
            if "out of bounds" in msg and "link_index [6]" in msg:
                found_error = True
                break
        assert found_error, "Failed to log out-of-bounds link_index error"

        # Case 2: Skeleton skeleton_indices out of bounds
        config_2 = {
            "file_path": "dummy.xyz",
            "label": "skel2",
            "link_index": [1],  # Valid
            "skeleton_indices": [1, 2, 8],  # 8 > 5
        }
        res, label = runner._load_molecule(config_2, "skeleton", 1)

        assert res is None
        assert label == "skel2"

        found_error = False
        for call_args in mock_logger.error.call_args_list:
            msg = call_args[0][0]
            if "out of bounds" in msg and "skeleton_indices [8]" in msg:
                found_error = True
                break
        assert (
            found_error
        ), "Failed to log out-of-bounds skeleton_indices error"

        # Case 3: Valid input
        config_3 = {
            "file_path": "dummy.xyz",
            "label": "skel3",
            "link_index": [1],
            "skeleton_indices": [1, 2, 3],
        }
        res, label = runner._load_molecule(config_3, "skeleton", 2)
        assert res is not None
        assert res == mock_mol


def test_iterate_cli_pipeline_success(
    iterate_configs_directory,
    iterate_input_directory,
    iterate_expected_output_directory,
    tmpdir,
):
    """
    Test the full Iterate pipeline via the CLI:
    1. Run 'chemsmart iterate yaml -f config.yaml'
    2. Verify success exit code
    3. Verify output file exists and matches expected content.
    This ensures that the CLI entry point
    correctly orchestrates the job runner.
    """

    from click.testing import CliRunner

    from chemsmart.cli.iterate.yaml_cmd import yaml_cmd

    # Use the renamed config file which represents a valid CLI happy path
    config_file = os.path.join(
        iterate_configs_directory, "cli_happy_path.yaml"
    )
    expected_output_file = os.path.join(
        iterate_expected_output_directory, "cli_happy_path.xyz"
    )

    # Define output path in tmp directory (without extension for -o argument)
    output_base_path = str(tmpdir / "cli_happy_path_out")
    output_xyz_path = output_base_path + ".xyz"

    # Change CWD to input directory so relative paths in configuration work
    original_cwd = os.getcwd()
    os.chdir(iterate_input_directory)

    runner = CliRunner()

    try:
        # Pass obj={} to initialize context object
        result = runner.invoke(
            yaml_cmd,
            [
                "-f",
                config_file,
                "-o",
                output_base_path,
            ],
            obj={},
        )

        assert (
            result.exit_code == 0
        ), f"CLI execution failed. Output:\n{result.output}"
        assert os.path.exists(
            output_xyz_path
        ), "Output XYZ file was not generated."

        # Verify content matches
        with open(expected_output_file, "r") as f_exp:
            expected_content = f_exp.read().strip()

        with open(output_xyz_path, "r") as f_out:
            generated_content = f_out.read().strip()

        # Basic check: Ensuring specific label
        # presence which confirms combination logic ran
        # The runner combines skeleton label,
        # link index, sub label, and link index.
        # e.g. Carbene1 + 8 + OTf -> Carbene1_8OTf
        assert "Carbene1_8OTf" in generated_content

        # Verify line count matches (structure completeness check)
        exp_lines = expected_content.splitlines()
        gen_lines = generated_content.splitlines()
        assert len(exp_lines) == len(
            gen_lines
        ), "Generated XYZ line count differs from expected."

    finally:
        os.chdir(original_cwd)


# ---------------------------------------------------------------------------
# Algorithm layer: config resolution, YAML validation, and CLI priority
# ---------------------------------------------------------------------------

_ALGO_SKELETONS_SUBS = """\
skeletons:
  - file_path: "skel.xyz"
    label: "s1"
    link_index: "1"
substituents:
  - file_path: "sub.xyz"
    label: "sub1"
    link_index: "1"
    groups: [1]
"""


def _write_iterate_config(tmpdir, name, algorithm_block=None):
    """Write a minimal iterate YAML config (optionally with algorithm block).

    Molecule file paths are placeholders; tests that use this helper patch
    IterateJob so molecules are never actually loaded.
    """
    content = ""
    if algorithm_block is not None:
        content += algorithm_block + "\n"
    content += _ALGO_SKELETONS_SUBS
    config_file = tmpdir / name
    with open(config_file, "w") as f:
        f.write(content)
    return str(config_file)


def _run_yaml_cmd_capture_settings(args, config_file):
    """Invoke yaml_cmd with IterateJob patched; return (result, settings).

    ``settings`` is the IterateJobSettings passed to IterateJob, or None if
    the job was never constructed.
    """
    from unittest.mock import patch

    from click.testing import CliRunner

    from chemsmart.cli.iterate.yaml_cmd import yaml_cmd

    runner = CliRunner()
    with patch("chemsmart.cli.iterate.yaml_cmd.IterateJob") as mock_job:
        from chemsmart.jobs.iterate.runner import IterateRunSummary

        mock_job.return_value.run.return_value = IterateRunSummary(
            total=1,
            succeeded=1,
            failed=0,
            timed_out=0,
            structures_written=1,
            error_codes=[],
            output_paths=["out.xyz"],
            exit_code=0,
        )
        result = runner.invoke(yaml_cmd, ["-f", config_file, *args], obj={})
    settings = None
    if mock_job.call_args is not None:
        settings = mock_job.call_args.kwargs.get("settings")
    return result, settings


def test_resolve_algorithm_config_priority():
    """Unit test the default < YAML < CLI merge rules."""
    from chemsmart.jobs.iterate.settings import resolve_algorithm_config

    # Default: lagrange_multipliers with spec defaults.
    cfg = resolve_algorithm_config()
    assert cfg.name == "lagrange_multipliers"
    assert cfg.options == {
        "sphere_direction_samples_num": 96,
        "axial_rotations_sample_num": 6,
    }

    yaml_algo = {
        "name": "lagrange_multipliers",
        "options": {
            "sphere_direction_samples_num": 200,
            "axial_rotations_sample_num": 10,
        },
    }

    # YAML overrides defaults.
    cfg = resolve_algorithm_config(yaml_algorithm=yaml_algo)
    assert cfg.options["sphere_direction_samples_num"] == 200
    assert cfg.options["axial_rotations_sample_num"] == 10

    # CLI options override YAML; unset CLI options keep YAML value.
    cfg = resolve_algorithm_config(
        yaml_algorithm=yaml_algo,
        cli_algorithm_name="lagrange_multipliers",
        cli_options={"sphere_direction_samples_num": 128},
    )
    assert cfg.options["sphere_direction_samples_num"] == 128
    assert cfg.options["axial_rotations_sample_num"] == 10

    # CLI algorithm name overrides YAML name.
    cfg = resolve_algorithm_config(
        yaml_algorithm=yaml_algo, cli_algorithm_name="etkdg"
    )
    assert cfg.name == "etkdg"

    # Switching algorithm discards the other algorithm's YAML options and
    # falls back to the target algorithm's own defaults.
    assert "sphere_direction_samples_num" not in cfg.options
    assert cfg.options["num_conformers"] == 10
    assert cfg.options["use_global_optimization"] is False

    # Alias normalization.
    cfg = resolve_algorithm_config(cli_algorithm_name="lagrange")
    assert cfg.name == "lagrange_multipliers"


@pytest.mark.parametrize(
    "algorithm,options",
    [
        ("lagrange_multipliers", {"sphere_direction_samples_num": 0}),
        ("lagrange_multipliers", {"axial_rotations_sample_num": -1}),
        ("lagrange_multipliers", {"sphere_direction_samples_num": 1.5}),
        ("etkdg", {"num_conformers": 0}),
        ("etkdg", {"max_iterations": -5}),
        ("etkdg", {"random_seed": "seed"}),
        ("etkdg", {"use_global_optimization": "yes"}),
        ("etkdg", {"enforce_chirality": 1}),
        ("etkdg", {"force_field": "xyz"}),
    ],
)
def test_validate_algorithm_options_rejects_bad_values(algorithm, options):
    """Out-of-range / wrong-type / unknown-enum option values are rejected."""
    from chemsmart.jobs.iterate.settings import validate_algorithm_options

    with pytest.raises(ValueError):
        validate_algorithm_options(algorithm, options)


def test_validate_algorithm_options_accepts_valid_values():
    """Valid option values pass validation unchanged."""
    from chemsmart.jobs.iterate.settings import validate_algorithm_options

    result = validate_algorithm_options(
        "etkdg", {"num_conformers": 5, "random_seed": -1, "force_field": "uff"}
    )
    assert result == {
        "num_conformers": 5,
        "random_seed": -1,
        "force_field": "uff",
    }


def test_settings_invalid_combination_mode_rejected():
    """IterateJobSettings rejects an unknown combination_mode."""
    with pytest.raises(ValueError, match="combination_mode"):
        IterateJobSettings(combination_mode="bogus")


def test_settings_normalizes_algorithm_and_fills_defaults():
    """A directly-built config is normalized and completed on construction."""
    from chemsmart.jobs.iterate.settings import IterateAlgorithmConfig

    # Alias name + no options -> canonical name + full defaults.
    settings = IterateJobSettings(
        algorithm_config=IterateAlgorithmConfig(name="lagrange")
    )
    assert settings.algorithm_config.name == "lagrange_multipliers"
    assert settings.algorithm_config.options == {
        "sphere_direction_samples_num": 96,
        "axial_rotations_sample_num": 6,
    }

    # Partial options are completed with the algorithm defaults.
    settings = IterateJobSettings(
        algorithm_config=IterateAlgorithmConfig(
            name="etkdg", options={"num_conformers": 3}
        )
    )
    options = settings.algorithm_config.options
    assert options["num_conformers"] == 3
    assert options["max_iterations"] == 2000
    assert options["force_field"] == "none"


def test_settings_rejects_unknown_algorithm():
    """A directly-built config with an unknown algorithm name is rejected."""
    from chemsmart.jobs.iterate.settings import IterateAlgorithmConfig

    with pytest.raises(ValueError, match="Unknown algorithm"):
        IterateJobSettings(
            algorithm_config=IterateAlgorithmConfig(name="bogus")
        )


def test_settings_rejects_invalid_option_value():
    """A directly-built config with an invalid option value is rejected."""
    from chemsmart.jobs.iterate.settings import IterateAlgorithmConfig

    with pytest.raises(ValueError, match="Invalid value"):
        IterateJobSettings(
            algorithm_config=IterateAlgorithmConfig(
                name="etkdg", options={"num_conformers": -1}
            )
        )


def test_algorithm_config_copy_is_isolated():
    """Mutating a copied config must not affect the original."""
    from chemsmart.jobs.iterate.settings import resolve_algorithm_config

    cfg = resolve_algorithm_config()
    clone = cfg.copy()
    clone.options["sphere_direction_samples_num"] = 999
    clone.name = "etkdg"
    assert cfg.options["sphere_direction_samples_num"] == 96
    assert cfg.name == "lagrange_multipliers"


def test_algorithm_spec_default_options_immutable():
    """A registered spec's default_options cannot be mutated in place."""
    from chemsmart.jobs.iterate.settings import get_algorithm_spec

    spec = get_algorithm_spec("etkdg")
    with pytest.raises(TypeError):
        spec.default_options["num_conformers"] = 1


def test_build_analyzer_dispatch_selects_analyzer():
    """The registry dispatch routes each algorithm to its own analyzer."""
    from chemsmart.io.molecules.structure import Molecule
    from chemsmart.jobs.iterate.iterate import (
        IterateAnalyzer,
        IterateETKDGAnalyzer,
    )
    from chemsmart.jobs.iterate.runner import build_analyzer
    from chemsmart.jobs.iterate.settings import resolve_algorithm_config

    skeleton = Molecule(symbols=["C", "H"], positions=[[0, 0, 0], [1, 0, 0]])
    substituent = Molecule(
        symbols=["O", "H"], positions=[[0, 0, 0], [1, 0, 0]]
    )

    lagrange = build_analyzer(
        resolve_algorithm_config(), skeleton, [(substituent, 1, 1)]
    )
    etkdg = build_analyzer(
        resolve_algorithm_config(cli_algorithm_name="etkdg"),
        skeleton,
        [(substituent, 1, 1)],
    )
    assert isinstance(lagrange, IterateAnalyzer)
    assert isinstance(etkdg, IterateETKDGAnalyzer)


def test_algorithm_spec_requires_builder():
    """Every AlgorithmSpec must ship an analyzer builder."""
    from chemsmart.jobs.iterate.settings import AlgorithmSpec

    # Omitting the builder is a construction-time TypeError (required field).
    with pytest.raises(TypeError):
        AlgorithmSpec(
            canonical_name="dummy",
            aliases=("dummy",),
            default_options={},
        )
    # Explicitly passing None is rejected by the consistency check.
    with pytest.raises(ValueError, match="analyzer_builder is required"):
        AlgorithmSpec(
            canonical_name="dummy",
            aliases=("dummy",),
            default_options={},
            analyzer_builder=None,
        )


def test_algorithm_spec_rejects_canonical_not_in_aliases():
    """The canonical name must be one of the spec's own aliases."""
    from chemsmart.jobs.iterate.settings import AlgorithmSpec

    with pytest.raises(ValueError, match="canonical name must be included"):
        AlgorithmSpec(
            canonical_name="dummy",
            aliases=("other",),
            default_options={},
            analyzer_builder=lambda *a: None,
        )


def test_algorithm_spec_rejects_default_without_validator():
    """Every default option must have a matching validator."""
    from chemsmart.jobs.iterate.settings import AlgorithmSpec

    with pytest.raises(ValueError, match="have no validator"):
        AlgorithmSpec(
            canonical_name="dummy",
            aliases=("dummy",),
            default_options={"foo": 1},
            analyzer_builder=lambda *a: None,
        )


def test_build_alias_map_rejects_duplicates():
    """The registry rejects duplicate canonical names and duplicate aliases."""
    from chemsmart.jobs.iterate.settings import (
        AlgorithmSpec,
        _build_alias_map,
    )

    def _spec(canonical, aliases):
        return AlgorithmSpec(
            canonical_name=canonical,
            aliases=aliases,
            default_options={},
            analyzer_builder=lambda *a: None,
        )

    with pytest.raises(ValueError, match="Duplicate canonical"):
        _build_alias_map([_spec("a", ("a",)), _spec("a", ("a", "a2"))])

    with pytest.raises(ValueError, match="Duplicate algorithm alias"):
        _build_alias_map(
            [_spec("a", ("a", "shared")), _spec("b", ("b", "shared"))]
        )


def test_build_analyzer_rejects_multi_substituent_for_single_only():
    """A single-substituent-only algorithm rejects multi-substituent input."""
    from chemsmart.io.molecules.structure import Molecule
    from chemsmart.jobs.iterate.runner import build_analyzer
    from chemsmart.jobs.iterate.settings import resolve_algorithm_config

    skeleton = Molecule(symbols=["C", "H"], positions=[[0, 0, 0], [1, 0, 0]])
    substituent = Molecule(
        symbols=["O", "H"], positions=[[0, 0, 0], [1, 0, 0]]
    )

    # Lagrange advertises supports_multi_substituent=False.
    with pytest.raises(NotImplementedError, match="multi-substituent"):
        build_analyzer(
            resolve_algorithm_config(),
            skeleton,
            [(substituent, 1, 1), (substituent, 2, 1)],
        )


def test_iterate_cli_invalid_yaml_option_value_errors(tmpdir):
    """An out-of-range option value in the YAML block is rejected by the CLI."""
    from click.testing import CliRunner

    from chemsmart.cli.iterate.yaml_cmd import yaml_cmd

    algo = "algorithm:\n  name: etkdg\n  options:\n    num_conformers: -1\n"
    config_file = _write_iterate_config(tmpdir, "bad_value.yaml", algo)
    result = CliRunner().invoke(yaml_cmd, ["-f", config_file], obj={})
    assert result.exit_code != 0
    assert "Invalid value" in result.output
    assert "num_conformers" in result.output


def test_validate_yaml_config_algorithm_block(tmpdir):
    """YAML algorithm block is validated and normalized."""
    # No algorithm block -> None.
    cfg = validate_yaml_config(yaml.safe_load(_ALGO_SKELETONS_SUBS), "c.yaml")
    assert cfg["algorithm"] is None

    # Alias name is normalized to canonical form.
    raw = yaml.safe_load(
        "algorithm:\n  name: lagrange\n" + _ALGO_SKELETONS_SUBS
    )
    cfg = validate_yaml_config(raw, "c.yaml")
    assert cfg["algorithm"]["name"] == "lagrange_multipliers"

    # Unknown algorithm raises.
    import click

    raw = yaml.safe_load("algorithm:\n  name: bogus\n" + _ALGO_SKELETONS_SUBS)
    with pytest.raises(click.BadParameter, match="Unknown algorithm"):
        validate_yaml_config(raw, "c.yaml")

    # Unknown option raises.
    raw = yaml.safe_load(
        "algorithm:\n  name: lagrange_multipliers\n  options:\n"
        "    bad_option: 1\n" + _ALGO_SKELETONS_SUBS
    )
    with pytest.raises(click.BadParameter, match="Unknown option"):
        validate_yaml_config(raw, "c.yaml")


def test_iterate_cli_default_algorithm_when_none_specified(tmpdir):
    """No YAML algorithm and no CLI subcommand -> default Lagrange."""
    config_file = _write_iterate_config(tmpdir, "default.yaml")
    result, settings = _run_yaml_cmd_capture_settings([], config_file)
    assert result.exit_code == 0, result.output
    assert settings is not None
    assert settings.algorithm_config.name == "lagrange_multipliers"
    assert settings.algorithm_config.options == {
        "sphere_direction_samples_num": 96,
        "axial_rotations_sample_num": 6,
    }


def test_iterate_cli_yaml_algorithm_applied(tmpdir):
    """YAML algorithm block is used when no CLI subcommand is given."""
    algo = (
        "algorithm:\n"
        "  name: lagrange_multipliers\n"
        "  options:\n"
        "    sphere_direction_samples_num: 200\n"
        "    axial_rotations_sample_num: 10\n"
    )
    config_file = _write_iterate_config(tmpdir, "yaml_algo.yaml", algo)
    result, settings = _run_yaml_cmd_capture_settings([], config_file)
    assert result.exit_code == 0, result.output
    assert settings.algorithm_config.name == "lagrange_multipliers"
    assert (
        settings.algorithm_config.options["sphere_direction_samples_num"]
        == 200
    )
    assert (
        settings.algorithm_config.options["axial_rotations_sample_num"] == 10
    )


def test_iterate_cli_lagrange_overrides_yaml(tmpdir):
    """CLI lagrange options override the matching YAML values."""
    algo = (
        "algorithm:\n"
        "  name: lagrange_multipliers\n"
        "  options:\n"
        "    sphere_direction_samples_num: 200\n"
        "    axial_rotations_sample_num: 10\n"
    )
    config_file = _write_iterate_config(tmpdir, "override.yaml", algo)
    result, settings = _run_yaml_cmd_capture_settings(
        [
            "lagrange",
            "--sphere-direction-samples-number",
            "128",
            "--axial-rotations-sample-number",
            "8",
        ],
        config_file,
    )
    assert result.exit_code == 0, result.output
    assert (
        settings.algorithm_config.options["sphere_direction_samples_num"]
        == 128
    )
    assert settings.algorithm_config.options["axial_rotations_sample_num"] == 8


def test_iterate_cli_unset_option_does_not_override_yaml(tmpdir):
    """CLI options not explicitly passed keep their YAML value."""
    algo = (
        "algorithm:\n"
        "  name: lagrange_multipliers\n"
        "  options:\n"
        "    sphere_direction_samples_num: 200\n"
        "    axial_rotations_sample_num: 10\n"
    )
    config_file = _write_iterate_config(tmpdir, "partial.yaml", algo)
    # Only sphere is passed explicitly; axial must remain the YAML value (10),
    # not the click default (6).
    result, settings = _run_yaml_cmd_capture_settings(
        ["lagrange", "--sphere-direction-samples-number", "128"],
        config_file,
    )
    assert result.exit_code == 0, result.output
    assert (
        settings.algorithm_config.options["sphere_direction_samples_num"]
        == 128
    )
    assert (
        settings.algorithm_config.options["axial_rotations_sample_num"] == 10
    )


def test_iterate_cli_algorithm_name_overrides_yaml(tmpdir):
    """CLI algorithm subcommand overrides the YAML algorithm name."""
    algo = "algorithm:\n  name: lagrange_multipliers\n"
    config_file = _write_iterate_config(tmpdir, "switch.yaml", algo)
    result, settings = _run_yaml_cmd_capture_settings(
        ["etkdg", "--num-conformers", "50"], config_file
    )
    # The resolved settings must reflect the CLI-selected algorithm.
    assert settings is not None
    assert settings.algorithm_config.name == "etkdg"
    assert settings.algorithm_config.options["num_conformers"] == 50
    # The YAML lagrange options must not leak into the etkdg config.
    assert (
        "sphere_direction_samples_num" not in settings.algorithm_config.options
    )


def test_iterate_cli_etkdg_default_is_local(tmpdir):
    """The etkdg subcommand defaults to local embedding."""
    config_file = _write_iterate_config(tmpdir, "etkdg_local.yaml")
    result, settings = _run_yaml_cmd_capture_settings(["etkdg"], config_file)
    assert result.exit_code == 0, result.output
    assert settings.algorithm_config.name == "etkdg"
    assert (
        settings.algorithm_config.options["use_global_optimization"] is False
    )


def test_iterate_cli_etkdg_global_flag(tmpdir):
    """The etkdg '--global' flag enables global embedding."""
    config_file = _write_iterate_config(tmpdir, "etkdg_global.yaml")
    result, settings = _run_yaml_cmd_capture_settings(
        ["etkdg", "--global"], config_file
    )
    assert result.exit_code == 0, result.output
    assert settings.algorithm_config.options["use_global_optimization"] is True


def test_iterate_cli_etkdg_all_options_exposed(tmpdir):
    """Every etkdg option is settable from the CLI and flows into settings."""
    config_file = _write_iterate_config(tmpdir, "etkdg_opts.yaml")
    result, settings = _run_yaml_cmd_capture_settings(
        [
            "etkdg",
            "--num-conformers",
            "4",
            "--random-seed",
            "5",
            "--max-iterations",
            "300",
            "--no-random-coords",
            "--enforce-chirality",
            "--force-field",
            "mmff94",
        ],
        config_file,
    )
    assert result.exit_code == 0, result.output
    options = settings.algorithm_config.options
    assert options["num_conformers"] == 4
    assert options["random_seed"] == 5
    assert options["max_iterations"] == 300
    assert options["use_random_coordinates"] is False
    assert options["enforce_chirality"] is True
    assert options["force_field"] == "mmff94"


def test_iterate_cli_etkdg_yaml_options_applied(tmpdir):
    """ETKDG options set in the YAML algorithm block are honored."""
    algo = (
        "algorithm:\n"
        "  name: etkdg\n"
        "  options:\n"
        "    use_global_optimization: true\n"
        "    max_iterations: 500\n"
        "    force_field: uff\n"
    )
    config_file = _write_iterate_config(tmpdir, "etkdg_yaml.yaml", algo)
    # Only --num-conformers is passed explicitly; the YAML options must be
    # preserved (not overridden by Click defaults).
    result, settings = _run_yaml_cmd_capture_settings(
        ["etkdg", "--num-conformers", "3"], config_file
    )
    assert result.exit_code == 0, result.output
    options = settings.algorithm_config.options
    assert options["use_global_optimization"] is True
    assert options["max_iterations"] == 500
    assert options["force_field"] == "uff"
    assert options["num_conformers"] == 3


@pytest.mark.parametrize(
    "removed_option",
    [
        ["-m", "lagrange_multipliers"],
        ["-s", "96"],
        ["-a", "6"],
    ],
)
def test_iterate_cli_legacy_options_removed(tmpdir, removed_option):
    """The legacy -m/-s/-a options are no longer accepted."""
    from click.testing import CliRunner

    from chemsmart.cli.iterate.yaml_cmd import yaml_cmd

    config_file = _write_iterate_config(tmpdir, "legacy.yaml")
    runner = CliRunner()
    result = runner.invoke(
        yaml_cmd, ["-f", config_file, *removed_option], obj={}
    )
    assert result.exit_code != 0
    assert "No such option" in result.output


def test_iterate_cli_unknown_algorithm_errors(tmpdir):
    """An unknown YAML algorithm name is rejected by the CLI."""
    from click.testing import CliRunner

    from chemsmart.cli.iterate.yaml_cmd import yaml_cmd

    config_file = _write_iterate_config(
        tmpdir, "unknown_algo.yaml", "algorithm:\n  name: nonexistent\n"
    )
    runner = CliRunner()
    result = runner.invoke(yaml_cmd, ["-f", config_file], obj={})
    assert result.exit_code != 0
    assert "Unknown algorithm" in result.output


def test_iterate_cli_unknown_option_errors(tmpdir):
    """An unknown YAML algorithm option is rejected by the CLI."""
    from click.testing import CliRunner

    from chemsmart.cli.iterate.yaml_cmd import yaml_cmd

    algo = (
        "algorithm:\n"
        "  name: lagrange_multipliers\n"
        "  options:\n"
        "    bad_option: 5\n"
    )
    config_file = _write_iterate_config(tmpdir, "unknown_opt.yaml", algo)
    runner = CliRunner()
    result = runner.invoke(yaml_cmd, ["-f", config_file], obj={})
    assert result.exit_code != 0
    assert "Unknown option" in result.output


def test_iterate_etkdg_local_preserves_skeleton(iterate_input_directory):
    """ETKDG local mode joins the substituent and keeps the skeleton fixed."""
    from chemsmart.io.molecules.structure import Molecule
    from chemsmart.jobs.iterate.iterate import (
        IterateETKDGAnalyzer,
        SkeletonPreprocessor,
        SubstituentPreprocessor,
    )

    original_cwd = os.getcwd()
    os.chdir(iterate_input_directory)
    try:
        skeleton = Molecule.from_filepath("carbene1_opt.com")
        substituent = Molecule.from_filepath("OTf_opt.xyz")

        # Preprocess exactly as the runner does (remove old group at link 8).
        skel_prep = SkeletonPreprocessor(skeleton.copy(), 8, None)
        processed_skeleton = skel_prep.run()
        new_skel_link = skel_prep.get_new_link_index()

        sub_prep = SubstituentPreprocessor(substituent.copy(), 8)
        processed_sub = sub_prep.run()
        new_sub_link = sub_prep.get_new_link_index()

        n_skeleton = len(processed_skeleton)

        analyzer = IterateETKDGAnalyzer(
            skeleton=processed_skeleton,
            substituents=[(processed_sub, new_skel_link, new_sub_link)],
            options={
                "num_conformers": 2,
                "random_seed": 7,
                "use_global_optimization": False,
            },
        )
        result = analyzer.run()

        assert result is not None, "ETKDG local embedding returned no molecule"
        # Combined molecule = skeleton atoms followed by substituent atoms.
        assert len(result) == n_skeleton + len(processed_sub)
        # Local mode must keep the skeleton atoms exactly fixed.
        skeleton_displacement = np.linalg.norm(
            np.asarray(result.positions[:n_skeleton])
            - np.asarray(processed_skeleton.positions),
            axis=1,
        ).max()
        assert skeleton_displacement < 1e-4
    finally:
        os.chdir(original_cwd)


def _prepare_etkdg_multi(input_dir, skel_links, sub_link=8):
    """Build ``(processed_skeleton, substituents)`` exactly as the runner does.

    Batch-preprocesses the carbene1 skeleton at every ``skel_links`` position
    and preprocesses one OTf substituent per position, returning the
    ``(substituent, skeleton_link, substituent_link)`` tuples that
    :class:`IterateETKDGAnalyzer` consumes. Link indices are 1-based.
    """
    from chemsmart.io.molecules.structure import Molecule
    from chemsmart.jobs.iterate.iterate import SubstituentPreprocessor
    from chemsmart.jobs.iterate.runner import _batch_preprocess_skeleton

    skeleton = Molecule.from_filepath(
        os.path.join(input_dir, "carbene1_opt.com")
    )
    processed_skeleton, index_map = _batch_preprocess_skeleton(
        skeleton.copy(), list(skel_links), None
    )
    substituents = []
    for link in skel_links:
        sub = Molecule.from_filepath(os.path.join(input_dir, "OTf_opt.xyz"))
        prep = SubstituentPreprocessor(sub.copy(), sub_link)
        processed_sub = prep.run()
        substituents.append(
            (processed_sub, index_map[link], prep.get_new_link_index())
        )
    return processed_skeleton, substituents


def test_iterate_etkdg_multi_single_embed_call(
    iterate_input_directory, monkeypatch
):
    """A multi-substituent combination triggers exactly one ETKDG embedding."""
    from chemsmart.jobs.iterate import iterate as iterate_mod
    from chemsmart.jobs.iterate.iterate import IterateETKDGAnalyzer

    processed_skeleton, substituents = _prepare_etkdg_multi(
        iterate_input_directory, [4, 6, 7]
    )

    calls = {"n": 0}
    real_embed = iterate_mod.AllChem.EmbedMultipleConfs

    def counting_embed(mol, **kwargs):
        calls["n"] += 1
        return real_embed(mol, **kwargs)

    monkeypatch.setattr(
        iterate_mod.AllChem, "EmbedMultipleConfs", counting_embed
    )

    analyzer = IterateETKDGAnalyzer(
        skeleton=processed_skeleton,
        substituents=substituents,
        options={"num_conformers": 2, "random_seed": 42},
    )
    result = analyzer.run()

    assert result is not None
    # Three substituents, but ETKDG must run only once for the whole assembly.
    assert calls["n"] == 1
    expected = len(processed_skeleton) + sum(
        len(s) for s, _, _ in substituents
    )
    assert len(result) == expected


def test_iterate_etkdg_multi_preserves_connection_bonds(
    iterate_input_directory,
):
    """Every skeleton-substituent connection bond exists before embedding."""
    from rdkit import Chem

    from chemsmart.jobs.iterate.iterate import IterateETKDGAnalyzer

    processed_skeleton, substituents = _prepare_etkdg_multi(
        iterate_input_directory, [4, 6, 7, 8]
    )
    analyzer = IterateETKDGAnalyzer(
        skeleton=processed_skeleton, substituents=substituents
    )
    combined, n_skeleton = analyzer._build_combined_rdmol()

    offset = n_skeleton
    for substituent, skel_link, sub_link in analyzer.substituents:
        bond = combined.GetBondBetweenAtoms(skel_link, offset + sub_link)
        assert bond is not None, "connection bond missing in combined molecule"
        assert bond.GetBondType() == Chem.BondType.SINGLE
        offset += len(substituent)


def test_iterate_etkdg_multi_local_fixes_only_skeleton(
    iterate_input_directory,
):
    """Local mode pins only the original skeleton atoms, not substituents."""
    from chemsmart.jobs.iterate.iterate import IterateETKDGAnalyzer

    processed_skeleton, substituents = _prepare_etkdg_multi(
        iterate_input_directory, [4, 7]
    )
    analyzer = IterateETKDGAnalyzer(
        skeleton=processed_skeleton,
        substituents=substituents,
        options={"use_global_optimization": False},
    )
    combined, n_skeleton = analyzer._build_combined_rdmol()
    reference = analyzer._fixed_reference_positions(combined, n_skeleton)

    assert n_skeleton == len(processed_skeleton)
    # Only the original skeleton indices are fixed; substituents stay free.
    assert sorted(reference.keys()) == list(range(n_skeleton))
    assert combined.GetNumAtoms() > n_skeleton


def test_iterate_etkdg_multi_stays_single_fragment(
    iterate_input_directory,
):
    """The combined multi-substituent molecule is one connected fragment."""
    from rdkit import Chem

    from chemsmart.jobs.iterate.iterate import IterateETKDGAnalyzer

    processed_skeleton, substituents = _prepare_etkdg_multi(
        iterate_input_directory, [4, 6, 7, 8]
    )
    analyzer = IterateETKDGAnalyzer(
        skeleton=processed_skeleton, substituents=substituents
    )
    combined, _ = analyzer._build_combined_rdmol()
    assert len(Chem.GetMolFrags(combined)) == 1


def test_iterate_etkdg_multi_high_substitution_completes(
    iterate_input_directory,
):
    """A four-substituent combination completes in one ETKDG run and keeps
    the original skeleton fixed in local mode."""
    from chemsmart.jobs.iterate.iterate import IterateETKDGAnalyzer

    processed_skeleton, substituents = _prepare_etkdg_multi(
        iterate_input_directory, [4, 6, 7, 8]
    )
    n_skeleton = len(processed_skeleton)
    analyzer = IterateETKDGAnalyzer(
        skeleton=processed_skeleton,
        substituents=substituents,
        options={
            "num_conformers": 3,
            "random_seed": 42,
            "use_global_optimization": False,
        },
    )
    result = analyzer.run()

    assert result is not None, "four-substituent ETKDG returned no molecule"
    expected = n_skeleton + sum(len(s) for s, _, _ in substituents)
    assert len(result) == expected
    # Local mode keeps the original skeleton atoms exactly fixed.
    skeleton_displacement = np.linalg.norm(
        np.asarray(result.positions[:n_skeleton])
        - np.asarray(processed_skeleton.positions),
        axis=1,
    ).max()
    assert skeleton_displacement < 1e-4


# ---------------------------------------------------------------------------
# Fix 3: ETKDG random_seed int32 bounds
# ---------------------------------------------------------------------------


def test_validate_random_seed_bounds():
    """random_seed must stay within RDKit's int32 range and reject bools."""
    from chemsmart.jobs.iterate.settings import validate_algorithm_options

    # Boundary and typical values are accepted.
    for seed in (-(2**31), 2**31 - 1, -1, 0, 42):
        validate_algorithm_options("etkdg", {"random_seed": seed})

    # Values outside the int32 range are rejected (RDKit would overflow).
    for seed in (2**31, -(2**31) - 1, 10**18):
        with pytest.raises(ValueError):
            validate_algorithm_options("etkdg", {"random_seed": seed})

    # A bool must not be accepted as an integer seed.
    with pytest.raises(ValueError):
        validate_algorithm_options("etkdg", {"random_seed": True})


# ---------------------------------------------------------------------------
# Fix 1: run-result contract (success/failure counts, real paths)
# ---------------------------------------------------------------------------


def test_write_outputs_merged_all_failed_writes_nothing(
    tmpdir, fake_iterate_jobrunner
):
    """A fully-failed merged run must not create an empty output file."""
    from types import SimpleNamespace

    from chemsmart.jobs.iterate.report import STATUS_FAILED, CombinationResult

    runner = fake_iterate_jobrunner
    outfile = str(tmpdir / "out.xyz")
    job = SimpleNamespace(
        separate_outputs=False, output_directory=None, outputfile=outfile
    )
    results = [
        CombinationResult(1, "a", STATUS_FAILED),
        CombinationResult(2, "b", STATUS_FAILED),
    ]
    paths = runner._write_outputs(results, job)
    assert paths == []
    assert not os.path.exists(outfile)


def test_write_outputs_separate_returns_written_paths(
    tmpdir, fake_iterate_jobrunner
):
    """Separate-outputs mode returns one path per successful structure."""
    from types import SimpleNamespace

    from chemsmart.io.molecules.structure import Molecule
    from chemsmart.jobs.iterate.report import (
        STATUS_FAILED,
        STATUS_SUCCESS,
        CombinationResult,
    )

    mol = Molecule(symbols=["C", "H"], positions=[[0, 0, 0], [1, 0, 0]])
    runner = fake_iterate_jobrunner
    job = SimpleNamespace(
        separate_outputs=True,
        output_directory=str(tmpdir),
        outputfile=str(tmpdir / "unused.xyz"),
    )
    results = [
        CombinationResult(1, "m1", STATUS_SUCCESS, molecule=mol),
        CombinationResult(2, "bad", STATUS_FAILED),
        CombinationResult(3, "m2", STATUS_SUCCESS, molecule=mol),
    ]
    paths = runner._write_outputs(results, job)
    assert sorted(os.path.basename(p) for p in paths) == ["m1.xyz", "m2.xyz"]
    for p in paths:
        assert os.path.exists(p)


def test_write_outputs_separate_rejects_duplicate_filename(
    tmpdir, fake_iterate_jobrunner
):
    """Separate-outputs mode must not silently overwrite a same-named file."""
    from types import SimpleNamespace

    from chemsmart.io.molecules.structure import Molecule
    from chemsmart.jobs.iterate.report import (
        STATUS_SUCCESS,
        CombinationResult,
    )

    mol = Molecule(symbols=["C", "H"], positions=[[0, 0, 0], [1, 0, 0]])
    runner = fake_iterate_jobrunner
    job = SimpleNamespace(
        separate_outputs=True,
        output_directory=str(tmpdir),
        outputfile=str(tmpdir / "unused.xyz"),
    )
    results = [
        CombinationResult(1, "dup", STATUS_SUCCESS, molecule=mol),
        CombinationResult(2, "dup", STATUS_SUCCESS, molecule=mol),
    ]
    with pytest.raises(ValueError, match="Duplicate output filename"):
        runner._write_outputs(results, job)


def test_iterate_cli_all_failed_exits_nonzero(tmpdir):
    """When no structure is produced the CLI fails with a non-zero exit."""
    from click.testing import CliRunner

    from chemsmart.cli.iterate.yaml_cmd import yaml_cmd

    content = (
        "skeletons:\n"
        '  - file_path: "./missing_skeleton.xyz"\n'
        '    label: "s1"\n'
        '    link_index: "1"\n'
        "substituents:\n"
        '  - file_path: "./missing_sub.xyz"\n'
        '    label: "sub1"\n'
        '    link_index: "1"\n'
        "    groups: [1]\n"
    )
    config_file = tmpdir / "all_failed.yaml"
    with open(config_file, "w") as f:
        f.write(content)

    runner = CliRunner()
    with runner.isolated_filesystem():
        result = runner.invoke(yaml_cmd, ["-f", str(config_file)], obj={})
    assert result.exit_code == 1
    assert "ITR-INPUT-001" in result.output


# ---------------------------------------------------------------------------
# Fix 2: relative molecule paths resolve against the config directory
# ---------------------------------------------------------------------------


def test_validate_yaml_config_resolves_relative_paths(tmpdir):
    """Relative file_path values resolve against the config file's directory."""
    from chemsmart.utils.iterate import validate_yaml_config

    config_dir = tmpdir.mkdir("cfg")
    config_path = str(config_dir / "c.yaml")
    raw = {
        "skeletons": [
            {"file_path": "./mol.xyz", "label": "s1", "link_index": "1"}
        ],
        "substituents": [
            {
                "file_path": "sub.xyz",
                "label": "sub1",
                "link_index": "1",
                "groups": [1],
            }
        ],
    }
    cfg = validate_yaml_config(raw, config_path)
    assert cfg["skeletons"][0]["file_path"] == os.path.join(
        str(config_dir), "mol.xyz"
    )
    assert cfg["skeletons"][0]["file_path_raw"] == "./mol.xyz"
    assert cfg["substituents"][0]["file_path"] == os.path.join(
        str(config_dir), "sub.xyz"
    )


def test_validate_yaml_config_preserves_absolute_paths(tmpdir):
    """Absolute file_path values are left unchanged."""
    from chemsmart.utils.iterate import validate_yaml_config

    config_dir = tmpdir.mkdir("cfg2")
    abs_path = os.path.join(str(tmpdir), "abs_mol.xyz")
    raw = {
        "skeletons": [
            {"file_path": abs_path, "label": "s1", "link_index": "1"}
        ],
        "substituents": [
            {
                "file_path": "sub.xyz",
                "label": "sub1",
                "link_index": "1",
                "groups": [1],
            }
        ],
    }
    cfg = validate_yaml_config(raw, str(config_dir / "c.yaml"))
    assert cfg["skeletons"][0]["file_path"] == abs_path


# ---------------------------------------------------------------------------
# Fix 4: label / link-site collision rejection
# ---------------------------------------------------------------------------


def test_iterate_rejects_duplicate_skeleton_label():
    """Two skeletons sharing a label are rejected."""
    import click

    from chemsmart.utils.iterate import validate_yaml_config

    raw = {
        "skeletons": [
            {"file_path": "a.xyz", "label": "dup", "link_index": "1"},
            {"file_path": "b.xyz", "label": "dup", "link_index": "2"},
        ],
        "substituents": [
            {
                "file_path": "s.xyz",
                "label": "sub1",
                "link_index": "1",
                "groups": [1, 2],
            }
        ],
    }
    with pytest.raises(click.BadParameter, match="duplicate label"):
        validate_yaml_config(raw, "c.yaml")


def test_iterate_rejects_duplicate_substituent_label():
    """Two substituents sharing a label are rejected."""
    import click

    from chemsmart.utils.iterate import validate_yaml_config

    raw = {
        "skeletons": [
            {"file_path": "a.xyz", "label": "s1", "link_index": "1"}
        ],
        "substituents": [
            {
                "file_path": "x.xyz",
                "label": "dup",
                "link_index": "1",
                "groups": [1],
            },
            {
                "file_path": "y.xyz",
                "label": "dup",
                "link_index": "1",
                "groups": [1],
            },
        ],
    }
    with pytest.raises(click.BadParameter, match="duplicate label"):
        validate_yaml_config(raw, "c.yaml")


def test_iterate_rejects_duplicate_link_index():
    """A skeleton listing the same link_index twice is rejected."""
    import click

    from chemsmart.utils.iterate import validate_yaml_config

    raw = {
        "skeletons": [
            {"file_path": "a.xyz", "label": "s1", "link_index": [1, 1]}
        ],
        "substituents": [
            {
                "file_path": "s.xyz",
                "label": "sub1",
                "link_index": "1",
                "groups": [1],
            }
        ],
    }
    with pytest.raises(click.BadParameter, match="Duplicate link_index"):
        validate_yaml_config(raw, "c.yaml")


def test_iterate_rejects_overlapping_slot_links():
    """The same connection site used by two slots is rejected."""
    import click

    from chemsmart.utils.iterate import validate_yaml_config

    raw = {
        "skeletons": [
            {
                "file_path": "a.xyz",
                "label": "s1",
                "slots": [
                    {"group": 1, "link_indices": "5"},
                    {"group": 2, "link_indices": "5"},
                ],
            }
        ],
        "substituents": [
            {
                "file_path": "s.xyz",
                "label": "sub1",
                "link_index": "1",
                "groups": [1, 2],
            }
        ],
    }
    with pytest.raises(click.BadParameter, match="multiple slots"):
        validate_yaml_config(raw, "c.yaml")


def test_iterate_rejects_duplicate_link_within_slot():
    """A slot listing the same link index twice is rejected."""
    import click

    from chemsmart.utils.iterate import validate_yaml_config

    raw = {
        "skeletons": [
            {
                "file_path": "a.xyz",
                "label": "s1",
                "slots": [{"group": 1, "link_indices": [5, 5]}],
            }
        ],
        "substituents": [
            {
                "file_path": "s.xyz",
                "label": "sub1",
                "link_index": "1",
                "groups": [1],
            }
        ],
    }
    with pytest.raises(click.BadParameter, match="duplicate value"):
        validate_yaml_config(raw, "c.yaml")


# ---------------------------------------------------------------------------
# Run report (plain-text .out)
# ---------------------------------------------------------------------------


def _minimal_report(**overrides):
    from datetime import datetime

    from chemsmart.jobs.iterate.report import IterateReport

    now = datetime(2024, 1, 1, 12, 0, 0)
    data = dict(
        run_id="abc123def456",
        chemsmart_version="9.9.9",
        rdkit_version="2024.03.1",
        started_at=now,
        finished_at=now,
        duration_seconds=1.5,
        working_directory="/tmp/work",
        command_line="chemsmart run iterate yaml -f x.yaml",
        config_file="/path/x.yaml",
        config_sha256="deadbeef",
    )
    data.update(overrides)
    return IterateReport(**data)


def test_report_render_contains_all_sections():
    from chemsmart.jobs.iterate.report import (
        STATUS_SUCCESS,
        CombinationResult,
    )

    report = _minimal_report(
        algorithm_name="etkdg",
        algorithm_options={"num_conformers": 10, "random_seed": 42},
        combination_mode="independent",
        nprocs=1,
        timeout_seconds=120,
        output_mode="merged",
        total_combinations=1,
        per_skeleton_counts=[("sk", 1)],
        results=[
            CombinationResult(
                combination_number=1,
                label="sk_5me",
                execution_status=STATUS_SUCCESS,
                duration_seconds=0.5,
                output_path="/out.xyz",
                structure_index=1,
            )
        ],
        output_location="/out.xyz",
    )
    text = report.render()
    for header in (
        "CHEMSMART ITERATE JOB REPORT",
        "GENERAL INFORMATION",
        "INPUT STRUCTURES",
        "INPUT ERRORS",
        "ALGORITHM",
        "COMBINATION GENERATION",
        "EXECUTION RESULTS",
        "FAILED COMBINATIONS",
        "TIMED-OUT COMBINATIONS",
        "OUTPUT WRITE FAILURES",
        "OUTPUT STRUCTURES",
        "FINAL SUMMARY",
    ):
        assert header in text
    # Empty error sections are shown explicitly, never omitted.
    assert "None." in text
    assert "Normal termination of ChemSmart Iterate" in text
    assert "num_conformers" in text and "random_seed" in text


def test_report_error_termination_and_failed_section():
    from chemsmart.jobs.iterate.report import (
        STATUS_FAILED,
        CombinationResult,
    )

    report = _minimal_report(
        algorithm_name="etkdg",
        output_mode="merged",
        total_combinations=1,
        results=[
            CombinationResult(
                combination_number=1,
                label="sk_5me",
                execution_status=STATUS_FAILED,
                duration_seconds=0.1,
                failure_stage="ALGORITHM",
                error_type="NoSolution",
                error_message="Algorithm produced no structure.",
            )
        ],
    )
    text = report.render()
    assert "Error termination of ChemSmart Iterate" in text
    assert "NoSolution" in text
    assert "ALGORITHM" in text


def test_summarize_results_counts():
    from chemsmart.jobs.iterate.report import (
        ERROR_CODE_EXEC,
        ERROR_CODE_INPUT,
        ERROR_CODE_TIMEOUT,
        ERROR_CODE_WRITE,
        STATUS_FAILED,
        STATUS_SUCCESS,
        STATUS_TIMED_OUT,
        STATUS_WRITE_FAILED,
        CombinationResult,
        summarize_results,
    )

    results = [
        CombinationResult(1, "a", STATUS_SUCCESS),
        CombinationResult(2, "b", STATUS_WRITE_FAILED),
        CombinationResult(3, "c", STATUS_FAILED),
        CombinationResult(4, "d", STATUS_TIMED_OUT),
    ]
    stats = summarize_results(results)
    assert stats["total"] == 4
    assert stats["generated"] == 2
    assert stats["write_succeeded"] == 1
    assert stats["write_failed"] == 1
    assert stats["failed"] == 1
    assert stats["timed_out"] == 1
    # Failure, timeout and write failure each contribute an error code and
    # force an error termination (exit 1); no priority/ordering is applied.
    assert stats["exit_code"] == 1
    assert set(stats["error_codes"]) == {
        ERROR_CODE_EXEC,
        ERROR_CODE_TIMEOUT,
        ERROR_CODE_WRITE,
    }

    ok = summarize_results([CombinationResult(1, "a", STATUS_SUCCESS)])
    assert ok["exit_code"] == 0
    assert ok["error_codes"] == []

    # A declared input file failure yields ITR-INPUT-001 (exit 1), even when
    # some structures were generated.
    with_input_err = summarize_results(
        [CombinationResult(1, "a", STATUS_SUCCESS)], input_error_count=1
    )
    assert with_input_err["error_codes"] == [ERROR_CODE_INPUT]
    assert with_input_err["exit_code"] == 1

    # Some succeed, some fail: there is no PARTIAL SUCCESS exit-0 special case
    # anymore -> the run is an error termination (exit 1) carrying
    # ITR-EXEC-001, while the successful structure is still delivered.
    partial = summarize_results(
        [
            CombinationResult(1, "a", STATUS_SUCCESS),
            CombinationResult(2, "b", STATUS_FAILED),
        ]
    )
    assert partial["error_codes"] == [ERROR_CODE_EXEC]
    assert partial["exit_code"] == 1


def test_write_report_atomically(tmpdir):
    import glob

    from chemsmart.jobs.iterate.report import write_report_atomically

    path = str(tmpdir / "r_iterate.out")
    write_report_atomically(path, "hello")
    with open(path) as f:
        assert f.read() == "hello"
    # Overwriting a pre-existing report is allowed.
    write_report_atomically(path, "world")
    with open(path) as f:
        assert f.read() == "world"
    # No half-written temp file is left behind.
    assert glob.glob(str(tmpdir / "*.tmp*")) == []


def test_write_report_uses_config_stem_merged(tmpdir, fake_iterate_jobrunner):
    from types import SimpleNamespace

    runner = fake_iterate_jobrunner
    settings = SimpleNamespace(config_file=str(tmpdir / "screening.v2.yaml"))
    job = SimpleNamespace(
        settings=settings,
        separate_outputs=False,
        outputfile=str(tmpdir / "out.xyz"),
        output_directory=None,
    )
    report = _minimal_report(output_mode="merged", output_location="/out.xyz")
    path, err = runner._write_report(job, report)
    assert err is None
    assert os.path.basename(path) == "screening.v2_iterate.out"
    assert os.path.dirname(path) == str(tmpdir)
    assert os.path.exists(path)


def test_write_report_uses_directory_separate(tmpdir, fake_iterate_jobrunner):
    from types import SimpleNamespace

    runner = fake_iterate_jobrunner
    out_dir = str(tmpdir / "outdir")
    settings = SimpleNamespace(config_file=str(tmpdir / "123.yaml"))
    job = SimpleNamespace(
        settings=settings,
        separate_outputs=True,
        outputfile=str(tmpdir / "unused.xyz"),
        output_directory=out_dir,
    )
    report = _minimal_report(output_mode="separate", output_location=out_dir)
    path, err = runner._write_report(job, report)
    assert err is None
    assert os.path.basename(path) == "123_iterate.out"
    assert os.path.dirname(path) == out_dir
    assert os.path.exists(path)


def test_iterate_cli_config_error_exit_code_2(tmpdir):
    """An invalid algorithm option value is a usage error (Click exit 2)."""
    from click.testing import CliRunner

    from chemsmart.cli.iterate.yaml_cmd import yaml_cmd

    algo = "algorithm:\n  name: etkdg\n  options:\n    num_conformers: -1\n"
    config_file = _write_iterate_config(tmpdir, "bad_exit2.yaml", algo)
    result = CliRunner().invoke(yaml_cmd, ["-f", config_file], obj={})
    assert result.exit_code == 2


def test_iterate_cli_yaml_parse_error_exit_code_2(tmpdir):
    """A malformed YAML file is a usage error (Click exit 2), not a crash."""
    from click.testing import CliRunner

    from chemsmart.cli.iterate.yaml_cmd import yaml_cmd

    config_file = tmpdir / "broken.yaml"
    # Unclosed flow sequence -> yaml.YAMLError on safe_load.
    with open(config_file, "w") as f:
        f.write("skeletons: [unclosed\n")

    result = CliRunner().invoke(yaml_cmd, ["-f", str(config_file)], obj={})
    assert result.exit_code == 2
    assert "not valid YAML" in result.output


def test_iterate_cli_non_mapping_yaml_exit_code_2(tmpdir):
    """A top-level YAML list/scalar is a usage error (Click exit 2)."""
    from click.testing import CliRunner

    from chemsmart.cli.iterate.yaml_cmd import yaml_cmd

    config_file = tmpdir / "list_top.yaml"
    with open(config_file, "w") as f:
        f.write("- a\n- b\n")

    result = CliRunner().invoke(yaml_cmd, ["-f", str(config_file)], obj={})
    assert result.exit_code == 2
    assert "must contain a YAML mapping" in result.output


def test_run_internal_error_writes_report(
    tmpdir, iterate_jobrunner, monkeypatch
):
    """An unexpected runtime failure yields a best-effort INTERNAL ERROR
    report and exit code 1."""
    from chemsmart.jobs.iterate.job import IterateJob
    from chemsmart.jobs.iterate.settings import IterateJobSettings

    settings = IterateJobSettings(config_file=str(tmpdir / "cfg.yaml"))
    settings.skeleton_list = []
    settings.substituent_list = []
    job = IterateJob(
        settings=settings,
        jobrunner=iterate_jobrunner,
        outputfile=str(tmpdir / "out"),
    )

    def _boom(_job):
        raise RuntimeError("kaboom")

    monkeypatch.setattr(iterate_jobrunner, "_generate_combinations", _boom)

    summary = iterate_jobrunner.run(job)
    assert summary.exit_code == 1
    assert summary.error_codes == ["ITR-INTERNAL-001"]
    assert summary.summary_path is not None
    with open(summary.summary_path) as handle:
        text = handle.read()
    assert "ITR-INTERNAL-001" in text
    assert "kaboom" in text
    assert "Error termination of ChemSmart Iterate" in text


def test_run_summary_write_failure_exit_1(
    tmpdir, iterate_jobrunner, monkeypatch
):
    """A summary-write failure is surfaced and forces exit code 1."""
    from chemsmart.jobs.iterate.job import IterateJob
    from chemsmart.jobs.iterate.settings import IterateJobSettings

    settings = IterateJobSettings(config_file=str(tmpdir / "cfg.yaml"))
    settings.skeleton_list = []
    settings.substituent_list = []
    job = IterateJob(
        settings=settings,
        jobrunner=iterate_jobrunner,
        outputfile=str(tmpdir / "out"),
    )

    def _boom_write(path, text):
        raise OSError("disk full")

    monkeypatch.setattr(
        "chemsmart.jobs.iterate.runner.write_report_atomically", _boom_write
    )

    summary = iterate_jobrunner.run(job)
    assert summary.exit_code == 1
    assert summary.summary_path is None
    assert summary.summary_write_error is not None
    assert "disk full" in summary.summary_write_error


def test_run_interrupt_writes_report_and_reraises(
    tmpdir, iterate_jobrunner, monkeypatch
):
    """SIGINT writes a best-effort INTERRUPTED report and re-raises."""
    from chemsmart.jobs.iterate.job import IterateJob
    from chemsmart.jobs.iterate.settings import IterateJobSettings

    settings = IterateJobSettings(config_file=str(tmpdir / "cfg.yaml"))
    settings.skeleton_list = []
    settings.substituent_list = []
    job = IterateJob(
        settings=settings,
        jobrunner=iterate_jobrunner,
        outputfile=str(tmpdir / "out"),
    )

    def _interrupt(_job):
        raise KeyboardInterrupt()

    monkeypatch.setattr(
        iterate_jobrunner, "_generate_combinations", _interrupt
    )

    with pytest.raises(KeyboardInterrupt):
        iterate_jobrunner.run(job)

    report_path = str(tmpdir / "cfg_iterate.out")
    assert os.path.exists(report_path)
    with open(report_path) as handle:
        text = handle.read()
    assert "ITR-INTERRUPTED-001" in text
    assert "Error termination of ChemSmart Iterate" in text


def test_iterate_cli_sigint_exit_code_130(tmpdir, monkeypatch):
    """A user interrupt during the run makes the CLI exit with code 130."""
    from click.testing import CliRunner

    from chemsmart.cli.iterate.yaml_cmd import yaml_cmd

    config_file = _write_iterate_config(tmpdir, "sigint.yaml")

    def _raise_interrupt(self, progress_callback=None):
        raise KeyboardInterrupt()

    monkeypatch.setattr(
        "chemsmart.jobs.iterate.job.IterateJob.run", _raise_interrupt
    )
    result = CliRunner().invoke(yaml_cmd, ["-f", config_file], obj={})
    assert result.exit_code == 130
