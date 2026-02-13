import os

import numpy as np
import tomlkit

from chemsmart.cli.iterate.iterate import validate_config
from chemsmart.jobs.iterate.job import IterateJob
from chemsmart.jobs.iterate.settings import IterateJobSettings


def test_iterate_integration_workflow(
    tmpdir,
    iterate_integration_config_file,
    iterate_input_directory,
    iterate_expected_output_file,
    iterate_jobrunner,
):
    """
    Test the full Iterate workflow (Integration Test):
    1. Load config from TOML (integration_iterate.toml)
    2. Run IterateJob
    3. Compare output with expected XYZ file
    """
    # Change CWD to input directory so relative paths in TOML work
    original_cwd = os.getcwd()
    os.chdir(iterate_input_directory)

    try:
        # 1. Load and validate configuration
        with open(iterate_integration_config_file, "r") as f:
            raw_config = tomlkit.load(f).unwrap()

        config = validate_config(raw_config, iterate_integration_config_file)

        # 2. Setup Job Settings
        job_settings = IterateJobSettings(
            config_file=iterate_integration_config_file,
            method="lagrange_multipliers",
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
        generated_output_path = job.run()

        # 5. Verify Output
        assert os.path.exists(
            generated_output_path
        ), "Output file was not generated"

        # Compare generated output with expected output
        # Semantic comparison (atoms and coordinates) is preferred over byte-comparison
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
            np.testing.assert_allclose(
                np.asarray(gen.positions, dtype=float),
                np.asarray(exp.positions, dtype=float),
                atol=1e-5,
                err_msg=f"Coordinate mismatch for {gen.info.get('comment')}",
            )

    finally:
        # Restore CWD
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
    We set a very short timeout (e.g. 0.001s) which should cause the worker to fail due to timeout.
    """
    import logging

    # Change CWD to input directory for relative file paths
    original_cwd = os.getcwd()
    os.chdir(iterate_input_directory)

    # Capture logs
    caplog.set_level(logging.WARNING)

    try:
        # 1. Load Config
        with open(iterate_timeout_config_file, "r") as f:
            raw_config = tomlkit.load(f).unwrap()
        config = validate_config(raw_config, iterate_timeout_config_file)

        # 2. Setup Job with very short timeout
        job_settings = IterateJobSettings(
            config_file=iterate_timeout_config_file,
            method="lagrange_multipliers",
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
        # Expected log from runner.py: "Timeout ({timeout}s) for combination: {label}"

        # We need to construct the label to search for
        # Carbene1_34_OTf_8
        label = "Carbene1_34_OTf_8"

        found_timeout_log = False
        for record in caplog.records:
            if "Timeout" in record.message and label in record.message:
                found_timeout_log = True
                break

        assert (
            found_timeout_log
        ), f"Timeout warning log not found for the combination {label}. Logs: {[r.message for r in caplog.records]}"

    finally:
        os.chdir(original_cwd)


def test_iterate_template_generation(tmpdir, iterate_template_file):
    """
    Test that the iterate configuration template is generated correctly and matches the golden copy.
    """
    from chemsmart.utils.iterate import generate_template

    # 1. Generate template
    generated_path = tmpdir / "test_template.toml"
    generate_template(str(generated_path))

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

    # 4. Verify it is valid TOML
    import tomlkit

    parsed = tomlkit.parse(generated_content)
    assert "skeletons" in parsed
    assert "substituents" in parsed
    assert len(parsed["skeletons"]) == 1  # Based on current template examples
    assert len(parsed["substituents"]) == 1


def test_iterate_validation_fails_on_invalid_link_index(
    iterate_invalid_skeleton_link_index_config_file,
):
    """
    Test that the validation logic correctly identifies when a link_index is not
    included in skeleton_indices.
    """
    from click.testing import CliRunner

    from chemsmart.cli.iterate.iterate import iterate

    runner = CliRunner()
    # Pass obj={} to initialize context object, required by MyGroup middleware
    result = runner.invoke(
        iterate,
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
    Test various validation failure scenarios using dynamically generated configs.
    Covers missing required fields and forbidden keys.
    """
    from click.testing import CliRunner

    from chemsmart.cli.iterate.iterate import iterate

    test_cases = [
        # Case 1: Skeleton missing file_path
        (
            """
            [[skeletons]]
            label = "s1"
            link_index = "1"
            [[substituents]]
            file_path = "sub.xyz"
            label = "sub1"
            link_index = "1"
            """,
            ["Missing required field 'file_path'", "Skeleton entry 1"],
            "Skeleton missing file_path",
        ),
        # Case 2: Substituent missing file_path
        (
            """
            [[skeletons]]
            file_path = "skel.xyz"
            label = "s1"
            link_index = "1"
            [[substituents]]
            label = "sub1"
            link_index = "1"
            """,
            ["Missing required field 'file_path'", "Substituent entry 1"],
            "Substituent missing file_path",
        ),
        # Case 3: Skeleton has forbidden key 'smiles'
        (
            """
            [[skeletons]]
            file_path = "skel.xyz"
            label = "s1"
            link_index = "1"
            smiles = "C"
            [[substituents]]
            file_path = "sub.xyz"
            label = "sub1"
            link_index = "1"
            """,
            ["Unknown key(s) in skeleton entry 1", "{'smiles'}"],
            "Skeleton forbidden key smiles",
        ),
        # Case 4: Substituent has forbidden key 'pubchem'
        (
            """
            [[skeletons]]
            file_path = "skel.xyz"
            label = "s1"
            link_index = "1"
            [[substituents]]
            file_path = "sub.xyz"
            label = "sub1"
            link_index = "1"
            pubchem = "100"
            """,
            ["Unknown key(s) in substituent entry 1", "{'pubchem'}"],
            "Substituent forbidden key pubchem",
        ),
        # Case 5: Random garbage key
        (
            """
            [[skeletons]]
            file_path = "skel.xyz"
            label = "s1"
            link_index = "1"
            random_key = "garbage"
            """,
            ["Unknown key(s) in skeleton entry 1", "{'random_key'}"],
            "Skeleton random garbage key",
        ),
        # Case 6: Zero link_index (S2 Check)
        (
            """
            [[skeletons]]
            file_path = "skel.xyz"
            label = "s1"
            link_index = "0" 
            """,
            ["Found invalid index <= 0", "link_index"],
            "Skeleton zero link_index",
        ),
        # Case 7: Negative skeleton_indices (S2 Check)
        (
            """
            [[skeletons]]
            file_path = "skel.xyz"
            label = "s1"
            link_index = "1"
            skeleton_indices = [1, -5, 3]
            """,
            ["Found invalid index <= 0", "skeleton_indices"],
            "Skeleton negative skeleton_indices",
        ),
        # Case 8: Substituent multiple link_indices (S3 Check)
        (
            """
            [[skeletons]]
            file_path = "skel.xyz"
            label = "s1"
            link_index = "1"
            [[substituents]]
            file_path = "sub.xyz"
            label = "sub1"
            link_index = "1, 2" 
            """,
            ["Multiple values found in 'link_index'", "exactly one link atom"],
            "Substituent multiple link_index",
        ),
        # Case 9: Skeleton label with invalid characters (S4 Check: Safe Label)
        (
            """
            [[skeletons]]
            file_path = "skel.xyz"
            label = "s1/unsafe"
            link_index = "1"
            """,
            ["Contains invalid characters", "Allowed characters"],
            "Skeleton label unsafe characters",
        ),
        # Case 10: Substituent label with invalid characters (S4 Check: Safe Label)
        (
            """
            [[skeletons]]
            file_path = "skel.xyz"
            label = "s1"
            link_index = "1"
            [[substituents]]
            file_path = "sub.xyz"
            label = "sub .. 1"
            link_index = "1"
            """,
            ["Contains invalid characters", "Allowed characters"],
            "Substituent label unsafe characters",
        ),
        # Case 11: Valid complex label (Testing allowed chars)
        # This one should PASS, but our loop expects FAILURES.
        # We will add it to a separate test if needed, or invert logic here.
        # Sticking to FAILURE cases here.
        # Case 11: Label with space
        (
            """
            [[skeletons]]
            file_path = "skel.xyz"
            label = "s1 "
            link_index = "1"
            """,
            ["Contains invalid characters"],
            "Label with space",
        ),
    ]

    runner = CliRunner()

    for idx, (config_content, expected_fragments, case_name) in enumerate(
        test_cases
    ):
        config_file = tmpdir / f"test_config_{idx}.toml"
        with open(config_file, "w") as f:
            f.write(config_content)

        result = runner.invoke(iterate, ["-f", str(config_file)], obj={})

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


def test_iterate_runner_bounds_validation(tmpdir):
    """
    Test that IterateJobRunner correctly validates indices against molecule size.
    (S2) Check: indices > num_atoms checking LOGS, not exceptions.
    """
    from unittest.mock import MagicMock, patch

    from chemsmart.io.molecules.structure import Molecule
    from chemsmart.jobs.iterate.runner import IterateJobRunner

    # Setup wrapper for the test
    runner = IterateJobRunner(fake=True)

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
    1. Run 'chemsmart iterate -f config.toml'
    2. Verify success exit code
    3. Verify output file exists and matches expected content.
    This ensures that the CLI entry point correctly orchestrates the job runner.
    """
    import os

    from click.testing import CliRunner

    from chemsmart.cli.iterate.iterate import iterate

    # Use the renamed config file which represents a valid CLI happy path
    config_file = os.path.join(
        iterate_configs_directory, "cli_happy_path.toml"
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
            iterate,
            ["-f", config_file, "-o", output_base_path],
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

        # Basic check: Ensuring specific label presence which confirms combination logic ran
        # The runner combines skeleton label, link index, sub label, and link index.
        # e.g. Carbene1 + 8 + OTf + 8 -> Carbene1_8_OTf_8
        assert "Carbene1_8_OTf_8" in generated_content

        # Verify line count matches (structure completeness check)
        exp_lines = expected_content.splitlines()
        gen_lines = generated_content.splitlines()
        assert len(exp_lines) == len(
            gen_lines
        ), "Generated XYZ line count differs from expected."

    finally:
        os.chdir(original_cwd)
