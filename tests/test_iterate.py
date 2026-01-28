import os

import numpy as np
import pytest
import tomlkit

from chemsmart.cli.iterate.iterate import validate_config
from chemsmart.jobs.iterate.job import IterateJob
from chemsmart.jobs.iterate.runner import IterateJobRunner
from chemsmart.jobs.iterate.settings import IterateJobSettings


def test_iterate_integration_workflow(
    iterate_integration_config_file,
    iterate_input_directory,
    iterate_expected_output_file,
    tmp_path,
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
            config_file=iterate_integration_config_file, method="lagrange_multipliers"
        )
        job_settings.skeleton_list = config["skeletons"]
        job_settings.substituent_list = config["substituents"]

        # 3. Setup Job
        # Use a temporary file for output
        output_file = tmp_path / "test_output"

        jobrunner = IterateJobRunner()
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

        # Helper to parse multi-structure XYZ
        def parse_multi_xyz(filepath):
            structures = []
            with open(filepath, "r") as f:
                lines = f.readlines()

            i = 0
            while i < len(lines):
                try:
                    num_atoms = int(lines[i].strip())
                    label = lines[i + 1].strip()
                    atoms = []
                    coords = []
                    for j in range(num_atoms):
                        parts = lines[i + 2 + j].split()
                        atoms.append(parts[0])
                        coords.append([float(x) for x in parts[1:4]])

                    structures.append(
                        {
                            "label": label,
                            "atoms": atoms,
                            "coords": np.array(coords),
                        }
                    )
                    i += 2 + num_atoms
                except (ValueError, IndexError):
                    break
            return structures

        generated_structures = parse_multi_xyz(generated_output_path)
        expected_structures = parse_multi_xyz(iterate_expected_output_file)

        # Sort structures by label to ensure order-independent comparison
        # (multiprocessing or config order might vary generation sequence)
        generated_structures.sort(key=lambda x: x["label"])
        expected_structures.sort(key=lambda x: x["label"])

        assert len(generated_structures) == len(expected_structures), (
            f"Number of generated structures ({len(generated_structures)}) "
            f"does not match expected ({len(expected_structures)})"
        )

        for gen, exp in zip(generated_structures, expected_structures):
            # Compare labels
            assert (
                gen["label"] == exp["label"]
            ), f"Label mismatch: {gen['label']} != {exp['label']}"

            # Compare atom count
            assert len(gen["atoms"]) == len(
                exp["atoms"]
            ), f"Atom count mismatch for {gen['label']}: {len(gen['atoms'])} != {len(exp['atoms'])}"

            # Compare atom symbols (ensure composition is correct)
            assert (
                gen["atoms"] == exp["atoms"]
            ), f"Atom symbols mismatch for {gen['label']}"

            # Check coordinates
            np.testing.assert_allclose(
                gen["coords"],
                exp["coords"],
                atol=1e-5,
                err_msg=f"Coordinate mismatch for {gen['label']}",
            )

    finally:
        # Restore CWD
        os.chdir(original_cwd)


def test_iterate_timeout(
    iterate_timeout_config_file,
    iterate_input_directory,
    tmp_path,
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
            config_file=iterate_timeout_config_file, method="lagrange_multipliers"
        )
        job_settings.skeleton_list = config["skeletons"]
        job_settings.substituent_list = config["substituents"]

        output_file = tmp_path / "timeout_output"

        jobrunner = IterateJobRunner()
        job = IterateJob(
            settings=job_settings,
            jobrunner=jobrunner,
            nprocs=1, # Use 1 proc to ensure we hit it
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
        
        assert found_timeout_log, f"Timeout warning log not found for the combination {label}. Logs: {[r.message for r in caplog.records]}"

    finally:
        os.chdir(original_cwd)


def test_iterate_template_generation(tmp_path, iterate_template_file):
    """
    Test that the iterate configuration template is generated correctly and matches the golden copy.
    """
    from chemsmart.utils.iterate import generate_template

    # 1. Generate template
    generated_path = tmp_path / "test_template.toml"
    generate_template(str(generated_path))

    # 2. Assert file exists
    assert generated_path.exists()

    # 3. Compare content with expected template
    with open(generated_path, "r") as f:
        generated_content = f.read()

    with open(iterate_template_file, "r") as f:
        expected_content = f.read()

    # Normalize newlines and strip whitespace for robust comparison
    assert generated_content.strip() == expected_content.strip(), (
        "Generated template does not match expected template content."
    )

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
        iterate, ["-f", iterate_invalid_skeleton_link_index_config_file], obj={}
    )

    assert result.exit_code != 0
    # Click 8.x format for BadParameter: "Error: Invalid value for ..."
    # We check for key parts of the error message
    assert "Invalid value" in result.output
    assert (
        "The link_index [6] is not included in 'skeleton_indices'" in result.output
    )


def test_iterate_validation_failures_comprehensive(tmp_path):
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
    ]

    runner = CliRunner()

    for idx, (config_content, expected_fragments, case_name) in enumerate(test_cases):
        config_file = tmp_path / f"test_config_{idx}.toml"
        with open(config_file, "w") as f:
            f.write(config_content)

        result = runner.invoke(iterate, ["-f", str(config_file)], obj={})

        assert result.exit_code != 0, f"Case '{case_name}' failed to raise error"
        assert (
            "Invalid value" in result.output
        ), f"Case '{case_name}' missing generic invalid value message"
        for fragment in expected_fragments:
            assert fragment in result.output, (
                f"Case '{case_name}' failed. "
                f"Expected fragment '{fragment}' not found in output:\n{result.output}"
            )


