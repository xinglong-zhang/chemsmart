import os

import numpy as np
import yaml

from chemsmart.cli.iterate.iterate import validate_config
from chemsmart.jobs.iterate.job import IterateJob
from chemsmart.jobs.iterate.runner import IterateJobRunner
from chemsmart.jobs.iterate.settings import IterateJobSettings


def test_iterate_workflow(
    iterate_config_file,
    iterate_input_directory,
    iterate_expected_output_file,
    tmp_path,
):
    """
    Test the full Iterate workflow:
    1. Load config from YAML
    2. Run IterateJob
    3. Compare output with expected XYZ file
    """
    # Change CWD to input directory so relative paths in YAML work
    original_cwd = os.getcwd()
    os.chdir(iterate_input_directory)

    try:
        # 1. Load and validate configuration
        with open(iterate_config_file, "r") as f:
            raw_config = yaml.safe_load(f)

        config = validate_config(raw_config, iterate_config_file)

        # 2. Setup Job Settings
        job_settings = IterateJobSettings(
            config_file=iterate_config_file, method="lagrange_multipliers"
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
