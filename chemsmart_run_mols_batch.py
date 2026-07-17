#!/usr/bin/env python
import os

os.environ["OMP_NUM_THREADS"] = "1"

from chemsmart.cli.run import run


def run_job():
    run(
        [
            "--server",
            "PBS",
            "--num-nodes",
            "1",
            "--no-run-in-parallel",
            "--no-fake",
            "--no-delete-scratch",
            "--no-debug",
            "--stream",
            "orca",
            "--project",
            "gas_solv",
            "--filename",
            "/Users/taipanlan/PycharmProjects/chemsmart/tests/data/StructuresTests/xyz/crest_conformers.xyz",
            "--index",
            "1,2",
            "--charge",
            "0",
            "--multiplicity",
            "1",
            "--defgrid",
            "defgrid2",
            "--no-forces",
            "--no-remove-solvent",
            "opt",
            "--skip-completed",
            "--no-remove-solvent",
            "--no-invert-constraints",
        ]
    )


if __name__ == "__main__":
    run_job()
