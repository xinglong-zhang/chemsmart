import logging
import os
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import TimeoutError as FuturesTimeoutError
from concurrent.futures import as_completed
from dataclasses import dataclass
from typing import TYPE_CHECKING, Optional

from chemsmart.analysis.iterate import (
    IterateAnalyzer,
    SkeletonPreprocessor,
    SubstituentPreprocessor,
)
from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.runner import JobRunner

if TYPE_CHECKING:
    from chemsmart.jobs.iterate.job import IterateJob

logger = logging.getLogger(__name__)

# Default timeout for each combination worker (in seconds)
DEFAULT_WORKER_TIMEOUT = 120  # 2 minutes


@dataclass
class IterateCombination:
    """
    Represents a single combination of skeleton, link_index, and substituent.
    """

    skeleton: Molecule
    skeleton_label: str
    skeleton_link_index: int  # 1-based
    skeleton_indices: Optional[list[int]]  # 1-based, or None
    substituent: Molecule
    substituent_label: str
    substituent_link_index: int  # 1-based
    algorithm: str = "lagrange_multipliers"
    sphere_direction_samples_num: int = 96
    axial_rotations_sample_num: int = 6

    @property
    def label(self) -> str:
        """
        Generate a unique label for this combination.

        Format: {skeleton_label}_{skeleton_link_index}_{substituent_label}_{substituent_link_index}
        Example: benzene_5_methyl_1
        """
        return f"{self.skeleton_label}_{self.skeleton_link_index}_{self.substituent_label}_{self.substituent_link_index}"


def _run_combination_worker(
    combination: IterateCombination,
) -> tuple[str, Optional[Molecule]]:
    """
    Worker function for multiprocessing Pool.

    This function is module-level to be picklable for multiprocessing.

    Parameters
    ----------
    combination : IterateCombination
        The combination to process

    Returns
    -------
    tuple[str, Molecule | None]
        (label, result_molecule) or (label, None) if failed
    """
    label = combination.label

    try:
        # Preprocess skeleton if needed
        skeleton_preprocessor = SkeletonPreprocessor(
            molecule=combination.skeleton,
            link_index=combination.skeleton_link_index,
            skeleton_indices=combination.skeleton_indices,
        )
        processed_skeleton = skeleton_preprocessor.run()
        new_skeleton_link_index = skeleton_preprocessor.get_new_link_index()

        # Preprocess substituent if needed
        substituent_preprocessor = SubstituentPreprocessor(
            molecule=combination.substituent,
            link_index=combination.substituent_link_index,
        )
        processed_substituent = substituent_preprocessor.run()
        new_substituent_link_index = (
            substituent_preprocessor.get_new_link_index()
        )

        # Run analysis to generate combined molecule
        analyzer = IterateAnalyzer(
            skeleton=processed_skeleton,
            substituent=processed_substituent,
            skeleton_link_index=new_skeleton_link_index,
            substituent_link_index=new_substituent_link_index,
            algorithm=combination.algorithm,
            sphere_direction_samples_num=combination.sphere_direction_samples_num,
            axial_rotations_sample_num=combination.axial_rotations_sample_num,
        )

        result = analyzer.run()

        if result is not None:
            logger.info(f"Generated molecule for {label}")
        else:
            logger.warning(f"Failed to generate molecule for {label}")

        return (label, result)

    except Exception as e:
        logger.error(f"Error in worker for {label}: {e}")
        return (label, None)


class IterateJobRunner(JobRunner):
    """
    Job runner for Iterate jobs.

    Iterate jobs are special in that they don't call external programs.
    They run purely in Python to generate molecular structures.

    This runner handles the execution of IterateCombination tasks,
    including multiprocessing support via run_combinations().
    """

    JOBTYPES = ["iterate"]
    PROGRAM = "Iterate"
    FAKE = False
    SCRATCH = False

    def __init__(
        self, server=None, scratch=None, fake=False, scratch_dir=None, **kwargs
    ):
        if scratch is None:
            scratch = self.SCRATCH
        super().__init__(
            server=server,
            scratch=scratch,
            scratch_dir=scratch_dir,
            fake=fake,
            **kwargs,
        )
        logger.debug("IterateJobRunner initialized")
        logger.debug(f"Jobrunner server: {self.server}")
        logger.debug(f"Jobrunner scratch: {self.scratch}")
        logger.debug(f"Jobrunner fake mode: {self.fake}")

    @property
    def executable(self):
        """Iterate jobs don't use an external executable."""
        return None

    def _get_command(self, job):
        """Iterate jobs don't need a command - they run in Python."""
        return None

    def run_single(
        self, combination: IterateCombination
    ) -> tuple[str, Optional[Molecule]]:
        """
        Run a single combination to generate a combined molecule.

        Parameters
        ----------
        combination : IterateCombination
            The combination to process

        Returns
        -------
        tuple[str, Molecule | None]
            (label, result_molecule) or (label, None) if failed
        """
        if self.fake:
            logger.info(
                f"Fake mode enabled, not actually running {combination.label}."
            )
            return (combination.label, None)

        return _run_combination_worker(combination)

    def run_combinations(
        self,
        combinations: list[IterateCombination],
        nprocs: int = 1,
        timeout: float = DEFAULT_WORKER_TIMEOUT,
    ) -> list[tuple[str, Optional[Molecule]]]:
        """
        Run multiple combinations, optionally using multiprocessing.

        Each worker has a timeout limit. Timeout or error in one worker
        does not affect other workers.

        Parameters
        ----------
        combinations : list[IterateCombination]
            List of combinations to process
        nprocs : int
            Number of processes for parallel execution. Default 1.
        timeout : float
            Timeout in seconds for each worker. Default 120 (2 minutes).

        Returns
        -------
        list[tuple[str, Molecule | None]]
            List of (label, result_molecule) tuples
        """
        if self.fake:
            logger.info(
                "Fake mode enabled, not actually running combinations."
            )
            return [(c.label, None) for c in combinations]

        if not combinations:
            logger.warning("No combinations to process.")
            return []

        logger.info(
            f"Running {len(combinations)} combination(s) with {nprocs} process(es), timeout={timeout}s"
        )

        # Track results by label to maintain order
        results_dict: dict[str, Optional[Molecule]] = {}
        failed_labels: list[str] = []
        timed_out_labels: list[str] = []

        max_workers = 1 if nprocs == 1 else nprocs

        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            # Submit all tasks
            future_to_comb = {
                executor.submit(_run_combination_worker, comb): comb
                for comb in combinations
            }

            # Use as_completed to process results as they finish
            # This allows other tasks to continue even if one times out
            for future in as_completed(future_to_comb, timeout=None):
                comb = future_to_comb[future]
                label = comb.label

                try:
                    # Get result with timeout for this specific future
                    result = future.result(timeout=timeout)
                    results_dict[label] = result[
                        1
                    ]  # result is (label, molecule)
                    if result[1] is not None:
                        logger.info(f"Completed: {label}")
                    else:
                        logger.warning(
                            f"Failed to generate molecule for: {label}"
                        )
                        failed_labels.append(label)
                except FuturesTimeoutError:
                    logger.warning(
                        f"Timeout ({timeout}s) for combination: {label}"
                    )
                    results_dict[label] = None
                    timed_out_labels.append(label)
                except Exception as e:
                    logger.warning(f"Error for combination {label}: {e}")
                    results_dict[label] = None
                    failed_labels.append(label)

        # Build results list in original order
        results = [
            (comb.label, results_dict.get(comb.label)) for comb in combinations
        ]

        # Log results summary
        successful_labels = [
            label for label, mol in results if mol is not None
        ]
        logger.info(
            f"Completed: {len(successful_labels)}/{len(results)} molecules generated successfully"
        )

        # Print summary of results
        if successful_labels or timed_out_labels or failed_labels:
            logger.info("=" * 40)
            logger.info("       SUMMARY OF RESULTS")
            logger.info("=" * 40)

            if successful_labels:
                logger.info(f"Successful ({len(successful_labels)}):")
                for label in successful_labels:
                    logger.info(f"  - {label}")

            if timed_out_labels:
                logger.warning(f"Timed out ({len(timed_out_labels)}):")
                for label in timed_out_labels:
                    logger.warning(f"  - {label}")

            if failed_labels:
                logger.warning(
                    f"Failed to find solution ({len(failed_labels)}):"
                )
                for label in failed_labels:
                    logger.warning(f"  - {label}")

            logger.info("=" * 40)

        return results

    def _load_molecule(
        self, mol_config: dict, mol_type: str, idx: int
    ) -> Optional[Molecule]:
        """
        Load a Molecule from configuration dict.

        Parameters
        ----------
        mol_config : dict
            Configuration containing file_path, smiles, or pubchem
        mol_type : str
            "skeleton" or "substituent" for logging
        idx : int
            Index for logging

        Returns
        -------
        Molecule or None
            Loaded molecule or None if loading fails
        """
        label = mol_config.get("label") or f"{mol_type}_{idx + 1}"
        file_path = mol_config.get("file_path")
        smiles = mol_config.get("smiles")
        pubchem = mol_config.get("pubchem")

        try:
            if file_path:
                molecule = Molecule.from_filepath(file_path)
                logger.debug(
                    f"Loaded {mol_type} '{label}' from file: {file_path}"
                )
            elif smiles:
                molecule = Molecule.from_smiles(smiles)
                logger.debug(
                    f"Loaded {mol_type} '{label}' from SMILES: {smiles}"
                )
            elif pubchem:
                molecule = Molecule.from_pubchem(pubchem)
                logger.debug(
                    f"Loaded {mol_type} '{label}' from PubChem: {pubchem}"
                )
            else:
                logger.warning(
                    f"{mol_type.capitalize()} '{label}' has no valid source "
                    f"(file_path, smiles, or pubchem), skipping."
                )
                return None
            return molecule
        except Exception as e:
            logger.error(f"Failed to load {mol_type} '{label}': {e}")
            return None

    def _generate_combinations(
        self, job: "IterateJob"
    ) -> list[IterateCombination]:
        """
        Generate all combinations of (skeleton, skeleton_link_index, substituent).

        Parameters
        ----------
        job : IterateJob
            The job containing settings with skeleton_list and substituent_list

        Returns
        -------
        list[IterateCombination]
            List of all valid combinations
        """
        combinations = []

        skeleton_list = job.settings.skeleton_list or []
        substituent_list = job.settings.substituent_list or []
        algorithm = job.settings.algorithm
        sphere_direction_samples_num = (
            job.settings.sphere_direction_samples_num
        )
        axial_rotations_sample_num = job.settings.axial_rotations_sample_num

        for skel_idx, skel_config in enumerate(skeleton_list):
            skeleton = self._load_molecule(skel_config, "skeleton", skel_idx)
            if skeleton is None:
                continue

            skel_label = skel_config.get("label") or f"skeleton_{skel_idx + 1}"
            skel_link_indices = skel_config.get(
                "link_index"
            )  # list[int], 1-based
            skeleton_indices = skel_config.get(
                "skeleton_indices"
            )  # list[int] or None

            if not skel_link_indices:
                logger.warning(
                    f"Skeleton '{skel_label}' has no link_index, skipping."
                )
                continue

            for sub_idx, sub_config in enumerate(substituent_list):
                substituent = self._load_molecule(
                    sub_config, "substituent", sub_idx
                )
                if substituent is None:
                    continue

                sub_label = (
                    sub_config.get("label") or f"substituent_{sub_idx + 1}"
                )
                sub_link_index = sub_config.get(
                    "link_index"
                )  # list[int], 1-based

                if not sub_link_index:
                    logger.warning(
                        f"Substituent '{sub_label}' has no link_index, skipping."
                    )
                    continue

                # Use the first substituent link index
                sub_link_idx = sub_link_index[0]

                # For each skeleton link position, create a combination
                for skel_link_idx in skel_link_indices:
                    combination = IterateCombination(
                        skeleton=skeleton.copy(),
                        skeleton_label=skel_label,
                        skeleton_link_index=skel_link_idx,
                        skeleton_indices=skeleton_indices,
                        substituent=substituent.copy(),
                        substituent_label=sub_label,
                        sphere_direction_samples_num=sphere_direction_samples_num,
                        axial_rotations_sample_num=axial_rotations_sample_num,
                        substituent_link_index=sub_link_idx,
                        algorithm=algorithm,
                    )
                    combinations.append(combination)
                    logger.info(f"Created combination: {combination.label}")

        return combinations

    def _write_output_xyz(
        self, results: list[tuple[str, Optional[Molecule]]], outputfile: str
    ) -> None:
        """
        Write all successful molecules to a single xyz file.

        Format:
        - First line: number of atoms
        - Second line: label of the combination
        - Following lines: atom coordinates

        Parameters
        ----------
        results : list[tuple[str, Molecule | None]]
            List of (label, molecule) tuples from run_combinations
        outputfile : str
            Path to the output xyz file
        """
        # Ensure output directory exists
        output_dir = os.path.dirname(outputfile)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)

        successful_count = 0
        with open(outputfile, "w") as f:
            for label, mol in results:
                if mol is not None:
                    # Build xyz string manually
                    # First line: number of atoms
                    f.write(f"{mol.num_atoms}\n")
                    # Second line: label as comment
                    f.write(f"       {label}\n")
                    # Following lines: atom symbol and coordinates
                    for symbol, pos in zip(
                        mol.chemical_symbols, mol.positions
                    ):
                        f.write(
                            f"{symbol:2s}  {pos[0]:15.10f}  {pos[1]:15.10f}  {pos[2]:15.10f}\n"
                        )
                    successful_count += 1

        logger.info(f"Wrote {successful_count} molecule(s) to {outputfile}")

    def run(self, job: "IterateJob", **kwargs) -> None:
        """
        Run an IterateJob.

        This method handles all execution logic:
        1. Load molecules from job settings
        2. Generate all combinations
        3. Run combinations with multiprocessing
        4. Write results to output xyz file

        Parameters
        ----------
        job : IterateJob
            The iterate job to run
        **kwargs
            Additional keyword arguments
        """
        logger.info(f"IterateJobRunner.run() called for job: {job}")

        if self.fake:
            logger.info("Fake mode enabled, not actually running job.")
            return

        # Generate combinations from job settings
        combinations = self._generate_combinations(job)
        logger.info(f"Generated {len(combinations)} combination(s)")

        if not combinations:
            logger.warning("No valid combinations to process.")
            return

        # Run combinations with multiprocessing
        results = self.run_combinations(
            combinations, nprocs=job.nprocs, timeout=job.timeout
        )

        # Write output xyz file
        self._write_output_xyz(results, job.outputfile)

        logger.info("IterateJobRunner completed job")
