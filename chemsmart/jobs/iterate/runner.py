import logging
import multiprocessing
import os
import time
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
    method: str = "lagrange_multipliers"
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


def _run_combination_task(
    combination: IterateCombination,
) -> tuple[str, Optional[Molecule]]:
    """
    Core logic to run a single combination.

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
            method=combination.method,
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
        logger.error(f"Error in task for {label}: {e}")
        return (label, None)


def _run_combination_worker(
    combination: IterateCombination,
    result_queue: "multiprocessing.Queue",
) -> None:
    """
    Worker function for multiprocessing.Process.

    Parameters
    ----------
    combination : IterateCombination
        The combination to process
    result_queue : multiprocessing.Queue
        Queue to put the result tuple
    """
    try:
        result_pair = _run_combination_task(combination)
        result_queue.put(result_pair)
    except Exception as e:
        logger.error(f"Worker process panic for {combination.label}: {e}")
        result_queue.put((combination.label, None))


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

        return _run_combination_task(combination)

    def run_combinations(
        self,
        combinations: list[IterateCombination],
        nprocs: int = 1,
        timeout: float = DEFAULT_WORKER_TIMEOUT,
    ) -> list[tuple[str, Optional[Molecule]]]:
        """
        Run multiple combinations using multiprocessing.Process with watchdog.

        This implementation manually manages processes to ensure they can be
        forcefully killed (terminate/kill) if they exceed the timeout.
        ProcessPoolExecutor does not support killing individual tasks/processes
        mid-execution easily.

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
            f"Running {len(combinations)} combination(s) with {nprocs} process(es) (manual management), timeout={timeout}s"
        )

        results_dict: dict[str, Optional[Molecule]] = {}
        failed_labels: list[str] = []
        timed_out_labels: list[str] = []

        max_workers = 1 if nprocs == 1 else nprocs

        # Use Manager Queue for IPC
        manager = multiprocessing.Manager()
        result_queue = manager.Queue()

        try:
            # State tracking
            pending_combinations = combinations.copy()
            # process_id -> (process, combination, start_time)
            # Use id(p) instead of p.pid as key to avoid potential None pid issues
            active_processes: dict[
                int, tuple[multiprocessing.Process, IterateCombination, float]
            ] = {}

            while pending_combinations or active_processes:
                # 1. Fill up empty slots
                while (
                    pending_combinations
                    and len(active_processes) < max_workers
                ):
                    comb = pending_combinations.pop(0)
                    p = multiprocessing.Process(
                        target=_run_combination_worker,
                        args=(comb, result_queue),
                        daemon=True,
                    )
                    p.start()
                    active_processes[id(p)] = (p, comb, time.time())

                # 2. Check for results in queue
                while not result_queue.empty():
                    try:
                        # Grab all available data
                        lbl, mol = result_queue.get_nowait()
                        results_dict[lbl] = mol
                    except Exception:  # Empty or other error
                        break

                # 3. Monitor running processes
                current_time = time.time()
                pids_to_remove = []

                for proc_id, (p, comb, start_time) in active_processes.items():
                    if not p.is_alive():
                        # Process finished naturally (or crashed).
                        # Result should be in queue (handled above or next loop before exit)
                        p.join()
                        pids_to_remove.append(proc_id)
                    else:
                        # Check timeout
                        if (current_time - start_time) > timeout:
                            logger.warning(
                                f"Timeout ({timeout}s) for {comb.label} (pid {p.pid}). Terminating..."
                            )
                            p.terminate()
                            # Give it a tiny bit to terminate gracefully-ish
                            p.join(timeout=0.5)
                            if p.is_alive():
                                logger.warning(
                                    f"Process {p.pid} stuck, killing..."
                                )
                                p.kill()  # SIGKILL
                                # Must join after kill to reap zombie process
                                p.join(timeout=1.0)

                            timed_out_labels.append(comb.label)
                            results_dict[comb.label] = (
                                None  # Explicitly mark as failed/timeout
                            )
                            pids_to_remove.append(proc_id)

                # 4. Cleanup removed processes
                for proc_id in pids_to_remove:
                    del active_processes[proc_id]

                # 5. Sleep briefly to yield CPU
                if active_processes:
                    time.sleep(0.1)

            # Double check queue one last time just in case
            while not result_queue.empty():
                try:
                    lbl, mol = result_queue.get_nowait()
                    results_dict[lbl] = mol
                except Exception:
                    break

        finally:
            # Properly shutdown the Manager to avoid resource leaks
            manager.shutdown()

        # Check for missing results (crashes that didn't write to queue)
        for comb in combinations:
            if comb.label not in results_dict:
                logger.error(
                    f"No result found for {comb.label} - assuming worker crash."
                )
                failed_labels.append(comb.label)
                results_dict[comb.label] = None
            elif (
                results_dict[comb.label] is None
                and comb.label not in timed_out_labels
            ):
                # It might have returned None explicitly (handled in worker)
                if comb.label not in failed_labels:
                    failed_labels.append(comb.label)

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
    ) -> tuple[Optional[Molecule], str]:
        """
        Load a Molecule from configuration dict.

        Parameters
        ----------
        mol_config : dict
            Configuration containing file_path
        mol_type : str
            "skeleton" or "substituent" for logging
        idx : int
            Index for logging

        Returns
        -------
        tuple[Molecule | None, str]
            (Loaded molecule or None if loading fails, label)
        """
        label = mol_config.get("label") or f"{mol_type}{idx + 1}"
        file_path = mol_config.get("file_path")

        try:
            if file_path:
                molecule = Molecule.from_filepath(file_path)
                logger.debug(
                    f"Loaded {mol_type} '{label}' from file: {file_path}"
                )
            else:
                logger.warning(
                    f"{mol_type.capitalize()} '{label}' has no valid source "
                    f"(file_path), skipping."
                )
                return None, label

            # S2 Check: Validate indices against atom count
            num_atoms = molecule.num_atoms

            # Validate link_index
            link_indices = mol_config.get("link_index")
            # Ensure it is a list if not None (normalized in CLI, but for safety in runner)
            if link_indices is not None and not isinstance(link_indices, list):
                if isinstance(link_indices, int):
                    link_indices = [link_indices]
                # If string or other, it should have been caught by CLI, but we assume it might be raw here if skipped CLI

            if link_indices:
                invalid_links = [i for i in link_indices if i > num_atoms]
                if invalid_links:
                    logger.error(
                        f"{mol_type.capitalize()} '{label}': link_index {invalid_links} "
                        f"out of bounds. Molecule has {num_atoms} atoms."
                    )
                    return None, label

            # Validate skeleton_indices (only for skeletons)
            if mol_type == "skeleton":
                skel_indices = mol_config.get("skeleton_indices")
                # Ensure list
                if skel_indices is not None and not isinstance(
                    skel_indices, list
                ):
                    if isinstance(skel_indices, int):
                        skel_indices = [skel_indices]

                if skel_indices:
                    invalid_skels = [i for i in skel_indices if i > num_atoms]
                    if invalid_skels:
                        logger.error(
                            f"{mol_type.capitalize()} '{label}': skeleton_indices {invalid_skels} "
                            f"out of bounds. Molecule has {num_atoms} atoms."
                        )
                        return None, label

            return molecule, label
        except Exception as e:
            logger.error(f"Failed to load {mol_type} '{label}': {e}")
            return None, label

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
        method = job.settings.method
        sphere_direction_samples_num = (
            job.settings.sphere_direction_samples_num
        )
        axial_rotations_sample_num = job.settings.axial_rotations_sample_num

        valid_skeletons = []
        for skel_idx, skel_config in enumerate(skeleton_list):
            skeleton, skel_label = self._load_molecule(
                skel_config, "skeleton", skel_idx
            )
            if skeleton is None:
                continue

            # skel_label is returned by _load_molecule
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

            valid_skeletons.append((skeleton, skel_label, skel_config))

        valid_substituents = []
        for sub_idx, sub_config in enumerate(substituent_list):
            substituent, sub_label = self._load_molecule(
                sub_config, "substituent", sub_idx
            )
            if substituent is None:
                continue

            # sub_label is returned by _load_molecule
            sub_link_index = sub_config.get("link_index")  # list[int], 1-based

            if not sub_link_index:
                logger.warning(
                    f"Substituent '{sub_label}' has no link_index, skipping."
                )
                continue

            valid_substituents.append((substituent, sub_label, sub_config))

        for skeleton, skel_label, skel_config in valid_skeletons:
            # skel_label is returned by _load_molecule
            skel_link_indices = skel_config.get(
                "link_index"
            )  # list[int], 1-based
            skeleton_indices = skel_config.get(
                "skeleton_indices"
            )  # list[int] or None

            for substituent, sub_label, sub_config in valid_substituents:
                # sub_label is returned by _load_molecule
                sub_link_index = sub_config.get(
                    "link_index"
                )  # list[int], 1-based

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
                        method=method,
                    )
                    combinations.append(combination)
                    logger.info(f"Created combination: {combination.label}")

        return combinations

    def _write_outputs(
        self, results: list[tuple[str, Optional[Molecule]]], job: "IterateJob"
    ) -> None:
        """
        Write execution results to output file(s).

        Parameters
        ----------
        results : list[tuple[str, Molecule | None]]
            List of (label, molecule) tuples from run_combinations
        job : IterateJob
             The job instance with configuration options
        """
        successful_count = 0

        if job.separate_outputs:
            # Separate files mode
            output_dir = job.output_directory
            if output_dir:
                os.makedirs(output_dir, exist_ok=True)
            else:
                # Should have been validated by CLI, but fallback
                output_dir = "."

            logger.info(
                f"Writing separate output files to directory: {output_dir}"
            )

            for label, mol in results:
                if mol is not None:
                    # Construct filename: directory + label + .xyz
                    filename = os.path.join(output_dir, f"{label}.xyz")
                    try:
                        with open(filename, "w") as f:
                            f.write(f"{mol.num_atoms}\n")
                            f.write(f"       {label}\n")
                            for symbol, pos in zip(
                                mol.chemical_symbols, mol.positions
                            ):
                                f.write(
                                    f"{symbol:2s}  {pos[0]:15.10f}  {pos[1]:15.10f}  {pos[2]:15.10f}\n"
                                )
                        successful_count += 1
                        logger.debug(f"Wrote {filename}")
                    except Exception as e:
                        logger.error(f"Failed to write {filename}: {e}")

            logger.info(
                f"Wrote {successful_count} separate molecule files to {output_dir}"
            )

        else:
            # Single merged output file mode
            outputfile = job.outputfile
            # Ensure output directory exists
            output_dir = os.path.dirname(outputfile)
            if output_dir:
                os.makedirs(output_dir, exist_ok=True)

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

            logger.info(
                f"Wrote {successful_count} molecule(s) to {outputfile}"
            )

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

        # Write output results
        self._write_outputs(results, job)

        logger.info("IterateJobRunner completed job")
