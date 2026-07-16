import hashlib
import logging
import multiprocessing
import os
import queue
import sys
import time
import uuid
from collections import deque
from dataclasses import dataclass, field
from datetime import datetime
from itertools import product
from pathlib import Path
from typing import TYPE_CHECKING, Optional

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.iterate.iterate import (
    SkeletonPreprocessor,
    SubstituentPreprocessor,
)
from chemsmart.jobs.iterate.report import (
    ERROR_CODE_INTERNAL,
    ERROR_CODE_INTERRUPTED,
    STAGE_ALGORITHM,
    STAGE_PREPROCESSING,
    STAGE_WRITE,
    STATUS_FAILED,
    STATUS_SUCCESS,
    STATUS_TIMED_OUT,
    STATUS_WRITE_FAILED,
    CombinationResult,
    IterateReport,
    summarize_results,
    write_report_atomically,
)
from chemsmart.jobs.iterate.settings import (
    IterateAlgorithmConfig,
    get_algorithm_spec,
    resolve_algorithm_config,
)
from chemsmart.jobs.runner import JobRunner

if TYPE_CHECKING:
    from chemsmart.jobs.iterate.job import IterateJob

logger = logging.getLogger(__name__)

# Default timeout for each combination worker (in seconds)
DEFAULT_WORKER_TIMEOUT = 120  # 2 minutes


def _silence_rdkit_warnings() -> None:
    """Suppress RDKit's benign C++ warnings in the calling (worker) process.

    ETKDG embedding and force-field typing emit noisy ``rdApp.warning``
    messages (e.g. "explicit Hs", "UFFTYPER ...") that do not affect the run
    outcome. Only warnings are disabled -- errors and Python exceptions are
    untouched and no computation parameters change.
    """
    try:
        from rdkit import RDLogger

        RDLogger.DisableLog("rdApp.warning")
    except Exception:
        pass


@dataclass
class IterateRunSummary:
    """Outcome of an iterate run.

    Attributes
    ----------
    total : int
        Number of combinations processed.
    succeeded : int
        Number of combinations that generated a structure (whether or not it
        was subsequently written).
    failed : int
        Number of combinations that failed (excluding timeouts).
    timed_out : int
        Number of combinations killed by the worker timeout.
    write_failed : int
        Number of generated structures that could not be written.
    structures_written : int
        Number of structures actually written to disk.
    input_error_count : int
        Number of declared input files that failed to load.
    error_codes : list[str]
        Gaussian-style error codes that apply to the run (empty on a normal
        termination). Several may be present at once.
    output_paths : list[str]
        Paths of the files actually written.
    summary_path : str or None
        Path of the run report, or ``None`` if it could not be written.
    summary_write_error : str or None
        Error message if the report failed to write.
    exit_code : int
        0 for a completely clean run, 1 for any error, 130 for a user
        interrupt.
    """

    total: int = 0
    succeeded: int = 0
    failed: int = 0
    timed_out: int = 0
    write_failed: int = 0
    structures_written: int = 0
    input_error_count: int = 0
    error_codes: list = field(default_factory=list)
    output_paths: list = field(default_factory=list)
    summary_path: Optional[str] = None
    summary_write_error: Optional[str] = None
    exit_code: int = 0


@dataclass
class IterateMoleculePool:
    """Shared pool of unique Molecule objects for iterate combinations.

    Avoids duplicating Molecule objects across combinations. Each combination
    references molecules by index; actual copies are made only at worker
    dispatch time.

    Attributes
    ----------
    skeletons : list[Molecule]
        Pool of unique skeleton molecules.
    substituents : list[Molecule]
        Pool of unique substituent molecules.
    """

    skeletons: list[Molecule] = field(default_factory=list)
    substituents: list[Molecule] = field(default_factory=list)


@dataclass
class IterateAssignment:
    """One substituent assignment at one skeleton position.

    Substituents are referenced by index into IterateMoleculePool.substituents.
    Actual Molecule copies are created only at worker dispatch time.

    Attributes
    ----------
    substituent_idx : int
        Index into IterateMoleculePool.substituents.
    substituent_label : str
        Human-readable label for the substituent (e.g. "methyl").
    substituent_link_index : int
        1-based atom index on the substituent that bonds to the skeleton.
    skeleton_link_index : int
        1-based atom index on the **original** skeleton where this
        substituent will be attached.
    """

    substituent_idx: int
    substituent_label: str
    substituent_link_index: int  # 1-based
    skeleton_link_index: int  # 1-based, on original skeleton


@dataclass
class IterateCombination:
    """A combination of one skeleton with one or more substituent assignments.

    Molecules are stored in a shared IterateMoleculePool and referenced
    by index. Actual Molecule copies are created only at worker dispatch
    time.

    Each combination has one or more assignments — attaching a substituent
    at a link position. Single-assignment combinations come from
    independent-mode expansion (one slot at a time); multi-assignment
    combinations come from global-mode Cartesian product.

    Attributes
    ----------
    skeleton_idx : int
        Index into IterateMoleculePool.skeletons.
    skeleton_label : str
        Human-readable label for the skeleton.
    skeleton_indices : list[int] or None
        1-based atom indices defining the skeleton core. None means all
        atoms are skeleton atoms.
    assignments : list[IterateAssignment]
        One or more assignments. len=1 for independent mode,
        >=1 for global mode.
    algorithm_config : IterateAlgorithmConfig
        Resolved algorithm configuration (name + options) used to build the
        per-assignment analyzer (Lagrange multipliers or ETKDG).
    """

    skeleton_idx: int
    skeleton_label: str
    skeleton_indices: Optional[list[int]]  # 1-based, or None
    assignments: list[IterateAssignment] = field(default_factory=list)
    algorithm_config: IterateAlgorithmConfig = field(
        default_factory=resolve_algorithm_config
    )

    @property
    def label(self) -> str:
        """Generate a unique label for this combination.

        Format: {skeleton}_{link1}{sub1}_{link2}{sub2}_...
        Examples:
            Single:  benzene_5methyl
            Multi:   benzene_5methyl_8ethyl
        """
        parts = [self.skeleton_label]
        for a in sorted(self.assignments, key=lambda x: x.skeleton_link_index):
            parts.append(f"{a.skeleton_link_index}{a.substituent_label}")
        return "_".join(parts)


def _batch_preprocess_skeleton(
    skeleton: Molecule,
    link_indices: list[int],
    skeleton_indices: Optional[list[int]],
) -> tuple[Molecule, dict[int, int]]:
    """Remove old substituent groups at one or more link positions in one pass.

    For each link position, checks whether there is already an available
    bonding position. If not, detects and removes the existing substituent
    group. All removals are batched so that atom index remapping happens
    only once.

    Parameters
    ----------
    skeleton : Molecule
        Original skeleton molecule.
    link_indices : list[int]
        All 1-based link indices to preprocess.
    skeleton_indices : list[int] or None
        1-based skeleton atom indices (for SkeletonPreprocessor).

    Returns
    -------
    tuple[Molecule, dict[int, int]]
        (processed skeleton, mapping from original 1-based index
        to new 1-based index).
    """
    all_remove_indices = set()  # 0-based

    for link_idx in link_indices:
        prep = SkeletonPreprocessor(skeleton, link_idx, skeleton_indices)
        if not prep._has_available_bonding_position():
            if skeleton_indices is not None:
                # Precise path: DFS + skeleton_indices
                branches = prep._find_non_skeleton_branches()
                for branch in branches:
                    all_remove_indices.update(branch)
            else:
                # Fallback: no skeleton_indices, use smallest-branch heuristic
                removed = prep.detect_substituent()  # 0-based indices
                all_remove_indices.update(removed)

    keep_indices = sorted(set(range(len(skeleton))) - all_remove_indices)

    if not keep_indices:
        raise ValueError(
            "All skeleton atoms would be removed during preprocessing."
        )

    # Verify all link atoms are kept
    for link_idx in link_indices:
        if (link_idx - 1) not in keep_indices:
            raise ValueError(
                f"Link atom at index {link_idx} (1-based) was removed during "
                f"batch preprocessing."
            )

    # Build processed skeleton
    symbols = [skeleton.chemical_symbols[i] for i in keep_indices]
    positions = skeleton.positions[keep_indices]
    frozen = None
    if skeleton.frozen_atoms is not None:
        frozen = [skeleton.frozen_atoms[i] for i in keep_indices]

    processed = Molecule(
        symbols=symbols,
        positions=positions,
        charge=skeleton.charge,
        multiplicity=skeleton.multiplicity,
        frozen_atoms=frozen,
    )

    # Build index map: original 1-based -> new 1-based
    index_map = {}
    for new_0, orig_0 in enumerate(keep_indices):
        index_map[orig_0 + 1] = new_0 + 1

    return processed, index_map


def build_analyzer(
    algorithm_config: IterateAlgorithmConfig,
    skeleton: Molecule,
    substituents: list,
):
    """Build the analyzer for a combination's algorithm via the registry.

    Dispatch is driven by the algorithm registry (:mod:`settings`): the
    resolved :class:`AlgorithmSpec` supplies the builder, so this function
    stays algorithm-agnostic. The analyzer joins the skeleton and every
    substituent into one molecule and produces the combined structure in a
    single ``run()``.

    The algorithm's ``supports_multi_substituent`` capability is enforced
    here, before the analyzer is built, so a single-substituent-only
    algorithm fails with a clear error instead of hiding the limit inside its
    builder.

    Parameters
    ----------
    algorithm_config : IterateAlgorithmConfig
        Resolved algorithm configuration (canonical name + validated options).
    skeleton : Molecule
        Batch-preprocessed skeleton shared by every attachment.
    substituents : list of tuple
        One ``(substituent, skeleton_link_index, substituent_link_index)``
        tuple per attachment (1-based link indices).

    Returns
    -------
    object
        An analyzer exposing a ``run() -> Molecule | None`` method.

    Raises
    ------
    NotImplementedError
        If more than one substituent is requested for an algorithm whose
        spec does not advertise ``supports_multi_substituent``.
    """
    spec = get_algorithm_spec(algorithm_config.name)
    if len(substituents) > 1 and not spec.supports_multi_substituent:
        raise NotImplementedError(
            f"Algorithm '{spec.canonical_name}' does not support "
            f"multi-substituent combinations (received {len(substituents)} "
            f"substituents); it attaches a single substituent per run."
        )
    return spec.analyzer_builder(
        algorithm_config.options,
        skeleton,
        substituents,
    )


def _run_combination_task(
    combination: IterateCombination,
    pool: IterateMoleculePool,
    number: int,
) -> CombinationResult:
    """Execute a combination and return a structured :class:`CombinationResult`.

    The skeleton and all assigned substituents are joined into one molecule
    and optimized together (e.g. ETKDG embeds the whole assembly once), so
    explicit connection bonds are preserved and no intermediate ``Molecule``
    round-trip fragments the topology. The failure stage
    (``PREPROCESSING`` / ``ALGORITHM``) and a concise error reason are
    recorded; full tracebacks go to the debug log only.

    Parameters
    ----------
    combination : IterateCombination
        The combination to process (lightweight, references pool by index).
    pool : IterateMoleculePool
        Shared molecule pool; molecules are copied here at execution time.
    number : int
        Stable 1-based combination number for the report.

    Returns
    -------
    CombinationResult
        Structured result carrying status, timing and failure details.
    """
    label = combination.label
    start = time.perf_counter()
    stage = STAGE_PREPROCESSING

    try:
        # Copy skeleton from pool.
        skeleton = pool.skeletons[combination.skeleton_idx].copy()

        # Batch preprocess: remove old groups at every link position in one
        # pass so the atom indices remap only once.
        all_link_indices = [
            a.skeleton_link_index for a in combination.assignments
        ]
        processed_skeleton, index_map = _batch_preprocess_skeleton(
            skeleton, all_link_indices, combination.skeleton_indices
        )

        # Preprocess every substituent, mapping each link index onto the
        # preprocessed skeleton. Descending link_index keeps the attachment
        # order deterministic.
        substituents = []
        for assignment in sorted(
            combination.assignments,
            key=lambda a: a.skeleton_link_index,
            reverse=True,
        ):
            substituent = pool.substituents[assignment.substituent_idx].copy()
            sub_prep = SubstituentPreprocessor(
                molecule=substituent,
                link_index=assignment.substituent_link_index,
            )
            processed_sub = sub_prep.run()
            substituents.append(
                (
                    processed_sub,
                    index_map[assignment.skeleton_link_index],
                    sub_prep.get_new_link_index(),
                )
            )

        # Build the combined structure in a single run.
        stage = STAGE_ALGORITHM
        analyzer = build_analyzer(
            combination.algorithm_config, processed_skeleton, substituents
        )
        result = analyzer.run()
        duration = time.perf_counter() - start

        if result is None:
            logger.debug(f"Failed to generate molecule for {label}")
            return CombinationResult(
                combination_number=number,
                label=label,
                execution_status=STATUS_FAILED,
                duration_seconds=duration,
                failure_stage=STAGE_ALGORITHM,
                error_type="NoSolution",
                error_message="Algorithm produced no structure.",
            )

        logger.debug(f"Generated molecule for {label}")
        return CombinationResult(
            combination_number=number,
            label=label,
            execution_status=STATUS_SUCCESS,
            molecule=result,
            duration_seconds=duration,
        )

    except Exception as e:
        duration = time.perf_counter() - start
        logger.debug(f"Error in task for {label} during {stage}: {e}")
        logger.debug("Traceback for %s", label, exc_info=True)
        return CombinationResult(
            combination_number=number,
            label=label,
            execution_status=STATUS_FAILED,
            duration_seconds=duration,
            failure_stage=stage,
            error_type=type(e).__name__,
            error_message=str(e),
        )


def _run_combination_worker(
    combination: IterateCombination,
    pool: IterateMoleculePool,
    result_queue: "multiprocessing.Queue",
    number: int,
) -> None:
    """Worker function for multiprocessing.Process.

    Parameters
    ----------
    combination : IterateCombination
        The combination to process.
    pool : IterateMoleculePool
        Shared molecule pool.
    result_queue : multiprocessing.Queue
        Queue to put the :class:`CombinationResult`.
    number : int
        Stable 1-based combination number for the report.
    """
    _silence_rdkit_warnings()
    try:
        result = _run_combination_task(combination, pool, number)
        result_queue.put(result)
    except Exception as e:
        logger.debug(f"Worker process panic for {combination.label}: {e}")
        result_queue.put(
            CombinationResult(
                combination_number=number,
                label=combination.label,
                execution_status=STATUS_FAILED,
                failure_stage=STAGE_ALGORITHM,
                error_type=type(e).__name__,
                error_message=str(e),
            )
        )


class IterateJobRunner(JobRunner):
    """Job runner for Iterate jobs.

    Iterate jobs run purely in Python to generate molecular structures
    by attaching substituents to skeletons. Supports single and
    multi-position substitution via a unified slot-based pipeline.

    All combinations share a single IterateMoleculePool; molecules are
    copied only at worker dispatch time.
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

    def run_combinations(
        self,
        pool: IterateMoleculePool,
        combinations: list[IterateCombination],
        nprocs: int = 1,
        timeout: float = DEFAULT_WORKER_TIMEOUT,
        progress_callback=None,
    ) -> list[CombinationResult]:
        """Run multiple combinations using multiprocessing.Process with watchdog.

        Manually manages processes to ensure they can be forcefully killed
        (terminate/kill) if they exceed the timeout.

        Parameters
        ----------
        pool : IterateMoleculePool
            Shared molecule pool; workers copy molecules at execution time.
        combinations : list[IterateCombination]
            List of combinations to process.
        nprocs : int
            Number of processes for parallel execution. Default 1.
        timeout : float
            Timeout in seconds for each worker. Default 120 (2 minutes).
        progress_callback : callable, optional
            ``callback(completed, total)`` invoked once when execution starts
            (``completed == 0``) and once each time a combination is finalized
            (SUCCESS / FAILED / TIMED OUT / crash). Each combination is
            counted exactly once; a late worker result for an already-finalized
            (e.g. timed-out) combination is ignored. Display-only: the runner
            never depends on the callback.

        Returns
        -------
        list[CombinationResult]
            Structured result for every combination, in original order.
        """
        if not combinations:
            logger.debug("No combinations to process.")
            return []

        total = len(combinations)
        logger.debug(
            f"Running {total} combination(s) with {nprocs} "
            f"process(es), timeout={timeout}s"
        )

        number_by_label = {
            comb.label: i for i, comb in enumerate(combinations, start=1)
        }
        results_by_label: dict[str, CombinationResult] = {}
        # Labels whose result has been finalized (and counted for progress).
        # A finalized combination is never re-counted or overwritten, even if
        # a late worker result arrives after a timeout was already recorded.
        finalized: set[str] = set()

        def _finalize(result: CombinationResult) -> None:
            if result.label in finalized:
                return
            results_by_label[result.label] = result
            finalized.add(result.label)
            if progress_callback is not None:
                progress_callback(len(finalized), total)

        if progress_callback is not None:
            progress_callback(0, total)

        max_workers = 1 if nprocs == 1 else nprocs

        manager = multiprocessing.Manager()
        result_queue = manager.Queue()

        try:
            pending_combinations = deque(combinations)
            active_processes: dict[
                int,
                tuple[multiprocessing.Process, IterateCombination, float],
            ] = {}

            while pending_combinations or active_processes:
                # 1. Fill up empty slots
                while (
                    pending_combinations
                    and len(active_processes) < max_workers
                ):
                    comb = pending_combinations.popleft()
                    p = multiprocessing.Process(
                        target=_run_combination_worker,
                        args=(
                            comb,
                            pool,
                            result_queue,
                            number_by_label[comb.label],
                        ),
                        daemon=True,
                    )
                    p.start()
                    active_processes[id(p)] = (p, comb, time.time())

                # 2. Drain the queue (late results for finalized labels are
                #    ignored by _finalize).
                while True:
                    try:
                        res = result_queue.get_nowait()
                        _finalize(res)
                    except queue.Empty:
                        break
                    except Exception:
                        break

                # 3. Monitor running processes
                current_time = time.time()
                pids_to_remove = []

                for proc_id, (
                    p,
                    comb,
                    start_time,
                ) in active_processes.items():
                    if not p.is_alive():
                        p.join()
                        pids_to_remove.append(proc_id)
                    else:
                        if (current_time - start_time) > timeout:
                            logger.debug(
                                f"Timeout ({timeout}s) for {comb.label} "
                                f"(pid {p.pid}). Terminating..."
                            )
                            p.terminate()
                            p.join(timeout=0.5)
                            if p.is_alive():
                                logger.debug(
                                    f"Process {p.pid} stuck, killing..."
                                )
                                p.kill()  # SIGKILL
                                p.join(timeout=1.0)

                            _finalize(
                                CombinationResult(
                                    combination_number=number_by_label[
                                        comb.label
                                    ],
                                    label=comb.label,
                                    execution_status=STATUS_TIMED_OUT,
                                    duration_seconds=current_time - start_time,
                                    error_type="Timeout",
                                    error_message=(
                                        f"Exceeded worker timeout of "
                                        f"{timeout} seconds."
                                    ),
                                )
                            )
                            pids_to_remove.append(proc_id)

                # 4. Cleanup removed processes
                for proc_id in pids_to_remove:
                    del active_processes[proc_id]

                # 5. Sleep briefly to yield CPU
                if active_processes:
                    time.sleep(0.1)

            # Final queue drain
            while True:
                try:
                    res = result_queue.get_nowait()
                    _finalize(res)
                except queue.Empty:
                    break
                except Exception:
                    break

        finally:
            manager.shutdown()

        # Fill in missing results (worker crashed without writing to queue).
        for comb in combinations:
            if comb.label not in finalized:
                logger.debug(
                    f"No result found for {comb.label} - assuming worker "
                    f"crash."
                )
                _finalize(
                    CombinationResult(
                        combination_number=number_by_label[comb.label],
                        label=comb.label,
                        execution_status=STATUS_FAILED,
                        failure_stage=STAGE_ALGORITHM,
                        error_type="WorkerCrash",
                        error_message="Worker process produced no result.",
                    )
                )

        # Build the results list in the original combination order.
        results = [results_by_label[comb.label] for comb in combinations]

        stats = summarize_results(results)
        logger.debug(
            f"Completed: {stats['generated']}/{stats['total']} "
            f"combinations generated a structure "
            f"({stats['failed']} failed, {stats['timed_out']} timed out)."
        )
        return results

    def _load_molecule(
        self, mol_config: dict, mol_type: str, idx: int, errors: list = None
    ) -> tuple[Optional[Molecule], str]:
        """Load a Molecule from configuration dict.

        Parameters
        ----------
        mol_config : dict
            Configuration containing file_path.
        mol_type : str
            "skeleton" or "substituent" for logging.
        idx : int
            Index for logging.
        errors : list, optional
            When provided, structured input-error records are appended here
            for the run report (never raises; failures still return None).

        Returns
        -------
        tuple[Molecule | None, str]
            (Loaded molecule or None if loading fails, label).
        """
        label = mol_config.get("label") or f"{mol_type}{idx + 1}"
        file_path = mol_config.get("file_path")
        raw_path = mol_config.get("file_path_raw", file_path)

        def _record(error_type: str, message: str) -> None:
            if errors is not None:
                errors.append(
                    {
                        "type": mol_type,
                        "label": label,
                        "raw_path": raw_path,
                        "resolved_path": file_path,
                        "error_type": error_type,
                        "error_message": message,
                    }
                )

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
                _record("MissingFilePath", "No file_path provided.")
                return None, label

            # S2 Check: Validate indices against atom count
            num_atoms = molecule.num_atoms

            # Validate link_index
            link_indices = mol_config.get("link_index")
            if link_indices is not None and not isinstance(link_indices, list):
                if isinstance(link_indices, int):
                    link_indices = [link_indices]

            if link_indices:
                invalid_links = [i for i in link_indices if i > num_atoms]
                if invalid_links:
                    logger.error(
                        f"{mol_type.capitalize()} '{label}': link_index "
                        f"{invalid_links} out of bounds. Molecule has "
                        f"{num_atoms} atoms."
                    )
                    _record(
                        "IndexError",
                        f"link_index {invalid_links} out of bounds "
                        f"(molecule has {num_atoms} atoms).",
                    )
                    return None, label

            # Validate skeleton_indices (only for skeletons)
            if mol_type == "skeleton":
                skel_indices = mol_config.get("skeleton_indices")
                if skel_indices is not None and not isinstance(
                    skel_indices, list
                ):
                    if isinstance(skel_indices, int):
                        skel_indices = [skel_indices]

                if skel_indices:
                    invalid_skels = [i for i in skel_indices if i > num_atoms]
                    if invalid_skels:
                        logger.error(
                            f"{mol_type.capitalize()} '{label}': "
                            f"skeleton_indices {invalid_skels} out of bounds. "
                            f"Molecule has {num_atoms} atoms."
                        )
                        _record(
                            "IndexError",
                            f"skeleton_indices {invalid_skels} out of bounds "
                            f"(molecule has {num_atoms} atoms).",
                        )
                        return None, label

            return molecule, label
        except Exception as e:
            logger.error(
                f"Failed to load {mol_type} '{label}' from '{raw_path}' "
                f"(resolved to '{file_path}'): {e}"
            )
            logger.debug("Load traceback for %s", label, exc_info=True)
            _record(type(e).__name__, str(e))
            return None, label

    def _generate_combinations(
        self, job: "IterateJob"
    ) -> tuple[IterateMoleculePool, list[IterateCombination], list, int]:
        """Generate all combinations from job settings.

        Unified path for all skeletons. Skeletons without explicit slots
        are converted to virtual single-position slots (one per
        link_index), with an implicit group number assigned from the
        global contiguous sequence.

        combination_mode controls how multi-group slots are processed:
        - 'independent': each group generates combinations independently,
          results are merged (union).
        - 'global': all groups merge position options first, then a single
          Cartesian product is taken across all positions.

        Parameters
        ----------
        job : IterateJob
            The job containing settings.

        Returns
        -------
        tuple[IterateMoleculePool, list[IterateCombination], list, int]
            Shared molecule pool, list of all valid combinations, structured
            input-error records, and the total attachment-site count.
        """
        pool = IterateMoleculePool()
        combinations: list[IterateCombination] = []
        input_errors: list = []
        attachment_site_count = 0

        skeleton_list = job.settings.skeleton_list or []
        substituent_list = job.settings.substituent_list or []
        algorithm_config = job.settings.algorithm_config
        combination_mode = job.settings.combination_mode

        # --- Load substituents into pool ---
        valid_substituents: list[tuple[int, str, dict]] = []
        # group_number -> [(pool_idx, label, link_idx)]
        sub_by_group: dict[int, list[tuple[int, str, int]]] = {}

        for sub_idx, sub_config in enumerate(substituent_list):
            mol, label = self._load_molecule(
                sub_config, "substituent", sub_idx, input_errors
            )
            if mol is None:
                continue

            sub_link_index = sub_config.get("link_index")
            if not sub_link_index:
                logger.warning(
                    f"Substituent '{label}' has no link_index, skipping."
                )
                continue

            pool_idx = len(pool.substituents)
            pool.substituents.append(mol)
            valid_substituents.append((pool_idx, label, sub_config))

            sub_link_idx = sub_link_index[0]

            # Index by group
            groups = sub_config.get("groups")
            if groups:
                for g in groups:
                    if g not in sub_by_group:
                        sub_by_group[g] = []
                    sub_by_group[g].append((pool_idx, label, sub_link_idx))

        # --- Compute global group map for all skeletons ---
        # Every skeleton participates: no-slots → 1 implicit group,
        # slots → N groups.
        implicit_group_map: dict[int, int] = {}  # skel_idx -> group
        _next_group = 1
        for _i, _sc in enumerate(skeleton_list):
            if _sc.get("slots"):
                _next_group += len(_sc["slots"])
            else:
                implicit_group_map[_i] = _next_group
                _next_group += 1
        total_groups = _next_group - 1

        # For substituents without 'groups' (e.g. from TOML), assign to
        # all groups so they are available everywhere.
        for sp_idx, sp_label, sp_config in valid_substituents:
            if not sp_config.get("groups"):
                sub_link_idx = sp_config.get("link_index")[0]
                for g in range(1, total_groups + 1):
                    if g not in sub_by_group:
                        sub_by_group[g] = []
                    sub_by_group[g].append((sp_idx, sp_label, sub_link_idx))

        # --- Load skeletons and generate combinations ---
        for skel_idx, skel_config in enumerate(skeleton_list):
            skel_mol, skel_label = self._load_molecule(
                skel_config, "skeleton", skel_idx, input_errors
            )
            if skel_mol is None:
                continue

            skel_pool_idx = len(pool.skeletons)
            pool.skeletons.append(skel_mol)
            skeleton_indices = skel_config.get("skeleton_indices")

            # Count attachment sites for this loaded skeleton.
            _slots = skel_config.get("slots")
            if _slots:
                attachment_site_count += sum(
                    len(s["link_indices"]) for s in _slots
                )
            else:
                attachment_site_count += len(
                    skel_config.get("link_index") or []
                )

            slots = skel_config.get("slots")

            if slots:
                # --- Slots mode ---
                if combination_mode == "global":
                    # Merge all slots into one position_options, single product
                    position_options = self._build_position_options(
                        skel_label, slots, sub_by_group
                    )
                    self._expand_position_options(
                        skel_pool_idx,
                        skel_label,
                        skeleton_indices,
                        position_options,
                        algorithm_config,
                        combinations,
                    )
                else:
                    # independent: each slot group independently, union results
                    for slot in slots:
                        position_options = self._build_position_options(
                            skel_label, [slot], sub_by_group
                        )
                        self._expand_position_options(
                            skel_pool_idx,
                            skel_label,
                            skeleton_indices,
                            position_options,
                            algorithm_config,
                            combinations,
                        )
            else:
                # --- No explicit slots: virtual single-position slots ---
                skel_link_indices = skel_config.get("link_index")
                if not skel_link_indices:
                    logger.warning(
                        f"Skeleton '{skel_label}' has no link_index "
                        f"or slots, skipping."
                    )
                    continue

                assigned_group = implicit_group_map[skel_idx]
                subs_for_group = sub_by_group.get(assigned_group, [])
                if not subs_for_group:
                    logger.warning(
                        f"Skeleton '{skel_label}' "
                        f"(implicit group {assigned_group}): "
                        f"no substituents assigned to this group, "
                        f"skipping."
                    )
                    continue
                virtual_sub_by_group: dict[int, list[tuple[int, str, int]]] = {
                    assigned_group: subs_for_group
                }
                for link_idx in skel_link_indices:
                    virtual_slot = {
                        "group": assigned_group,
                        "link_indices": [link_idx],
                    }
                    position_options = self._build_position_options(
                        skel_label,
                        [virtual_slot],
                        virtual_sub_by_group,
                    )
                    self._expand_position_options(
                        skel_pool_idx,
                        skel_label,
                        skeleton_indices,
                        position_options,
                        algorithm_config,
                        combinations,
                    )

        # Guard: combination labels must be unique so results are not
        # silently overwritten and separate-output files never collide.
        seen_labels: dict = {}
        for combo in combinations:
            if combo.label in seen_labels:
                raise ValueError(
                    f"Duplicate combination label '{combo.label}'. This "
                    f"usually means duplicate skeleton/substituent labels or "
                    f"overlapping link positions in the configuration."
                )
            seen_labels[combo.label] = combo

        return pool, combinations, input_errors, attachment_site_count

    def _build_position_options(
        self,
        skel_label: str,
        slots: list[dict],
        sub_by_group: dict[int, list[tuple[int, str, int]]],
    ) -> dict[int, list[IterateAssignment]]:
        """Build per-position candidate lists from slots and groups.

        Parameters
        ----------
        skel_label : str
            Skeleton label (for logging).
        slots : list[dict]
            Slot definitions with 'group' and 'link_indices'.
        sub_by_group : dict
            Group number -> list of (pool_idx, label, link_idx).

        Returns
        -------
        dict[int, list[IterateAssignment]]
            Mapping from skeleton link_index to list of candidate
            IterateAssignment objects.
        """
        position_options: dict[int, list[IterateAssignment]] = {}

        for slot in slots:
            group = slot["group"]
            link_indices = slot["link_indices"]

            subs_for_group = sub_by_group.get(group, [])
            if not subs_for_group:
                logger.warning(
                    f"Skeleton '{skel_label}': no substituents "
                    f"for group {group}, skipping slot."
                )
                continue

            for sub_pool_idx, sub_label, sub_link_idx in subs_for_group:
                for link_idx in link_indices:
                    if link_idx not in position_options:
                        position_options[link_idx] = []
                    # Deduplicate: same substituent at same position
                    is_dup = any(
                        e.substituent_idx == sub_pool_idx
                        and e.substituent_link_index == sub_link_idx
                        for e in position_options[link_idx]
                    )
                    if not is_dup:
                        position_options[link_idx].append(
                            IterateAssignment(
                                substituent_idx=sub_pool_idx,
                                substituent_label=sub_label,
                                substituent_link_index=sub_link_idx,
                                skeleton_link_index=link_idx,
                            )
                        )

        return position_options

    def _expand_position_options(
        self,
        skel_pool_idx: int,
        skel_label: str,
        skeleton_indices: Optional[list[int]],
        position_options: dict[int, list[IterateAssignment]],
        algorithm_config: IterateAlgorithmConfig,
        combinations: list[IterateCombination],
    ) -> None:
        """Take the Cartesian product of position options and create combinations.

        Parameters
        ----------
        skel_pool_idx : int
            Index of this skeleton in the pool.
        skel_label : str
            Skeleton label.
        skeleton_indices : list[int] or None
            Skeleton core indices.
        position_options : dict[int, list[IterateAssignment]]
            Per-position candidate lists.
        algorithm_config : IterateAlgorithmConfig
            Resolved algorithm configuration for the generated combinations.
        combinations : list[IterateCombination]
            Output list to append to.
        """
        if not position_options:
            logger.warning(
                f"Skeleton '{skel_label}': no valid position options."
            )
            return

        # For each position: [None (don't fill), candidate1, candidate2, ...]
        sorted_positions = sorted(position_options.keys())
        per_position_choices = []
        for pos in sorted_positions:
            per_position_choices.append([None] + position_options[pos])

        # Cartesian product across all positions
        for combo in product(*per_position_choices):
            assignments = [a for a in combo if a is not None]
            if not assignments:
                continue  # skip the all-empty case

            combination = IterateCombination(
                skeleton_idx=skel_pool_idx,
                skeleton_label=skel_label,
                skeleton_indices=skeleton_indices,
                assignments=list(assignments),
                algorithm_config=algorithm_config,
            )
            combinations.append(combination)
            logger.debug(f"Created combination: {combination.label}")

    @staticmethod
    def _write_xyz_file(
        filename: str, items: list[tuple[str, Molecule]]
    ) -> None:
        """Write one or more labelled molecules to an XYZ file."""
        with open(filename, "w") as f:
            for label, mol in items:
                f.write(f"{mol.num_atoms}\n")
                f.write(f"       {label}\n")
                for symbol, pos in zip(mol.chemical_symbols, mol.positions):
                    f.write(
                        f"{symbol:2s}  {pos[0]:15.10f}  "
                        f"{pos[1]:15.10f}  {pos[2]:15.10f}\n"
                    )

    def _write_outputs(
        self, results: list[CombinationResult], job: "IterateJob"
    ) -> list[str]:
        """Write successful structures and annotate each result.

        Mutates each successful :class:`CombinationResult` with its
        ``output_path`` and ``structure_index``; a structure that fails to
        write is flipped to ``WRITE FAILED`` with the failure recorded, so the
        report can distinguish generation success from delivery success.

        Returns
        -------
        list[str]
            Paths of the files actually written (empty when nothing was
            delivered; no empty placeholder file is ever created).
        """
        written_paths: list[str] = []
        successful = [
            r
            for r in results
            if r.execution_status == STATUS_SUCCESS and r.molecule is not None
        ]

        if job.separate_outputs:
            output_dir = job.output_directory or "."
            os.makedirs(output_dir, exist_ok=True)
            logger.debug(
                f"Writing separate output files to directory: {output_dir}"
            )
            used_filenames: set[str] = set()
            for r in successful:
                filename = os.path.join(output_dir, f"{r.label}.xyz")
                if filename in used_filenames:
                    raise ValueError(
                        f"Duplicate output filename '{filename}' for label "
                        f"'{r.label}'; refusing to overwrite a previously "
                        f"written structure."
                    )
                used_filenames.add(filename)
                try:
                    self._write_xyz_file(filename, [(r.label, r.molecule)])
                    r.output_path = filename
                    r.structure_index = len(written_paths) + 1
                    written_paths.append(filename)
                    logger.debug(f"Wrote {filename}")
                except Exception as e:
                    logger.debug(f"Failed to write {filename}: {e}")
                    r.execution_status = STATUS_WRITE_FAILED
                    r.failure_stage = STAGE_WRITE
                    r.error_type = type(e).__name__
                    r.error_message = str(e)
                    r.output_path = filename
            logger.debug(
                f"Wrote {len(written_paths)} separate molecule file(s) to "
                f"{output_dir}"
            )
            return written_paths

        # Merged mode: only write when at least one structure succeeded, so a
        # fully-failed run does not leave an empty file that looks like success.
        if not successful:
            logger.debug(
                "No structures were generated; not creating an empty merged "
                "output file."
            )
            return written_paths

        outputfile = job.outputfile
        output_dir = os.path.dirname(outputfile)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)

        try:
            self._write_xyz_file(
                outputfile, [(r.label, r.molecule) for r in successful]
            )
        except Exception as e:
            logger.debug(f"Failed to write merged output {outputfile}: {e}")
            for r in successful:
                r.execution_status = STATUS_WRITE_FAILED
                r.failure_stage = STAGE_WRITE
                r.error_type = type(e).__name__
                r.error_message = str(e)
                r.output_path = outputfile
            return written_paths

        for index, r in enumerate(successful, start=1):
            r.output_path = outputfile
            r.structure_index = index
        written_paths.append(outputfile)
        logger.debug(f"Wrote {len(successful)} molecule(s) to {outputfile}")
        return written_paths

    def run(
        self,
        job: "IterateJob",
        progress_callback=None,
        **kwargs,
    ) -> "IterateRunSummary":
        """Run an IterateJob and emit a plain-text run report.

        Generates combinations, executes them, writes the output XYZ(s), and
        finally renders a Gaussian-log-style ``<config_stem>_iterate.out``
        report next to the output. The report is always attempted once the
        run phase is entered, even on full failure or timeout.

        Parameters
        ----------
        job : IterateJob
            The iterate job to run.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        IterateRunSummary
            Execution/delivery counts, the paths written, and the report path.
        """
        logger.debug(f"IterateJobRunner.run() called for job: {job}")

        if self.fake:
            logger.debug("Fake mode enabled, not actually running job.")
            return IterateRunSummary()

        run_id = uuid.uuid4().hex[:12]
        started_at = datetime.now()
        start_perf = time.perf_counter()

        try:
            pool, combinations, input_errors, attachment_site_count = (
                self._generate_combinations(job)
            )
            logger.debug(f"Generated {len(combinations)} combination(s)")

            if combinations:
                results = self.run_combinations(
                    pool,
                    combinations,
                    nprocs=job.nprocs,
                    timeout=job.timeout,
                    progress_callback=progress_callback,
                )
                output_paths = self._write_outputs(results, job)
            else:
                logger.debug("No valid combinations to process.")
                results = []
                output_paths = []

            finished_at = datetime.now()
            duration = time.perf_counter() - start_perf
            stats = summarize_results(results, len(input_errors))

            report = self._build_report(
                job=job,
                run_id=run_id,
                started_at=started_at,
                finished_at=finished_at,
                duration=duration,
                pool=pool,
                combinations=combinations,
                input_errors=input_errors,
                attachment_site_count=attachment_site_count,
                results=results,
            )
            summary_path, summary_error = self._write_report(job, report)

            # A summary-write failure is itself an error termination.
            exit_code = stats["exit_code"] if summary_error is None else 1

            logger.debug(
                f"IterateJobRunner completed job: {stats['generated']}/"
                f"{stats['total']} generated, "
                f"{stats['write_succeeded']} written, "
                f"{stats['failed']} failed, {stats['timed_out']} timed out "
                f"(error codes {stats['error_codes'] or 'none'}, "
                f"exit {exit_code})."
            )
            return IterateRunSummary(
                total=stats["total"],
                succeeded=stats["generated"],
                failed=stats["failed"],
                timed_out=stats["timed_out"],
                write_failed=stats["write_failed"],
                structures_written=stats["write_succeeded"],
                input_error_count=len(input_errors),
                error_codes=stats["error_codes"],
                output_paths=output_paths,
                summary_path=summary_path,
                summary_write_error=summary_error,
                exit_code=exit_code,
            )
        except KeyboardInterrupt:
            # Best-effort interrupt (SIGINT) summary; keep any structures
            # already written and re-raise so the CLI exits with 130.
            logger.debug("Iterate interrupted by user (SIGINT).")
            finished_at = datetime.now()
            duration = time.perf_counter() - start_perf
            try:
                report = self._build_error_report(
                    job,
                    run_id,
                    started_at,
                    finished_at,
                    duration,
                    "Interrupted by user (SIGINT).",
                    error_codes=[ERROR_CODE_INTERRUPTED],
                    exit_code=130,
                )
                self._write_report(job, report)
            except Exception as report_error:
                logger.error(
                    f"Failed to write interrupt report: {report_error}"
                )
            raise
        except Exception as error:
            # Unexpected runtime failure: best-effort INTERNAL ERROR report.
            logger.error(f"Unexpected error during iterate run: {error}")
            logger.debug("Iterate run traceback", exc_info=True)
            finished_at = datetime.now()
            duration = time.perf_counter() - start_perf
            summary_path, summary_error = None, None
            try:
                report = self._build_error_report(
                    job, run_id, started_at, finished_at, duration, str(error)
                )
                summary_path, summary_error = self._write_report(job, report)
            except Exception as report_error:
                summary_error = str(report_error)
                logger.error(
                    f"Failed to write internal-error report: {report_error}"
                )
            return IterateRunSummary(
                error_codes=[ERROR_CODE_INTERNAL],
                summary_path=summary_path,
                summary_write_error=summary_error,
                exit_code=1,
            )

    @staticmethod
    def _sha256_of_file(path: Optional[str]) -> str:
        """Return the SHA256 hex digest of a file, or ``N/A`` if unavailable."""
        if not path or not os.path.isfile(path):
            return "N/A"
        try:
            digest = hashlib.sha256()
            with open(path, "rb") as handle:
                for chunk in iter(lambda: handle.read(65536), b""):
                    digest.update(chunk)
            return digest.hexdigest()
        except OSError:
            return "N/A"

    def _build_report(
        self,
        job: "IterateJob",
        run_id: str,
        started_at: datetime,
        finished_at: datetime,
        duration: float,
        pool: IterateMoleculePool,
        combinations: list,
        input_errors: list,
        attachment_site_count: int,
        results: list,
    ) -> IterateReport:
        """Assemble the :class:`IterateReport` from the run state."""
        settings = job.settings
        config_file = settings.config_file

        # Per-skeleton combination counts, preserving first-appearance order.
        per_counts: dict = {}
        for combo in combinations:
            per_counts[combo.skeleton_label] = (
                per_counts.get(combo.skeleton_label, 0) + 1
            )

        if job.separate_outputs:
            output_mode = "separate"
            output_location = os.path.abspath(
                job.output_directory or os.getcwd()
            )
        else:
            output_mode = "merged"
            output_location = os.path.abspath(job.outputfile)

        command_line = getattr(job, "command_line", None) or " ".join(sys.argv)

        import rdkit

        from chemsmart import __version__ as chemsmart_version

        return IterateReport(
            run_id=run_id,
            chemsmart_version=chemsmart_version,
            rdkit_version=getattr(rdkit, "__version__", "unknown"),
            started_at=started_at,
            finished_at=finished_at,
            duration_seconds=duration,
            working_directory=os.getcwd(),
            command_line=command_line,
            config_file=config_file,
            config_sha256=self._sha256_of_file(config_file),
            skeleton_entries=list(settings.skeleton_list or []),
            substituent_entries=list(settings.substituent_list or []),
            input_errors=input_errors,
            algorithm_name=settings.algorithm_config.name,
            algorithm_options=dict(settings.algorithm_config.options),
            combination_mode=settings.combination_mode,
            nprocs=job.nprocs,
            timeout_seconds=job.timeout,
            output_mode=output_mode,
            loaded_skeletons=len(pool.skeletons),
            loaded_substituents=len(pool.substituents),
            attachment_site_count=attachment_site_count,
            total_combinations=len(combinations),
            per_skeleton_counts=list(per_counts.items()),
            results=results,
            output_location=output_location,
        )

    def _build_error_report(
        self,
        job: "IterateJob",
        run_id: str,
        started_at: datetime,
        finished_at: datetime,
        duration: float,
        error_message: str,
        error_codes: Optional[list] = None,
        exit_code: int = 1,
    ) -> IterateReport:
        """Assemble a minimal report for an abnormal termination.

        Used for both an unexpected runtime failure (INTERNAL ERROR) and a
        user interrupt (INTERRUPTED). Uses whatever job/settings context is
        available; the results section is empty because the run did not
        complete normally.
        """
        import rdkit

        from chemsmart import __version__ as chemsmart_version

        settings = job.settings
        config_file = settings.config_file
        command_line = getattr(job, "command_line", None) or " ".join(sys.argv)
        return IterateReport(
            run_id=run_id,
            chemsmart_version=chemsmart_version,
            rdkit_version=getattr(rdkit, "__version__", "unknown"),
            started_at=started_at,
            finished_at=finished_at,
            duration_seconds=duration,
            working_directory=os.getcwd(),
            command_line=command_line,
            config_file=config_file,
            config_sha256=self._sha256_of_file(config_file),
            skeleton_entries=list(settings.skeleton_list or []),
            substituent_entries=list(settings.substituent_list or []),
            algorithm_name=settings.algorithm_config.name,
            algorithm_options=dict(settings.algorithm_config.options),
            combination_mode=settings.combination_mode,
            nprocs=job.nprocs,
            timeout_seconds=job.timeout,
            output_mode="separate" if job.separate_outputs else "merged",
            extra_error_codes=(
                error_codes
                if error_codes is not None
                else [ERROR_CODE_INTERNAL]
            ),
            exit_code_override=exit_code,
            error_message=error_message,
        )

    def _write_report(
        self, job: "IterateJob", report: IterateReport
    ) -> tuple[Optional[str], Optional[str]]:
        """Render and atomically write the run report.

        The filename is derived from the YAML config stem
        (``<stem>_iterate.out``) and placed next to the output (merged file's
        directory, the separate-outputs directory, or the current working
        directory). Returns ``(summary_path, error_message)``; on failure the
        path is ``None`` and a clear error is logged to the terminal.
        """
        config_file = job.settings.config_file
        stem = Path(config_file).stem if config_file else "iterate"
        filename = f"{stem}_iterate.out"
        if job.separate_outputs:
            summary_dir = job.output_directory or os.getcwd()
        else:
            out_dir = os.path.dirname(job.outputfile)
            summary_dir = out_dir if out_dir else os.getcwd()
        summary_path = os.path.join(summary_dir, filename)

        try:
            write_report_atomically(summary_path, report.render())
            logger.debug(f"Wrote run report to {summary_path}")
            return summary_path, None
        except Exception as e:
            message = f"could not write report to {summary_path}: {e}"
            logger.debug(f"Failed to write run report: {message}")
            return None, message
