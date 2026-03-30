import logging
import multiprocessing
import os
import queue
import time
from collections import deque
from dataclasses import dataclass, field
from itertools import product
from typing import TYPE_CHECKING, Optional

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.iterate.iterate import (
    IterateAnalyzer,
    SkeletonPreprocessor,
    SubstituentPreprocessor,
)
from chemsmart.jobs.runner import JobRunner

if TYPE_CHECKING:
    from chemsmart.jobs.iterate.job import IterateJob

logger = logging.getLogger(__name__)

# Default timeout for each combination worker (in seconds)
DEFAULT_WORKER_TIMEOUT = 120  # 2 minutes


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
    method : str
        Optimization method for orientation search.
    sphere_direction_samples_num : int
        Number of sphere direction samples for orientation search (default 96).
    axial_rotations_sample_num : int
        Number of axial rotation samples per direction (default 6).
    """

    skeleton_idx: int
    skeleton_label: str
    skeleton_indices: Optional[list[int]]  # 1-based, or None
    assignments: list[IterateAssignment] = field(default_factory=list)
    method: str = "lagrange_multipliers"
    sphere_direction_samples_num: int = 96
    axial_rotations_sample_num: int = 6

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


def _run_combination_task(
    combination: IterateCombination,
    pool: IterateMoleculePool,
) -> tuple[str, Optional[Molecule]]:
    """Execute a combination: batch-preprocess skeleton, then sequentially
    insert all assigned substituents.

    Works for single-assignment and multi-assignment combinations.

    Parameters
    ----------
    combination : IterateCombination
        The combination to process (lightweight, references pool by index).
    pool : IterateMoleculePool
        Shared molecule pool; molecules are copied here at execution time.

    Returns
    -------
    tuple[str, Molecule | None]
        (label, result_molecule) or (label, None) if failed.
    """
    label = combination.label

    try:
        # Copy skeleton from pool
        skeleton = pool.skeletons[combination.skeleton_idx].copy()

        # Collect all skeleton link indices from assignments
        all_link_indices = [
            a.skeleton_link_index for a in combination.assignments
        ]

        # Batch preprocess: remove old groups at all link positions
        processed_skeleton, index_map = _batch_preprocess_skeleton(
            skeleton, all_link_indices, combination.skeleton_indices
        )

        # Sequential insertion (descending link_index for determinism)
        sorted_assignments = sorted(
            combination.assignments,
            key=lambda a: a.skeleton_link_index,
            reverse=True,
        )

        current_skeleton = processed_skeleton

        for assignment in sorted_assignments:
            # Copy substituent from pool
            substituent = pool.substituents[assignment.substituent_idx].copy()

            # Map original link index to preprocessed skeleton index
            new_skel_link = index_map[assignment.skeleton_link_index]

            # Preprocess substituent
            sub_prep = SubstituentPreprocessor(
                molecule=substituent,
                link_index=assignment.substituent_link_index,
            )
            processed_sub = sub_prep.run()
            new_sub_link = sub_prep.get_new_link_index()

            # Run IterateAnalyzer to combine skeleton + substituent
            analyzer = IterateAnalyzer(
                skeleton=current_skeleton,
                substituent=processed_sub,
                skeleton_link_index=new_skel_link,
                substituent_link_index=new_sub_link,
                method=combination.method,
                sphere_direction_samples_num=combination.sphere_direction_samples_num,
                axial_rotations_sample_num=combination.axial_rotations_sample_num,
            )
            result = analyzer.run()

            if result is None:
                logger.warning(
                    f"Failed at {assignment.substituent_label}"
                    f"@{assignment.skeleton_link_index} in {label}"
                )
                return (label, None)

            # Use result as new skeleton for next insertion.
            # Skeleton atoms are at positions 0..n_skel-1 in the combined
            # result, so index_map values remain valid for remaining
            # assignments.
            current_skeleton = result

        logger.info(f"Generated molecule for {label}")
        return (label, current_skeleton)

    except Exception as e:
        logger.error(f"Error in task for {label}: {e}")
        return (label, None)


def _run_combination_worker(
    combination: IterateCombination,
    pool: IterateMoleculePool,
    result_queue: "multiprocessing.Queue",
) -> None:
    """Worker function for multiprocessing.Process.

    Parameters
    ----------
    combination : IterateCombination
        The combination to process.
    pool : IterateMoleculePool
        Shared molecule pool.
    result_queue : multiprocessing.Queue
        Queue to put the result tuple.
    """
    try:
        result_pair = _run_combination_task(combination, pool)
        result_queue.put(result_pair)
    except Exception as e:
        logger.error(f"Worker process panic for {combination.label}: {e}")
        result_queue.put((combination.label, None))


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
    ) -> list[tuple[str, Optional[Molecule]]]:
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

        Returns
        -------
        list[tuple[str, Molecule | None]]
            List of (label, result_molecule) tuples.
        """
        if not combinations:
            logger.warning("No combinations to process.")
            return []

        logger.info(
            f"Running {len(combinations)} combination(s) with {nprocs} "
            f"process(es), timeout={timeout}s"
        )

        results_dict: dict[str, Optional[Molecule]] = {}
        failed_labels: list[str] = []
        timed_out_labels: list[str] = []

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
                        args=(comb, pool, result_queue),
                        daemon=True,
                    )
                    p.start()
                    active_processes[id(p)] = (p, comb, time.time())

                # 2. Check for results in queue
                while True:
                    try:
                        lbl, mol = result_queue.get_nowait()
                        results_dict[lbl] = mol
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
                            logger.warning(
                                f"Timeout ({timeout}s) for {comb.label} "
                                f"(pid {p.pid}). Terminating..."
                            )
                            p.terminate()
                            p.join(timeout=0.5)
                            if p.is_alive():
                                logger.warning(
                                    f"Process {p.pid} stuck, killing..."
                                )
                                p.kill()  # SIGKILL
                                p.join(timeout=1.0)

                            timed_out_labels.append(comb.label)
                            results_dict[comb.label] = None
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
                    lbl, mol = result_queue.get_nowait()
                    results_dict[lbl] = mol
                except queue.Empty:
                    break
                except Exception:
                    break

        finally:
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
            f"Completed: {len(successful_labels)}/{len(results)} "
            f"molecules generated successfully"
        )

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
        """Load a Molecule from configuration dict.

        Parameters
        ----------
        mol_config : dict
            Configuration containing file_path.
        mol_type : str
            "skeleton" or "substituent" for logging.
        idx : int
            Index for logging.

        Returns
        -------
        tuple[Molecule | None, str]
            (Loaded molecule or None if loading fails, label).
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
                        return None, label

            return molecule, label
        except Exception as e:
            logger.error(f"Failed to load {mol_type} '{label}': {e}")
            return None, label

    def _generate_combinations(
        self, job: "IterateJob"
    ) -> tuple[IterateMoleculePool, list[IterateCombination]]:
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
        tuple[IterateMoleculePool, list[IterateCombination]]
            Shared molecule pool and list of all valid combinations.
        """
        pool = IterateMoleculePool()
        combinations: list[IterateCombination] = []

        skeleton_list = job.settings.skeleton_list or []
        substituent_list = job.settings.substituent_list or []
        method = job.settings.method
        sphere_samples = job.settings.sphere_direction_samples_num
        axial_samples = job.settings.axial_rotations_sample_num
        combination_mode = job.settings.combination_mode

        # --- Load substituents into pool ---
        valid_substituents: list[tuple[int, str, dict]] = []
        # group_number -> [(pool_idx, label, link_idx)]
        sub_by_group: dict[int, list[tuple[int, str, int]]] = {}

        for sub_idx, sub_config in enumerate(substituent_list):
            mol, label = self._load_molecule(
                sub_config, "substituent", sub_idx
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
                    sub_by_group[g].append(
                        (sp_idx, sp_label, sub_link_idx)
                    )

        # --- Load skeletons and generate combinations ---
        for skel_idx, skel_config in enumerate(skeleton_list):
            skel_mol, skel_label = self._load_molecule(
                skel_config, "skeleton", skel_idx
            )
            if skel_mol is None:
                continue

            skel_pool_idx = len(pool.skeletons)
            pool.skeletons.append(skel_mol)
            skeleton_indices = skel_config.get("skeleton_indices")

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
                        method,
                        sphere_samples,
                        axial_samples,
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
                            method,
                            sphere_samples,
                            axial_samples,
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
                virtual_sub_by_group: dict[
                    int, list[tuple[int, str, int]]
                ] = {assigned_group: subs_for_group}
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
                        method,
                        sphere_samples,
                        axial_samples,
                        combinations,
                    )

        return pool, combinations

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
        method: str,
        sphere_samples: int,
        axial_samples: int,
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
        method, sphere_samples, axial_samples
            Algorithm parameters.
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
                method=method,
                sphere_direction_samples_num=sphere_samples,
                axial_rotations_sample_num=axial_samples,
            )
            combinations.append(combination)
            logger.info(f"Created combination: {combination.label}")

    def _write_outputs(
        self, results: list[tuple[str, Optional[Molecule]]], job: "IterateJob"
    ) -> None:
        """Write execution results to output file(s).

        Parameters
        ----------
        results : list[tuple[str, Molecule | None]]
            List of (label, molecule) tuples from run_combinations.
        job : IterateJob
             The job instance with configuration options.
        """
        successful_count = 0

        if job.separate_outputs:
            output_dir = job.output_directory or "."
            os.makedirs(output_dir, exist_ok=True)

            logger.info(
                f"Writing separate output files to directory: {output_dir}"
            )

            for label, mol in results:
                if mol is not None:
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
                f"Wrote {successful_count} separate molecule files to "
                f"{output_dir}"
            )

        else:
            outputfile = job.outputfile
            output_dir = os.path.dirname(outputfile)
            if output_dir:
                os.makedirs(output_dir, exist_ok=True)

            with open(outputfile, "w") as f:
                for label, mol in results:
                    if mol is not None:
                        f.write(f"{mol.num_atoms}\n")
                        f.write(f"       {label}\n")
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
        """Run an IterateJob.

        1. Generate all combinations
        2. Run combinations with multiprocessing
        3. Write results to output xyz file

        Parameters
        ----------
        job : IterateJob
            The iterate job to run.
        **kwargs
            Additional keyword arguments.
        """
        logger.info(f"IterateJobRunner.run() called for job: {job}")

        if self.fake:
            logger.info("Fake mode enabled, not actually running job.")
            return

        pool, combinations = self._generate_combinations(job)
        logger.info(f"Generated {len(combinations)} combination(s)")

        if not combinations:
            logger.warning("No valid combinations to process.")
            return

        results = self.run_combinations(
            pool, combinations, nprocs=job.nprocs, timeout=job.timeout
        )

        self._write_outputs(results, job)

        logger.info("IterateJobRunner completed job")
