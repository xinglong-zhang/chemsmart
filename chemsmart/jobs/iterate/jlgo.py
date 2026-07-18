"""Configurable Version 2 Joint-Lagrange substituent placement.

The optimized variable is a 6K vector.  Each substituent contributes the
Cartesian coordinates of its linking atom and three intrinsic xyz Euler
angles.  The skeleton is fixed, every linking atom is constrained to its bond
sphere, and all substituents are optimized together by SLSQP.

The default strategy is the best Test7.2 fusion, while the atomic radii,
steric scales, WCA parameters, and quality gate use the RDKit calibration from
Test8:

region exclusion + adaptive sampling + feasible-domain pruning + greedy
diversity selection + callback-based SLSQP convergence detection + structural
quality early stop + quality repair + ordinary-SLSQP fallback.

The numerical core is self-contained apart from NumPy, SciPy, and RDKit. All
public numerical behavior is controlled through ``JointLagrangeConfig``.

The Version 2 numerical core is ported from the standalone
``joint_lagrange_optimizer`` (only the test-only ``write_xyz`` helper is
dropped, and storage-only internal records are compacted without changing
the numerical path). ``IterateJointLagrangeAnalyzer`` provides the CHEMSMART
integration at the end of this algorithm module. The registry imports this
module lazily, so the ETKDG path never loads the Joint-Lagrange SciPy core.
"""

from __future__ import annotations

import itertools
import logging
import time
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from rdkit import Chem
from scipy.optimize import minimize
from scipy.spatial.transform import Rotation
from scipy.special import logsumexp

from chemsmart.io.molecules.structure import Molecule

logger = logging.getLogger(__name__)

_PERIODIC_TABLE = Chem.GetPeriodicTable()
_ELEMENT_SYMBOLS = tuple(
    _PERIODIC_TABLE.GetElementSymbol(atomic_number)
    for atomic_number in range(1, 119)
)
COVALENT_RADII = {
    element: float(
        _PERIODIC_TABLE.GetRcovalent(_PERIODIC_TABLE.GetAtomicNumber(element))
    )
    for element in _ELEMENT_SYMBOLS
}

VDW_RADII = {
    element: float(
        _PERIODIC_TABLE.GetRvdw(_PERIODIC_TABLE.GetAtomicNumber(element))
    )
    for element in _ELEMENT_SYMBOLS
}

_TWO_POW_ONE_SIXTH = 2.0 ** (1.0 / 6.0)


@dataclass
class Mol:
    """Lightweight molecular coordinate container used by Joint-Lagrange.

    Parameters
    ----------
    positions:
        Cartesian coordinates of shape ``(n_atoms, 3)`` in Angstrom.
    elements:
        Element symbols aligned row-by-row with ``positions``, e.g.
        ``["C", "H", "Ni"]``.
    """

    positions: np.ndarray
    elements: list[str]

    def __post_init__(self) -> None:
        self.positions = np.asarray(self.positions, dtype=float)
        self.elements = list(self.elements)
        if self.positions.ndim != 2 or self.positions.shape[1] != 3:
            raise ValueError("positions must have shape (n_atoms, 3)")
        if self.positions.shape[0] != len(self.elements):
            raise ValueError(
                "positions and elements must contain the same atoms"
            )

    @property
    def n(self) -> int:
        return self.positions.shape[0]


@dataclass(frozen=True)
class JointLagrangeConfig:
    """Sampling, filtering, SLSQP, and quality-gate parameters.

    In region-exclusion mode, start generation proceeds in this order:

    ``sphere/orientation sampling`` -> ``candidate_pool_size single-
    substituent poses`` -> ``joint combinations bounded by beam_width`` ->
    ``preselect top-scoring combinations`` -> ``at most max_starts SLSQP
    starts``.

    ``candidate_pool_size``, ``beam_width``, ``preselect`` and ``max_starts``
    therefore sit at four distinct levels and cannot substitute for one
    another.

    Parameters
    ----------
    use_region_exclusion:
        Whether to check both the linking atom and the full substituent
        against the VDW exclusion boundary while generating starts. When
        disabled, a simpler baseline start scheme is used and adaptive
        sampling and feasible-domain pruning have no effect.
    use_adaptive_sampling:
        Whether to first run the ``coarse_*`` sampling; if the coarse stage
        does not produce an acceptable optimized structure, it then expands
        to the full sampling parameters. Usually avoids full start generation
        for simple systems.
    use_feasible_pruning:
        Whether to immediately drop partial combinations that incur
        substituent-substituent collisions as each substituent is added,
        using ``beam_width`` to bound how many are kept. When disabled, a
        ``u^K`` Cartesian product is used instead, where
        ``u = unpruned_candidates_per_sub``.
    use_greedy_selection:
        Whether to preferentially select geometrically diverse starts from
        the ``preselect`` top-scoring joint combinations. When disabled, the
        highest-scoring ``max_starts`` combinations are taken directly.
    use_convergence_detection:
        Whether to attach a callback to SLSQP that ends the current start
        early once the objective reaches a plateau and the constraints are
        already satisfied.
    use_early_stop:
        Whether to stop the remaining multi-start computation after the first
        acceptable result is obtained.
    use_quality_early_stop:
        Whether a result must pass the structural quality gate before it is
        accepted and triggers early stop. When disabled, any numerically
        feasible SLSQP result can trigger early stop.
    use_slsqp_fallback:
        When none of the callback-SLSQP starts are accepted, whether to rerun
        ordinary SLSQP (without the callback) on the same starts. Only
        meaningful when convergence detection is enabled.
    n_link_sphere:
        Full-stage number of linking-atom positions sampled on each
        skeleton connection atom's bond sphere.
    n_orientation_sphere:
        Full-stage number of target directions for the substituent principal
        axis on the sphere.
    n_axial:
        Number of uniform rotations added about each principal-axis
        direction. The per-substituent raw-pose upper bound is about
        ``n_link_sphere * n_orientation_sphere * n_axial``.
    candidate_pool_size:
        Maximum number ``M`` of high-scoring, geometrically diverse poses
        kept per substituent after skeleton region exclusion. This is the
        "single-substituent candidate pool size", not the number of SLSQP
        starts. Increasing ``M`` improves the chance of finding a joint
        feasible combination, but the pairwise substituent check grows
        roughly as ``O(K²M²)`` and the subsequent combination search also
        slows; too small may discard an orientation ultimately needed.
    preselect:
        The top ``P`` combinations, after the full joint combinations are
        sorted by score, fed into greedy diversity selection. At most
        ``max_starts`` SLSQP starts are still returned. Increasing it lets
        greedy selection see more candidate geometries; decreasing it speeds
        up slightly, but ``P < max_starts`` directly caps the number of
        starts.
    beam_width:
        Maximum number ``B`` of high-scoring "partial joint combinations"
        kept after each substituent is added during feasible-domain pruning.
        It truncates the theoretical ``M^K`` combinatorial explosion into a
        bounded beam search. Increasing it lowers the risk of pruning a
        potentially good combination too early, but time and memory grow
        roughly with ``B``; too small may leave no complete joint start.
        Effective only when ``use_feasible_pruning=True``.
    max_starts:
        Maximum number of full ``6K``-dimensional joint starts fed into SLSQP
        after start selection. Increasing it usually improves robustness, but
        the worst-case solve time increases roughly linearly; not necessarily
        all are executed when early stop is enabled.
    coarse_n_link_sphere:
        Linking-atom sphere sampling count for the adaptive coarse stage.
    coarse_n_orientation_sphere:
        Principal-axis direction count for the adaptive coarse stage.
    coarse_n_axial:
        Axial rotations per principal-axis direction for the adaptive coarse
        stage.
    coarse_candidate_pool_size:
        Per-substituent candidate pool cap for the adaptive coarse stage;
        same meaning as ``candidate_pool_size``.
    coarse_preselect:
        Combination cap entering joint greedy selection for the adaptive
        coarse stage; same meaning as ``preselect``.
    coarse_beam_width:
        Beam width for the adaptive coarse stage; same meaning as
        ``beam_width``.
    baseline_top_m:
        Effective only when ``use_region_exclusion=False``. The baseline
        generates and keeps the top ``m`` orientations per substituent, then
        builds ``m^K`` joint combinations. The combination count grows
        exponentially as it increases; the default region-exclusion path
        never uses it.
    unpruned_candidates_per_sub:
        Effective only when region exclusion is on but
        ``use_feasible_pruning=False``. Takes the top ``u`` candidates from
        each substituent pool and fully enumerates ``u^K`` combinations, then
        checks feasibility one by one. Used for an unpruned control or
        small-K systems; e.g. ``u=4, K=6`` already gives 4096 combinations.
    candidate_rms_radius:
        Pose RMS-distance threshold within a single-substituent candidate
        pool, in Angstrom. Larger values make candidates more spread out but
        may not fill the pool; smaller values keep more mutually similar
        high-scoring poses.
    greedy_rms_radius:
        RMS-distance threshold between full joint starts, in Angstrom, used
        for greedy diversity selection.
    region_margin:
        Extra minimum clearance required by region exclusion, in Angstrom; 0
        requires only not crossing the current VDW constraint boundary, a
        positive value is stricter.
    vdw_constraint_scale:
        Dimensionless scale in the non-bonded minimum distance
        ``scale * (r_i + r_j)``, used for both the initial region exclusion
        and the SLSQP inequality constraints.
    softmin_alpha:
        Smoothing coefficient of the soft-min distance objective, in units of
        about Å⁻¹. Larger is closer to the hard minimum but gives a sharper
        objective surface; smaller is smoother but more averaged.
    wca_epsilon:
        Energy scale of the WCA repulsive potential; together with
        ``wca_weight`` it determines the repulsion contribution.
    wca_sigma_scale:
        The ``sigma_ij = scale * (r_i + r_j)`` scale for WCA.
    wca_weight:
        Multiplicative weight when adding the WCA energy into the soft-min
        objective. Too large may make the optimization care only about
        repulsion; too small makes it hard to keep structures away from the
        constraint boundary.
    slsqp_maxiter:
        Maximum number of iterations per SLSQP start.
    slsqp_ftol:
        SciPy SLSQP objective-convergence tolerance; smaller is usually more
        accurate but slower.
    equality_tolerance:
        Maximum squared-distance residual allowed when accepting the
        linking-bond sphere equality constraints, in Å².
    inequality_tolerance:
        Maximum absolute value of a negative squared-distance residual
        allowed when accepting the non-bonded inequality constraints, in Å².
    accept_feasible_slsqp:
        When SciPy reports ``success=False``, whether to still accept a
        candidate whose constraints satisfy the tolerances above.
    detection_patience:
        Number of most-recent objective values the callback plateau
        detection observes.
    detection_min_iterations:
        Minimum number of callback iterations that must run before a plateau
        convergence may be declared.
    detection_tolerance:
        Maximum allowed ratio of the recent objective fluctuation relative to
        the current objective scale.
    quality_max_close_20:
        Number of atom pairs with non-bonded distance below 2.00 Å that the
        quality gate allows.
    quality_max_vdw75_overlap:
        Upper bound on the total VDW75 overlap the quality gate allows, in
        Angstrom.
    quality_min_nonbond:
        Minimum non-bonded atom distance the quality gate requires, in
        Angstrom.
    """

    use_region_exclusion: bool = True
    use_adaptive_sampling: bool = True
    use_feasible_pruning: bool = True
    use_greedy_selection: bool = True
    use_convergence_detection: bool = True
    use_early_stop: bool = True
    use_quality_early_stop: bool = True
    use_slsqp_fallback: bool = True

    n_link_sphere: int = 48
    n_orientation_sphere: int = 24
    n_axial: int = 4
    candidate_pool_size: int = 20
    preselect: int = 48
    beam_width: int = 4096
    max_starts: int = 8

    coarse_n_link_sphere: int = 24
    coarse_n_orientation_sphere: int = 12
    coarse_n_axial: int = 4
    coarse_candidate_pool_size: int = 12
    coarse_preselect: int = 32
    coarse_beam_width: int = 2048

    baseline_top_m: int = 2
    unpruned_candidates_per_sub: int = 4
    candidate_rms_radius: float = 0.15
    greedy_rms_radius: float = 0.35
    region_margin: float = 0.0

    vdw_constraint_scale: float = 0.48640625
    softmin_alpha: float = 5.0
    wca_epsilon: float = 1.0
    wca_sigma_scale: float = 0.884375
    wca_weight: float = 4.897903992977276e-4

    slsqp_maxiter: int = 200
    slsqp_ftol: float = 1e-6
    equality_tolerance: float = 1e-4
    inequality_tolerance: float = 1e-6
    accept_feasible_slsqp: bool = True

    detection_patience: int = 5
    detection_min_iterations: int = 9
    detection_tolerance: float = 2e-4

    quality_max_close_20: int = 0
    quality_max_vdw75_overlap: float = 0.793427816544153
    quality_min_nonbond: float = 2.00

    def validate(self) -> None:
        integer_fields = {
            "n_link_sphere": self.n_link_sphere,
            "n_orientation_sphere": self.n_orientation_sphere,
            "n_axial": self.n_axial,
            "candidate_pool_size": self.candidate_pool_size,
            "preselect": self.preselect,
            "beam_width": self.beam_width,
            "max_starts": self.max_starts,
            "coarse_n_link_sphere": self.coarse_n_link_sphere,
            "coarse_n_orientation_sphere": self.coarse_n_orientation_sphere,
            "coarse_n_axial": self.coarse_n_axial,
            "coarse_candidate_pool_size": self.coarse_candidate_pool_size,
            "coarse_preselect": self.coarse_preselect,
            "coarse_beam_width": self.coarse_beam_width,
            "baseline_top_m": self.baseline_top_m,
            "unpruned_candidates_per_sub": self.unpruned_candidates_per_sub,
            "slsqp_maxiter": self.slsqp_maxiter,
            "detection_patience": self.detection_patience,
            "detection_min_iterations": self.detection_min_iterations,
        }
        invalid = [
            name for name, value in integer_fields.items() if value <= 0
        ]
        if invalid:
            raise ValueError(
                f"positive integer parameters required: {invalid}"
            )
        if self.vdw_constraint_scale <= 0.0:
            raise ValueError("vdw_constraint_scale must be positive")
        if self.softmin_alpha <= 0.0:
            raise ValueError("softmin_alpha must be positive")
        if self.wca_epsilon < 0.0 or self.wca_weight < 0.0:
            raise ValueError("WCA epsilon and weight cannot be negative")
        if self.slsqp_ftol <= 0.0:
            raise ValueError("slsqp_ftol must be positive")
        if self.equality_tolerance <= 0.0 or self.inequality_tolerance < 0.0:
            raise ValueError("constraint tolerances are invalid")


@dataclass
class StartBuildStats:
    """Diagnostic statistics for the start-generation stage.

    ``per_sub_*`` records, per substituent, the sampling counts, region-
    exclusion pass counts, and final candidate pool size;
    ``joint_combinations_*`` records the joint-search size;
    ``selected_starts`` is the number finally handed to SLSQP, while
    ``initial_feasible_starts`` is how many of those already satisfy all
    equality and inequality tolerances at the initial point.
    """

    seed_time_s: float = 0.0
    adaptive_levels: int = 1
    adaptive_expanded: bool = False
    per_sub_link_samples: list[int] = field(default_factory=list)
    per_sub_link_pass: list[int] = field(default_factory=list)
    per_sub_pose_samples: list[int] = field(default_factory=list)
    per_sub_pose_pass: list[int] = field(default_factory=list)
    per_sub_pool_size: list[int] = field(default_factory=list)
    joint_combinations_examined: int = 0
    joint_combinations_feasible: int = 0
    generated_starts: int = 0
    selected_starts: int = 0
    initial_feasible_starts: int = 0


@dataclass
class OptimizationStats:
    """Start, solve, repair, and fallback statistics for one run."""

    build: StartBuildStats
    solve_time_s: float
    attempted_starts: int
    successful_starts: int
    accepted_successes: int
    detected_early_starts: int
    stopped_early: bool
    fallback_triggered: bool
    fallback_attempted_starts: int
    fallback_successful_starts: int
    fallback_accepted_successes: int
    quality_repair_triggered: bool
    quality_repair_attempted_starts: int
    quality_repair_accepted_successes: int
    total_iterations: int
    total_function_evaluations: int
    best_equality_error: Optional[float]
    best_inequality_slack: Optional[float]
    messages: dict[str, int] = field(default_factory=dict)

    @property
    def total_time_s(self) -> float:
        return self.build.seed_time_s + self.solve_time_s

    @property
    def raw_failure_rate(self) -> float:
        if self.attempted_starts == 0:
            return 1.0
        return 1.0 - self.successful_starts / self.attempted_starts


@dataclass
class JointLagrangeResult:
    """Final return value of Joint-Lagrange.

    ``success`` means at least one acceptable structure is returned;
    ``quality_ok`` indicates whether the chosen structure passes the quality
    gate; ``final`` orders atoms with the skeleton first and each substituent
    appended in input order; ``substituent_atom_ranges`` gives each
    substituent's half-open index range within ``final``.
    """

    success: bool
    final: Optional[Mol]
    n_skeleton_atoms: int
    substituent_atom_ranges: dict[int, tuple[int, int]]
    objective_function_value: float
    quality_ok: bool
    quality_metrics: dict[str, float]
    stats: OptimizationStats
    message: str


@dataclass(slots=True)
class _SamplingLevel:
    n_link_sphere: int
    n_orientation_sphere: int
    n_axial: int
    candidate_pool_size: int
    preselect: int
    beam_width: int


@dataclass(slots=True)
class _SeedCandidate:
    linking_position: np.ndarray
    euler: np.ndarray
    positions: np.ndarray
    score: float
    skeleton_clearance: float
    link_clearance: float


@dataclass(slots=True)
class _StartRecord:
    x0: np.ndarray
    score: float
    minimum_inequality: float
    feature: np.ndarray


@dataclass(slots=True)
class _BeamState:
    candidates: tuple[_SeedCandidate, ...]
    score: float
    minimum_substituent_clearance: float
    combination_indices: tuple[int, ...]


@dataclass(slots=True)
class _CandidateSolution:
    x: np.ndarray
    objective: float
    raw_ok: bool
    equality_error: float
    inequality_slack: float
    message: str
    iterations: int
    function_evaluations: int
    detected_early: bool


class _DetectedConvergence(Exception):
    pass


class _ConvergenceDetector:
    def __init__(self, optimizer: "JointLagrangeOptimizer"):
        self.optimizer = optimizer
        self.values: list[float] = []
        self.last_x: Optional[np.ndarray] = None

    def __call__(self, x: np.ndarray) -> None:
        config = self.optimizer.config
        self.last_x = np.asarray(x, dtype=float).copy()
        self.values.append(float(self.optimizer._objective(self.last_x)))
        required = max(
            config.detection_patience,
            config.detection_min_iterations,
        )
        if len(self.values) < required:
            return
        recent = np.asarray(self.values[-config.detection_patience :])
        scale = 1.0 + abs(float(recent[-1]))
        if float(np.ptp(recent)) > config.detection_tolerance * scale:
            return
        equality_error = float(
            np.max(np.abs(self.optimizer._equality_constraints(self.last_x)))
        )
        inequality_slack = float(
            self.optimizer._inequality_constraints(self.last_x).min()
        )
        if (
            equality_error <= config.equality_tolerance
            and inequality_slack >= -config.inequality_tolerance
        ):
            raise _DetectedConvergence


def fibonacci_sphere(count: int) -> np.ndarray:
    """Return ``count`` near-uniform Fibonacci directions on the sphere."""

    indices = np.arange(count, dtype=float)
    golden_ratio = (1.0 + np.sqrt(5.0)) / 2.0
    golden_angle = 2.0 * np.pi * (1.0 - 1.0 / golden_ratio)
    z = 1.0 - 2.0 * (indices + 0.5) / count
    radius = np.sqrt(np.maximum(0.0, 1.0 - z * z))
    theta = golden_angle * indices
    return np.stack(
        [radius * np.cos(theta), radius * np.sin(theta), z],
        axis=1,
    )


class JointLagrangeOptimizer:
    """Joint-Lagrange solver that optimizes K rigid substituents at once.

    Each substituent contributes its linking atom's three Cartesian
    coordinates and three intrinsic-xyz Euler angles, so the total SLSQP
    variable dimension is ``6K``. The linking atom is held on the covalent
    bond-length sphere by an equality constraint; substituent-skeleton and
    substituent-substituent collisions are handled by VDW inequality
    constraints.

    Parameters
    ----------
    skeleton:
        The skeleton, whose coordinates are held fixed.
    substituents:
        The K rigid substituents; the optimization only translates and
        rotates them, without changing their internal bond lengths or angles.
    skeleton_link_indices:
        The 0-based index of the skeleton connection atom for each
        substituent; length must equal K.
    substituent_link_indices:
        The 0-based index of the linking atom inside each substituent; length
        must equal K.
    config:
        Sampling, filtering, SLSQP, early-stop, and quality-gate
        configuration; ``JointLagrangeConfig()`` is used when omitted.
    covalent_radii:
        Optional per-element covalent-radius override dict, in Angstrom;
        elements not provided use the full RDKit table. It only affects the
        linking-bond equality-constraint target bond length.
    vdw_radii:
        Optional per-element VDW-radius override dict, in Angstrom; elements
        not provided use the full RDKit table. It affects region exclusion,
        SLSQP non-bonded constraints, WCA, and the VDW quality metrics.
    """

    def __init__(
        self,
        skeleton: Mol,
        substituents: list[Mol],
        skeleton_link_indices: list[int],
        substituent_link_indices: list[int],
        config: Optional[JointLagrangeConfig] = None,
        covalent_radii: Optional[dict[str, float]] = None,
        vdw_radii: Optional[dict[str, float]] = None,
    ):
        self.config = config or JointLagrangeConfig()
        self.config.validate()
        if not (
            len(substituents)
            == len(skeleton_link_indices)
            == len(substituent_link_indices)
        ):
            raise ValueError("substituents and linking-index lists must match")
        if not substituents:
            raise ValueError("at least one substituent is required")

        self.skeleton = skeleton
        self.substituents = substituents
        self.skeleton_link_indices = list(skeleton_link_indices)
        self.substituent_link_indices = list(substituent_link_indices)
        self.k_substituents = len(substituents)
        self.covalent_radii = dict(COVALENT_RADII)
        self.vdw_radii = dict(VDW_RADII)
        if covalent_radii:
            self.covalent_radii.update(covalent_radii)
        if vdw_radii:
            self.vdw_radii.update(vdw_radii)
        self._validate_inputs()

        self.skeleton_positions = np.asarray(skeleton.positions, dtype=float)
        self.skeleton_covalent = self._radii(
            skeleton.elements, self.covalent_radii
        )
        self.skeleton_vdw = self._radii(skeleton.elements, self.vdw_radii)
        self.substituent_covalent = [
            self._radii(substituent.elements, self.covalent_radii)
            for substituent in substituents
        ]
        self.substituent_vdw = [
            self._radii(substituent.elements, self.vdw_radii)
            for substituent in substituents
        ]
        self.relative_positions = [
            substituent.positions - substituent.positions[link_index]
            for substituent, link_index in zip(
                substituents,
                self.substituent_link_indices,
            )
        ]
        self.bond_distances = [
            self.substituent_covalent[k][self.substituent_link_indices[k]]
            + self.skeleton_covalent[self.skeleton_link_indices[k]]
            for k in range(self.k_substituents)
        ]
        self.minimum_squared_sub_skeleton = [
            (
                (self.substituent_vdw[k][:, None] + self.skeleton_vdw[None, :])
                * self.config.vdw_constraint_scale
            )
            ** 2
            for k in range(self.k_substituents)
        ]
        self.sub_skeleton_masks = []
        for k, substituent in enumerate(substituents):
            mask = np.ones(
                (substituent.n, skeleton.n),
                dtype=bool,
            )
            mask[
                self.substituent_link_indices[k],
                self.skeleton_link_indices[k],
            ] = False
            self.sub_skeleton_masks.append(mask)
        self.minimum_squared_sub_sub = {
            (left, right): (
                (
                    self.substituent_vdw[left][:, None]
                    + self.substituent_vdw[right][None, :]
                )
                * self.config.vdw_constraint_scale
            )
            ** 2
            for left, right in itertools.combinations(
                range(self.k_substituents),
                2,
            )
        }
        nonbond_vdw_chunks = [
            (self.substituent_vdw[k][:, None] + self.skeleton_vdw[None, :])[
                self.sub_skeleton_masks[k]
            ]
            for k in range(self.k_substituents)
        ]
        nonbond_vdw_chunks.extend(
            (
                self.substituent_vdw[left][:, None]
                + self.substituent_vdw[right][None, :]
            ).reshape(-1)
            for left, right in itertools.combinations(
                range(self.k_substituents),
                2,
            )
        )
        self._nonbond_vdw_sums = np.concatenate(nonbond_vdw_chunks)
        self.wca_sigma_sub_skeleton = [
            (self.substituent_vdw[k][:, None] + self.skeleton_vdw[None, :])
            * self.config.wca_sigma_scale
            for k in range(self.k_substituents)
        ]
        self.wca_sigma_sub_sub = {
            (left, right): (
                self.substituent_vdw[left][:, None]
                + self.substituent_vdw[right][None, :]
            )
            * self.config.wca_sigma_scale
            for left, right in itertools.combinations(
                range(self.k_substituents),
                2,
            )
        }
        self.principal_axes = []
        for relative in self.relative_positions:
            center = relative.mean(axis=0)
            norm = np.linalg.norm(center)
            self.principal_axes.append(
                center / norm if norm > 1e-3 else np.array([0.0, 0.0, 1.0])
            )

    @staticmethod
    def _radii(elements: list[str], table: dict[str, float]) -> np.ndarray:
        missing = sorted(set(elements) - set(table))
        if missing:
            raise ValueError(f"missing atomic radii for elements: {missing}")
        return np.asarray(
            [table[element] for element in elements], dtype=float
        )

    def _validate_inputs(self) -> None:
        for index in self.skeleton_link_indices:
            if index < 0 or index >= self.skeleton.n:
                raise IndexError(
                    f"skeleton linking index out of range: {index}"
                )
        for substituent, index in zip(
            self.substituents,
            self.substituent_link_indices,
        ):
            if index < 0 or index >= substituent.n:
                raise IndexError(
                    f"substituent linking index out of range: {index}"
                )
        all_elements = list(self.skeleton.elements)
        for substituent in self.substituents:
            all_elements.extend(substituent.elements)
        missing_covalent = sorted(set(all_elements) - set(self.covalent_radii))
        missing_vdw = sorted(set(all_elements) - set(self.vdw_radii))
        if missing_covalent or missing_vdw:
            raise ValueError(
                "missing radii: "
                f"covalent={missing_covalent}, vdw={missing_vdw}"
            )

    def _pack(self, positions: np.ndarray, eulers: np.ndarray) -> np.ndarray:
        return np.concatenate([positions.reshape(-1), eulers.reshape(-1)])

    def _unpack(self, x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        split = 3 * self.k_substituents
        positions = x[:split].reshape(self.k_substituents, 3)
        eulers = x[split : 2 * split].reshape(self.k_substituents, 3)
        return positions, eulers

    def _substituent_positions(self, x: np.ndarray) -> list[np.ndarray]:
        linking_positions, eulers = self._unpack(x)
        positions = []
        for k in range(self.k_substituents):
            rotation = Rotation.from_euler("xyz", eulers[k]).as_matrix()
            positions.append(
                linking_positions[k] + self.relative_positions[k] @ rotation.T
            )
        return positions

    def _all_nonbond_distances(self, x: np.ndarray) -> np.ndarray:
        positions = self._substituent_positions(x)
        chunks = []
        for k in range(self.k_substituents):
            distances = np.linalg.norm(
                positions[k][:, None, :] - self.skeleton_positions[None, :, :],
                axis=-1,
            )
            chunks.append(distances[self.sub_skeleton_masks[k]])
        for left, right in itertools.combinations(
            range(self.k_substituents),
            2,
        ):
            distances = np.linalg.norm(
                positions[left][:, None, :] - positions[right][None, :, :],
                axis=-1,
            )
            chunks.append(distances.reshape(-1))
        return np.concatenate(chunks)

    def _wca_energy(self, x: np.ndarray) -> float:
        positions = self._substituent_positions(x)
        energy = 0.0
        epsilon = self.config.wca_epsilon
        for k in range(self.k_substituents):
            distances = np.linalg.norm(
                positions[k][:, None, :] - self.skeleton_positions[None, :, :],
                axis=-1,
            )
            sigma = self.wca_sigma_sub_skeleton[k]
            active = self.sub_skeleton_masks[k] & (
                distances <= sigma * _TWO_POW_ONE_SIXTH
            )
            if np.any(active):
                safe_distances = np.maximum(distances[active], 1e-8)
                ratio = sigma[active] / safe_distances
                ratio_six = ratio**6
                energy += float(
                    np.sum(
                        4.0 * epsilon * (ratio_six**2 - ratio_six) + epsilon
                    )
                )
        for left, right in itertools.combinations(
            range(self.k_substituents),
            2,
        ):
            distances = np.linalg.norm(
                positions[left][:, None, :] - positions[right][None, :, :],
                axis=-1,
            )
            sigma = self.wca_sigma_sub_sub[(left, right)]
            active = distances <= sigma * _TWO_POW_ONE_SIXTH
            if np.any(active):
                safe_distances = np.maximum(distances[active], 1e-8)
                ratio = sigma[active] / safe_distances
                ratio_six = ratio**6
                energy += float(
                    np.sum(
                        4.0 * epsilon * (ratio_six**2 - ratio_six) + epsilon
                    )
                )
        return energy

    def _objective(self, x: np.ndarray) -> float:
        distances = self._all_nonbond_distances(x)
        soft_minimum = float(
            -logsumexp(-self.config.softmin_alpha * distances)
            / self.config.softmin_alpha
        )
        return -soft_minimum + self.config.wca_weight * self._wca_energy(x)

    def _equality_constraints(self, x: np.ndarray) -> np.ndarray:
        linking_positions, _ = self._unpack(x)
        constraints = np.empty(self.k_substituents)
        for k in range(self.k_substituents):
            displacement = (
                linking_positions[k]
                - self.skeleton_positions[self.skeleton_link_indices[k]]
            )
            constraints[k] = (
                displacement @ displacement - self.bond_distances[k] ** 2
            )
        return constraints

    def _inequality_constraints(self, x: np.ndarray) -> np.ndarray:
        positions = self._substituent_positions(x)
        chunks = []
        for k in range(self.k_substituents):
            squared_distances = np.sum(
                (
                    positions[k][:, None, :]
                    - self.skeleton_positions[None, :, :]
                )
                ** 2,
                axis=-1,
            )
            chunks.append(
                (squared_distances - self.minimum_squared_sub_skeleton[k])[
                    self.sub_skeleton_masks[k]
                ]
            )
        for left, right in itertools.combinations(
            range(self.k_substituents),
            2,
        ):
            squared_distances = np.sum(
                (positions[left][:, None, :] - positions[right][None, :, :])
                ** 2,
                axis=-1,
            )
            chunks.append(
                (
                    squared_distances
                    - self.minimum_squared_sub_sub[(left, right)]
                ).reshape(-1)
            )
        return np.concatenate(chunks)

    @staticmethod
    def _rms_feature_distance(left: np.ndarray, right: np.ndarray) -> float:
        return float(np.sqrt(np.mean((left - right) ** 2)))

    def _link_clearance(self, k: int, position: np.ndarray) -> float:
        link_index = self.substituent_link_indices[k]
        distances = np.linalg.norm(
            position[None, :] - self.skeleton_positions,
            axis=1,
        )
        required = np.sqrt(self.minimum_squared_sub_skeleton[k][link_index])
        mask = self.sub_skeleton_masks[k][link_index]
        return float((distances - required)[mask].min())

    def _skeleton_clearance(self, k: int, positions: np.ndarray) -> float:
        distances = np.linalg.norm(
            positions[:, None, :] - self.skeleton_positions[None, :, :],
            axis=-1,
        )
        required = np.sqrt(self.minimum_squared_sub_skeleton[k])
        return float((distances - required)[self.sub_skeleton_masks[k]].min())

    def _substituent_clearance(
        self,
        left_k: int,
        left_positions: np.ndarray,
        right_k: int,
        right_positions: np.ndarray,
    ) -> float:
        key = (left_k, right_k) if left_k < right_k else (right_k, left_k)
        required = np.sqrt(self.minimum_squared_sub_sub[key])
        if left_k > right_k:
            required = required.T
        distances = np.linalg.norm(
            left_positions[:, None, :] - right_positions[None, :, :],
            axis=-1,
        )
        return float((distances - required).min())

    def _select_diverse_candidates(
        self,
        candidates: list[_SeedCandidate],
        count: int,
    ) -> list[_SeedCandidate]:
        ordered = sorted(
            candidates,
            key=lambda candidate: (
                candidate.score,
                candidate.skeleton_clearance,
                candidate.link_clearance,
            ),
            reverse=True,
        )
        selected: list[_SeedCandidate] = []
        selected_features: list[np.ndarray] = []
        for candidate in ordered:
            feature = candidate.positions.reshape(-1)
            if all(
                self._rms_feature_distance(feature, previous)
                >= self.config.candidate_rms_radius
                for previous in selected_features
            ):
                selected.append(candidate)
                selected_features.append(feature)
            if len(selected) >= count:
                return selected
        selected_ids = {id(candidate) for candidate in selected}
        for candidate in ordered:
            if id(candidate) not in selected_ids:
                selected.append(candidate)
            if len(selected) >= count:
                break
        return selected

    def _candidate_pool(
        self,
        k: int,
        level: _SamplingLevel,
    ) -> tuple[list[_SeedCandidate], dict[str, int]]:
        link_directions = fibonacci_sphere(level.n_link_sphere)
        orientation_directions = fibonacci_sphere(level.n_orientation_sphere)
        skeleton_link = self.skeleton_positions[self.skeleton_link_indices[k]]
        axis = self.principal_axes[k]
        relative = self.relative_positions[k]
        rotations: list[Rotation] = []
        for direction in orientation_directions:
            alignment, _ = Rotation.align_vectors(
                direction.reshape(1, 3),
                axis.reshape(1, 3),
            )
            for axial_index in range(level.n_axial):
                angle = axial_index * 2.0 * np.pi / level.n_axial
                spin = Rotation.from_rotvec(angle * axis)
                rotations.append(alignment * spin)

        candidates = []
        link_pass = 0
        pose_samples = 0
        pose_pass = 0
        for direction in link_directions:
            linking_position = (
                skeleton_link + direction * self.bond_distances[k]
            )
            link_clearance = self._link_clearance(k, linking_position)
            if link_clearance < self.config.region_margin:
                continue
            link_pass += 1
            for rotation in rotations:
                pose_samples += 1
                positions = (
                    linking_position + relative @ rotation.as_matrix().T
                )
                skeleton_clearance = self._skeleton_clearance(k, positions)
                if skeleton_clearance < self.config.region_margin:
                    continue
                pose_pass += 1
                distance_sum = float(
                    np.linalg.norm(
                        positions[:, None, :]
                        - self.skeleton_positions[None, :, :],
                        axis=-1,
                    ).sum()
                )
                score = (
                    30.0 * skeleton_clearance
                    + 5.0 * link_clearance
                    + 0.002 * distance_sum
                )
                candidates.append(
                    _SeedCandidate(
                        linking_position=linking_position.copy(),
                        euler=rotation.as_euler("xyz"),
                        positions=positions,
                        score=score,
                        skeleton_clearance=skeleton_clearance,
                        link_clearance=link_clearance,
                    )
                )
        selected = self._select_diverse_candidates(
            candidates,
            level.candidate_pool_size,
        )
        return selected, {
            "link_samples": level.n_link_sphere,
            "link_pass": link_pass,
            "pose_samples": pose_samples,
            "pose_pass": pose_pass,
            "pool_size": len(selected),
        }

    def _build_pruned_records(
        self,
        pools: list[list[_SeedCandidate]],
        level: _SamplingLevel,
    ) -> tuple[list[_StartRecord], int, int]:
        pair_clearances = {}
        for left, right in itertools.combinations(
            range(self.k_substituents),
            2,
        ):
            pair_clearances[(left, right)] = np.asarray(
                [
                    [
                        self._substituent_clearance(
                            left,
                            left_candidate.positions,
                            right,
                            right_candidate.positions,
                        )
                        for right_candidate in pools[right]
                    ]
                    for left_candidate in pools[left]
                ]
            )

        states = [
            _BeamState(
                candidates=(),
                score=0.0,
                minimum_substituent_clearance=float("inf"),
                combination_indices=(),
            )
        ]
        examined = 0
        for k, pool in enumerate(pools):
            next_states = []
            for state in states:
                for candidate_index, candidate in enumerate(pool):
                    examined += 1
                    minimum_clearance = state.minimum_substituent_clearance
                    feasible = True
                    for previous_k, previous_index in enumerate(
                        state.combination_indices
                    ):
                        clearance = pair_clearances[(previous_k, k)][
                            previous_index,
                            candidate_index,
                        ]
                        minimum_clearance = min(
                            minimum_clearance,
                            clearance,
                        )
                        if clearance < self.config.region_margin:
                            feasible = False
                            break
                    if not feasible:
                        continue
                    score = state.score + candidate.score
                    if np.isfinite(minimum_clearance):
                        score += 20.0 * minimum_clearance
                    next_states.append(
                        _BeamState(
                            candidates=state.candidates + (candidate,),
                            score=score,
                            minimum_substituent_clearance=minimum_clearance,
                            combination_indices=(
                                state.combination_indices + (candidate_index,)
                            ),
                        )
                    )
            next_states.sort(
                key=lambda state: (
                    state.score,
                    state.minimum_substituent_clearance,
                ),
                reverse=True,
            )
            states = next_states[: level.beam_width]
            if not states:
                break

        records = []
        for state in states:
            if len(state.candidates) != self.k_substituents:
                continue
            positions = np.asarray(
                [candidate.linking_position for candidate in state.candidates]
            )
            eulers = np.asarray(
                [candidate.euler for candidate in state.candidates]
            )
            x0 = self._pack(positions, eulers)
            minimum_inequality = float(self._inequality_constraints(x0).min())
            if minimum_inequality < -self.config.inequality_tolerance:
                continue
            feature = np.concatenate(
                [
                    candidate.positions.reshape(-1)
                    for candidate in state.candidates
                ]
            )
            records.append(
                _StartRecord(
                    x0=x0,
                    score=state.score,
                    minimum_inequality=minimum_inequality,
                    feature=feature,
                )
            )
        records.sort(
            key=lambda record: (record.score, record.minimum_inequality),
            reverse=True,
        )
        return records, examined, len(records)

    def _build_unpruned_records(
        self,
        pools: list[list[_SeedCandidate]],
    ) -> tuple[list[_StartRecord], int, int]:
        limit = self.config.unpruned_candidates_per_sub
        limited = [pool[: min(limit, len(pool))] for pool in pools]
        cartesian_size = int(
            np.prod([len(pool) for pool in limited], dtype=int)
        )
        records = []
        for combination_indices in itertools.product(
            *[range(len(pool)) for pool in limited]
        ):
            candidates = [
                limited[k][candidate_index]
                for k, candidate_index in enumerate(combination_indices)
            ]
            positions = np.asarray(
                [candidate.linking_position for candidate in candidates]
            )
            eulers = np.asarray([candidate.euler for candidate in candidates])
            x0 = self._pack(positions, eulers)
            minimum_inequality = float(self._inequality_constraints(x0).min())
            if minimum_inequality < -self.config.inequality_tolerance:
                continue
            records.append(
                _StartRecord(
                    x0=x0,
                    score=float(
                        sum(candidate.score for candidate in candidates)
                    ),
                    minimum_inequality=minimum_inequality,
                    feature=np.concatenate(
                        [
                            candidate.positions.reshape(-1)
                            for candidate in candidates
                        ]
                    ),
                )
            )
        records.sort(
            key=lambda record: (record.score, record.minimum_inequality),
            reverse=True,
        )
        return records, cartesian_size, len(records)

    def _select_start_records(
        self,
        records: list[_StartRecord],
        preselect: int,
    ) -> tuple[list[_StartRecord], int]:
        records = records[: min(preselect, len(records))]
        if not self.config.use_greedy_selection:
            selected = records[: self.config.max_starts]
        else:
            selected = []
            for record in records:
                if all(
                    self._rms_feature_distance(
                        record.feature,
                        previous.feature,
                    )
                    >= self.config.greedy_rms_radius
                    for previous in selected
                ):
                    selected.append(record)
                if len(selected) >= self.config.max_starts:
                    break
            selected_ids = {id(record) for record in selected}
            for record in records:
                if len(selected) >= self.config.max_starts:
                    break
                if id(record) not in selected_ids:
                    selected.append(record)
        initial_feasible = sum(
            record.minimum_inequality >= -self.config.inequality_tolerance
            and float(np.max(np.abs(self._equality_constraints(record.x0))))
            <= self.config.equality_tolerance
            for record in selected
        )
        return selected, initial_feasible

    def _build_region_level(
        self,
        level: _SamplingLevel,
    ) -> tuple[list[_StartRecord], StartBuildStats]:
        stats = StartBuildStats()
        pools = []
        for k in range(self.k_substituents):
            pool, sub_stats = self._candidate_pool(k, level)
            pools.append(pool)
            stats.per_sub_link_samples.append(sub_stats["link_samples"])
            stats.per_sub_link_pass.append(sub_stats["link_pass"])
            stats.per_sub_pose_samples.append(sub_stats["pose_samples"])
            stats.per_sub_pose_pass.append(sub_stats["pose_pass"])
            stats.per_sub_pool_size.append(sub_stats["pool_size"])
        if any(not pool for pool in pools):
            return [], stats
        if self.config.use_feasible_pruning:
            records, examined, feasible = self._build_pruned_records(
                pools,
                level,
            )
        else:
            records, examined, feasible = self._build_unpruned_records(pools)
        stats.joint_combinations_examined = examined
        stats.joint_combinations_feasible = feasible
        stats.generated_starts = len(records)
        selected, initial_feasible = self._select_start_records(
            records,
            level.preselect,
        )
        stats.selected_starts = len(selected)
        stats.initial_feasible_starts = initial_feasible
        return selected, stats

    def _baseline_seeds_for_sub(
        self,
        k: int,
    ) -> list[tuple[np.ndarray, np.ndarray]]:
        skeleton_link = self.skeleton_positions[self.skeleton_link_indices[k]]
        skeleton_center = self.skeleton_positions.mean(axis=0)
        direction = skeleton_link - skeleton_center
        norm = np.linalg.norm(direction)
        direction = (
            direction / norm if norm > 1e-6 else np.array([1.0, 0.0, 0.0])
        )
        linking_position = skeleton_link + direction * self.bond_distances[k]
        other_link_positions = []
        for other_k in range(self.k_substituents):
            if other_k == k:
                continue
            other_skeleton_link = self.skeleton_positions[
                self.skeleton_link_indices[other_k]
            ]
            other_direction = other_skeleton_link - skeleton_center
            other_norm = np.linalg.norm(other_direction)
            other_direction = (
                other_direction / other_norm
                if other_norm > 1e-6
                else np.array([1.0, 0.0, 0.0])
            )
            other_link_positions.append(
                other_skeleton_link
                + other_direction * self.bond_distances[other_k]
            )
        other_links = (
            np.asarray(other_link_positions)
            if other_link_positions
            else np.zeros((0, 3))
        )
        sphere = fibonacci_sphere(self.config.n_orientation_sphere)
        axis = self.principal_axes[k]
        candidates = []
        for direction in sphere:
            alignment, _ = Rotation.align_vectors(
                direction.reshape(1, 3),
                axis.reshape(1, 3),
            )
            for axial_index in range(self.config.n_axial):
                angle = axial_index * 2.0 * np.pi / self.config.n_axial
                rotation = alignment * Rotation.from_rotvec(angle * axis)
                positions = (
                    linking_position
                    + self.relative_positions[k] @ rotation.as_matrix().T
                )
                squared_distances = np.sum(
                    (
                        positions[:, None, :]
                        - self.skeleton_positions[None, :, :]
                    )
                    ** 2,
                    axis=-1,
                )
                slack = float(
                    (squared_distances - self.minimum_squared_sub_skeleton[k])[
                        self.sub_skeleton_masks[k]
                    ].min()
                )
                distance_sum = float(np.sqrt(squared_distances).sum())
                feasibility_bonus = 0.0 if slack < 0.0 else 10.0 + slack
                anchor_penalty = 0.0
                if other_links.shape[0] > 0:
                    anchor_distances = np.linalg.norm(
                        positions[:, None, :] - other_links[None, :, :],
                        axis=-1,
                    )
                    anchor_penalty = float(
                        np.maximum(0.0, 1.5 - anchor_distances).sum()
                    )
                candidates.append(
                    (
                        distance_sum
                        + feasibility_bonus
                        - 5.0 * anchor_penalty,
                        linking_position.copy(),
                        rotation.as_euler("xyz"),
                    )
                )
        candidates.sort(key=lambda item: item[0], reverse=True)
        return [
            (position, euler)
            for _, position, euler in candidates[: self.config.baseline_top_m]
        ]

    def _build_baseline_starts(
        self,
    ) -> tuple[list[_StartRecord], StartBuildStats]:
        started = time.perf_counter()
        seeds = [
            self._baseline_seeds_for_sub(k) for k in range(self.k_substituents)
        ]
        records = []
        for combination in itertools.product(*seeds):
            positions = np.asarray([item[0] for item in combination])
            eulers = np.asarray([item[1] for item in combination])
            x0 = self._pack(positions, eulers)
            minimum_inequality = float(self._inequality_constraints(x0).min())
            score = -float(self._objective(x0)) + 100.0 * min(
                0.0, minimum_inequality
            )
            records.append(
                _StartRecord(
                    x0=x0,
                    score=score,
                    minimum_inequality=minimum_inequality,
                    feature=np.concatenate(
                        [
                            positions.reshape(-1)
                            for positions in self._substituent_positions(x0)
                        ]
                    ),
                )
            )
        records.sort(
            key=lambda record: (record.score, record.minimum_inequality),
            reverse=True,
        )
        selected, initial_feasible = self._select_start_records(
            records,
            len(records),
        )
        stats = StartBuildStats(
            seed_time_s=time.perf_counter() - started,
            per_sub_pool_size=[self.config.baseline_top_m]
            * self.k_substituents,
            joint_combinations_examined=len(records),
            joint_combinations_feasible=sum(
                record.minimum_inequality >= -self.config.inequality_tolerance
                for record in records
            ),
            generated_starts=len(records),
            selected_starts=len(selected),
            initial_feasible_starts=initial_feasible,
        )
        return selected, stats

    def _sampling_levels(self) -> list[_SamplingLevel]:
        full = _SamplingLevel(
            n_link_sphere=self.config.n_link_sphere,
            n_orientation_sphere=self.config.n_orientation_sphere,
            n_axial=self.config.n_axial,
            candidate_pool_size=self.config.candidate_pool_size,
            preselect=self.config.preselect,
            beam_width=self.config.beam_width,
        )
        if not self.config.use_adaptive_sampling:
            return [full]
        coarse = _SamplingLevel(
            n_link_sphere=self.config.coarse_n_link_sphere,
            n_orientation_sphere=self.config.coarse_n_orientation_sphere,
            n_axial=self.config.coarse_n_axial,
            candidate_pool_size=self.config.coarse_candidate_pool_size,
            preselect=self.config.coarse_preselect,
            beam_width=self.config.coarse_beam_width,
        )
        return [coarse, full]

    def _build_starts(
        self,
        level: Optional[_SamplingLevel] = None,
    ) -> tuple[list[_StartRecord], StartBuildStats]:
        """Build starts for one sampling level.

        Adaptive coarse-to-full orchestration belongs to :meth:`optimize`,
        because expansion depends on whether optimization produces an
        acceptable structure rather than on the number of starts alone.
        """
        if not self.config.use_region_exclusion:
            return self._build_baseline_starts()
        started = time.perf_counter()
        selected_level = (
            level if level is not None else self._sampling_levels()[0]
        )
        records, stats = self._build_region_level(selected_level)
        stats.seed_time_s = time.perf_counter() - started
        return records, stats

    def _make_result_molecule(
        self,
        x: np.ndarray,
    ) -> tuple[Mol, dict[int, tuple[int, int]]]:
        substituent_positions = self._substituent_positions(x)
        position_blocks = [self.skeleton_positions]
        elements = list(self.skeleton.elements)
        ranges = {}
        cursor = self.skeleton.n
        for k, positions in enumerate(substituent_positions):
            ranges[k] = (cursor, cursor + self.substituents[k].n)
            cursor += self.substituents[k].n
            position_blocks.append(positions)
            elements.extend(self.substituents[k].elements)
        return Mol(np.vstack(position_blocks), elements), ranges

    def _quality_metrics(
        self,
        molecule: Mol,
        ranges: dict[int, tuple[int, int]],
    ) -> dict[str, float]:
        metrics = {
            "distance_sum_nonbond": 0.0,
            "distance_mean_nonbond": 0.0,
            "min_sub_skeleton_nonbond": float("inf"),
            "min_sub_sub": float("inf"),
            "min_nonbond": float("inf"),
            "n_close_20": 0,
            "n_close_25": 0,
            "nonbond_pairs": 0,
            "vdw055_overlap_count": 0,
            "vdw055_overlap_sum": 0.0,
            "vdw055_max_overlap": 0.0,
            "vdw075_overlap_count": 0,
            "vdw075_overlap_sum": 0.0,
            "vdw075_max_overlap": 0.0,
            "vdw100_overlap_count": 0,
            "vdw100_overlap_sum": 0.0,
            "vdw100_max_overlap": 0.0,
            "link_bond_max_abs_error": 0.0,
            "link_bond_mean_abs_error": 0.0,
        }
        skeleton_indices = np.arange(self.skeleton.n)
        link_errors = []

        def accumulate(
            distances: np.ndarray,
            required_by_scale: dict[float, np.ndarray],
        ) -> None:
            metrics["distance_sum_nonbond"] += float(distances.sum())
            metrics["nonbond_pairs"] += int(distances.size)
            metrics["n_close_20"] += int((distances < 2.0).sum())
            metrics["n_close_25"] += int((distances < 2.5).sum())
            for scale, required in required_by_scale.items():
                overlap = np.maximum(0.0, required - distances)
                positive = overlap[overlap > 1e-9]
                prefix = f"vdw{int(scale * 100):03d}"
                metrics[f"{prefix}_overlap_count"] += int(positive.size)
                if positive.size:
                    metrics[f"{prefix}_overlap_sum"] += float(positive.sum())
                    metrics[f"{prefix}_max_overlap"] = max(
                        metrics[f"{prefix}_max_overlap"],
                        float(positive.max()),
                    )

        for k in range(self.k_substituents):
            start, end = ranges[k]
            sub_indices = np.arange(start, end)
            distances = np.linalg.norm(
                molecule.positions[sub_indices][:, None, :]
                - molecule.positions[skeleton_indices][None, :, :],
                axis=-1,
            )
            sub_vdw = self.substituent_vdw[k]
            required = {
                scale: (sub_vdw[:, None] + self.skeleton_vdw[None, :]) * scale
                for scale in (0.55, 0.75, 1.00)
            }
            mask = self.sub_skeleton_masks[k]
            nonbond = distances[mask]
            accumulate(
                nonbond,
                {scale: values[mask] for scale, values in required.items()},
            )
            metrics["min_sub_skeleton_nonbond"] = min(
                metrics["min_sub_skeleton_nonbond"],
                float(nonbond.min()),
            )
            link_distance = float(
                distances[
                    self.substituent_link_indices[k],
                    self.skeleton_link_indices[k],
                ]
            )
            link_errors.append(abs(link_distance - self.bond_distances[k]))

        for left, right in itertools.combinations(
            range(self.k_substituents),
            2,
        ):
            left_start, left_end = ranges[left]
            right_start, right_end = ranges[right]
            distances = np.linalg.norm(
                molecule.positions[left_start:left_end][:, None, :]
                - molecule.positions[right_start:right_end][None, :, :],
                axis=-1,
            )
            required = {
                scale: (
                    self.substituent_vdw[left][:, None]
                    + self.substituent_vdw[right][None, :]
                )
                * scale
                for scale in (0.55, 0.75, 1.00)
            }
            flat_distances = distances.reshape(-1)
            accumulate(
                flat_distances,
                {
                    scale: values.reshape(-1)
                    for scale, values in required.items()
                },
            )
            metrics["min_sub_sub"] = min(
                metrics["min_sub_sub"],
                float(flat_distances.min()),
            )

        metrics["distance_mean_nonbond"] = (
            metrics["distance_sum_nonbond"] / metrics["nonbond_pairs"]
        )
        metrics["min_nonbond"] = min(
            metrics["min_sub_skeleton_nonbond"],
            metrics["min_sub_sub"],
        )
        metrics["link_bond_max_abs_error"] = max(link_errors)
        metrics["link_bond_mean_abs_error"] = float(np.mean(link_errors))
        return metrics

    def _passes_quality(self, metrics: dict[str, float]) -> bool:
        return (
            int(metrics["n_close_20"]) <= self.config.quality_max_close_20
            and float(metrics["vdw075_overlap_sum"])
            <= self.config.quality_max_vdw75_overlap
            and float(metrics["min_nonbond"])
            >= self.config.quality_min_nonbond
        )

    def _quality_constraints(self, x: np.ndarray) -> np.ndarray:
        """Return constraints matching the exact structural quality gate.

        The per-pair constraints target the configured minimum non-bonded
        distance. The final scalar constrains the same aggregate VDW75
        overlap used by :meth:`_passes_quality`. With the default 2.00
        Angstrom minimum, the ``n_close_20 == 0`` requirement is covered by
        the distance constraints and is still verified by the exact gate.
        """
        distances = self._all_nonbond_distances(x)
        vdw75_overlap = float(
            np.maximum(
                0.0,
                0.75 * self._nonbond_vdw_sums - distances,
            ).sum()
        )
        return np.concatenate(
            [
                distances - self.config.quality_min_nonbond,
                np.asarray(
                    [self.config.quality_max_vdw75_overlap - vdw75_overlap]
                ),
            ]
        )

    def _quality_violation(self, solution: _CandidateSolution) -> float:
        """Rank rejected candidates by squared quality-constraint deficit."""
        deficits = np.minimum(0.0, self._quality_constraints(solution.x))
        return float(deficits @ deficits)

    def _solve_one(
        self,
        record: _StartRecord,
        use_detection: bool,
        enforce_quality: bool = False,
    ) -> _CandidateSolution:
        constraints = [
            {"type": "eq", "fun": self._equality_constraints},
            {"type": "ineq", "fun": self._inequality_constraints},
        ]
        if enforce_quality:
            constraints.append(
                {"type": "ineq", "fun": self._quality_constraints}
            )
        detector = _ConvergenceDetector(self) if use_detection else None
        detected = False
        try:
            result = minimize(
                self._objective,
                record.x0.copy(),
                method="SLSQP",
                constraints=constraints,
                callback=detector,
                options={
                    "ftol": self.config.slsqp_ftol,
                    "maxiter": self.config.slsqp_maxiter,
                    "disp": False,
                },
            )
            x = np.asarray(result.x, dtype=float)
            scipy_success = bool(result.success)
            message = str(result.message)
            iterations = int(getattr(result, "nit", 0))
            function_evaluations = int(getattr(result, "nfev", 0))
        except _DetectedConvergence:
            detected = True
            x = detector.last_x if detector is not None else record.x0.copy()
            scipy_success = False
            message = "convergence detected"
            iterations = len(detector.values) if detector is not None else 0
            function_evaluations = 0
        except Exception as error:
            return _CandidateSolution(
                x=record.x0.copy(),
                objective=float("inf"),
                raw_ok=False,
                equality_error=float("inf"),
                inequality_slack=-float("inf"),
                message=f"exception: {type(error).__name__}",
                iterations=0,
                function_evaluations=0,
                detected_early=False,
            )

        equality_error = float(np.max(np.abs(self._equality_constraints(x))))
        inequality_slack = float(self._inequality_constraints(x).min())
        feasible = (
            equality_error <= self.config.equality_tolerance
            and inequality_slack >= -self.config.inequality_tolerance
        )
        raw_ok = scipy_success or (
            self.config.accept_feasible_slsqp and feasible
        )
        return _CandidateSolution(
            x=x.copy(),
            objective=float(self._objective(x)),
            raw_ok=bool(raw_ok),
            equality_error=equality_error,
            inequality_slack=inequality_slack,
            message=message,
            iterations=iterations,
            function_evaluations=function_evaluations,
            detected_early=detected,
        )

    def optimize(self) -> JointLagrangeResult:
        """Generate starts, solve them, and repair the best quality rejects.

        Each adaptive sampling level receives at most one quality-repair
        attempt. A repair starts from the numerically feasible candidate
        closest to the exact quality constraints; repaired structures are
        accepted only through the unchanged exact quality gate.
        """

        attempted = 0
        successful = 0
        accepted = 0
        detected = 0
        stopped_early = False
        total_iterations = 0
        total_function_evaluations = 0
        messages: dict[str, int] = {}
        best_solution: Optional[_CandidateSolution] = None
        best_molecule: Optional[Mol] = None
        best_ranges: dict[int, tuple[int, int]] = {}
        best_metrics: dict[str, float] = {}
        best_quality_ok = False
        quality_rejected: list[_CandidateSolution] = []
        build_stats = StartBuildStats()
        build_time = 0.0
        total_examined = 0
        used_levels = 0
        solve_time = 0.0

        def consume(
            solution: _CandidateSolution,
            message_prefix: str = "",
            track_quality_rejection: bool = True,
        ) -> bool:
            nonlocal attempted, successful, accepted, detected
            nonlocal total_iterations, total_function_evaluations
            nonlocal best_solution, best_molecule, best_ranges
            nonlocal best_metrics, best_quality_ok
            attempted += 1
            message = f"{message_prefix}{solution.message}"
            messages[message] = messages.get(message, 0) + 1
            total_iterations += solution.iterations
            total_function_evaluations += solution.function_evaluations
            detected += int(solution.detected_early)
            if not solution.raw_ok:
                return False
            successful += 1
            molecule, ranges = self._make_result_molecule(solution.x)
            metrics = self._quality_metrics(molecule, ranges)
            quality_ok = self._passes_quality(metrics)
            if self.config.use_quality_early_stop and not quality_ok:
                if track_quality_rejection:
                    quality_rejected.append(solution)
                return False
            accepted += 1
            if (
                best_solution is None
                or solution.objective < best_solution.objective
            ):
                best_solution = solution
                best_molecule = molecule
                best_ranges = ranges
                best_metrics = metrics
                best_quality_ok = quality_ok
            return True

        def execute_phase(
            starts: list[_StartRecord],
            use_detection: bool,
            message_prefix: str = "",
        ) -> None:
            nonlocal solve_time, stopped_early
            started = time.perf_counter()
            try:
                for record in starts:
                    accepted_now = consume(
                        self._solve_one(record, use_detection),
                        message_prefix,
                    )
                    should_stop = (
                        self.config.use_early_stop
                        or self.config.use_quality_early_stop
                    )
                    if should_stop and accepted_now:
                        stopped_early = True
                        break
            finally:
                solve_time += time.perf_counter() - started

        fallback_triggered = False
        fallback_attempted = 0
        fallback_successful = 0
        fallback_accepted = 0
        quality_repair_triggered = False
        quality_repair_attempted = 0
        quality_repair_accepted = 0

        def execute_quality_repair(
            candidates: list[_CandidateSolution],
        ) -> bool:
            """Repair the closest quality reject; return whether attempted."""
            nonlocal quality_repair_triggered, quality_repair_attempted
            nonlocal quality_repair_accepted, solve_time, stopped_early
            if not candidates:
                return False
            candidate = min(candidates, key=self._quality_violation)
            quality_repair_triggered = True
            quality_repair_attempted += 1
            accepted_before = accepted
            repair_record = _StartRecord(
                x0=candidate.x,
                score=-candidate.objective,
                minimum_inequality=candidate.inequality_slack,
                feature=np.empty(0),
            )
            repair_started = time.perf_counter()
            try:
                repaired = self._solve_one(
                    repair_record,
                    use_detection=False,
                    enforce_quality=True,
                )
                accepted_now = consume(
                    repaired,
                    "quality repair: ",
                    track_quality_rejection=False,
                )
            finally:
                solve_time += time.perf_counter() - repair_started
            quality_repair_accepted += accepted - accepted_before
            if accepted_now:
                stopped_early = True
            return True

        levels: list[Optional[_SamplingLevel]] = (
            self._sampling_levels()
            if self.config.use_region_exclusion
            else [None]
        )
        for level in levels:
            used_levels += 1
            rejected_before_stage = len(quality_rejected)
            starts, stage_stats = self._build_starts(level)
            build_time += stage_stats.seed_time_s
            total_examined += stage_stats.joint_combinations_examined
            build_stats = stage_stats

            execute_phase(starts, self.config.use_convergence_detection)

            stage_rejected = quality_rejected[rejected_before_stage:]
            repair_attempted_for_stage = False
            if best_solution is None and stage_rejected:
                repair_attempted_for_stage = execute_quality_repair(
                    stage_rejected
                )

            stage_fallback = bool(
                self.config.use_convergence_detection
                and self.config.use_slsqp_fallback
                and starts
                and best_solution is None
            )
            if stage_fallback:
                fallback_triggered = True
                attempted_before = attempted
                successful_before = successful
                accepted_before = accepted
                execute_phase(starts, False, "fallback: ")
                fallback_attempted += attempted - attempted_before
                fallback_successful += successful - successful_before
                fallback_accepted += accepted - accepted_before

            if (
                best_solution is None
                and not repair_attempted_for_stage
                and len(quality_rejected) > rejected_before_stage
            ):
                execute_quality_repair(
                    quality_rejected[rejected_before_stage:]
                )

            if best_solution is not None:
                break

        build_stats.seed_time_s = build_time
        build_stats.adaptive_levels = used_levels
        build_stats.adaptive_expanded = bool(
            self.config.use_region_exclusion
            and self.config.use_adaptive_sampling
            and used_levels > 1
        )
        build_stats.joint_combinations_examined = total_examined
        stats = OptimizationStats(
            build=build_stats,
            solve_time_s=solve_time,
            attempted_starts=attempted,
            successful_starts=successful,
            accepted_successes=accepted,
            detected_early_starts=detected,
            stopped_early=stopped_early,
            fallback_triggered=fallback_triggered,
            fallback_attempted_starts=fallback_attempted,
            fallback_successful_starts=fallback_successful,
            fallback_accepted_successes=fallback_accepted,
            quality_repair_triggered=quality_repair_triggered,
            quality_repair_attempted_starts=quality_repair_attempted,
            quality_repair_accepted_successes=quality_repair_accepted,
            total_iterations=total_iterations,
            total_function_evaluations=total_function_evaluations,
            best_equality_error=(
                best_solution.equality_error
                if best_solution is not None
                else None
            ),
            best_inequality_slack=(
                best_solution.inequality_slack
                if best_solution is not None
                else None
            ),
            messages=messages,
        )
        if best_solution is None or best_molecule is None:
            return JointLagrangeResult(
                success=False,
                final=None,
                n_skeleton_atoms=self.skeleton.n,
                substituent_atom_ranges={},
                objective_function_value=float("inf"),
                quality_ok=False,
                quality_metrics={},
                stats=stats,
                message=f"failed: {accepted}/{attempted}",
            )
        return JointLagrangeResult(
            success=True,
            final=best_molecule,
            n_skeleton_atoms=self.skeleton.n,
            substituent_atom_ranges=best_ranges,
            objective_function_value=best_solution.objective,
            quality_ok=best_quality_ok,
            quality_metrics=best_metrics,
            stats=stats,
            message=f"success: {accepted}/{attempted}",
        )


def optimize_joint_lagrange(
    skeleton: Mol,
    substituents: list[Mol],
    skeleton_link_indices: list[int],
    substituent_link_indices: list[int],
    config: Optional[JointLagrangeConfig] = None,
    covalent_radii: Optional[dict[str, float]] = None,
    vdw_radii: Optional[dict[str, float]] = None,
) -> JointLagrangeResult:
    """Construct a ``JointLagrangeOptimizer`` and run it immediately.

    The parameters have exactly the same meaning as the
    ``JointLagrangeOptimizer`` constructor. This function is convenient for a
    one-off call; to inspect the precomputed data or repeat experiments,
    create an optimizer instance explicitly.
    """

    optimizer = JointLagrangeOptimizer(
        skeleton=skeleton,
        substituents=substituents,
        skeleton_link_indices=skeleton_link_indices,
        substituent_link_indices=substituent_link_indices,
        config=config,
        covalent_radii=covalent_radii,
        vdw_radii=vdw_radii,
    )
    return optimizer.optimize()


class IterateJointLagrangeAnalyzer:
    """CHEMSMART adapter for the Version 2 Joint-Lagrange optimizer.

    Joins a skeleton and one *or more* substituents into a single molecule via
    one 6K-dimensional joint optimization. The same path is used for a single
    attachment (K=1) and for multi-substituent combinations (K>1). The
    analyzer satisfies the iterate ``AnalyzerProtocol``: :meth:`run` returns
    the combined :class:`~chemsmart.io.molecules.structure.Molecule`, or
    ``None`` when the optimizer reports no acceptable structure (never a
    partial insertion).

    The heavy algorithm lives in :class:`JointLagrangeOptimizer` in this
    module, which the registry imports lazily so the ETKDG path never loads the
    Joint-Lagrange SciPy core. This class only performs the CHEMSMART <->
    Version 2 translation: it converts each :class:`Molecule` into the
    lightweight :class:`Mol` container, converts the 1-based iterate link
    indices to the 0-based indices the optimizer expects (exactly once), runs
    the optimizer a single time on the full substituent list, and converts the
    result back into a :class:`Molecule`.
    """

    #: Lagrange options mapped onto ``JointLagrangeConfig`` fields. Every other
    #: config field keeps its Version 2 (Test7.2/Test8) default.
    _CONFIG_OPTION_KEYS = (
        "use_adaptive_sampling",
        "n_link_sphere",
        "n_orientation_sphere",
        "n_axial",
        "candidate_pool_size",
        "preselect",
        "beam_width",
        "max_starts",
        "slsqp_maxiter",
    )

    def __init__(
        self,
        skeleton: Molecule,
        substituents: list[tuple[Molecule, int, int]],
        options: Optional[dict[str, object]] = None,
    ):
        """
        Parameters
        ----------
        skeleton : Molecule
            Preprocessed skeleton molecule shared by every attachment.
        substituents : list of (Molecule, int, int)
            One ``(substituent, skeleton_link_index, substituent_link_index)``
            tuple per attachment, where ``substituent`` is a preprocessed
            :class:`Molecule` and both link indices are 1-based. A single
            attachment is simply a list holding one tuple; the whole list is
            optimized jointly.
        options : dict of str to bool or int, optional
            Resolved Lagrange options (see the ``lagrange_multipliers`` entry
            in the algorithm registry). Only the nine keys in
            :attr:`_CONFIG_OPTION_KEYS` are consumed; anything absent keeps the
            Version 2 default.
        """
        self.skeleton = skeleton
        # Keep the received (Molecule, 1-based skeleton link, 1-based
        # substituent link) tuples verbatim; the 1-based -> 0-based conversion
        # happens exactly once inside run().
        self.substituents = list(substituents)
        self.options = dict(options or {})

    def _build_config(self) -> JointLagrangeConfig:
        """Map the exposed Lagrange options onto a ``JointLagrangeConfig``.

        Only options actually present override the Version 2 defaults, so the
        Test7.2/Test8 calibration is preserved for every field not exposed to
        the user.
        """
        config_kwargs = {
            key: self.options[key]
            for key in self._CONFIG_OPTION_KEYS
            if key in self.options
        }
        return JointLagrangeConfig(**config_kwargs)

    def _merge_frozen(self) -> Optional[list]:
        """Merge skeleton and substituent frozen-atom flags (skeleton first).

        The skeleton flags come first, each substituent's flags follow in
        input order, missing flags are padded with zeros of the matching
        length, and ``None`` is returned only when no molecule carries any
        frozen flags (the same convention as the ETKDG analyzer).
        """
        skeleton_frozen = self.skeleton.frozen_atoms
        substituent_frozens = [
            substituent.frozen_atoms for substituent, _, _ in self.substituents
        ]
        if skeleton_frozen is None and all(
            frozen is None for frozen in substituent_frozens
        ):
            return None
        merged = list(
            skeleton_frozen
            if skeleton_frozen is not None
            else [0] * len(self.skeleton)
        )
        for (substituent, _, _), frozen in zip(
            self.substituents, substituent_frozens
        ):
            merged.extend(
                frozen if frozen is not None else [0] * len(substituent)
            )
        return merged

    @staticmethod
    def _to_joint_mol(molecule: Molecule) -> Mol:
        """Convert a CHEMSMART Molecule into the Version 2 ``Mol``."""
        return Mol(
            positions=molecule.positions,
            elements=list(molecule.chemical_symbols),
        )

    def run(self) -> Optional[Molecule]:
        """Run one joint optimization and return the combined molecule.

        Returns
        -------
        Molecule or None
            Combined skeleton + substituents molecule (skeleton atoms first,
            then each substituent in input order), inheriting the skeleton's
            charge/multiplicity and the merged frozen-atom flags; ``None`` if
            the optimizer reports no acceptable structure.
        """
        skeleton_mol = self._to_joint_mol(self.skeleton)
        substituent_mols = []
        skeleton_link_indices = []
        substituent_link_indices = []
        for (
            substituent,
            skeleton_link_index,
            substituent_link_index,
        ) in self.substituents:
            substituent_mols.append(self._to_joint_mol(substituent))
            # Convert 1-based iterate indices to 0-based joint indices once.
            skeleton_link_indices.append(skeleton_link_index - 1)
            substituent_link_indices.append(substituent_link_index - 1)

        result = optimize_joint_lagrange(
            skeleton=skeleton_mol,
            substituents=substituent_mols,
            skeleton_link_indices=skeleton_link_indices,
            substituent_link_indices=substituent_link_indices,
            config=self._build_config(),
        )

        logger.debug(
            "Joint Lagrange %s (objective=%.6g, quality_ok=%s)",
            result.message,
            result.objective_function_value,
            result.quality_ok,
        )
        logger.debug(
            "Joint Lagrange quality metrics: %s", result.quality_metrics
        )

        if not result.success or result.final is None:
            return None

        # result.final already orders atoms as skeleton + each substituent in
        # input order; merge frozen flags in the same order.
        return Molecule(
            symbols=list(result.final.elements),
            positions=result.final.positions,
            charge=self.skeleton.charge,
            multiplicity=self.skeleton.multiplicity,
            frozen_atoms=self._merge_frozen(),
        )
