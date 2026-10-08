import os
import re
from io import StringIO
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from chemsmart.analysis.thermochemistry import (
    MECPProjectedFrequencyOutput,
    MECPThermochemistry,
    Thermochemistry,
    thermochemistry_from_file,
)
from chemsmart.cli.gaussian.mecp_options import add_mecp_method_suffix
from chemsmart.io.gaussian.output import Gaussian16Output
from chemsmart.jobs.gaussian.mecp import GaussianMECPJob
from chemsmart.jobs.gaussian.runner import GaussianJobRunner
from chemsmart.jobs.gaussian.settings import (
    GaussianJobSettings,
    GaussianLinkJobSettings,
    GaussianMECPJobSettings,
)
from chemsmart.jobs.gaussian.writer import GaussianInputWriter
from chemsmart.jobs.thermochemistry.job import ThermochemistryJob
from chemsmart.utils.mixins import FileMixin
from chemsmart.utils.utils import remove_word_from_parenthesized_option


def test_remove_arbitrary_word_from_parenthesized_option():
    match = re.fullmatch(r"(\w+)=\(([^)]*)\)", "opt=(calcfc,modredundant)")

    assert (
        remove_word_from_parenthesized_option(match, "calcfc")
        == "opt=(modredundant)"
    )
    assert remove_word_from_parenthesized_option(match, "MODREDUNDANT") == (
        "opt=(calcfc)"
    )


def test_first_mecp_step_removes_unavailable_guess_read():
    remove_read = GaussianMECPJob._route_without_guess_read

    assert remove_read("scf=xqc guess=read nosymm") == "scf=xqc nosymm"
    assert remove_read("guess=(mix,read) nosymm") == "guess=(mix) nosymm"
    assert remove_read("guess=(read) nosymm") == "nosymm"
    assert remove_read("scf=xqc nosymm") == "scf=xqc nosymm"


@pytest.mark.parametrize("keyword", ["nosymm", "NoSymm", "nosymmetry"])
def test_mecp_accepts_disabled_symmetry_keywords(keyword):
    ensure_nosymm = GaussianMECPJob._with_required_nosymm

    route = f"{keyword} scf=xqc"
    assert ensure_nosymm(route) == route


@pytest.mark.parametrize("keyword", ["symm", "symmetry", "symmetry=loose"])
def test_mecp_rejects_enabled_symmetry_keywords(keyword):
    ensure_nosymm = GaussianMECPJob._with_required_nosymm

    with pytest.raises(ValueError, match="requires nosymm"):
        ensure_nosymm(keyword)


def test_mecp_adds_nosymm_when_symmetry_is_unspecified():
    assert GaussianMECPJob._with_required_nosymm("scf=xqc") == (
        "scf=xqc nosymm"
    )


def test_gaussian_spin_squared_parser_returns_final_values():
    output = object.__new__(Gaussian16Output)
    output.__dict__["contents"] = [
        "S**2 before annihilation     2.1040, after     2.0060",
        "S**2 before annihilation     2.0550, after     2.0015",
    ]

    assert output.spin_squared_before_annihilation == pytest.approx(2.0550)
    assert output.spin_squared_after_annihilation == pytest.approx(2.0015)

    output.__dict__["contents"] = ["no spin contamination data"]
    output.__dict__.pop("_final_spin_squared_values", None)
    assert output.spin_squared_before_annihilation is None
    assert output.spin_squared_after_annihilation is None


def test_gaussian_header_reuses_rolling_checkpoint():
    job = SimpleNamespace(
        label="mecp_step2_A",
        chkfile="steps/mecp_A.chk",
        oldchkfile="steps/mecp_A.chk",
        settings=SimpleNamespace(chk=True),
        jobrunner=SimpleNamespace(num_cores=4, mem_gb=8),
    )
    output = StringIO()

    GaussianInputWriter(job)._write_gaussian_header(output)

    assert output.getvalue().splitlines()[0] == "%chk=mecp_A.chk"
    assert "%oldchk" not in output.getvalue().lower()


def test_seam_checkpoint_cleanup_preserves_final_mecp_checkpoints(
    monkeypatch,
):
    job = object.__new__(GaussianMECPJob)
    job._state_checkpoint_files = {
        (None, "A"): "mecp_A.chk",
        (None, "B"): "mecp_B.chk",
        ("check", "A"): "mecp_check_A.chk",
        ("check", "B"): "mecp_check_B.chk",
    }
    removed = []
    monkeypatch.setattr("os.path.isfile", lambda path: True)
    monkeypatch.setattr("os.remove", removed.append)

    job._remove_checkpoint_set("check")

    assert removed == ["mecp_check_A.chk", "mecp_check_B.chk"]
    assert job._state_checkpoint_files == {
        (None, "A"): "mecp_A.chk",
        (None, "B"): "mecp_B.chk",
    }


def test_mecp_scratch_jobs_are_grouped_in_steps_folder():
    runner = object.__new__(GaussianJobRunner)
    runner._scratch_dir = "/scratch/project"
    job = SimpleNamespace(
        label="mecp_step2_A", scratch_parent_folder="mecp_steps"
    )

    scratch_directory = runner._scratch_job_directory(job)

    assert scratch_directory == os.path.join(
        "/scratch/project", "mecp_steps", "mecp_step2_A"
    )


def test_inverse_bfgs_update_satisfies_secant_condition():
    inv_hessian = np.diag([0.7, 1.1, 1.6])
    delta_x = np.array([0.12, -0.04, 0.08])
    delta_g = np.array([0.20, -0.03, 0.10])

    updated = GaussianMECPJob._update_inverse_hessian(
        inv_hessian, delta_x, delta_g
    )

    np.testing.assert_allclose(updated @ delta_g, delta_x, atol=1.0e-12)
    np.testing.assert_allclose(updated, updated.T, atol=1.0e-14)
    assert np.all(np.linalg.eigvalsh(updated) > 0.0)


def test_inverse_bfgs_update_accepts_negative_curvature_like_easymecp():
    inv_hessian = np.eye(3)
    delta_x = np.array([1.0, 0.0, 0.0])
    delta_g = np.array([-1.0, 0.0, 0.0])

    updated = GaussianMECPJob._update_inverse_hessian(
        inv_hessian, delta_x, delta_g
    )

    np.testing.assert_allclose(updated @ delta_g, delta_x)
    assert np.linalg.eigvalsh(updated).min() < 0.0


def test_inverse_bfgs_update_skips_singular_secant_pair():
    inv_hessian = np.eye(3)
    delta_x = np.zeros(3)
    delta_g = np.ones(3)

    updated = GaussianMECPJob._update_inverse_hessian(
        inv_hessian, delta_x, delta_g
    )

    np.testing.assert_array_equal(updated, inv_hessian)


def test_removed_harvey_bfgs_setting_is_rejected():
    with pytest.raises(ValueError, match="Unknown MECP step_size_method"):
        GaussianMECPJobSettings(step_size_method="harvey_bfgs")


@pytest.mark.parametrize(
    "kwargs",
    [
        {"max_steps": 0},
        {"hess_step_size": 0.0},
        {"step_size_min": 2.0, "step_size_max": 1.0},
        {"multiplicity_a": 0},
        {"multiplicity_a": 1, "multiplicity_b": 1},
    ],
)
def test_invalid_mecp_settings_are_rejected(kwargs):
    with pytest.raises(ValueError):
        GaussianMECPJobSettings(**kwargs)


def test_reduced_hessian_preserves_negative_seam_curvature():
    hessian = np.diag([-2.0, 0.0, 3.0])
    projected = [np.array([0.0, 1.0, 0.0])]

    reduced = GaussianMECPJob._reduced_hessian(hessian, projected)
    eigenvalues = np.linalg.eigvalsh(reduced)

    np.testing.assert_allclose(eigenvalues, [-2.0, 3.0], atol=1.0e-12)


def test_projected_mecp_frequencies_have_3n_minus_7_modes():
    job = object.__new__(GaussianMECPJob)
    job.molecule = SimpleNamespace(
        symbols=["H", "H", "H"],
        most_abundant_masses=np.ones(3),
    )
    positions = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    diff_grad = np.array([[1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]])

    frequencies, modes, eigenvalues, n_projected = (
        job._projected_frequencies_and_modes(np.eye(9), positions, diff_grad)
    )

    assert n_projected == 7
    assert frequencies.shape == (2,)
    assert modes.shape == (2, 9)
    assert eigenvalues.shape == (2,)
    assert np.all(frequencies > 0.0)
    np.testing.assert_allclose(np.linalg.norm(modes, axis=1), 1.0)


def test_mass_weighted_frequency_preserves_imaginary_mode_sign():
    positive = GaussianMECPJob._frequency_from_mass_weighted_eigenvalue(1.0)
    negative = GaussianMECPJob._frequency_from_mass_weighted_eigenvalue(-1.0)

    assert positive > 0.0
    assert negative == pytest.approx(-positive)


def test_mecp_frequency_log_is_accepted_by_thermochemistry(tmp_path):
    job = object.__new__(GaussianMECPJob)
    job.folder = str(tmp_path)
    job.label = "crossing"
    job.molecule = SimpleNamespace(symbols=["H", "H"])
    result = {
        "positions_angstrom": np.array([[0.0, 0.0, 0.0], [0.7, 0.0, 0.0]]),
        "atomic_masses": np.array([1.007825, 1.007825]),
        "energy_a": -1.0,
        "energy_b": -0.99998,
        "mecp_energy": -0.99999,
        "energy_diff": -2.0e-5,
        "lagrange_multiplier": 0.5,
        "n_projected": 6,
        "n_negative": 0,
        "frequencies": np.array([1234.5]),
        "modes": np.array([[0.5, 0.0, 0.0, -0.5, 0.0, 0.0]]),
        "multiplicity_a": 1,
        "multiplicity_b": 3,
        "rotational_symmetry_number": 1,
    }

    job._write_mecp_frequency_log(result, 1.0e-3)
    frequency_file = tmp_path / "crossing_mecp_freq.log"
    frequency_output = MECPProjectedFrequencyOutput(str(frequency_file))
    thermochemistry = thermochemistry_from_file(
        str(frequency_file),
        temperature=298.15,
        electronic_degeneracy=4,
    )
    thermochemistry_job = ThermochemistryJob.from_filename(
        str(frequency_file),
        jobrunner=object(),
    )

    assert isinstance(frequency_output, FileMixin)
    assert frequency_output.normal_termination
    assert frequency_output.freq
    assert frequency_output.symbols == ["H", "H"]
    assert frequency_output.vibrational_frequencies == pytest.approx([1234.5])
    assert isinstance(thermochemistry, MECPThermochemistry)
    assert isinstance(thermochemistry, Thermochemistry)
    assert "_load_molecule" not in MECPThermochemistry.__dict__
    assert thermochemistry.jobtype == "mecp"
    assert thermochemistry.vibrational_frequencies == pytest.approx([1234.5])
    assert thermochemistry.file_object.energies == pytest.approx([-0.99999])
    assert thermochemistry.multiplicity == 4
    assert thermochemistry.rotational_symmetry_number == 1
    assert thermochemistry_job.molecule.symbols == ["H", "H"]
    assert thermochemistry_job.label == "crossing_mecp_freq"

    result["frequencies"] = np.array([-1.0])
    job._write_mecp_frequency_log(result, 1.0e-3)
    contents = frequency_file.read_text(encoding="utf-8")
    assert "n_imaginary=1\n" in contents
    assert "n_significant_imaginary=0\n" in contents
    assert contents.rstrip().endswith(
        MECPProjectedFrequencyOutput.TERMINATION_MARKER
    )


def test_mecp_frequency_status_is_determined_from_file_content(tmp_path):
    frequency_file = tmp_path / "incomplete_mecp_freq.log"
    frequency_file.write_text(
        f"{MECPProjectedFrequencyOutput.HEADER}\n",
        encoding="utf-8",
    )

    output = MECPProjectedFrequencyOutput(str(frequency_file))

    assert not output.normal_termination
    assert not output.freq
    with pytest.raises(ValueError, match="did not terminate normally"):
        thermochemistry_from_file(str(frequency_file), temperature=298.15)

    frequency_file.write_text(
        f"{MECPProjectedFrequencyOutput.HEADER}\n"
        f"{MECPProjectedFrequencyOutput.FREQUENCIES_HEADER}\n"
        f"{MECPProjectedFrequencyOutput.TERMINATION_MARKER}\n",
        encoding="utf-8",
    )
    output = MECPProjectedFrequencyOutput(str(frequency_file))

    assert output.normal_termination
    assert output.freq


def test_lagrangian_hessian_weight_is_not_forced_to_average():
    grad_a = np.array([2.0, 0.0])
    grad_b = np.array([-1.0, 0.0])
    hessian_a = np.diag([1.0, 2.0])
    hessian_b = np.diag([4.0, 8.0])
    hessian, multiplier = GaussianMECPJob._lagrangian_hessian(
        hessian_a, hessian_b, grad_a, grad_b
    )

    assert multiplier == pytest.approx(2.0 / 3.0)
    np.testing.assert_allclose(
        hessian, (1.0 / 3.0) * hessian_a + (2.0 / 3.0) * hessian_b
    )


def test_optimizer_state_round_trip(tmp_path):
    job = object.__new__(GaussianMECPJob)
    job.folder = str(tmp_path)
    job.label = "restartable"
    job.molecule = SimpleNamespace(
        symbols=["H", "H"], positions=np.zeros((2, 3))
    )
    job.settings = SimpleNamespace(step_size_method="harvey")
    positions = np.arange(6, dtype=float).reshape(2, 3)
    inverse_hessian = np.eye(6)

    job._save_optimizer_state(
        next_step=12,
        positions_bohr=positions,
        current_step_size=0.1,
        prev_merit=None,
        prev_positions=None,
        prev_proj_grad=None,
        inv_hessian=inverse_hessian,
        prev_eff_grad=np.ones(6),
        prev_positions_bfgs=positions - 0.1,
    )
    restored = job._load_optimizer_state()

    assert restored["next_step"] == 12
    assert restored["prev_merit"] is None
    assert restored["prev_positions"] is None
    np.testing.assert_array_equal(restored["positions_bohr"], positions)
    np.testing.assert_array_equal(restored["inv_hessian"], inverse_hessian)


def test_optimizer_state_rejects_method_change(tmp_path):
    job = object.__new__(GaussianMECPJob)
    job.folder = str(tmp_path)
    job.label = "restartable"
    job.molecule = SimpleNamespace(symbols=["H"], positions=np.zeros((1, 3)))
    job.settings = SimpleNamespace(step_size_method="harvey")
    job._save_optimizer_state(
        next_step=2,
        positions_bohr=np.zeros((1, 3)),
        current_step_size=0.1,
        prev_merit=None,
        prev_positions=None,
        prev_proj_grad=None,
        inv_hessian=np.eye(3),
        prev_eff_grad=np.zeros(3),
        prev_positions_bfgs=np.zeros((1, 3)),
    )
    job.settings.step_size_method = "bb"

    with pytest.raises(RuntimeError, match="optimizer method differs"):
        job._load_optimizer_state()


def test_job_is_complete_requires_post_verification_marker(tmp_path):
    job = object.__new__(GaussianMECPJob)
    job.settings = GaussianMECPJobSettings()
    job.folder = str(tmp_path)
    job.label = "verified"
    report = Path(job.report_file)
    report.parent.mkdir()
    report.write_text("Optimization converged at step 8.\n", encoding="utf-8")

    assert job._job_is_complete() is False

    report.write_text("Converged at step 8.\n", encoding="utf-8")
    assert job._job_is_complete() is True


def test_mecp_output_paths_follow_functional_layout(tmp_path):
    job = object.__new__(GaussianMECPJob)
    job.folder = str(tmp_path)
    job.label = "crossing"

    assert Path(job.final_report_file) == tmp_path / "crossing_final_report.log"
    assert Path(job.report_file) == (
        tmp_path / "crossing_optimization" / "crossing_report.log"
    )
    assert Path(job.trajectory_file) == (
        tmp_path / "crossing_optimization" / "crossing_traj.xyz"
    )
    assert Path(job.numfreq_folder) == tmp_path / "crossing_numfreq"


def test_state_subjobs_are_routed_by_calculation_phase(tmp_path, monkeypatch):
    import chemsmart.jobs.gaussian.mecp as mecp_module

    calls = []

    class Molecule:
        symbols = ["H"]

        def __init__(self):
            self.positions = np.zeros((1, 3))

        def copy(self):
            return Molecule()

    class FakeSubjob:
        def __init__(self, **kwargs):
            self.label = kwargs["label"]
            self.scratch_parent_folder = kwargs["scratch_parent_folder"]
            self.chkfile = str(tmp_path / "missing.chk")

        def set_folder(self, folder):
            self.folder = folder

        def run(self):
            calls.append(
                (self.label, Path(self.folder), self.scratch_parent_folder)
            )

        def _output(self):
            return SimpleNamespace(
                energies=[-1.0],
                forces=[np.zeros((1, 3))],
                spin_squared_after_annihilation=None,
            )

    monkeypatch.setattr(mecp_module, "GaussianGeneralJob", FakeSubjob)
    job = object.__new__(GaussianMECPJob)
    job.folder = str(tmp_path)
    job.label = "crossing"
    job.jobrunner = object()
    job.molecule = Molecule()
    job.settings = SimpleNamespace(
        charge_a=0,
        charge_b=0,
        multiplicity_a=1,
        multiplicity_b=3,
        title_a="A",
        title_b="B",
        use_link=False,
    )
    job.steps_folder = job.optimization_folder
    job._state_checkpoint_files = {}
    job._state_settings = lambda **kwargs: SimpleNamespace(
        additional_route_parameters=None
    )

    job._run_state(np.zeros((1, 3)), 1, "A")
    job._run_state(
        np.zeros((1, 3)), "macro07_check_coord01_plus", "B", "check"
    )

    assert calls == [
        (
            "crossing_step1_A",
            Path(job.optimization_folder),
            "crossing_optimization",
        ),
        (
            "crossing_macro07_check_coord01_plus_B",
            Path(job.numfreq_folder),
            "crossing_numfreq",
        ),
    ]
    assert Path(job.optimization_folder).is_dir()
    assert Path(job.numfreq_folder).is_dir()


def test_completed_mecp_requires_requested_frequency_outputs(tmp_path):
    job = object.__new__(GaussianMECPJob)
    job.settings = GaussianMECPJobSettings(mecp_numfreq=True)
    job.folder = str(tmp_path)
    job.label = "crossing"
    Path(job.optimization_folder).mkdir()
    Path(job.report_file).write_text("Converged at step 3.\n")
    assert not job._job_is_complete()
    Path(job.numfreq_folder).mkdir()
    (Path(job.numfreq_folder) / "crossing_seam_check.log").write_text(
        "checked\n"
    )
    assert not job._job_is_complete()
    (tmp_path / "crossing_mecp_freq.log").write_text("frequencies\n")
    assert job._job_is_complete()


def test_seam_check_uses_same_tolerance_with_and_without_frequencies():
    job = object.__new__(GaussianMECPJob)
    job.label = "crossing"
    job.settings = GaussianMECPJobSettings()
    job.molecule = SimpleNamespace(
        symbols=["H"] * 3,
        positions=np.zeros((3, 3)),
        most_abundant_masses=np.ones(3),
    )
    job._run_state = lambda *args, **kwargs: (0.0, np.ones((3, 3)))
    job._compute_numerical_hessian = lambda *args, **kwargs: (None, None)
    job._lagrangian_hessian = lambda *args: (np.eye(9), 0.5)
    job._build_projection_vectors = lambda *args: []
    job._projected_frequencies_and_modes = lambda *args: (
        np.array([-10.0]),
        np.zeros((1, 9)),
        np.array([-4e-6]),
        7,
    )
    job._write_seam_check_log = lambda *args: None
    job._write_mecp_frequency_log = lambda *args: None
    job._remove_checkpoint_set = lambda *args: None
    for write_frequencies in (False, True):
        result = job.verify_seam_minimum(write_frequencies=write_frequencies)
        assert result["n_negative"] == 1
        assert not result["is_minimum"]


@pytest.mark.parametrize("swap_states", [False, True])
def test_harvey_optimizer_converges_on_analytic_crossing(swap_states):
    job = object.__new__(GaussianMECPJob)
    job.settings = GaussianMECPJobSettings()
    positions = np.array([[0.7, 1.0, 0.0]])
    prev_positions = None
    prev_gradient = None
    inverse_hessian = None

    for _ in range(100):
        x, y, _ = positions[0]
        energy_a = 0.5 * ((x - 1.0) ** 2 + y**2)
        energy_b = 0.5 * ((x + 1.0) ** 2 + y**2)
        gradient_a = np.array([[x - 1.0, y, 0.0]])
        gradient_b = np.array([[x + 1.0, y, 0.0]])
        if swap_states:
            energy_a, energy_b = energy_b, energy_a
            gradient_a, gradient_b = gradient_b, gradient_a

        (
            displacement,
            effective_gradient,
            _,
            inverse_hessian,
            _,
            _,
            _,
            optimizer_gradient,
        ) = job._bfgs_displacement(
            energy_a,
            energy_b,
            gradient_a,
            gradient_b,
            prev_positions,
            positions,
            prev_gradient,
            inverse_hessian,
        )
        if job._is_converged(
            energy_a - energy_b, effective_gradient, displacement
        ):
            break
        prev_positions = positions.copy()
        prev_gradient = optimizer_gradient.ravel().copy()
        positions = positions + displacement
    else:
        pytest.fail("Harvey optimizer did not converge on analytic surfaces")

    np.testing.assert_allclose(positions, np.zeros((1, 3)), atol=1.0e-5)


def _adaptive_step_job():
    job = object.__new__(GaussianMECPJob)
    job.settings = SimpleNamespace(
        step_size=0.1,
        step_size_grow=1.2,
        step_size_shrink=0.5,
        step_size_min=1.0e-4,
        step_size_max=1.0,
    )
    return job


@pytest.mark.parametrize(
    ("previous_merit", "current_merit", "expected"),
    [
        (100.0, 80.0, 0.12),
        (100.0, 95.0, 0.10),
        (100.0, 101.0, 0.10),
        (100.0, 103.0, 0.05),
    ],
)
def test_grow_shrink_uses_relative_dead_band(
    previous_merit, current_merit, expected
):
    job = _adaptive_step_job()

    step = job._adapt_step_size(0.1, previous_merit, current_merit)

    assert step == pytest.approx(expected)


def test_bb_step_limits_abrupt_growth():
    job = _adaptive_step_job()

    step = job._bb_step_size(
        0.1,
        np.array([[0.0, 0.0, 0.0]]),
        np.array([[1.0, 0.0, 0.0]]),
        np.array([[0.0, 0.0, 0.0]]),
        np.array([[0.1, 0.0, 0.0]]),
    )

    assert step == pytest.approx(0.2)


def test_bb_step_damps_unreliable_curvature():
    job = _adaptive_step_job()

    step = job._bb_step_size(
        0.1,
        np.array([[0.0, 0.0, 0.0]]),
        np.array([[1.0, 0.0, 0.0]]),
        np.array([[0.0, 0.0, 0.0]]),
        np.array([[-1.0, 0.0, 0.0]]),
    )

    assert step == pytest.approx(0.05)


def test_bb_step_uses_only_seam_tangent_displacement():
    job = _adaptive_step_job()

    step = job._bb_step_size(
        0.4,
        np.array([[0.0, 0.0, 0.0]]),
        np.array([[10.0, 1.0, 0.0]]),
        np.array([[0.0, 0.0, 0.0]]),
        np.array([[0.0, 2.0, 0.0]]),
        np.array([[1.0, 0.0, 0.0]]),
    )

    assert step == pytest.approx(0.5)


@pytest.mark.parametrize(
    ("label", "method", "expected"),
    [
        ("sample", "bb", "sample_bb_mecp"),
        ("sample_mecp", "grow_shrink", "sample_grow_shrink_mecp"),
        ("sample_harvey_mecp", "harvey", "sample_harvey_mecp"),
    ],
)
def test_add_mecp_method_suffix(label, method, expected):
    assert add_mecp_method_suffix(label, method) == expected


def test_mecp_settings_presets_and_explicit_overrides():
    tight = GaussianMECPJobSettings(convergence_preset="tight")
    overridden = GaussianMECPJobSettings(
        convergence_preset="tight",
        energy_diff_tol=9.0e-5,
        trust_radius=0.25,
    )

    assert tight.energy_diff_tol == pytest.approx(1.0e-5)
    assert tight.force_rms_tol == pytest.approx(1.0e-4)
    assert tight.trust_radius == pytest.approx(0.1)
    assert overridden.energy_diff_tol == pytest.approx(9.0e-5)
    assert overridden.trust_radius == pytest.approx(0.25)

    with pytest.raises(ValueError, match="Unknown convergence_preset"):
        GaussianMECPJobSettings(convergence_preset="extreme")


def test_mecp_settings_are_copied_from_generic_settings():
    generic = GaussianJobSettings(
        functional="pbe0",
        basis="def2svp",
        charge=1,
        multiplicity=2,
    )

    converted = GaussianMECPJobSettings.from_settings(generic)
    copied = GaussianMECPJobSettings.from_settings(converted)

    assert converted.functional == "pbe0"
    assert converted.basis == "def2svp"
    assert copied == converted
    assert copied is not converted


@pytest.mark.parametrize("use_link", [False, True])
def test_state_settings_build_force_jobs_with_nosymm(use_link):
    job = object.__new__(GaussianMECPJob)
    job.settings = GaussianMECPJobSettings(
        use_link=use_link,
        additional_route_parameters="scf=xqc",
        stable="qrhf",
        guess="mix",
    )

    state_settings = job._state_settings(1, 2, "state A")

    expected_type = (
        GaussianLinkJobSettings if use_link else GaussianJobSettings
    )
    assert isinstance(state_settings, expected_type)
    assert state_settings.jobtype == "sp"
    assert state_settings.forces is True
    assert state_settings.freq is False
    assert state_settings.numfreq is False
    assert state_settings.charge == 1
    assert state_settings.multiplicity == 2
    assert "nosymm" in state_settings.additional_route_parameters.lower()
    if use_link:
        assert state_settings.link is True
        assert state_settings.stable == "qrhf"
        assert state_settings.guess == "mix"


def test_optimizer_state_rejects_symbols_and_coordinate_shape(tmp_path):
    job = object.__new__(GaussianMECPJob)
    job.folder = str(tmp_path)
    job.label = "restartable"
    job.molecule = SimpleNamespace(
        symbols=["H", "H"], positions=np.zeros((2, 3))
    )
    job.settings = SimpleNamespace(step_size_method="bb")
    state = dict(
        next_step=2,
        positions_bohr=np.zeros((2, 3)),
        current_step_size=0.1,
        prev_merit=1.0,
        prev_positions=np.zeros((2, 3)),
        prev_proj_grad=np.ones((2, 3)),
        inv_hessian=None,
        prev_eff_grad=None,
        prev_positions_bfgs=None,
    )
    job._save_optimizer_state(**state)

    job.molecule.symbols = ["H", "He"]
    with pytest.raises(RuntimeError, match="atom symbols differ"):
        job._load_optimizer_state()

    job.molecule.symbols = ["H", "H"]
    job.molecule.positions = np.zeros((1, 3))
    with pytest.raises(RuntimeError, match="coordinate shape differs"):
        job._load_optimizer_state()


def test_trust_radius_scales_only_atoms_that_exceed_limit():
    job = object.__new__(GaussianMECPJob)
    job.settings = SimpleNamespace(trust_radius=0.5)

    limited = job._apply_trust_radius(
        np.array([[3.0, 4.0, 0.0], [0.1, 0.2, 0.0]])
    )

    np.testing.assert_allclose(limited[0], [0.3, 0.4, 0.0])
    np.testing.assert_allclose(limited[1], [0.1, 0.2, 0.0])


def test_convergence_requires_every_threshold():
    job = object.__new__(GaussianMECPJob)
    job.settings = SimpleNamespace(
        energy_diff_tol=0.1,
        force_max_tol=0.2,
        force_rms_tol=0.2,
        disp_max_tol=0.3,
        disp_rms_tol=0.3,
    )
    gradient = np.array([[0.1, -0.1, 0.0]])
    displacement = np.array([[0.2, 0.0, 0.0]])

    assert job._is_converged(0.05, gradient, displacement)
    assert not job._is_converged(0.11, gradient, displacement)
    assert not job._is_converged(0.05, gradient * 3, displacement)
    assert not job._is_converged(0.05, gradient, displacement * 2)


def test_constrained_seam_step_restores_both_linear_constraints():
    grad_a = np.array([[2.0, 1.0, 3.0]])
    grad_b = np.array([[1.0, 1.0, 3.0]])
    mode = np.array([0.0, 1.0, 0.0])

    displacement, tangent_gradient, correction = (
        GaussianMECPJob._constrained_seam_displacement(
            energy_diff=0.2,
            progress_error=0.3,
            grad_a=grad_a,
            grad_b=grad_b,
            progress_mode=mode,
            step_size=0.1,
        )
    )

    diff_grad = (grad_a - grad_b).ravel()
    np.testing.assert_allclose(np.dot(diff_grad, correction.ravel()), -0.2)
    np.testing.assert_allclose(np.dot(mode, correction.ravel()), -0.3)
    np.testing.assert_allclose(np.dot(diff_grad, tangent_gradient.ravel()), 0.0)
    np.testing.assert_allclose(np.dot(mode, tangent_gradient.ravel()), 0.0)
    np.testing.assert_allclose(displacement, [[-0.2, -0.3, -0.3]])


def test_negative_mode_tracking_uses_overlap_and_preserves_orientation():
    result = {
        "frequency_eigenvalues": np.array([-3.0, -2.0, 1.0]),
        "modes": np.array(
            [
                [1.0, 0.0, 0.0],
                [0.0, -1.0, 0.0],
                [0.0, 0.0, 1.0],
            ]
        ),
    }

    index, mode, overlap = GaussianMECPJob._overlap_tracked_negative_mode(
        result, np.array([0.1, 0.9, 0.0])
    )

    assert index == 1
    assert overlap > 0.99
    np.testing.assert_allclose(mode, [0.0, 1.0, 0.0])


def test_mode_displacement_reshapes_flat_mode_to_cartesian_geometry():
    positions = np.zeros((2, 3))
    mode = np.arange(6, dtype=float)

    displaced = GaussianMECPJob._displace_along_mode(
        positions, mode, distance=0.1
    )

    assert displaced.shape == (2, 3)
    np.testing.assert_allclose(displaced, 0.1 * mode.reshape(2, 3))


def test_numerical_hessian_uses_central_differences_and_symmetrizes():
    job = object.__new__(GaussianMECPJob)
    job.molecule = SimpleNamespace(symbols=["H"])
    calls = []
    hessian_a = np.diag([1.0, 2.0, 3.0])
    hessian_b = np.diag([4.0, 5.0, 6.0])

    def run_state(positions, step, state, checkpoint_tag=None):
        calls.append((step, state, checkpoint_tag))
        matrix = hessian_a if state == "A" else hessian_b
        gradient = (matrix @ np.asarray(positions).ravel()).reshape(1, 3)
        return 0.0, gradient

    job._run_state = run_state
    actual_a, actual_b = job._compute_numerical_hessian(
        np.zeros((1, 3)), h=1.0e-3, step_prefix=10
    )

    np.testing.assert_allclose(actual_a, hessian_a)
    np.testing.assert_allclose(actual_b, hessian_b)
    assert len(calls) == 12
    assert calls[0] == (10, "A", "check")
    assert calls[-1] == (15, "B", "check")


def test_follow_seam_hessian_uses_macro_coordinate_and_sign_labels():
    job = object.__new__(GaussianMECPJob)
    job.molecule = SimpleNamespace(symbols=["H"])
    calls = []

    def run_state(positions, step, state, checkpoint_tag=None):
        calls.append((step, state, checkpoint_tag))
        return 0.0, np.asarray(positions, dtype=float)

    job._run_state = run_state
    job._compute_numerical_hessian(
        np.zeros((1, 3)), h=1.0e-3, step_prefix=1, macro_step=7
    )

    assert calls[0] == ("macro07_check_coord01_plus", "A", "check")
    assert calls[3] == ("macro07_check_coord01_minus", "B", "check")
    assert calls[-1] == ("macro07_check_coord03_minus", "B", "check")


def test_follow_seam_inner_steps_use_macro_labels(tmp_path):
    job = _driver_job(tmp_path)
    job.settings.max_steps = 2
    seen = []

    def run_state(positions, step, state):
        seen.append((step, state))
        gradient = (
            np.array([[1.0, 0.0, 0.0]]) if state == "A" else np.zeros((1, 3))
        )
        return 0.0, gradient

    job._run_state = run_state
    job._optimize_on_seam_progress_plane(
        positions_bohr=np.zeros((1, 3)),
        plane_point_bohr=np.zeros((1, 3)),
        progress_mode=np.array([0.0, 1.0, 0.0]),
        trace=StringIO(),
        macro_step=7,
    )

    assert seen == [("macro07_inner001", "A"), ("macro07_inner001", "B")]


def test_lagrangian_hessian_rejects_identical_gradients():
    gradient = np.ones((1, 3))
    with pytest.raises(RuntimeError, match="Difference gradient is too small"):
        GaussianMECPJob._lagrangian_hessian(
            np.eye(3), np.eye(3), gradient, gradient
        )


def test_seam_check_runner_accepts_minimum_and_rejects_saddle():
    job = object.__new__(GaussianMECPJob)
    job.label = "crossing"
    job.settings = SimpleNamespace(mecp_numfreq=True)
    job.verify_seam_minimum = lambda **kwargs: {
        "is_minimum": True,
        "n_negative": 0,
    }
    job._run_seam_minimum_check(np.zeros((1, 3)))

    job.verify_seam_minimum = lambda **kwargs: {
        "is_minimum": False,
        "n_negative": 2,
    }
    with pytest.raises(RuntimeError, match="2 negative eigenvalue"):
        job._run_seam_minimum_check(np.zeros((1, 3)))


def test_seam_check_runner_follows_saddle_when_requested():
    job = object.__new__(GaussianMECPJob)
    job.label = "crossing"
    job.settings = SimpleNamespace(
        mecp_numfreq=True, follow_seam_imaginary_mode=True
    )
    saddle = {"is_minimum": False, "n_negative": 1}
    minimum = {"is_minimum": True, "n_negative": 0}
    job.verify_seam_minimum = lambda **kwargs: saddle
    seen = []
    job._follow_seam_imaginary_mode = lambda result: (
        seen.append(result),
        minimum,
    )[1]

    assert job._run_seam_minimum_check(np.zeros((1, 3))) is minimum
    assert seen == [saddle]
    assert job._initial_seam_result is saddle
    assert job._last_seam_result is minimum


def test_step_and_seam_logs_include_diagnostics(tmp_path):
    job = object.__new__(GaussianMECPJob)
    job.folder = str(tmp_path)
    job.label = "crossing"
    report = StringIO()
    job.log_step(
        report,
        3,
        -10.0,
        -10.1,
        np.array([[0.1, -0.2, 0.0]]),
        np.array([[0.01, 0.02, 0.0]]),
        np.array([[0.03, 0.0, 0.0]]),
        0.15,
    )
    assert "step=3" in report.getvalue()
    assert "dE=+1.000000e-01" in report.getvalue()
    assert "step_size=1.5000e-01" in report.getvalue()

    job._write_seam_check_log(
        {
            "energy_diff": 1.0e-5,
            "n_projected": 7,
            "lagrange_multiplier": 0.4,
            "n_negative": 1,
            "is_minimum": False,
            "eigenvalues": np.array([-0.2, 0.5]),
        },
        1.0e-3,
        20,
    )
    contents = (Path(job.numfreq_folder) / "crossing_seam_check.log").read_text()
    assert "SADDLE POINT ON SEAM" in contents
    assert "** NEGATIVE **" in contents


def _driver_job(tmp_path, max_steps=2):
    job = object.__new__(GaussianMECPJob)
    job.folder = str(tmp_path)
    job.label = "driver"
    job.molecule = SimpleNamespace(symbols=["H"], positions=np.zeros((1, 3)))
    job.settings = SimpleNamespace(
        step_size=0.1,
        max_steps=max_steps,
        trust_radius=0.3,
        adaptive_step_size=False,
        step_size_method="bb",
        restart=False,
        mecp_numfreq=False,
        energy_diff_tol=1.0e-4,
        force_max_tol=1.0e-3,
        force_rms_tol=1.0e-3,
        disp_max_tol=1.0e-3,
        disp_rms_tol=1.0e-3,
    )
    return job


def test_mecp_driver_writes_final_marker_on_convergence(tmp_path):
    job = _driver_job(tmp_path)
    states = []

    def run_state(positions, step, state):
        states.append((step, state))
        return 0.0, np.zeros((1, 3))

    job._run_state = run_state
    job._mecp_displacement = lambda **kwargs: (
        np.zeros((1, 3)),
        np.zeros((1, 3)),
        np.zeros((1, 3)),
    )
    job._run()

    assert states == [(1, "A"), (1, "B")]
    assert job.molecule.positions == pytest.approx(np.zeros((1, 3)))
    report = Path(job.report_file).read_text()
    assert "Optimization converged at step 1." in report
    assert report.endswith("Converged at step 1.\n")
    assert Path(job.trajectory_file).is_file()
    assert Path(job.final_report_file).is_file()
    assert "seam_minimum=NOT_CHECKED" in Path(job.final_report_file).read_text()
    assert not Path(job.numfreq_folder).exists()


@pytest.mark.parametrize("failed_branch", [None, "plus"])
def test_seam_follow_selects_lower_verified_branch(tmp_path, failed_branch):
    class Molecule:
        symbols = ["H"]

        def __init__(self):
            self.positions = np.zeros((1, 3))

        def copy(self):
            result = Molecule()
            result.positions = self.positions.copy()
            return result

    class FakeFollowingJob(GaussianMECPJob):
        def __init__(self, molecule, settings, label, jobrunner, **kwargs):
            self.molecule = molecule
            self.settings = settings
            self.label = label
            self.jobrunner = jobrunner

        def _optimize_on_seam_progress_plane(self, positions_bohr, **kwargs):
            return positions_bohr

        def verify_seam_minimum(self, **kwargs):
            return {
                "mecp_energy": -1.0,
                "n_negative": 0,
                "frequency_eigenvalues": np.array([1.0e-3]),
                "frequencies": np.array([100.0]),
                "modes": np.array([[1.0, 0.0, 0.0]]),
            }

        def run(self):
            if failed_branch and self.label.endswith(failed_branch):
                raise RuntimeError("Branch did not converge")
            energy = -2.0 if self.label.endswith("minus") else -1.0
            self._last_seam_result = {
                "is_minimum": True,
                "mecp_energy": energy,
                "positions_angstrom": self.molecule.positions.copy(),
            }
            self._final_convergence_metrics = {"energy_diff": 0.0}
            self._final_optimization_steps = 3

    job = FakeFollowingJob(
        molecule=Molecule(),
        settings=GaussianMECPJobSettings(seam_mode_max_steps=1),
        label="crossing",
        jobrunner=object(),
    )
    job.folder = str(tmp_path)
    written = []
    job._write_seam_check_log = lambda result, *args: written.append(
        ("check", result["mecp_energy"])
    )
    job._write_mecp_frequency_log = lambda result, *args: written.append(
        ("freq", result["mecp_energy"])
    )
    initial = {
        "frequencies": np.array([-40.0]),
        "modes": np.array([[1.0, 0.0, 0.0]]),
        "positions_angstrom": np.zeros((1, 3)),
    }

    selected = job._follow_seam_imaginary_mode(initial)

    assert selected["mecp_energy"] == -2.0
    assert selected["optimization_steps"] == 3
    assert selected["convergence_thresholds"]["energy_diff"] == pytest.approx(
        GaussianMECPJobSettings.CONVERGENCE_PRESETS["tight"]["energy_diff_tol"]
    )
    assert selected["convergence_thresholds"]["pgrad_max"] == pytest.approx(
        GaussianMECPJobSettings.CONVERGENCE_PRESETS["tight"]["force_max_tol"]
    )
    assert job._selected_seam_follow_branch == "minus"
    assert job._selected_seam_follow_macro_steps == 1
    assert written == [("check", -2.0), ("freq", -2.0)]
    follow_folder = tmp_path / "crossing_seam_follow"
    summary = (follow_folder / "crossing_seam_follow.log").read_text()
    assert "selected_branch=minus" in summary
    assert "selected_macro_steps=1" in summary
    assert not (tmp_path / "crossing_seam_follow.log").exists()
    for suffix in ("plus", "minus"):
        assert (
            follow_folder / f"crossing_seam_follow_{suffix}_mode_follow.log"
        ).is_file()


def test_seam_follow_without_imaginary_mode_creates_no_directory(tmp_path):
    job = object.__new__(GaussianMECPJob)
    job.folder = str(tmp_path)
    job.label = "crossing"
    result = {
        "frequencies": np.array([20.0]),
        "modes": np.array([[1.0, 0.0, 0.0]]),
    }

    assert job._follow_seam_imaginary_mode(result) is result
    assert not (tmp_path / "crossing_seam_follow").exists()


def test_intermediate_follow_branch_has_no_final_summary(tmp_path):
    job = _driver_job(tmp_path)
    job._is_seam_follow_branch = True
    job._mecp_displacement = lambda **kwargs: (
        np.zeros((1, 3)),
        np.zeros((1, 3)),
        np.zeros((1, 3)),
    )
    job._run_state = lambda positions, step, state: (
        0.0,
        np.zeros((1, 3)),
    )

    job._run()

    assert Path(job.report_file).is_file()
    assert not Path(job.final_report_file).exists()


def test_mecp_report_distinguishes_initial_saddle_and_final_minimum(tmp_path):
    job = _driver_job(tmp_path)
    job.settings.mecp_numfreq = True
    job._mecp_displacement = lambda **kwargs: (
        np.zeros((1, 3)),
        np.zeros((1, 3)),
        np.zeros((1, 3)),
    )
    job._run_state = lambda positions, step, state: (
        0.0,
        np.zeros((1, 3)),
    )
    initial = {
        "is_minimum": False,
        "n_negative": 1,
        "frequencies": np.array([-38.7, 250.0]),
    }
    final = {
        "energy_a": -505.732011947,
        "energy_b": -505.732002706,
        "mecp_energy": -505.7320073265,
        "energy_diff": -9.241e-6,
        "n_negative": 0,
        "is_minimum": True,
        "positions_angstrom": np.zeros((1, 3)),
    }

    def verify(_positions):
        job._initial_seam_result = initial
        job._selected_seam_follow_branch = "minus"
        job._selected_seam_follow_macro_steps = 11
        return final

    job._run_seam_minimum_check = verify
    job._run()

    report = Path(job.report_file).read_text(encoding="utf-8")
    assert "Initial MECP optimization converged at step 1." in report
    assert "Initial seam status: SADDLE" in report
    assert "Initial lowest projected frequency=-38.700000 cm^-1" in report
    assert "Seam following: completed; selected branch=minus" in report
    assert "Seam-following macro steps=11" in report
    assert "Final significant imaginary modes=0" in report
    assert "Final status: VERIFIED MECP MINIMUM" in report
    assert report.endswith("Converged at step 1.\n")
    final_report = Path(job.final_report_file).read_text(encoding="utf-8")
    assert "initial_optimization_steps=1" in final_report
    assert "seam_follow_macro_steps=11" in final_report
    assert "mecp_energy=-505.732007326500 Hartree" in final_report


def test_default_seam_mode_max_steps_is_thirty():
    settings = GaussianMECPJobSettings()
    assert settings.seam_mode_max_steps == 30


def test_final_report_summarizes_selected_structure(tmp_path):
    job = _driver_job(tmp_path)
    job.molecule = SimpleNamespace(
        symbols=["H"], positions=np.array([[0.0, 0.0, 0.0]])
    )
    job._final_convergence_metrics = {
        "energy_diff": 1.0e-6,
        "pgrad_max": 2.0e-4,
        "pgrad_rms": 1.0e-4,
        "disp_max": 3.0e-4,
        "disp_rms": 2.0e-4,
    }
    job._selected_seam_follow_macro_steps = 11
    job._selected_seam_follow_branch = "minus"
    result = {
        "energy_a": -10.0,
        "energy_b": -10.000001,
        "is_minimum": True,
        "n_negative": 0,
        "positions_angstrom": np.array([[1.0, 2.0, 3.0]]),
        "frequencies": np.array([100.0]),
        "convergence_metrics": {
            **job._final_convergence_metrics,
            "pgrad_max": 1.0e-4,
        },
        "convergence_thresholds": {
            "energy_diff": 1.0e-5,
            "pgrad_max": 5.0e-4,
            "pgrad_rms": 1.0e-4,
            "disp_max": 5.0e-4,
            "disp_rms": 3.0e-4,
        },
        "optimization_steps": 4,
    }

    job._write_final_report(15, -11.0, -11.000001, result)

    report = Path(job.final_report_file).read_text(encoding="utf-8")
    assert "initial_optimization_steps=15" in report
    assert "seam_follow_macro_steps=11" in report
    assert "seam_follow_selected_branch=minus" in report
    assert "final_branch_optimization_steps=4" in report
    assert "energy_A=-10.000000000000 Hartree" in report
    assert "pgrad_max: value=1.000000e-04" in report
    assert "pgrad_max: value=1.000000e-04 threshold=5.000000e-04" in report
    assert "seam_minimum=PASS" in report
    assert "H      +1.00000000" in report
    assert "mode    1:    +100.000000" in report


def test_mecp_report_records_failed_seam_verification(tmp_path):
    job = _driver_job(tmp_path)
    job.settings.mecp_numfreq = True
    job._mecp_displacement = lambda **kwargs: (
        np.zeros((1, 3)),
        np.zeros((1, 3)),
        np.zeros((1, 3)),
    )
    job._run_state = lambda positions, step, state: (
        0.0,
        np.zeros((1, 3)),
    )

    def fail(_positions):
        raise RuntimeError("No verified seam minimum")

    job._run_seam_minimum_check = fail
    with pytest.raises(RuntimeError, match="No verified seam minimum"):
        job._run()

    report = Path(job.report_file).read_text(encoding="utf-8")
    assert "Final status: FAILED SEAM VERIFICATION" in report
    assert "Converged at step 1.\n" not in report


def test_mecp_driver_raises_when_max_steps_are_exhausted(tmp_path):
    job = _driver_job(tmp_path, max_steps=1)
    job._run_state = lambda positions, step, state: (
        1.0 if state == "A" else 0.0,
        np.array([[1.0, 0.0, 0.0]])
        if state == "A"
        else np.array([[-1.0, 0.0, 0.0]]),
    )

    with pytest.raises(
        RuntimeError, match="did not converge within max_steps"
    ):
        job._run()

    assert Path(job.state_file).is_file()
