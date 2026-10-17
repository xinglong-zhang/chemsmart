"""Tests for native ORCA SurfCrossOpt support."""

import re
import shutil
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pytest
from click.testing import CliRunner

from chemsmart.cli.orca.orca import orca
from chemsmart.cli.sub import sub
from chemsmart.io.molecules.structure import Molecule
from chemsmart.io.orca.hessian import read_hessian
from chemsmart.io.orca.output import ORCAOutput
from chemsmart.jobs.mecp import MECPSeamMixin
from chemsmart.jobs.orca.job import ORCAGeneralJob
from chemsmart.jobs.orca.mecp import ORCAMECPJob
from chemsmart.jobs.orca.settings import ORCAMECPJobSettings
from chemsmart.jobs.orca.writer import ORCAInputWriter
from chemsmart.settings.orca import ORCAProjectSettings

MECP_CASES_FILE = Path(
    "tests/data/ORCATests/outputs/orca_mecp_cases.txt"
).resolve()


def _write_test_hessian(path, matrix):
    n = len(matrix)
    lines = ["$orca_hessian_file", "$hessian", str(n)]
    for first in range(0, n, 5):
        columns = list(range(first, min(first + 5, n)))
        lines.append(" ".join(map(str, columns)))
        for row in range(n):
            lines.append(
                str(row)
                + " "
                + " ".join(f"{matrix[row, col]:.12E}" for col in columns)
            )
    lines.append("$end")
    path.write_text("\n".join(lines))


def test_orca_blocked_hessian_reader(tmp_path):
    matrix = np.arange(81, dtype=float).reshape(9, 9)
    matrix = matrix + matrix.T
    path = tmp_path / "state.hess"
    _write_test_hessian(path, matrix)
    np.testing.assert_allclose(read_hessian(path, 3), matrix)
    with pytest.raises(ValueError, match="dimension"):
        read_hessian(path, 2)
    text = path.read_text()
    path.write_text(text.replace("$hessian", "$missing"))
    with pytest.raises(ValueError, match="section"):
        read_hessian(path, 3)
    path.write_text(text[: text.index("5 6 7 8")])
    with pytest.raises(ValueError, match="Incomplete"):
        read_hessian(path, 3)


@pytest.mark.parametrize("failure", ["nonfinite", "asymmetric", "duplicate"])
def test_orca_rejects_invalid_hessian(tmp_path, failure):
    matrix = np.eye(9)
    if failure == "nonfinite":
        matrix[0, 0] = np.nan
    if failure == "asymmetric":
        matrix[0, 1] = 1.0
    path = tmp_path / "bad.hess"
    _write_test_hessian(path, matrix)
    if failure == "duplicate":
        path.write_text(path.read_text().replace("5 6 7 8", "0 1 2 3"))
    with pytest.raises(ValueError):
        read_hessian(path, 3)


@pytest.mark.parametrize("method", ["analytic", "numerical"])
@pytest.mark.parametrize("broken_sym", [None, [1, 1]], ids=["ordinary", "bs"])
@pytest.mark.parametrize("energy_gap", [0.0, 2.1142623040759645e-5])
def test_orca_frequency_only_runs_shared_analysis(
    tmp_path, monkeypatch, method, broken_sym, energy_gap
):
    from ase import units

    molecule = Molecule(
        symbols=["O", "H", "H"],
        positions=np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], dtype=float),
    )
    origin = molecule.positions / units.Bohr
    normal = np.arange(1, 10, dtype=float).reshape(3, 3)
    normal /= np.linalg.norm(normal)
    calls = []

    def run_state(child, **kwargs):
        route = child.settings.route_string
        calls.append(route)
        assert "SurfCross" not in route and "Opt" not in route
        prefix = Path(child.folder) / child.label
        prefix.with_suffix(".out").write_text("ORCA TERMINATED NORMALLY\n")
        prefix.with_suffix(".gbw").write_text("mock orbitals")
        if "AnFreq" in route:
            assert "MORead" in route and "%moinp" in route
            if broken_sym and child.label.endswith("_B"):
                assert child.settings.scf_hftyp == "UHF"
                assert child.settings.scf_broken_sym is None
                assert child.settings.multiplicity == 1
                assert child.settings.moinp.endswith("hessian_reference_B.gbw")
            _write_test_hessian(prefix.with_suffix(".hess"), np.eye(9))
        else:
            sign = 1 if child.label.endswith("_A") else -1
            if broken_sym and child.label.endswith("_B"):
                assert child.settings.scf_hftyp == "UHF"
                initial = "hessian_reference_B" in child.label
                assert child.settings.scf_broken_sym == (
                    [1, 1] if initial else None
                )
                assert child.settings.multiplicity == (3 if initial else 1)
            delta = child.molecule.positions / units.Bohr - origin
            energy = sign * np.sum(normal * delta) + 0.5 * np.sum(delta**2)
            energy += sign * energy_gap / 2
            gradient = sign * normal + delta
            lines = [
                "# Number of atoms",
                "#",
                "3",
                "# The current total energy",
                "#",
                str(energy),
                "# The current gradient",
                "#",
            ]
            lines.extend(str(value) for value in gradient.ravel())
            lines += ["# The atomic numbers and current coordinates", "#"]
            for z, xyz in zip(
                [8, 1, 1], child.molecule.positions / units.Bohr
            ):
                lines.append(f"{z} " + " ".join(map(str, xyz)))
            prefix.with_suffix(".engrad").write_text("\n".join(lines))

    monkeypatch.setattr(ORCAGeneralJob, "_run", run_state)
    job = ORCAMECPJob(
        molecule=molecule,
        settings=mecp_settings(
            charge=0,
            multiplicity1=3 if broken_sym else 1,
            multiplicity2=1 if broken_sym else 3,
            broken_sym=broken_sym,
            mode="anafreq" if method == "analytic" else "numfreq",
            mecp_numfreq_only=True,
        ),
        label="water",
        jobrunner=MagicMock(),
    )
    job.folder = str(tmp_path)
    job.run()
    assert job.is_complete()
    assert job._last_seam_result["is_minimum"]
    assert len(job._last_seam_result["frequencies"]) == 2
    assert len(calls) == (4 if method == "analytic" else 38)
    assert (
        f"hessian_method={method}"
        in (tmp_path / "water_mecp_freq.log").read_text()
    )
    np.testing.assert_allclose(job.molecule.positions, molecule.positions)
    job.settings.hessian_method = (
        "numerical" if method == "analytic" else "analytic"
    )
    assert not job.is_complete()


@pytest.mark.parametrize("missing_reference", [False, True])
def test_orca_analytic_failure_does_not_reuse_old_hessian(
    tmp_path, monkeypatch, missing_reference
):
    job = ORCAMECPJob(
        molecule=Molecule(
            symbols=["H", "H"],
            positions=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]),
        ),
        settings=mecp_settings(
            charge=0, multiplicity1=1, multiplicity2=3, mode="anafreq"
        ),
        label="crossing",
        jobrunner=MagicMock(),
    )
    job.folder = str(tmp_path)
    folder = Path(job.frequency_type_folder)
    folder.mkdir()
    reference = folder / "reference_A.gbw"
    if not missing_reference:
        reference.write_text("reference orbitals")
    job._state_checkpoint_files[("check", "A")] = str(reference)
    stale_hessian = folder / "crossing_hessian_analytic_A.hess"
    _write_test_hessian(stale_hessian, np.eye(6))
    calls = []

    def incomplete_run(child, **kwargs):
        calls.append(child.label)
        Path(child.outputfile).write_text("ORCA TERMINATED NORMALLY\n")

    monkeypatch.setattr(ORCAGeneralJob, "_run", incomplete_run)
    with pytest.raises(
        RuntimeError,
        match="reference GBW" if missing_reference else "analytic Hessian",
    ):
        job._compute_analytic_hessians(np.asarray(job.molecule.positions))
    assert len(calls) == (0 if missing_reference else 1)
    assert "A" not in job._analytic_hessians
    if not missing_reference:
        assert not stale_hessian.exists()


@pytest.mark.parametrize("preset", ["standard", "tight"])
def test_orca_gap_default_is_independent_of_shared_presets(preset):
    from chemsmart.jobs.mecp import MECP_CONVERGENCE_PRESETS

    settings = mecp_settings(convergence_preset=preset)
    assert settings.energy_diff_tol == 1e-4
    assert (
        settings.force_max_tol
        == MECP_CONVERGENCE_PRESETS[preset]["force_max_tol"]
    )
    assert MECP_CONVERGENCE_PRESETS["tight"]["energy_diff_tol"] == 1e-5
    assert MECP_CONVERGENCE_PRESETS["standard"]["energy_diff_tol"] == 5e-5
    assert mecp_settings(energy_diff_tol=1e-6).energy_diff_tol == 1e-6


@pytest.mark.parametrize(
    "gap,tangent,tolerance",
    [
        (2e-4, 0.0, 1e-4),
        (2.114262e-5, 0.0, 1e-5),
        (0.0, 0.01, 1e-4),
    ],
)
def test_orca_stationarity_checks_still_reject_before_hessians(
    gap, tangent, tolerance
):
    from types import SimpleNamespace

    job = object.__new__(ORCAMECPJob)
    job.label = "nonstationary"
    job.settings = mecp_settings(
        mode="anafreq", mecp_numfreq_only=True, energy_diff_tol=tolerance
    )
    job.molecule = SimpleNamespace(positions=np.zeros((1, 3)))
    job._run_state = lambda positions, step, state, **kwargs: (
        (gap, np.array([[1.0, tangent, 0.0]]))
        if state == "A"
        else (0.0, np.array([[-1.0, tangent, 0.0]]))
    )

    def forbidden(*args, **kwargs):
        pytest.fail("Hessians must not run on a nonstationary crossing")

    job._get_state_hessians = forbidden
    with pytest.raises(RuntimeError, match="stationary crossing"):
        job.verify_seam_minimum(write_frequencies=True)


def test_both_backends_use_same_follow_algorithm():
    from chemsmart.jobs.gaussian.mecp import GaussianMECPJob

    assert (
        ORCAMECPJob._follow_seam_imaginary_mode
        is MECPSeamMixin._follow_seam_imaginary_mode
    )
    assert (
        GaussianMECPJob._follow_seam_imaginary_mode
        is MECPSeamMixin._follow_seam_imaginary_mode
    )
    assert (
        ORCAMECPJob._get_state_hessians is GaussianMECPJob._get_state_hessians
    )


@pytest.mark.parametrize("fail_plus", [False, True])
@pytest.mark.parametrize("broken_sym", [None, [1, 1]], ids=["ordinary", "bs"])
def test_orca_followseam_releases_and_selects_verified_branch(
    tmp_path, monkeypatch, fail_plus, broken_sym
):
    from types import SimpleNamespace

    from chemsmart.jobs.orca.job import ORCAJob

    molecule = Molecule(
        symbols=["O", "H", "H"],
        positions=np.array(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
        ),
    )
    normal = np.arange(1, 10, dtype=float).reshape(3, 3)
    normal /= np.linalg.norm(normal)
    native_calls = []
    hessian_calls = []

    def native_run(job, **kwargs):
        assert job.settings.broken_sym == broken_sym
        native_calls.append(job.label)
        if fail_plus and job.label.endswith("plus"):
            raise RuntimeError("Synthetic native optimization failure")

    def state_run(
        job,
        positions_bohr,
        step_idx,
        state,
        checkpoint_tag=None,
        analytic=False,
    ):
        is_branch = "seam_follow" in job.label
        if analytic:
            hessian_calls.append((job.label, step_idx, state))
            job._analytic_hessians[state] = np.eye(9) * (
                1 if is_branch else -1
            )
            return None, None
        energy = -2.0 if "minus" in job.label else -1.0
        return energy, normal if state == "A" else -normal

    monkeypatch.setattr(ORCAJob, "_run", native_run)
    monkeypatch.setattr(ORCAMECPJob, "_run_state", state_run)
    monkeypatch.setattr(
        ORCAMECPJob,
        "results",
        property(
            lambda job: SimpleNamespace(
                normal_termination=True,
                converged=True,
                final_structure=job.molecule.copy(),
                energy_gap_history=[0.0],
            )
        ),
    )
    job = ORCAMECPJob(
        molecule=molecule,
        settings=mecp_settings(
            charge=0,
            multiplicity1=3 if broken_sym else 1,
            multiplicity2=1 if broken_sym else 3,
            broken_sym=broken_sym,
            mode="followseam",
            hessian_method="analytic",
            seam_mode_max_steps=2,
        ),
        label="crossing",
        jobrunner=MagicMock(),
    )
    job.folder = str(tmp_path)
    job.run()
    assert job.is_complete()
    assert job._selected_seam_follow_branch == "minus"
    assert job._last_seam_result["is_minimum"]
    assert job._last_seam_result["mecp_energy"] == -2.0
    assert native_calls == [
        "crossing",
        "crossing_seam_follow_plus",
        "crossing_seam_follow_minus",
    ]
    assert any("macro01" in step for _, step, _ in hessian_calls)
    assert len(hessian_calls) == (8 if fail_plus else 10)
    assert (tmp_path / "crossing_mecp_freq.log").is_file()


@pytest.mark.parametrize(
    "options",
    [
        {"mode": "anafreq", "hessian_method": "numerical"},
        {"mode": "followseam", "mecp_numfreq_only": True},
        {"mode": "anafreq", "moinp": "pes2.gbw"},
        {"mode": "followseam", "casscf_nel": 2, "casscf_norb": 2},
    ],
)
def test_unsupported_seam_requests_fail_explicitly(options):
    with pytest.raises(ValueError):
        mecp_settings(**options).validate()


@pytest.fixture
def mecp_case(tmp_path):
    """Materialize one trimmed real ORCA output from the combined fixture."""
    contents = MECP_CASES_FILE.read_text(encoding="utf-8")

    def load(name):
        marker = f"===== {name} ====="
        section = contents.split(marker, 1)[1].split("===== END =====", 1)[0]
        output = tmp_path / f"{name}.out"
        output.write_text(section.strip() + "\n", encoding="utf-8")
        return output

    return load


def mecp_settings(**kwargs):
    values = {
        "charge": 1,
        "multiplicity1": 6,
        "multiplicity2": 4,
        "functional": "B3LYP",
        "basis": "TZVP",
    }
    values.update(kwargs)
    return ORCAMECPJobSettings(**values)


def test_route_and_validation():
    settings = mecp_settings()
    assert settings.route_string.startswith("! Opt SurfCrossOpt")
    assert settings.validate().multiplicity == 6
    assert settings.maxiter == 200
    with pytest.raises(ValueError, match="different"):
        mecp_settings(multiplicity1=4, multiplicity2=4).validate()


def test_numfreq_and_custom_route():
    settings = mecp_settings(
        mode="NumFreq", route_to_be_written="B3LYP TZVP TightSCF"
    )
    route = settings.route_string
    assert route.startswith("! Opt SurfCrossOpt SurfCrossNumFreq")
    assert "TightSCF" in route


def test_feo_official_example_input(
    tmpdir,
    orca_yaml_settings_orca_project_name,
    orca_jobrunner_no_scratch,
):
    project = ORCAProjectSettings.from_project(
        orca_yaml_settings_orca_project_name
    )
    settings = ORCAMECPJobSettings.from_settings(project.opt_settings())
    settings.charge = 1
    settings.multiplicity1 = 6
    settings.multiplicity2 = 4
    settings.maxiter = 200
    settings.validate()
    feo_xyz = Path("tests/data/ORCATests/inputs/xyz/feo_plus.xyz").resolve()
    job = ORCAMECPJob.from_filename(
        filename=str(feo_xyz),
        settings=settings,
        label="feo_mecp",
        jobrunner=orca_jobrunner_no_scratch,
    )
    ORCAInputWriter(job=job).write(target_directory=tmpdir)
    content = Path(str(tmpdir), "feo_mecp.inp").read_text()
    assert "Opt SurfCrossOpt" in content
    assert "%mecp\n  Mult 4\nend" in content
    assert "%geom\n  MaxIter 200\nend" in content
    assert "* xyz 1 6" in content


def test_casscf_options_are_written(
    tmpdir, single_molecule_xyz_file, orca_jobrunner_no_scratch
):
    settings = mecp_settings(
        casscf_nel=6,
        casscf_norb=6,
        casscf_mult=[1, 3],
        casscf_nroots=[4, 2],
        casscf_bweight=[2, 1],
    )
    job = ORCAMECPJob.from_filename(
        filename=single_molecule_xyz_file,
        settings=settings,
        label="casscf_mecp",
        jobrunner=orca_jobrunner_no_scratch,
    )
    ORCAInputWriter(job=job).write(target_directory=tmpdir)
    content = Path(str(tmpdir), "casscf_mecp.inp").read_text()
    assert "casscf_nel 6" in content
    assert "casscf_mult 1,3" in content
    assert "casscf_nroots 4,2" in content


def test_real_output_markers_are_parsed(tmp_path):
    output = tmp_path / "feo_mecp.out"
    output.write_text("""|  1> ! B3LYP TZVP Opt SurfCrossOpt SurfCrossNumFreq
|  2> %mecp Mult 4
FINAL SINGLE POINT ENERGY     -1340.123456789
FINAL SINGLE POINT ENERGY     -1340.123454391
Energy difference between both states      -0.000002398
THE OPTIMIZATION HAS CONVERGED
****ORCA TERMINATED NORMALLY****
""")
    parsed = ORCAOutput(str(output)).mecp_result
    assert parsed.converged
    assert parsed.normal_termination
    assert parsed.state_1_energy == pytest.approx(-1340.123456789)
    assert parsed.state_2_energy == pytest.approx(-1340.123454391)
    assert parsed.energy_gap == pytest.approx(-0.000002398)
    assert parsed.numfreq_requested


def test_real_feo_stationary_point_geometry(mecp_case):
    output = mecp_case("feo_plus")
    parsed = ORCAOutput(str(output)).mecp_result
    assert parsed.converged
    assert parsed.normal_termination
    assert parsed.final_structure.energy is None
    fe_o_distance = np.linalg.norm(
        parsed.final_structure.positions[0]
        - parsed.final_structure.positions[1]
    )
    assert fe_o_distance == pytest.approx(1.993648, abs=1e-6)


def test_numfreq_rejects_two_atom_system(orca_jobrunner_no_scratch):
    with pytest.raises(
        ValueError, match="SurfCrossNumFreq requires at least 3 atoms"
    ):
        ORCAMECPJob.from_filename(
            filename=str(
                Path("tests/data/ORCATests/inputs/xyz/feo_plus.xyz").resolve()
            ),
            settings=mecp_settings(mode="numfreq"),
            label="feo_numfreq",
            jobrunner=orca_jobrunner_no_scratch,
        )


def test_broken_symmetry_rejects_incompatible_pes2(
    orca_jobrunner_no_scratch,
):
    with pytest.raises(ValueError, match="generates PES2 multiplicity 1"):
        ORCAMECPJob.from_filename(
            filename=str(
                Path("tests/data/ORCATests/inputs/xyz/feo_plus.xyz").resolve()
            ),
            settings=mecp_settings(broken_sym=[1, 1]),
            label="invalid_feo_broken_symmetry",
            jobrunner=orca_jobrunner_no_scratch,
        )


def test_broken_symmetry_rejects_incompatible_electron_parity(
    orca_jobrunner_no_scratch,
):
    with pytest.raises(ValueError, match="incompatible.*33 electrons"):
        ORCAMECPJob.from_filename(
            filename=str(
                Path("tests/data/ORCATests/inputs/xyz/feo_plus.xyz").resolve()
            ),
            settings=mecp_settings(
                multiplicity2=1,
                broken_sym=[1, 1],
            ),
            label="invalid_feo_broken_symmetry_parity",
            jobrunner=orca_jobrunner_no_scratch,
        )


def test_twisted_ethylene_broken_symmetry_input(
    tmp_path, orca_jobrunner_no_scratch
):
    job = ORCAMECPJob.from_filename(
        filename=str(
            Path(
                "tests/data/ORCATests/inputs/xyz/ethylene_twisted_mecp.xyz"
            ).resolve()
        ),
        settings=mecp_settings(
            charge=0,
            multiplicity1=3,
            multiplicity2=1,
            functional="B3LYP",
            basis="def2-SVP",
            scf_tol="TightSCF",
            broken_sym=[1, 1],
            maxiter=50,
        ),
        label="ethylene_twisted_mecp",
        jobrunner=orca_jobrunner_no_scratch,
    )
    ORCAInputWriter(job=job).write(target_directory=tmp_path)
    input_text = (tmp_path / "ethylene_twisted_mecp.inp").read_text()
    assert "B3LYP def2-SVP TightSCF" in input_text
    assert "%mecp\n  Mult 1\n  brokenSym 1,1\nend" in input_text
    assert "* xyz 0 3" in input_text


def test_bs_state_initialization_and_hessian_inputs(
    tmp_path, orca_jobrunner_no_scratch
):
    molecule = Molecule(
        symbols=["O", "H", "H"],
        positions=np.array(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
        ),
    )
    job = ORCAMECPJob(
        molecule=molecule,
        settings=mecp_settings(
            charge=0,
            multiplicity1=3,
            multiplicity2=1,
            broken_sym=[1, 1],
            mode="anafreq",
            scf_maxiter=350,
            scf_convergence="tight",
            additional_route_parameters="RKS",
        ),
        label="crossing",
        jobrunner=orca_jobrunner_no_scratch,
    )
    job.folder = str(tmp_path)
    settings_a = job._state_settings("A", "check", False)
    assert settings_a.scf_broken_sym is None
    assert settings_a.scf_hftyp is None
    settings_b = job._state_settings("B", "check", False)
    first = job._state_job(
        np.zeros((3, 3)), "hessian_reference", "B", "check", settings_b
    )
    ORCAInputWriter(job=first).write(target_directory=first.folder)
    text = Path(first.inputfile).read_text()
    assert text.count("%scf") == 1
    assert "HFTyp UHF" in text and "BrokenSym 1,1" in text
    assert "maxiter 350" in text and "convergence tight" in text
    assert "* xyz 0 3" in text
    assert "%mecp" not in text and not re.search(r"\bRKS\b", text)

    Path(first.gbwfile).write_text("converged BS reference orbitals")
    job._state_checkpoint_files[("check", "B")] = first.gbwfile
    settings_b = job._state_settings("B", "check", True)
    second = job._state_job(
        np.zeros((3, 3)), "hessian_analytic", "B", "check", settings_b
    )
    ORCAInputWriter(job=second).write(target_directory=second.folder)
    text = Path(second.inputfile).read_text()
    assert "AnFreq" in text and "MORead" in text
    assert "* xyz 0 1" in text and "HFTyp UHF" in text
    assert 'moinp "crossing_hessian_reference_B.gbw"' in text
    assert "BrokenSym" not in text and "%mecp" not in text
    assert job.settings.multiplicity == 3
    assert job.settings.broken_sym == [1, 1]


def test_bs_state_requires_its_generated_orbitals(tmp_path, monkeypatch):
    job = ORCAMECPJob(
        molecule=Molecule(symbols=["H", "H"], positions=np.zeros((2, 3))),
        settings=mecp_settings(
            charge=0,
            multiplicity1=3,
            multiplicity2=1,
            broken_sym=[1, 1],
            mode="anafreq",
        ),
        label="crossing",
        jobrunner=MagicMock(),
    )
    job.folder = str(tmp_path)

    def no_orbitals(child, **kwargs):
        Path(child.outputfile).write_text("ORCA TERMINATED NORMALLY\n")

    monkeypatch.setattr(ORCAGeneralJob, "_run", no_orbitals)
    with pytest.raises(
        RuntimeError, match="Missing broken-symmetry state B orbitals"
    ):
        job._run_state(np.zeros((2, 3)), "hessian_reference", "B", "check")
    assert ("check", "B") not in job._state_checkpoint_files


def test_orca_official_ch3o_ch2oh_numfreq_input(
    tmpdir, orca_jobrunner_no_scratch
):
    xyz = Path("tests/data/ORCATests/inputs/xyz/ch3o_ch2oh_mecp.xyz").resolve()
    settings = mecp_settings(
        charge=1,
        multiplicity1=3,
        multiplicity2=1,
        mode="numfreq",
    )
    job = ORCAMECPJob.from_filename(
        filename=str(xyz),
        settings=settings,
        label="ch3o_ch2oh_mecp",
        jobrunner=orca_jobrunner_no_scratch,
    )
    ORCAInputWriter(job=job).write(target_directory=tmpdir)
    content = Path(str(tmpdir), "ch3o_ch2oh_mecp.inp").read_text()
    assert "Opt SurfCrossOpt SurfCrossNumFreq" in content
    assert "%mecp\n  Mult 1\nend" in content
    assert "* xyz 1 3" in content


def test_real_ch3o_ch2oh_numfreq_result(mecp_case):
    output = mecp_case("ch3o_ch2oh_numfreq")
    result = ORCAOutput(str(output)).mecp_result
    assert result.numfreq_requested
    assert result.numfreq_completed
    assert result.is_minimum is True
    assert len(result.state_1_frequencies) == 15
    assert len(result.state_2_frequencies) == 15
    assert result.state_1_frequencies[7] == pytest.approx(775.23)
    assert result.state_2_frequencies[7] == pytest.approx(601.21)
    assert result.state_1_imaginary_frequencies == ()
    assert result.state_2_imaginary_frequencies == ()


def test_orca_mecp_quality_report(
    tmp_path, orca_jobrunner_no_scratch, mecp_case
):
    source = mecp_case("ch3o_ch2oh_numfreq")
    xyz = Path("tests/data/ORCATests/inputs/xyz/ch3o_ch2oh_mecp.xyz").resolve()
    job = ORCAMECPJob.from_filename(
        filename=str(xyz),
        settings=mecp_settings(
            charge=1,
            multiplicity1=3,
            multiplicity2=1,
            mode="numfreq",
        ),
        label="ch3o_ch2oh_mecp",
        jobrunner=orca_jobrunner_no_scratch,
    )
    job.set_folder(str(tmp_path))
    shutil.copy(source, job.outputfile)

    report_path = job.log_result(energy_gap_tolerance=2.0e-4)
    report = Path(report_path).read_text(encoding="utf-8")

    assert Path(report_path).name == "ch3o_ch2oh_mecp_report.log"
    assert "CHEMSMART ORCA MECP quality report" in report
    assert "status=PASSED" in report
    assert "normal_termination=True" in report
    assert "optimization_converged=True" in report
    assert "energy_gap_accepted=True" in report
    assert "numfreq_completed=True" in report
    assert "is_minimum=True" in report
    assert "label=ch3o_ch2oh_mecp" in report
    assert "backend=orca" in report
    assert "energy_A=" in report and "energy_B=" in report
    assert "energy_difference=" in report and " Hartree" in report


def test_orca_mecp_quality_report_warns_for_large_final_gap(
    tmp_path, orca_jobrunner_no_scratch
):
    xyz = Path("tests/data/ORCATests/inputs/xyz/ch3o_ch2oh_mecp.xyz").resolve()
    job = ORCAMECPJob.from_filename(
        filename=str(xyz),
        settings=mecp_settings(
            charge=1,
            multiplicity1=3,
            multiplicity2=1,
        ),
        label="large_gap_mecp",
        jobrunner=orca_jobrunner_no_scratch,
    )
    job.set_folder(str(tmp_path))
    Path(job.outputfile).write_text(
        """|  1> ! B3LYP TZVP Opt SurfCrossOpt
|  2> %mecp Mult 1
FINAL SINGLE POINT ENERGY     -114.621194606632
FINAL SINGLE POINT ENERGY     -114.621307200948
Energy difference between both states        0.000112594
THE OPTIMIZATION HAS CONVERGED
****ORCA TERMINATED NORMALLY****
""",
        encoding="utf-8",
    )

    report = Path(job.log_result()).read_text(encoding="utf-8")

    assert "status=WARNING" in report
    assert "energy_gap_accepted=False" in report
    assert "final two-state energy gap exceeds" in report


def test_orca_runner_writes_mecp_report_after_run(
    orca_jobrunner_no_scratch,
):
    job = MagicMock()
    job.TYPE = "orcamecp"

    orca_jobrunner_no_scratch._postrun(job)

    job.log_result.assert_called_once_with()


def test_numfreq_imaginary_mode_is_not_minimum(tmp_path, mecp_case):
    source = mecp_case("ch3o_ch2oh_numfreq").read_text()
    output = tmp_path / "mecp_imaginary.out"
    output.write_text(source.replace("601.21 cm**-1", "-42.00 cm**-1"))
    result = ORCAOutput(str(output)).mecp_result
    assert result.numfreq_completed
    assert result.state_2_imaginary_frequencies == (-42.0,)
    assert result.is_minimum is False


@pytest.mark.parametrize(
    "surface", ["VIBRATIONAL FREQUENCIES", "VIBRATIONAL FREQUENCIES PES2"]
)
@pytest.mark.parametrize(
    "damage", ["missing", "duplicate", "empty", "truncated", "last_empty"]
)
def test_numfreq_incomplete_table_is_not_completed(
    tmp_path, surface, damage, mecp_case
):
    source = mecp_case("ch3o_ch2oh_numfreq").read_text(encoding="utf-8")
    lines = source.splitlines(keepends=True)
    start = next(i for i, line in enumerate(lines) if line.strip() == surface)
    rows = []
    for i in range(start + 1, len(lines)):
        if re.match(r"^\s*\d+:.*cm\*\*-1", lines[i]):
            rows.append(i)
        elif rows:
            break
    assert len(rows) == 15
    if damage == "missing":
        del lines[rows[7]]
    elif damage == "duplicate":
        # Preserve the count while replacing mode 7 with a second mode 6.
        lines[rows[7]] = lines[rows[6]]
    elif damage == "empty":
        del lines[rows[0] : rows[-1] + 1]
    elif damage == "truncated":
        lines = lines[: rows[7]]
    else:
        lines.append(f"\n{surface}\n-----------------------\n")
    output = tmp_path / "incomplete.out"
    output.write_text("".join(lines), encoding="utf-8")
    result = ORCAOutput(str(output)).mecp_result
    assert result.numfreq_requested
    assert result.numfreq_completed is False
    assert result.is_minimum is None
    if damage != "truncated":
        # Even a normal-termination marker must not hide missing modes.
        assert result.normal_termination


def test_real_twisted_ethylene_broken_sym_numfreq_result(mecp_case):
    """The saved ORCA 6.1.0 HPC run has a PES2 imaginary mode.

    Source: orca_mecp_test/orca_mecp/testnumfreq200/
    ethylene_twisted_mecp_max200.out. This checks a recorded result,
    not whether every new ethylene calculation reproduces this mode.
    """
    source = mecp_case("ethylene_twisted_numfreq")
    assert "brokenSym 1,1" in source.read_text(encoding="utf-8")
    result = ORCAOutput(str(source.resolve())).mecp_result
    assert result.normal_termination
    assert result.converged
    assert result.numfreq_requested
    assert result.numfreq_completed
    assert len(result.state_1_frequencies) == 18
    assert len(result.state_2_frequencies) == 18
    assert result.state_1_imaginary_frequencies == ()
    assert result.state_2_imaginary_frequencies == pytest.approx((-969.14,))
    assert result.is_minimum is False


def test_real_twisted_ethylene_numfreq_report_warns(
    tmp_path, orca_jobrunner_no_scratch, mecp_case
):
    source = mecp_case("ethylene_twisted_numfreq")
    xyz = Path(
        "tests/data/ORCATests/inputs/xyz/ethylene_twisted_mecp.xyz"
    ).resolve()
    job = ORCAMECPJob.from_filename(
        filename=str(xyz),
        settings=mecp_settings(
            charge=0,
            multiplicity1=3,
            multiplicity2=1,
            basis="def2-SVP",
            mode="numfreq",
            broken_sym=[1, 1],
            maxiter=200,
        ),
        label="ethylene_numfreq",
        jobrunner=orca_jobrunner_no_scratch,
    )
    job.set_folder(str(tmp_path))
    shutil.copy(source, job.outputfile)
    report = Path(job.log_result()).read_text(encoding="utf-8")
    assert "status=WARNING" in report
    assert "optimization_converged=True" in report
    assert "numfreq_completed=True" in report
    assert "is_minimum=False" in report
    assert "state_2_imaginary_frequencies_cm-1=[-969.14]" in report
    assert "imaginary mode detected on the crossing hyperline" in report


@pytest.mark.parametrize(
    "case,energies,gap,imaginary,is_minimum,status",
    [
        (
            "co",
            (-651.748192575514, -651.748053038959),
            -0.000139537,
            (-16.03,),
            False,
            "WARNING",
        ),
        (
            "ooh",
            (-651.746166826518, -651.746242182235),
            0.000075356,
            (),
            True,
            "PASSED",
        ),
    ],
    ids=["orcacase3-co-warning", "orcacase4-ooh-passed"],
)
def test_real_pathway_numfreq_results_and_reports(
    case,
    energies,
    gap,
    imaginary,
    is_minimum,
    status,
    tmp_path,
    orca_jobrunner_no_scratch,
    mecp_case,
):
    """Regress the recorded HPC results, including the default gap cutoff."""
    source = mecp_case(f"{case}_pathway_numfreq")
    text = source.read_text(encoding="utf-8")
    assert "brokenSym 1,1" in text
    assert "SMD(acetonitrile)" in text
    result = ORCAOutput(str(source)).mecp_result
    assert result.normal_termination
    assert result.converged
    assert result.numfreq_requested
    assert result.numfreq_completed
    assert result.state_1_energy == pytest.approx(energies[0], abs=1e-10)
    assert result.state_2_energy == pytest.approx(energies[1], abs=1e-10)
    assert result.energy_gap == pytest.approx(gap, abs=1e-12)
    assert len(result.final_structure) == 26
    assert len(result.state_1_frequencies) == 78
    assert len(result.state_2_frequencies) == 78
    assert result.state_1_imaginary_frequencies == ()
    assert result.state_2_imaginary_frequencies == pytest.approx(imaginary)
    assert result.is_minimum is is_minimum

    job = ORCAMECPJob(
        molecule=result.final_structure,
        settings=mecp_settings(
            charge=0,
            multiplicity1=3,
            multiplicity2=1,
            functional="M062X",
            basis="maug-cc-pV(D+d)Z",
            solvent_model="smd",
            solvent_id="acetonitrile",
            mode="numfreq",
            broken_sym=[1, 1],
            maxiter=200,
        ),
        label=f"{case}_pathway",
        jobrunner=orca_jobrunner_no_scratch,
    )
    job.set_folder(str(tmp_path))
    shutil.copy(source, job.outputfile)
    report = Path(job.log_result()).read_text(encoding="utf-8")
    assert f"status={status}" in report
    assert "normal_termination=True" in report
    assert "optimization_converged=True" in report
    assert "numfreq_completed=True" in report
    assert f"is_minimum={is_minimum}" in report
    assert f"energy_gap_accepted={case == 'ooh'}" in report
    if case == "co":
        assert "final two-state energy gap exceeds" in report
        assert "imaginary mode detected on the crossing hyperline" in report
    else:
        assert "All requested MECP quality checks passed." in report


def test_mecp_solvent_constraint_and_moinp_are_written(
    tmpdir, tmp_path, orca_jobrunner_no_scratch
):
    xyz = Path("tests/data/ORCATests/inputs/xyz/ch3o_ch2oh_mecp.xyz").resolve()
    gbw = tmp_path / "pes2.gbw"
    gbw.write_bytes(b"test orbitals")
    settings = mecp_settings(
        charge=1,
        multiplicity1=3,
        multiplicity2=1,
        solvent_model="cpcm",
        solvent_id="water",
        moinp=str(gbw),
        maxiter=80,
        invert_constraints=True,
    )
    job = ORCAMECPJob.from_filename(
        filename=str(xyz),
        settings=settings,
        label="advanced_mecp",
        jobrunner=orca_jobrunner_no_scratch,
    )
    job.molecule.frozen_atoms = [-1, 0, 0, 0, 0]
    ORCAInputWriter(job=job).write(target_directory=tmpdir)
    content = Path(str(tmpdir), "advanced_mecp.inp").read_text()
    assert 'moinp "pes2.gbw"' in content
    assert Path(str(tmpdir), "pes2.gbw").read_bytes() == b"test orbitals"
    assert "CPCM(water)" in content
    assert content.count("%geom") == 1
    assert (
        "%geom\n  MaxIter 80\n  { C 0 C }\n  InvertConstraints True\nend"
        in content
    )


def test_cli_rejects_equal_multiplicities(
    single_molecule_xyz_file, run_orca_and_capture_settings
):
    result, _ = run_orca_and_capture_settings(
        "chemsmart.jobs.orca.mecp.ORCAMECPJob",
        [
            "-p",
            "test",
            "-f",
            single_molecule_xyz_file,
            "mecp",
            "--multiplicity1",
            "4",
            "--multiplicity2",
            "4",
            "--charge",
            "1",
        ],
    )
    assert result.exit_code != 0


@pytest.mark.parametrize(
    "m1,m2,expected",
    [
        ("0", "4", "range x>=1"),
        ("-1", "4", "range x>=1"),
        ("6", "0", "range x>=1"),
        ("6", "-1", "range x>=1"),
    ],
)
def test_cli_rejects_nonpositive_multiplicities(
    single_molecule_xyz_file,
    run_orca_and_capture_settings,
    m1,
    m2,
    expected,
):
    result, _ = run_orca_and_capture_settings(
        "chemsmart.jobs.orca.mecp.ORCAMECPJob",
        [
            "-p",
            "test",
            "-f",
            single_molecule_xyz_file,
            "-c",
            "1",
            "mecp",
            "--multiplicity1",
            m1,
            "--multiplicity2",
            m2,
        ],
    )
    assert result.exit_code != 0
    assert expected in result.output


@pytest.mark.parametrize(
    "args", [["--multiplicity1", "6"], ["--multiplicity2", "4"], []]
)
def test_cli_requires_both_multiplicities(
    single_molecule_xyz_file, run_orca_and_capture_settings, args
):
    result, _ = run_orca_and_capture_settings(
        "chemsmart.jobs.orca.mecp.ORCAMECPJob",
        [
            "-p",
            "test",
            "-f",
            single_molecule_xyz_file,
            "-c",
            "1",
            "mecp",
            *args,
        ],
    )
    assert result.exit_code != 0
    assert "Missing option" in result.output


def test_orca_short_a_appends_label(
    single_molecule_xyz_file, orca_jobrunner_no_scratch
):
    with patch("chemsmart.jobs.orca.mecp.ORCAMECPJob") as job_class:
        job_class.return_value = object()
        result = CliRunner().invoke(
            orca,
            [
                "-p",
                "test",
                "-f",
                single_molecule_xyz_file,
                "-a",
                "testnumfreq",
                "-c",
                "1",
                "mecp",
                "--multiplicity1",
                "6",
                "--multiplicity2",
                "4",
                "--mode",
                "numfreq",
                "--freeze-atoms",
                "1",
            ],
            obj={"jobrunner": orca_jobrunner_no_scratch},
            catch_exceptions=False,
        )

    assert result.exit_code == 0, result.output
    call = job_class.call_args.kwargs
    assert call["label"].endswith("_testnumfreq")
    assert call["settings"].aux_basis is None
    assert call["molecule"].frozen_atoms[0] == -1


def test_sub_preserves_mecp_arguments(orca_jobrunner_no_scratch):
    xyz = Path("tests/data/ORCATests/inputs/xyz/ch3o_ch2oh_mecp.xyz").resolve()
    server = orca_jobrunner_no_scratch.server
    with (
        patch(
            "chemsmart.cli.sub.Server.from_servername",
            return_value=server,
        ),
        patch.object(server, "submit") as submit,
        patch("chemsmart.jobs.orca.mecp.ORCAMECPJob") as job_class,
    ):
        job_class.return_value = MagicMock()
        result = CliRunner().invoke(
            sub,
            [
                "--server",
                "cuhk",
                "--test",
                "orca",
                "-p",
                "test",
                "-f",
                str(xyz),
                "-a",
                "testnumfreq",
                "-c",
                "1",
                "-x",
                "B3LYP",
                "-b",
                "TZVP",
                "mecp",
                "--multiplicity1",
                "3",
                "--multiplicity2",
                "1",
                "--mode",
                "numfreq",
                "--freeze-atoms",
                "1",
            ],
            catch_exceptions=False,
        )

    assert result.exit_code == 0, result.output
    submitted_args = submit.call_args.kwargs["cli_args"]
    for option, value in (
        ("--append-label", "testnumfreq"),
        ("--multiplicity1", "3"),
        ("--multiplicity2", "1"),
        ("--mode", "numfreq"),
        ("--freeze-atoms", "1"),
    ):
        index = submitted_args.index(option)
        assert submitted_args[index + 1] == value
    assert "--aux-basis" not in submitted_args


def test_cli_numbered_state_options(
    single_molecule_xyz_file,
    run_orca_and_capture_settings,
    orca_jobrunner_no_scratch,
):
    result, settings = run_orca_and_capture_settings(
        "chemsmart.jobs.orca.mecp.ORCAMECPJob",
        [
            "-p",
            "test",
            "-f",
            single_molecule_xyz_file,
            "-c",
            "1",
            "mecp",
            "-m1",
            "6",
            "-m2",
            "4",
        ],
        ctx_obj={"jobrunner": orca_jobrunner_no_scratch},
    )
    assert result.exit_code == 0, result.output
    assert settings.multiplicity == 6
    assert settings.multiplicity1 == 6
    assert settings.multiplicity2 == 4


@pytest.mark.parametrize(
    "mode_args, expected",
    [
        ([], "opt"),
        (["opt"], "opt"),
        (["numfreq"], "numfreq"),
        (["--mode", "numfreq"], "numfreq"),
        (["NUMFREQ", "--mode", "numfreq"], "numfreq"),
    ],
)
def test_mecp_positional_mode(
    mode_args,
    expected,
    single_molecule_xyz_file,
    run_orca_and_capture_settings,
    orca_jobrunner_no_scratch,
):
    result, settings = run_orca_and_capture_settings(
        "chemsmart.jobs.orca.mecp.ORCAMECPJob",
        [
            "-p",
            "test",
            "-f",
            single_molecule_xyz_file,
            "-c",
            "1",
            "mecp",
            "-m1",
            "3",
            "-m2",
            "1",
            "-bs",
            "1,1",
        ]
        + mode_args,
        ctx_obj={"jobrunner": orca_jobrunner_no_scratch},
    )
    assert result.exit_code == 0, result.output
    assert settings.mode == expected
    assert settings.broken_sym == [1, 1]
    assert ("SurfCrossNumFreq" in settings.route_string) == (
        expected == "numfreq"
    )


@pytest.mark.parametrize("args", [["numfrq"], ["opt", "--mode", "numfreq"]])
def test_mecp_invalid_mode(args):
    from chemsmart.cli.orca.mecp import mecp

    result = CliRunner().invoke(mecp, ["-m1", "3", "-m2", "1"] + args)
    assert result.exit_code == 2


@pytest.mark.parametrize(
    "args,method,only",
    [
        (["anafreq"], "analytic", False),
        (["anafreq", "--only"], "analytic", True),
        (["anafreq", "--only", "--energy-diff-tol", "1e-5"], "analytic", True),
        (["numfreq", "--only"], "numerical", True),
        (["followseam"], "numerical", False),
        (["followseam", "--hessian-method", "analytic"], "analytic", False),
    ],
)
def test_mecp_shared_workflow_cli(
    args,
    method,
    only,
    single_molecule_xyz_file,
    run_orca_and_capture_settings,
    orca_jobrunner_no_scratch,
):
    result, settings = run_orca_and_capture_settings(
        "chemsmart.jobs.orca.mecp.ORCAMECPJob",
        [
            "-p",
            "test",
            "-f",
            single_molecule_xyz_file,
            "-c",
            "1",
            "mecp",
            "-m1",
            "3",
            "-m2",
            "1",
        ]
        + args,
        ctx_obj={"jobrunner": orca_jobrunner_no_scratch},
    )
    assert result.exit_code == 0, result.output
    assert settings.hessian_method == method
    assert settings.mecp_numfreq_only == only
    assert settings.energy_diff_tol == (
        1e-5 if "--energy-diff-tol" in args else 1e-4
    )
    assert settings.uses_seam_workflow
    if args[0] != "numfreq":
        assert "SurfCrossNumFreq" not in settings.route_string
    assert settings.follow_seam_imaginary_mode == (args[0] == "followseam")


def test_mecp_positional_mode_command_reconstruction():
    from chemsmart.cli.orca.mecp import mecp
    from chemsmart.utils.cli import (
        CtxObjArguments,
        _add_subcommand_info_to_ctx,
    )

    with mecp.make_context("mecp", ["-m1", "3", "-m2", "1", "numfreq"]) as ctx:
        _add_subcommand_info_to_ctx(ctx)
        arguments = CtxObjArguments(
            ctx.obj["subcommand"]
        ).reconstruct_command_line()
    assert "numfreq" in arguments
    assert "--mode-argument" not in arguments
    with mecp.make_context("mecp", arguments[1:]) as ctx:
        assert ctx.params["mode_argument"] == "numfreq"
