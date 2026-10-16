from pathlib import Path
from unittest.mock import MagicMock

import pytest
from click.testing import CliRunner

from chemsmart.analysis.redox import (
    RedoxReference,
    compute_redox_potential,
    format_redox_summary,
    get_redox_reference,
    list_redox_references,
    register_redox_reference,
)
from chemsmart.cli.run import run
from chemsmart.cli.sub import sub
from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.gaussian.job import GaussianJob
from chemsmart.jobs.gaussian.opt import GaussianOptJob
from chemsmart.jobs.gaussian.redox import (
    GaussianRedoxJob,
    GaussianRedoxJobSettings,
)
from chemsmart.jobs.gaussian.runner import GaussianJobRunner
from chemsmart.jobs.gaussian.settings import GaussianJobSettings
from chemsmart.jobs.gaussian.singlepoint import GaussianSinglePointJob
from chemsmart.jobs.redox import reduced_charge_and_multiplicity
from chemsmart.utils.constants import FARADAY, energy_conversion


def _write_h2_xyz(path, comment="h2"):
    Path(path).write_text(f"2\n{comment}\nH 0.0 0.0 0.0\nH 0.0 0.0 0.74\n")
    return str(path)


def _h2_molecule(charge=1, multiplicity=2):
    return Molecule(
        symbols=["H", "H"],
        positions=[[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]],
        charge=charge,
        multiplicity=multiplicity,
    )


@pytest.fixture
def isolated_redox_registry(monkeypatch):
    from chemsmart.analysis import redox as redox_analysis

    monkeypatch.setattr(
        redox_analysis, "_REGISTRY", dict(redox_analysis._REGISTRY)
    )
    return redox_analysis._REGISTRY


class TestRedoxReferenceRegistry:
    def test_builtin_fc_fc_plus(self):
        reference = get_redox_reference("fc_fc+")
        assert reference.name == "fc_fc+"
        assert reference.E_ref_V == 0.0
        assert reference.n_electrons == 1
        assert reference.scale == "Fc/Fc+"
        assert reference.couple_label == "Fc/Fc+"
        assert reference.ox_file is None
        assert reference.red_file is None
        names = [item.name for item in list_redox_references()]
        assert "fc_fc+" in names

    def test_unknown_reference_errors(self):
        with pytest.raises(ValueError, match="Unknown redox reference"):
            get_redox_reference("not_a_couple")

    def test_register_rejects_duplicates(self):
        with pytest.raises(ValueError, match="already registered"):
            register_redox_reference(get_redox_reference("fc_fc+"))

    def test_register_extension(self, isolated_redox_registry, tmp_path):
        ox = _write_h2_xyz(tmp_path / "ref_ox.xyz")
        red = _write_h2_xyz(tmp_path / "ref_red.xyz")
        custom = RedoxReference(
            name="custom_she",
            E_ref_V=0.40,
            n_electrons=1,
            scale="SHE",
            couple_label="Custom/Custom+",
            ox_file=ox,
            red_file=red,
            ox_charge=1,
            ox_multiplicity=2,
            red_charge=0,
            red_multiplicity=1,
        )
        register_redox_reference(custom)
        loaded = get_redox_reference("custom_she")
        assert loaded is custom
        assert loaded.E_ref_V == 0.40
        assert loaded.scale == "SHE"

    def test_invalid_n_electrons(self):
        with pytest.raises(ValueError, match="positive integer"):
            RedoxReference(
                name="bad",
                E_ref_V=0.0,
                n_electrons=0,
                scale="SHE",
                couple_label="X/X+",
            )

    def test_cli_reexports_registry_api(self):
        from chemsmart.analysis.redox import (
            RedoxReference as AnalysisRedoxReference,
        )
        from chemsmart.analysis.redox import (
            get_redox_reference as analysis_get,
        )
        from chemsmart.cli.redox import RedoxReference as CliRedoxReference
        from chemsmart.cli.redox import get_redox_reference as cli_get

        assert CliRedoxReference is AnalysisRedoxReference
        assert cli_get is analysis_get


class TestReducedChargeAndMultiplicity:
    def test_odd_n_singlet_becomes_doublet(self):
        assert reduced_charge_and_multiplicity(1, 1, 1) == (0, 2)

    def test_odd_n_doublet_becomes_singlet(self):
        assert reduced_charge_and_multiplicity(1, 2, 1) == (0, 1)

    def test_even_n_keeps_multiplicity(self):
        assert reduced_charge_and_multiplicity(2, 1, 2) == (0, 1)


class TestComputeRedoxPotential:
    def _install_fake_energies(self, monkeypatch, gas, solv):
        def fake_gas(filepath, **kwargs):
            return gas[Path(filepath).name]

        def fake_solv(filepath):
            return solv[Path(filepath).name]

        monkeypatch.setattr(
            "chemsmart.analysis.redox.gas_phase_data", fake_gas
        )
        monkeypatch.setattr(
            "chemsmart.analysis.redox.solvent_scf_energy", fake_solv
        )

    def _result_for_n(self, monkeypatch, tmp_path, n_electrons, e_ref):
        files = {
            name: _write_h2_xyz(tmp_path / name)
            for name in (
                "ox_gas.log",
                "red_gas.log",
                "ref_ox_gas.log",
                "ref_red_gas.log",
                "ox_solv.log",
                "red_solv.log",
                "ref_ox_solv.log",
                "ref_red_solv.log",
            )
        }
        # E_gas, G_corr
        gas = {
            "ox_gas.log": (1.00, 0.01),
            "red_gas.log": (1.10, 0.02),
            "ref_ox_gas.log": (2.00, 0.03),
            "ref_red_gas.log": (2.20, 0.04),
        }
        solv = {
            "ox_solv.log": 0.90,
            "red_solv.log": 1.00,
            "ref_ox_solv.log": 1.80,
            "ref_red_solv.log": 2.00,
        }
        self._install_fake_energies(monkeypatch, gas, solv)
        name = f"test_n{n_electrons}"
        register_redox_reference(
            RedoxReference(
                name=name,
                E_ref_V=e_ref,
                n_electrons=n_electrons,
                scale="Fc/Fc+",
                couple_label="Test/Test+",
            )
        )
        return compute_redox_potential(
            ox_gas_file=files["ox_gas.log"],
            red_gas_file=files["red_gas.log"],
            ref_ox_gas_file=files["ref_ox_gas.log"],
            ref_red_gas_file=files["ref_red_gas.log"],
            ox_solv_file=files["ox_solv.log"],
            red_solv_file=files["red_solv.log"],
            ref_ox_solv_file=files["ref_ox_solv.log"],
            ref_red_solv_file=files["ref_red_solv.log"],
            reference=name,
        )

    def test_n1_exchange_arithmetic(
        self, isolated_redox_registry, tmp_path, monkeypatch
    ):
        result = self._result_for_n(monkeypatch, tmp_path, 1, e_ref=0.0)
        # G_soln = E_solv + G_corr
        g_ox = 0.90 + 0.01
        g_red = 1.00 + 0.02
        g_ref_ox = 1.80 + 0.03
        g_ref_red = 2.00 + 0.04
        delta_g_au = g_red + g_ref_ox - g_ox - g_ref_red
        assert result["delta_G_exchange_au"] == pytest.approx(delta_g_au)
        assert result["delta_G_exchange_au"] == pytest.approx(-0.10)
        delta_g_j = energy_conversion("hartree", "j/mol", delta_g_au)
        expected_e = 0.0 - delta_g_j / (1 * FARADAY)
        assert result["E_target_V"] == pytest.approx(expected_e)
        assert result["E_ref_V"] == 0.0
        assert result["n_electrons"] == 1
        assert result["scale"] == "Fc/Fc+"
        assert result["couple_label"] == "Test/Test+"
        assert result["reference_name"] == "test_n1"
        assert result["G_soln_ox_au"] == pytest.approx(g_ox)
        assert result["G_soln_red_au"] == pytest.approx(g_red)

    def test_n2_halves_potential_contribution(
        self, isolated_redox_registry, tmp_path, monkeypatch
    ):
        result = self._result_for_n(monkeypatch, tmp_path, 2, e_ref=0.25)
        delta_g_j = energy_conversion("hartree", "j/mol", -0.10)
        expected_e = 0.25 - delta_g_j / (2 * FARADAY)
        assert result["n_electrons"] == 2
        assert result["E_ref_V"] == 0.25
        assert result["E_target_V"] == pytest.approx(expected_e)
        n1_shift = delta_g_j / FARADAY
        n2_shift = delta_g_j / (2 * FARADAY)
        assert n2_shift == pytest.approx(n1_shift / 2)

    def test_format_summary_includes_units(
        self, isolated_redox_registry, tmp_path, monkeypatch
    ):
        result = self._result_for_n(monkeypatch, tmp_path, 1, e_ref=0.0)
        text = format_redox_summary(result)
        assert "ΔG_exchange" in text
        assert "kcal/mol" in text
        assert "eV" in text
        assert "J/mol" in text
        assert "E_target" in text
        assert "Fc/Fc+" in text

    def test_missing_files_error(self):
        with pytest.raises(ValueError, match="Missing required files"):
            compute_redox_potential(
                ox_gas_file="ox.log",
                red_gas_file="red.log",
                ref_ox_gas_file=None,
                ref_red_gas_file="ref_red.log",
                ox_solv_file="oxs.log",
                red_solv_file="reds.log",
                ref_ox_solv_file="refoxs.log",
                ref_red_solv_file="refreds.log",
            )

    def test_e_ref_without_registry(self, tmp_path, monkeypatch):
        files = {
            name: _write_h2_xyz(tmp_path / name)
            for name in (
                "ox_gas.log",
                "red_gas.log",
                "ref_ox_gas.log",
                "ref_red_gas.log",
                "ox_solv.log",
                "red_solv.log",
                "ref_ox_solv.log",
                "ref_red_solv.log",
            )
        }
        self._install_fake_energies(
            monkeypatch,
            {
                "ox_gas.log": (1.00, 0.01),
                "red_gas.log": (1.10, 0.02),
                "ref_ox_gas.log": (2.00, 0.03),
                "ref_red_gas.log": (2.20, 0.04),
            },
            {
                "ox_solv.log": 0.90,
                "red_solv.log": 1.00,
                "ref_ox_solv.log": 1.80,
                "ref_red_solv.log": 2.00,
            },
        )
        result = compute_redox_potential(
            ox_gas_file=files["ox_gas.log"],
            red_gas_file=files["red_gas.log"],
            ref_ox_gas_file=files["ref_ox_gas.log"],
            ref_red_gas_file=files["ref_red_gas.log"],
            ox_solv_file=files["ox_solv.log"],
            red_solv_file=files["red_solv.log"],
            ref_ox_solv_file=files["ref_ox_solv.log"],
            ref_red_solv_file=files["ref_red_solv.log"],
            e_ref=0.2,
        )
        assert result["E_ref_V"] == pytest.approx(0.2)
        assert result["n_electrons"] == 1
        assert result["reference_name"] is None
        delta_g_j = energy_conversion("hartree", "j/mol", -0.10)
        assert result["E_target_V"] == pytest.approx(0.2 - delta_g_j / FARADAY)
        text = format_redox_summary(result)
        assert "E_ref = 0.2000 V" in text
        assert "Fc/Fc+" not in text

    def test_e_ref_overrides_registry_potential(
        self, isolated_redox_registry, tmp_path, monkeypatch
    ):
        files = {
            name: _write_h2_xyz(tmp_path / name)
            for name in (
                "ox_gas.log",
                "red_gas.log",
                "ref_ox_gas.log",
                "ref_red_gas.log",
                "ox_solv.log",
                "red_solv.log",
                "ref_ox_solv.log",
                "ref_red_solv.log",
            )
        }
        self._install_fake_energies(
            monkeypatch,
            {Path(p).name: (0.0, 0.0) for p in files.values()},
            {Path(p).name: 0.0 for p in files.values()},
        )
        result = compute_redox_potential(
            ox_gas_file=files["ox_gas.log"],
            red_gas_file=files["red_gas.log"],
            ref_ox_gas_file=files["ref_ox_gas.log"],
            ref_red_gas_file=files["ref_red_gas.log"],
            ox_solv_file=files["ox_solv.log"],
            red_solv_file=files["red_solv.log"],
            ref_ox_solv_file=files["ref_ox_solv.log"],
            ref_red_solv_file=files["ref_red_solv.log"],
            reference="fc_fc+",
            e_ref=0.2,
        )
        assert result["E_ref_V"] == pytest.approx(0.2)
        assert result["reference_name"] == "fc_fc+"

    def test_n_mismatch_errors(
        self, isolated_redox_registry, tmp_path, monkeypatch
    ):
        files = {
            name: _write_h2_xyz(tmp_path / name)
            for name in (
                "ox_gas.log",
                "red_gas.log",
                "ref_ox_gas.log",
                "ref_red_gas.log",
                "ox_solv.log",
                "red_solv.log",
                "ref_ox_solv.log",
                "ref_red_solv.log",
            )
        }
        self._install_fake_energies(
            monkeypatch,
            {Path(p).name: (0.0, 0.0) for p in files.values()},
            {Path(p).name: 0.0 for p in files.values()},
        )
        with pytest.raises(ValueError, match="must match the reference"):
            compute_redox_potential(
                ox_gas_file=files["ox_gas.log"],
                red_gas_file=files["red_gas.log"],
                ref_ox_gas_file=files["ref_ox_gas.log"],
                ref_red_gas_file=files["ref_red_gas.log"],
                ox_solv_file=files["ox_solv.log"],
                red_solv_file=files["red_solv.log"],
                ref_ox_solv_file=files["ref_ox_solv.log"],
                ref_red_solv_file=files["ref_red_solv.log"],
                reference="fc_fc+",
                n_electrons=2,
            )

    def test_output_class_wrappers(
        self, isolated_redox_registry, tmp_path, monkeypatch
    ):
        from chemsmart.io.gaussian.output import Gaussian16RedoxOutput
        from chemsmart.io.orca.output import ORCARedoxOutput

        result = self._result_for_n(monkeypatch, tmp_path, 1, e_ref=0.0)
        text = Gaussian16RedoxOutput.format_redox_summary(result)
        assert "E_target" in text
        assert "E_target" in ORCARedoxOutput.format_redox_summary(result)


def _redox_settings(tmp_path, **kwargs):
    ref_ox = kwargs.pop("ref_ox_file", None)
    if ref_ox is None:
        ref_ox = _write_h2_xyz(tmp_path / "ref_ox.xyz")
    include_ref_red = "ref_red_file" in kwargs
    ref_red = kwargs.pop("ref_red_file", None)
    if not include_ref_red:
        ref_red = _write_h2_xyz(tmp_path / "ref_red.xyz")
    defaults = dict(
        functional="B3LYP",
        basis="6-31G*",
        charge=1,
        multiplicity=2,
        solvent_model="SMD",
        solvent_id="acetonitrile",
        ref_ox_file=ref_ox,
        ref_ox_charge=1,
        ref_ox_multiplicity=2,
        ref_red_charge=0,
        ref_red_multiplicity=1,
    )
    if ref_red is not None:
        defaults["ref_red_file"] = ref_red
    defaults.update(kwargs)
    return GaussianRedoxJobSettings(**defaults)


class TestGaussianRedoxJob:
    def test_creates_opt_and_ref_jobs(
        self, tmp_path, gaussian_jobrunner_no_scratch
    ):
        mol = _h2_molecule()
        settings = _redox_settings(tmp_path)
        job = GaussianRedoxJob(
            molecule=mol,
            settings=settings,
            label="mol_redox",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert isinstance(job, GaussianJob)
        assert job.TYPE == "g16redox"
        assert [phase.name for phase in job.phases] == [
            "Opt",
            "Ref Opt",
            "SP",
            "Ref SP",
        ]
        assert job.ox_job.label == "mol_redox_ox_opt"
        assert job.red_job.label == "mol_redox_red_opt"
        assert job.ref_ox_job.label == "mol_redox_RefOx_opt"
        assert job.ref_red_job.label == "mol_redox_RefRed_opt"
        assert isinstance(job.ox_job, GaussianOptJob)
        assert job.ox_job.settings.charge == 1
        assert job.ox_job.settings.multiplicity == 2
        assert job.ox_job.settings.freq is True
        assert job.ox_job.settings.solvent_model is None
        assert job.red_job.settings.charge == 0
        assert job.red_job.settings.multiplicity == 1
        assert job.ref_ox_job.settings.charge == 1
        assert job.ref_red_job.settings.charge == 0

    def test_sp_jobs_use_solvent_and_labels(
        self, tmp_path, gaussian_jobrunner_no_scratch
    ):
        job = GaussianRedoxJob(
            molecule=_h2_molecule(),
            settings=_redox_settings(tmp_path),
            label="mol_redox",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        sp_jobs = job._sp_jobs_for_phase()
        ref_sp_jobs = job._ref_sp_jobs_for_phase()
        assert [child.label for child in sp_jobs] == [
            "mol_redox_ox_sp",
            "mol_redox_red_sp",
        ]
        assert [child.label for child in ref_sp_jobs] == [
            "mol_redox_RefOx_sp",
            "mol_redox_RefRed_sp",
        ]
        assert isinstance(job.ox_sp_job, GaussianSinglePointJob)
        assert job.ox_sp_job.settings.jobtype == "sp"
        assert job.ox_sp_job.settings.freq is False
        assert job.ox_sp_job.settings.solvent_model == "SMD"
        assert job.ox_sp_job.settings.solvent_id == "acetonitrile"

    def test_red_file_uses_separate_geometry(
        self, tmp_path, gaussian_jobrunner_no_scratch
    ):
        red_xyz = tmp_path / "red.xyz"
        red_xyz.write_text("3\nred\nH 0 0 0\nH 0 0 0.74\nH 0 0 1.5\n")
        settings = _redox_settings(tmp_path, red_file=str(red_xyz))
        job = GaussianRedoxJob(
            molecule=_h2_molecule(),
            settings=settings,
            label="mol_redox",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert len(job.red_job.molecule) == 3
        assert len(job.ox_job.molecule) == 2

    def test_ref_ox_only_derives_reduced_reference(
        self, tmp_path, gaussian_jobrunner_no_scratch
    ):
        ref_ox_xyz = tmp_path / "ref_ox.xyz"
        ref_ox_xyz.write_text("3\nref\nH 0 0 0\nH 0 0 0.74\nH 0 0 1.5\n")
        settings = _redox_settings(
            tmp_path, ref_ox_file=str(ref_ox_xyz), ref_red_file=None
        )
        job = GaussianRedoxJob(
            molecule=_h2_molecule(),
            settings=settings,
            label="mol_redox",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert len(job.ref_ox_job.molecule) == 3
        assert len(job.ref_red_job.molecule) == 3
        assert job.ref_ox_job.settings.charge == 1
        assert job.ref_red_job.settings.charge == 0
        assert job.ref_ox_job.settings.multiplicity == 2
        assert job.ref_red_job.settings.multiplicity == 1

    def test_ref_red_file_uses_separate_geometry(
        self, tmp_path, gaussian_jobrunner_no_scratch
    ):
        ref_red_xyz = tmp_path / "ref_red.xyz"
        ref_red_xyz.write_text("2\nref\nH 0 0 0\nH 0 0 0.74\n")
        settings = _redox_settings(tmp_path, ref_red_file=str(ref_red_xyz))
        job = GaussianRedoxJob(
            molecule=_h2_molecule(),
            settings=settings,
            label="mol_redox",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert len(job.ref_ox_job.molecule) == 2
        assert len(job.ref_red_job.molecule) == 2

    def test_missing_reference_geometries_error(
        self, tmp_path, gaussian_jobrunner_no_scratch
    ):
        settings = GaussianRedoxJobSettings(
            functional="B3LYP",
            basis="6-31G*",
            charge=1,
            multiplicity=2,
        )
        with pytest.raises(ValueError, match="oxidized reference"):
            GaussianRedoxJob(
                molecule=_h2_molecule(),
                settings=settings,
                label="mol_redox",
                jobrunner=gaussian_jobrunner_no_scratch,
            )

    def test_invalid_settings_type(
        self, tmp_path, gaussian_jobrunner_no_scratch
    ):
        settings = GaussianJobSettings(functional="B3LYP", basis="6-31G*")
        with pytest.raises(
            ValueError, match="must be instance of GaussianRedoxJobSettings"
        ):
            GaussianRedoxJob(
                molecule=_h2_molecule(),
                settings=settings,
                label="mol_redox",
                jobrunner=gaussian_jobrunner_no_scratch,
            )

    def test_custom_reference_without_editing_job_class(
        self,
        isolated_redox_registry,
        tmp_path,
        gaussian_jobrunner_no_scratch,
    ):
        ox = _write_h2_xyz(tmp_path / "custom_ox.xyz")
        red = _write_h2_xyz(tmp_path / "custom_red.xyz")
        register_redox_reference(
            RedoxReference(
                name="n2_couple",
                E_ref_V=0.55,
                n_electrons=2,
                scale="SHE",
                couple_label="N2/N2+",
                ox_file=ox,
                red_file=red,
                ox_charge=2,
                ox_multiplicity=1,
                red_charge=0,
                red_multiplicity=1,
            )
        )
        settings = GaussianRedoxJobSettings(
            functional="B3LYP",
            basis="6-31G*",
            charge=2,
            multiplicity=1,
            reference="n2_couple",
        )
        job = GaussianRedoxJob(
            molecule=_h2_molecule(charge=2, multiplicity=1),
            settings=settings,
            label="mol_redox",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert job.settings.reference.name == "n2_couple"
        assert job.settings.n_electrons == 2
        assert job.settings.reference.E_ref_V == 0.55
        assert job.ox_job.settings.charge == 2
        assert job.red_job.settings.charge == 0
        assert job.ref_ox_job.settings.charge == 2
        assert job.ref_red_job.settings.charge == 0

    def test_registry_ox_file_only_derives_reduced_reference(
        self,
        isolated_redox_registry,
        tmp_path,
        gaussian_jobrunner_no_scratch,
    ):
        ox = _write_h2_xyz(tmp_path / "registry_ox.xyz")
        register_redox_reference(
            RedoxReference(
                name="ox_only",
                E_ref_V=0.25,
                n_electrons=1,
                scale="SHE",
                couple_label="OxOnly/OxOnly+",
                ox_file=ox,
                ox_charge=1,
                ox_multiplicity=2,
            )
        )
        settings = GaussianRedoxJobSettings(
            functional="B3LYP",
            basis="6-31G*",
            charge=1,
            multiplicity=2,
            reference="ox_only",
        )
        job = GaussianRedoxJob(
            molecule=_h2_molecule(),
            settings=settings,
            label="mol_redox",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert len(job.ref_ox_job.molecule) == 2
        assert len(job.ref_red_job.molecule) == 2
        assert job.ref_ox_job.settings.charge == 1
        assert job.ref_red_job.settings.charge == 0

    def test_jobtypes_includes_redox(self):
        assert "g16redox" in GaussianJobRunner.JOBTYPES


class TestORCARedoxJob:
    def test_creates_opt_jobs(self, tmp_path, orca_jobrunner_no_scratch):
        from chemsmart.jobs.orca.job import ORCAJob
        from chemsmart.jobs.orca.opt import ORCAOptJob
        from chemsmart.jobs.orca.redox import (
            ORCARedoxJob,
            ORCARedoxJobSettings,
        )
        from chemsmart.jobs.orca.runner import ORCAJobRunner

        ref_ox = _write_h2_xyz(tmp_path / "ref_ox.xyz")
        ref_red = _write_h2_xyz(tmp_path / "ref_red.xyz")
        settings = ORCARedoxJobSettings(
            functional="B3LYP",
            basis="def2-SVP",
            charge=1,
            multiplicity=2,
            solvent_model="CPCM",
            solvent_id="acetonitrile",
            ref_ox_file=ref_ox,
            ref_red_file=ref_red,
            ref_ox_charge=1,
            ref_ox_multiplicity=2,
            ref_red_charge=0,
            ref_red_multiplicity=1,
        )
        job = ORCARedoxJob(
            molecule=_h2_molecule(),
            settings=settings,
            label="mol_redox",
            jobrunner=orca_jobrunner_no_scratch,
        )
        assert isinstance(job, ORCAJob)
        assert job.TYPE == "orcaredox"
        assert [phase.name for phase in job.phases] == [
            "Opt",
            "Ref Opt",
            "SP",
            "Ref SP",
        ]
        assert isinstance(job.ox_job, ORCAOptJob)
        assert job.ox_job.label == "mol_redox_ox_opt"
        assert job.red_job.label == "mol_redox_red_opt"
        assert job.ref_ox_job.label == "mol_redox_RefOx_opt"
        assert job.ref_red_job.label == "mol_redox_RefRed_opt"
        assert "orcaredox" in ORCAJobRunner.JOBTYPES


def _analyze_file_args(files):
    return [
        "--ox-gas",
        files["ox_gas.log"],
        "--red-gas",
        files["red_gas.log"],
        "--ref-ox-gas",
        files["ref_ox_gas.log"],
        "--ref-red-gas",
        files["ref_red_gas.log"],
        "--ox-solv",
        files["ox_solv.log"],
        "--red-solv",
        files["red_solv.log"],
        "--ref-ox-solv",
        files["ref_ox_solv.log"],
        "--ref-red-solv",
        files["ref_red_solv.log"],
    ]


def _analyze_cli_args(files, *extra):
    return ["analyze", "--e-ref", "0.0", *_analyze_file_args(files), *extra]


def _write_chem_smart_redox_outputs(tmp_path, basename="mol_redox"):
    """Create eight redox outputs using CHEMSMART job labels."""
    names = (
        f"{basename}_ox_opt.log",
        f"{basename}_red_opt.log",
        f"{basename}_RefOx_opt.log",
        f"{basename}_RefRed_opt.log",
        f"{basename}_ox_sp.log",
        f"{basename}_red_sp.log",
        f"{basename}_RefOx_sp.log",
        f"{basename}_RefRed_sp.log",
    )
    files = {}
    for name in names:
        path = tmp_path / name
        path.write_text("Gaussian, Inc.\n")
        files[name] = str(path)
    return files


class TestRedoxCLI:
    def test_run_help_lists_redox(self):
        result = CliRunner().invoke(run, ["--help"])
        assert result.exit_code == 0, result.output
        assert "\n  redox" in result.output
        result = CliRunner().invoke(run, ["chain", "--help"])
        assert result.exit_code == 0, result.output
        assert "\n  redox" in result.output

    def test_sub_help_does_not_list_analysis_redox(self):
        result = CliRunner().invoke(sub, ["--help"])
        assert result.exit_code == 0, result.output
        assert "\n  redox" not in result.output

    def test_gaussian_and_orca_help_list_redox_submit(self):
        runner = CliRunner()
        gaussian_help = runner.invoke(run, ["gaussian", "--help"])
        assert gaussian_help.exit_code == 0, gaussian_help.output
        assert "\n  redox" in gaussian_help.output
        orca_help = runner.invoke(run, ["orca", "--help"])
        assert orca_help.exit_code == 0, orca_help.output
        assert "\n  redox" in orca_help.output

    def test_gaussian_redox_help_is_submission_only(
        self, single_molecule_xyz_file
    ):
        result = CliRunner().invoke(
            run,
            [
                "--fake",
                "gaussian",
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "redox",
                "--help",
            ],
        )
        assert result.exit_code == 0, result.output
        assert "--ref-ox" in result.output
        assert "-ro" in result.output
        assert "--ref-red" in result.output
        assert "-rr" in result.output
        assert "--red" in result.output
        assert "-rd" in result.output
        assert "\n  analyze" not in result.output

    def test_orca_redox_help_is_submission_only(
        self, single_molecule_xyz_file
    ):
        result = CliRunner().invoke(
            run,
            [
                "--fake",
                "orca",
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "redox",
                "--help",
            ],
        )
        assert result.exit_code == 0, result.output
        assert "--ref-ox" in result.output
        assert "\n  analyze" not in result.output

    def test_run_redox_help_keeps_analyze(self):
        result = CliRunner().invoke(run, ["redox", "--help"])
        assert result.exit_code == 0, result.output
        assert "\n  analyze" in result.output
        assert "-T, --temperature" not in result.output
        assert "-csg" not in result.output
        assert "-r, --reference" not in result.output
        assert "-er, --e-ref" not in result.output

    def test_run_redox_analyze_help_lists_short_flags(self):
        result = CliRunner().invoke(run, ["redox", "analyze", "--help"])
        assert result.exit_code == 0, result.output
        assert "-oxg" in result.output
        assert "-rdg" in result.output
        assert "-rog" in result.output
        assert "-rrg" in result.output
        assert "-oxs" in result.output
        assert "-rds" in result.output
        assert "-ros" in result.output
        assert "-rrs" in result.output
        assert "-T, --temperature" in result.output
        assert "-c, --concentration" in result.output
        assert "-csg" in result.output
        assert "-ch" in result.output
        assert "-o, --output" in result.output
        assert "--e-ref" in result.output
        assert "-er" in result.output

    def test_gaussian_redox_builds_job(
        self, tmp_path, single_molecule_xyz_file
    ):
        from chemsmart.cli.gaussian.gaussian import gaussian
        from chemsmart.jobs.gaussian.redox import GaussianRedoxJob

        ref_ox = _write_h2_xyz(tmp_path / "ref_ox.xyz")
        ref_red = _write_h2_xyz(tmp_path / "ref_red.xyz")
        result = CliRunner().invoke(
            gaussian,
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "1",
                "-m",
                "2",
                "redox",
                "--ref-ox",
                ref_ox,
                "--ref-red",
                ref_red,
            ],
            obj={"jobrunner": MagicMock()},
            catch_exceptions=False,
            standalone_mode=False,
        )
        assert result.exit_code == 0, result.output
        job = result.return_value
        assert isinstance(job, GaussianRedoxJob)
        assert job.ox_job.settings.charge == 1
        assert job.red_job.settings.charge == 0
        assert job.ox_job.settings.solvent_model is None
        assert job.settings.solvent_model == "smd"
        assert job.settings.reference.name == "fc_fc+"

    def test_gaussian_redox_accepts_short_reference_flags(
        self, tmp_path, single_molecule_xyz_file
    ):
        from chemsmart.cli.gaussian.gaussian import gaussian
        from chemsmart.jobs.gaussian.redox import GaussianRedoxJob

        ref_ox = _write_h2_xyz(tmp_path / "ref_ox.xyz")
        ref_red = _write_h2_xyz(tmp_path / "ref_red.xyz")
        result = CliRunner().invoke(
            gaussian,
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "1",
                "-m",
                "2",
                "redox",
                "-ro",
                ref_ox,
                "-rr",
                ref_red,
            ],
            obj={"jobrunner": MagicMock()},
            catch_exceptions=False,
            standalone_mode=False,
        )
        assert result.exit_code == 0, result.output
        job = result.return_value
        assert isinstance(job, GaussianRedoxJob)
        assert job.settings.ref_ox_file == ref_ox
        assert job.settings.ref_red_file == ref_red

    def test_sub_reconstructs_ref_ox_not_dest_name(
        self, tmp_path, single_molecule_xyz_file, monkeypatch
    ):
        from chemsmart.settings.server import Server

        ref_ox = _write_h2_xyz(tmp_path / "ref_ox.xyz")
        ref_red = _write_h2_xyz(tmp_path / "ref_red.xyz")
        captured = {}
        fake_server = Server(name="dummy")
        fake_server.submit = (
            lambda job, test=False, cli_args=None, **kw: captured.update(
                cli_args=cli_args
            )
        )
        monkeypatch.setattr(
            "chemsmart.settings.server.Server.from_servername",
            lambda _name: fake_server,
        )

        result = CliRunner().invoke(
            sub,
            [
                "--test",
                "--server",
                "dummy",
                "gaussian",
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "1",
                "-m",
                "2",
                "redox",
                "--ref-ox",
                ref_ox,
                "--ref-red",
                ref_red,
            ],
        )
        assert result.exit_code == 0, result.output
        cli_args = captured["cli_args"]
        assert "--ref-ox" in cli_args
        assert "--ref-red" in cli_args
        assert "--ref-ox-file" not in cli_args
        assert "--ref-red-file" not in cli_args
        assert cli_args[cli_args.index("--ref-ox") + 1] == ref_ox
        assert cli_args[cli_args.index("--ref-red") + 1] == ref_red

    def test_orca_redox_builds_job(self, tmp_path, single_molecule_xyz_file):
        from chemsmart.cli.orca.orca import orca
        from chemsmart.jobs.orca.redox import ORCARedoxJob

        ref_ox = _write_h2_xyz(tmp_path / "ref_ox.xyz")
        ref_red = _write_h2_xyz(tmp_path / "ref_red.xyz")
        result = CliRunner().invoke(
            orca,
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "1",
                "-m",
                "2",
                "redox",
                "--ref-ox",
                ref_ox,
                "--ref-red",
                ref_red,
            ],
            obj={"jobrunner": MagicMock()},
            catch_exceptions=False,
            standalone_mode=False,
        )
        assert result.exit_code == 0, result.output
        assert isinstance(result.return_value, ORCARedoxJob)

    def test_gaussian_redox_requires_reference_geometries(
        self, single_molecule_xyz_file
    ):
        from chemsmart.cli.gaussian.gaussian import gaussian

        result = CliRunner().invoke(
            gaussian,
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "1",
                "-m",
                "2",
                "redox",
            ],
            obj={"jobrunner": MagicMock()},
        )
        assert result.exit_code == 2, result.output
        assert "oxidized reference geometry" in result.output

    def test_run_redox_analyze_requires_e_ref(
        self, isolated_redox_registry, tmp_path, monkeypatch
    ):
        files = {
            name: _write_h2_xyz(tmp_path / name)
            for name in (
                "ox_gas.log",
                "red_gas.log",
                "ref_ox_gas.log",
                "ref_red_gas.log",
                "ox_solv.log",
                "red_solv.log",
                "ref_ox_solv.log",
                "ref_red_solv.log",
            )
        }
        monkeypatch.setattr(
            "chemsmart.analysis.redox.gas_phase_data",
            lambda filepath, **kwargs: (0.0, 0.0),
        )
        monkeypatch.setattr(
            "chemsmart.analysis.redox.solvent_scf_energy",
            lambda filepath: 0.0,
        )
        result = CliRunner().invoke(
            run, ["redox", "analyze", *_analyze_file_args(files)]
        )
        assert result.exit_code == 2, result.output
        assert "--e-ref" in result.output
        assert "E_target" not in result.output

    def test_chain_redox_analyze_requires_e_ref(
        self, isolated_redox_registry, tmp_path, monkeypatch
    ):
        from chemsmart.cli.chain.chain import chain

        files = {
            name: _write_h2_xyz(tmp_path / name)
            for name in (
                "ox_gas.log",
                "red_gas.log",
                "ref_ox_gas.log",
                "ref_red_gas.log",
                "ox_solv.log",
                "red_solv.log",
                "ref_ox_solv.log",
                "ref_red_solv.log",
            )
        }
        monkeypatch.setattr(
            "chemsmart.analysis.redox.gas_phase_data",
            lambda filepath, **kwargs: (0.0, 0.0),
        )
        monkeypatch.setattr(
            "chemsmart.analysis.redox.solvent_scf_energy",
            lambda filepath: 0.0,
        )
        result = CliRunner().invoke(
            chain,
            ["redox", "analyze", *_analyze_file_args(files)],
            obj={"jobrunner": MagicMock()},
        )
        assert result.exit_code == 2, result.output
        assert "--e-ref" in result.output
        assert "E_target" not in result.output

    def test_run_redox_analyze_uses_e_ref(
        self, isolated_redox_registry, tmp_path, monkeypatch
    ):
        files = {
            name: _write_h2_xyz(tmp_path / name)
            for name in (
                "ox_gas.log",
                "red_gas.log",
                "ref_ox_gas.log",
                "ref_red_gas.log",
                "ox_solv.log",
                "red_solv.log",
                "ref_ox_solv.log",
                "ref_red_solv.log",
            )
        }
        monkeypatch.setattr(
            "chemsmart.analysis.redox.gas_phase_data",
            lambda filepath, **kwargs: (0.0, 0.0),
        )
        monkeypatch.setattr(
            "chemsmart.analysis.redox.solvent_scf_energy",
            lambda filepath: 0.0,
        )
        result = CliRunner().invoke(
            run,
            ["redox", "analyze", "--e-ref", "0.2", *_analyze_file_args(files)],
        )
        assert result.exit_code == 0, result.output
        assert "E_ref = 0.2000 V" in result.output
        assert "Fc/Fc+" not in result.output

    def test_run_redox_analyze_arithmetic(
        self, isolated_redox_registry, tmp_path, monkeypatch
    ):
        files = {
            name: _write_h2_xyz(tmp_path / name)
            for name in (
                "ox_gas.log",
                "red_gas.log",
                "ref_ox_gas.log",
                "ref_red_gas.log",
                "ox_solv.log",
                "red_solv.log",
                "ref_ox_solv.log",
                "ref_red_solv.log",
            )
        }
        gas = {
            "ox_gas.log": (1.00, 0.01),
            "red_gas.log": (1.10, 0.02),
            "ref_ox_gas.log": (2.00, 0.03),
            "ref_red_gas.log": (2.20, 0.04),
        }
        solv = {
            "ox_solv.log": 0.90,
            "red_solv.log": 1.00,
            "ref_ox_solv.log": 1.80,
            "ref_red_solv.log": 2.00,
        }

        def fake_gas(filepath, **kwargs):
            return gas[Path(filepath).name]

        def fake_solv(filepath):
            return solv[Path(filepath).name]

        monkeypatch.setattr(
            "chemsmart.analysis.redox.gas_phase_data", fake_gas
        )
        monkeypatch.setattr(
            "chemsmart.analysis.redox.solvent_scf_energy", fake_solv
        )
        result = CliRunner().invoke(run, ["redox", *_analyze_cli_args(files)])
        assert result.exit_code == 0, result.output
        assert "E_target" in result.output
        assert "E_ref = 0.0000 V" in result.output
        assert "kcal/mol" in result.output
        assert "J/mol" in result.output

    def test_run_redox_analyze_auto_discover_companion_outputs(
        self, isolated_redox_registry, tmp_path, monkeypatch
    ):
        files = _write_chem_smart_redox_outputs(tmp_path)
        gas = {
            "mol_redox_ox_opt.log": (1.00, 0.01),
            "mol_redox_red_opt.log": (1.10, 0.02),
            "mol_redox_RefOx_opt.log": (2.00, 0.03),
            "mol_redox_RefRed_opt.log": (2.20, 0.04),
        }
        solv = {
            "mol_redox_ox_sp.log": 0.90,
            "mol_redox_red_sp.log": 1.00,
            "mol_redox_RefOx_sp.log": 1.80,
            "mol_redox_RefRed_sp.log": 2.00,
        }

        def fake_gas(filepath, **kwargs):
            return gas[Path(filepath).name]

        def fake_solv(filepath):
            return solv[Path(filepath).name]

        monkeypatch.setattr(
            "chemsmart.analysis.redox.gas_phase_data", fake_gas
        )
        monkeypatch.setattr(
            "chemsmart.analysis.redox.solvent_scf_energy", fake_solv
        )
        result = CliRunner().invoke(
            run,
            [
                "redox",
                "analyze",
                "--e-ref",
                "0.0",
                "--ox-gas",
                files["mol_redox_ox_opt.log"],
                "--ref-ox-gas",
                files["mol_redox_RefOx_opt.log"],
            ],
        )
        assert result.exit_code == 0, result.output
        assert "E_target" in result.output

    def test_run_redox_analyze_requires_ox_and_ref_ox_gas(
        self, isolated_redox_registry, tmp_path, monkeypatch
    ):
        files = _write_chem_smart_redox_outputs(tmp_path)
        monkeypatch.setattr(
            "chemsmart.analysis.redox.gas_phase_data",
            lambda filepath, **kwargs: (0.0, 0.0),
        )
        monkeypatch.setattr(
            "chemsmart.analysis.redox.solvent_scf_energy",
            lambda filepath: 0.0,
        )
        result = CliRunner().invoke(
            run,
            [
                "redox",
                "analyze",
                "--e-ref",
                "0.0",
                "--ref-ox-gas",
                files["mol_redox_RefOx_opt.log"],
            ],
        )
        assert result.exit_code == 2, result.output
        assert "--ox-gas is required" in result.output

    def test_run_redox_analyze_accepts_thermochemistry_after_analyze(
        self, isolated_redox_registry, tmp_path, monkeypatch
    ):
        files = {
            name: _write_h2_xyz(tmp_path / name)
            for name in (
                "ox_gas.log",
                "red_gas.log",
                "ref_ox_gas.log",
                "ref_red_gas.log",
                "ox_solv.log",
                "red_solv.log",
                "ref_ox_solv.log",
                "ref_red_solv.log",
            )
        }
        captured = {}

        def fake_gas(filepath, **kwargs):
            captured.update(kwargs)
            return 0.0, 0.0

        monkeypatch.setattr(
            "chemsmart.analysis.redox.gas_phase_data", fake_gas
        )
        monkeypatch.setattr(
            "chemsmart.analysis.redox.solvent_scf_energy", lambda filepath: 0.0
        )
        result = CliRunner().invoke(
            run,
            [
                "redox",
                *_analyze_cli_args(
                    files,
                    "-T",
                    "333.15",
                    "-c",
                    "1.0",
                    "-csg",
                    "100",
                    "-ch",
                    "100",
                ),
            ],
        )
        assert result.exit_code == 0, result.output
        assert captured["temperature"] == 333.15
        assert captured["concentration"] == 1.0
        assert captured["cutoff_entropy_grimme"] == 100.0
        assert captured["cutoff_enthalpy"] == 100.0

    def test_run_redox_analyze_writes_output_file(
        self, isolated_redox_registry, tmp_path, monkeypatch
    ):
        files = {
            name: _write_h2_xyz(tmp_path / name)
            for name in (
                "ox_gas.log",
                "red_gas.log",
                "ref_ox_gas.log",
                "ref_red_gas.log",
                "ox_solv.log",
                "red_solv.log",
                "ref_ox_solv.log",
                "ref_red_solv.log",
            )
        }
        monkeypatch.setattr(
            "chemsmart.analysis.redox.gas_phase_data",
            lambda filepath, **kwargs: (0.0, 0.0),
        )
        monkeypatch.setattr(
            "chemsmart.analysis.redox.solvent_scf_energy",
            lambda filepath: 0.0,
        )
        output_path = tmp_path / "subdir" / "redox.dat"
        result = CliRunner().invoke(
            run,
            [
                "redox",
                *_analyze_cli_args(files, "-o", str(output_path)),
            ],
        )
        assert result.exit_code == 0, result.output
        assert "E_target" not in result.output
        text = output_path.read_text()
        assert "E_target" in text
        assert "E_ref = 0.0000 V" in text

    def test_chain_redox_analyze_arithmetic(
        self, isolated_redox_registry, tmp_path, monkeypatch
    ):
        from chemsmart.cli.chain.chain import chain

        files = {
            name: _write_h2_xyz(tmp_path / name)
            for name in (
                "ox_gas.log",
                "red_gas.log",
                "ref_ox_gas.log",
                "ref_red_gas.log",
                "ox_solv.log",
                "red_solv.log",
                "ref_ox_solv.log",
                "ref_red_solv.log",
            )
        }
        monkeypatch.setattr(
            "chemsmart.analysis.redox.gas_phase_data",
            lambda filepath, **kwargs: (0.0, 0.0),
        )
        monkeypatch.setattr(
            "chemsmart.analysis.redox.solvent_scf_energy",
            lambda filepath: 0.0,
        )
        result = CliRunner().invoke(
            chain,
            ["redox", *_analyze_cli_args(files)],
            obj={"jobrunner": MagicMock()},
        )
        assert result.exit_code == 0, result.output
        assert "E_target" in result.output
        assert "-p/--project" not in result.output

    def test_run_redox_analyze_n2_arithmetic(
        self, isolated_redox_registry, tmp_path, monkeypatch
    ):
        files = {
            name: _write_h2_xyz(tmp_path / name)
            for name in (
                "ox_gas.log",
                "red_gas.log",
                "ref_ox_gas.log",
                "ref_red_gas.log",
                "ox_solv.log",
                "red_solv.log",
                "ref_ox_solv.log",
                "ref_red_solv.log",
            )
        }
        monkeypatch.setattr(
            "chemsmart.analysis.redox.gas_phase_data",
            lambda filepath, **kwargs: {
                "ox_gas.log": (1.00, 0.01),
                "red_gas.log": (1.10, 0.02),
                "ref_ox_gas.log": (2.00, 0.03),
                "ref_red_gas.log": (2.20, 0.04),
            }[Path(filepath).name],
        )
        monkeypatch.setattr(
            "chemsmart.analysis.redox.solvent_scf_energy",
            lambda filepath: {
                "ox_solv.log": 0.90,
                "red_solv.log": 1.00,
                "ref_ox_solv.log": 1.80,
                "ref_red_solv.log": 2.00,
            }[Path(filepath).name],
        )
        result = CliRunner().invoke(
            run,
            [
                "redox",
                "analyze",
                "--e-ref",
                "0.25",
                "-n",
                "2",
                *_analyze_file_args(files),
            ],
        )
        assert result.exit_code == 0, result.output
        assert "n = 2" in result.output
        assert "E_ref = 0.2500 V" in result.output

    def test_chain_redox_analyze_n2_arithmetic(
        self, isolated_redox_registry, tmp_path, monkeypatch
    ):
        from chemsmart.cli.chain.chain import chain

        files = {
            name: _write_h2_xyz(tmp_path / name)
            for name in (
                "ox_gas.log",
                "red_gas.log",
                "ref_ox_gas.log",
                "ref_red_gas.log",
                "ox_solv.log",
                "red_solv.log",
                "ref_ox_solv.log",
                "ref_red_solv.log",
            )
        }
        monkeypatch.setattr(
            "chemsmart.analysis.redox.gas_phase_data",
            lambda filepath, **kwargs: (0.0, 0.0),
        )
        monkeypatch.setattr(
            "chemsmart.analysis.redox.solvent_scf_energy",
            lambda filepath: 0.0,
        )
        result = CliRunner().invoke(
            chain,
            [
                "redox",
                "analyze",
                "--e-ref",
                "0.25",
                "-n",
                "2",
                *_analyze_file_args(files),
            ],
            obj={"jobrunner": MagicMock()},
        )
        assert result.exit_code == 0, result.output
        assert "n = 2" in result.output
        assert "E_ref = 0.2500 V" in result.output
        assert "-p/--project" not in result.output
