from types import SimpleNamespace

import numpy as np
import pytest

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.chain_steps import (
    _CHAIN_PROGRAMS,
    CHAIN_NESTED_ONLY_JOBS,
    CHAIN_STEP_SPECS,
    ChainStep,
    build_chain_job,
    get_chain_step_spec,
    molecule_from_completed_job,
    parse_chain_step,
)
from chemsmart.jobs.crest.conformers import CRESTConformerSearchJob
from chemsmart.jobs.crest.job import CRESTJob
from chemsmart.jobs.crest.settings import CRESTJobSettings
from chemsmart.jobs.gaussian.opt import GaussianOptJob
from chemsmart.jobs.gaussian.singlepoint import GaussianSinglePointJob
from chemsmart.jobs.gaussian.ts import GaussianTSJob
from chemsmart.jobs.job import Job
from chemsmart.jobs.orca.opt import ORCAOptJob
from chemsmart.jobs.orca.singlepoint import ORCASinglePointJob
from chemsmart.jobs.orca.ts import ORCATSJob
from chemsmart.jobs.xtb.hess import XTBHessJob
from chemsmart.jobs.xtb.opt import XTBOptJob
from chemsmart.jobs.xtb.singlepoint import XTBSinglePointJob
from chemsmart.settings.chain import ChainProjectSettings
from chemsmart.settings.user import CHEMSMARTUserSettings


@pytest.fixture()
def isolated_config_dir(tmp_path, monkeypatch):
    config_dir = tmp_path / "chemsmart_config"
    config_dir.mkdir()
    monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_dir))
    monkeypatch.setattr(
        CHEMSMARTUserSettings,
        "USER_CONFIG_DIR",
        str(config_dir),
    )
    return config_dir


def _h2(z=0.74, charge=0, multiplicity=1):
    return Molecule(
        symbols=["H", "H"],
        positions=[[0.0, 0.0, 0.0], [0.0, 0.0, z]],
        charge=charge,
        multiplicity=multiplicity,
    )


def _chain_settings(**aliases):
    return ChainProjectSettings(aliases=aliases, project_name="combined")


class DummyHandoffJob(Job):
    TYPE = "dummy"

    def __init__(self, label, molecule, complete=False):
        super().__init__(molecule=molecule, label=label, jobrunner=None)
        self._complete = complete
        self._handoff = molecule

    def is_complete(self):
        return self._complete

    def optimized_structure(self):
        if not self._complete:
            return None
        return self._handoff

    def _run(self, **kwargs):
        self._complete = True


class DummySPJob(Job):
    TYPE = "dummy-sp"

    def __init__(self, label, molecule, output_molecule, complete=True):
        super().__init__(molecule=molecule, label=label, jobrunner=None)
        self._complete = complete
        self._output_molecule = output_molecule

    def is_complete(self):
        return self._complete

    def optimized_structure(self):
        return None

    def _output(self):
        return SimpleNamespace(molecule=self._output_molecule)

    def _run(self, **kwargs):
        self._complete = True


class DummyCRESTJob(CRESTJob):
    def __init__(self, molecule, best, complete=True):
        settings = CRESTJobSettings.default()
        settings.charge = 0
        settings.multiplicity = 1
        super().__init__(
            molecule=molecule,
            settings=settings,
            label="crest_dummy",
        )
        self._complete = complete
        self._best = best

    def is_complete(self):
        return self._complete

    def _output(self):
        return SimpleNamespace(best_conformer=self._best)


_EXPECTED_SPECS = {
    "crest": {"conformers": CRESTConformerSearchJob},
    "xtb": {
        "opt": XTBOptJob,
        "sp": XTBSinglePointJob,
        "hess": XTBHessJob,
    },
    "gaussian": {
        "opt": GaussianOptJob,
        "ts": GaussianTSJob,
        "sp": GaussianSinglePointJob,
    },
    "orca": {
        "opt": ORCAOptJob,
        "ts": ORCATSJob,
        "sp": ORCASinglePointJob,
    },
}


class TestParseChainStep:
    def test_parses_program_job(self):
        assert parse_chain_step("gaussian:opt") == ChainStep(
            program="gaussian", job="opt"
        )

    def test_strips_whitespace(self):
        assert parse_chain_step(" gaussian : opt ") == ChainStep(
            program="gaussian", job="opt"
        )

    def test_parses_options_before_job(self):
        step = parse_chain_step("gaussian: -o maxstep=8,maxsize=12 opt")
        assert step == ChainStep(
            program="gaussian",
            job="opt",
            extra_option_tokens=("-o", "maxstep=8,maxsize=12"),
        )

    def test_options_after_job_error(self):
        with pytest.raises(ValueError, match="before the job name"):
            parse_chain_step("gaussian:opt -o maxstep=8")

    def test_missing_colon_errors(self):
        with pytest.raises(ValueError, match="expected PROGRAM:JOB"):
            parse_chain_step("gaussian")

    def test_empty_parts_error(self):
        with pytest.raises(ValueError, match="expected PROGRAM:JOB"):
            parse_chain_step("gaussian:")


class TestChainStepRegistry:
    def test_chain_program_registry_matches_project_settings(self):
        assert tuple(_CHAIN_PROGRAMS) == ChainProjectSettings.PROGRAMS

    def test_supported_pairs_map_to_job_classes(self):
        assert CHAIN_STEP_SPECS == _EXPECTED_SPECS
        for program, jobs in _EXPECTED_SPECS.items():
            for job, job_class in jobs.items():
                spec = get_chain_step_spec(program, job)
                assert spec.job_class is job_class

    def test_nested_only_jobs_are_explicit_workflow_names(self):
        assert CHAIN_NESTED_ONLY_JOBS == frozenset(
            {"pka", "fukui", "redox", "reaction"}
        )
        for jobs in CHAIN_STEP_SPECS.values():
            assert CHAIN_NESTED_ONLY_JOBS.isdisjoint(jobs)

    def test_nested_only_pka_errors(self):
        with pytest.raises(
            ValueError,
            match="do not support gaussian 'pka'",
        ):
            get_chain_step_spec("gaussian", "pka")
        assert "pka" in CHAIN_NESTED_ONLY_JOBS
        assert "redox" in CHAIN_NESTED_ONLY_JOBS

    def test_nested_only_redox_errors(self):
        with pytest.raises(
            ValueError,
            match="do not support gaussian 'redox'",
        ):
            get_chain_step_spec("gaussian", "redox")

    def test_removed_yaml_job_errors_as_nested(self):
        with pytest.raises(
            ValueError,
            match="do not support gaussian 'fukui'",
        ):
            get_chain_step_spec("gaussian", "fukui")

    def test_non_workflow_cli_job_is_unsupported_not_nested(self):
        with pytest.raises(
            ValueError,
            match="Unsupported chain step job 'qmmm' for program 'gaussian'",
        ):
            get_chain_step_spec("gaussian", "qmmm")
        with pytest.raises(
            ValueError,
            match="Unsupported chain step job 'irc' for program 'gaussian'",
        ):
            get_chain_step_spec("gaussian", "irc")

    def test_unknown_job_lists_supported_jobs(self):
        with pytest.raises(
            ValueError,
            match="Unsupported chain step job 'foo' for program 'gaussian'",
        ):
            get_chain_step_spec("gaussian", "foo")

    def test_unsupported_program_errors(self):
        with pytest.raises(
            ValueError, match="Unsupported chain step program 'foo'"
        ):
            get_chain_step_spec("foo", "opt")


class TestMoleculeFromCompletedJob:
    def test_returns_none_when_incomplete(self):
        mol = _h2()
        job = DummyHandoffJob("opt", mol, complete=False)
        assert molecule_from_completed_job(job) is None

    def test_returns_fallback_when_incomplete(self):
        mol = _h2()
        fallback = _h2(z=0.80)
        job = DummyHandoffJob("opt", mol, complete=False)
        assert molecule_from_completed_job(job, fallback=fallback) is fallback

    def test_opt_uses_optimized_structure(self):
        mol = _h2(z=0.90)
        job = DummyHandoffJob("opt", mol, complete=True)
        assert molecule_from_completed_job(job) is mol

    def test_opt_ignores_fallback_when_complete(self):
        mol = _h2(z=0.90)
        fallback = _h2(z=0.80)
        job = DummyHandoffJob("opt", mol, complete=True)
        assert molecule_from_completed_job(job, fallback=fallback) is mol

    def test_sp_falls_back_to_output_molecule(self):
        input_mol = _h2(z=0.74)
        output_mol = _h2(z=0.90)
        job = DummySPJob("sp", input_mol, output_mol, complete=True)
        assert molecule_from_completed_job(job) is output_mol

    def test_crest_uses_best_conformer(self, temporary_working_dir):
        input_mol = _h2(z=0.74)
        best = _h2(z=0.82)
        job = DummyCRESTJob(input_mol, best, complete=True)
        assert molecule_from_completed_job(job) is best


class TestBuildChainJob:
    def test_empty_steps_error(self):
        settings = _chain_settings(gaussian="gas_solv")
        with pytest.raises(ValueError, match="requires at least one step"):
            build_chain_job(settings, steps=(), molecule=_h2(), label="mol")

    def test_nested_only_pka_step_errors(self, isolated_config_dir):
        settings = _chain_settings(gaussian="gas_solv")
        with pytest.raises(
            ValueError,
            match="do not support gaussian 'pka'",
        ):
            build_chain_job(
                settings,
                steps=(ChainStep(program="gaussian", job="pka"),),
                molecule=_h2(),
                label="mol",
            )

    def test_supported_pairs_build_typed_children(
        self, isolated_config_dir, temporary_working_dir
    ):
        aliases = {
            "crest": "test",
            "xtb": "test",
            "gaussian": "gas_solv",
            "orca": "gas_solv",
        }
        cases = [
            (program, job, job_class)
            for program, jobs in _EXPECTED_SPECS.items()
            for job, job_class in jobs.items()
        ]
        mol = _h2()
        for program, job, job_class in cases:
            settings = _chain_settings(**{program: aliases[program]})
            chain = build_chain_job(
                settings,
                steps=(ChainStep(program=program, job=job),),
                molecule=mol,
                label="mol",
            )
            assert chain.label == "mol"
            assert chain.PROGRAM == "chain"
            assert chain.programs == (program,)
            assert len(chain.phases) == 1
            assert chain.phases[0].name == f"00_{program}_{job}"
            children = chain.phases[0].resolve_jobs()
            assert len(children) == 1
            child = children[0]
            assert isinstance(child, job_class)
            assert child.label == f"mol_00_{program}_{job}"
            assert child.settings.charge == 0
            assert child.settings.multiplicity == 1

    def test_duplicate_steps_get_unique_labels(
        self, isolated_config_dir, temporary_working_dir
    ):
        settings = _chain_settings(gaussian="gas_solv")
        steps = (
            ChainStep(program="gaussian", job="opt"),
            ChainStep(program="gaussian", job="opt"),
        )
        chain = build_chain_job(
            settings, steps=steps, molecule=_h2(), label="mol"
        )
        assert len(chain.phases) == 2
        assert chain.programs == ("gaussian",)
        assert chain.phases[0].name == "00_gaussian_opt"
        assert chain.phases[1].name == "01_gaussian_opt"
        first = chain.phases[0].resolve_jobs()[0]
        assert first.label == "mol_00_gaussian_opt"
        dummy = DummyHandoffJob("opt", _h2(z=0.90), complete=True)
        chain.phases[0]._resolved_jobs = [dummy]
        second = chain.phases[1].resolve_jobs()[0]
        assert second.label == "mol_01_gaussian_opt"

    def test_cli_charge_multiplicity_override_molecule(
        self, isolated_config_dir, temporary_working_dir
    ):
        settings = _chain_settings(gaussian="gas_solv")
        chain = build_chain_job(
            settings,
            steps=(ChainStep(program="gaussian", job="opt"),),
            molecule=_h2(charge=0, multiplicity=1),
            label="mol",
            charge=1,
            multiplicity=2,
        )
        child = chain.phases[0].resolve_jobs()[0]
        assert child.settings.charge == 1
        assert child.settings.multiplicity == 2

    def test_incomplete_previous_phase_factory_returns_no_jobs(
        self, isolated_config_dir, temporary_working_dir
    ):
        settings = _chain_settings(xtb="test", gaussian="gas_solv")
        steps = (
            ChainStep(program="xtb", job="opt"),
            ChainStep(program="gaussian", job="opt"),
        )
        chain = build_chain_job(
            settings, steps=steps, molecule=_h2(), label="mol"
        )
        first = DummyHandoffJob("xtb", _h2(z=0.90), complete=False)
        chain.phases[0]._resolved_jobs = [first]
        assert chain.phases[1].resolve_jobs() == []

    def test_geometry_handoff_from_completed_dummy_child(
        self, isolated_config_dir, temporary_working_dir
    ):
        settings = _chain_settings(xtb="test", gaussian="gas_solv")
        steps = (
            ChainStep(program="xtb", job="opt"),
            ChainStep(program="gaussian", job="opt"),
        )
        chain = build_chain_job(
            settings, steps=steps, molecule=_h2(), label="mol"
        )
        handed_off = _h2(z=0.91)
        dummy = DummyHandoffJob("xtb", handed_off, complete=True)
        chain.phases[0]._resolved_jobs = [dummy]
        children = chain.phases[1].resolve_jobs()
        assert len(children) == 1
        child = children[0]
        assert isinstance(child, GaussianOptJob)
        assert child.label == "mol_01_gaussian_opt"
        assert np.allclose(child.molecule.positions, handed_off.positions)

    def test_geometry_handoff_from_completed_dummy_crest(
        self, isolated_config_dir, temporary_working_dir
    ):
        settings = _chain_settings(crest="test", xtb="test")
        steps = (
            ChainStep(program="crest", job="conformers"),
            ChainStep(program="xtb", job="opt"),
        )
        chain = build_chain_job(
            settings, steps=steps, molecule=_h2(), label="mol"
        )
        best = _h2(z=0.85)
        dummy = DummyCRESTJob(_h2(), best, complete=True)
        chain.phases[0]._resolved_jobs = [dummy]
        children = chain.phases[1].resolve_jobs()
        assert len(children) == 1
        assert isinstance(children[0], XTBOptJob)
        assert np.allclose(children[0].molecule.positions, best.positions)
