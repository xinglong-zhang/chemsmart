import os
import shutil
from types import SimpleNamespace

import numpy as np
import pytest

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.chain_steps import (
    CHAIN_NESTED_ONLY_JOBS,
    CHAIN_PROGRAM_SETTINGS,
    PROGRAM_ADAPTERS,
    ChainStep,
    ChainStepParser,
    build_chain_job,
    molecule_from_completed_job,
    parse_chain_step,
    validate_chain_step,
)
from chemsmart.jobs.crest.conformers import CRESTConformerSearchJob
from chemsmart.jobs.crest.job import CRESTJob
from chemsmart.jobs.crest.settings import CRESTJobSettings
from chemsmart.jobs.gaussian.crest import GaussianCrestJob
from chemsmart.jobs.gaussian.irc import GaussianIRCJob
from chemsmart.jobs.gaussian.link import GaussianLinkJob
from chemsmart.jobs.gaussian.opt import GaussianOptJob
from chemsmart.jobs.gaussian.qrc import GaussianQRCJob
from chemsmart.jobs.gaussian.settings import (
    GaussianIRCJobSettings,
    GaussianJobSettings,
    GaussianLinkJobSettings,
)
from chemsmart.jobs.gaussian.singlepoint import GaussianSinglePointJob
from chemsmart.jobs.gaussian.ts import GaussianTSJob
from chemsmart.jobs.grouper.job import GrouperJob
from chemsmart.jobs.job import Job
from chemsmart.jobs.orca.opt import ORCAOptJob
from chemsmart.jobs.orca.qrc import ORCAQRCJob
from chemsmart.jobs.orca.settings import ORCAJobSettings
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


def _pipeline_job_names(program):
    return ChainStepParser().pipeline_job_names(program)


def _program_commands():
    return PROGRAM_ADAPTERS


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
        "crest": GaussianCrestJob,
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

    def test_subcommand_options_stay_after_the_job(self):
        step = parse_chain_step("gaussian:opt -f 1-3")
        assert step == ChainStep(
            program="gaussian",
            job="opt",
            suffix_option_tokens=("-f", "1-3"),
        )

    def test_crest_num_conformers_after_job(self):
        step = parse_chain_step("crest:conformers -nc 10")
        assert step == ChainStep(
            program="crest", job="conformers", num_conformers=10
        )

    def test_crest_num_conformers_rejects_bad_values(self):
        with pytest.raises(ValueError, match="integer >= 1"):
            parse_chain_step("crest:conformers -nc 0")
        with pytest.raises(ValueError, match="use -nc N"):
            parse_chain_step("crest:conformers -nc")

    def test_gaussian_crest_keeps_options_after_job(self):
        step = parse_chain_step(
            "gaussian: -l structure_from_lowest -c 0 -m 1 crest -j opt -nc 10"
        )
        assert step == ChainStep(
            program="gaussian",
            job="crest",
            extra_option_tokens=(
                "-l",
                "structure_from_lowest",
                "-c",
                "0",
                "-m",
                "1",
            ),
            suffix_option_tokens=("-j", "opt", "-nc", "10"),
        )

    def test_grouper_keeps_strategy_options(self):
        step = parse_chain_step("grouper: -T 0.3 -N 10 irmsd --inversion on")
        assert step.program == "grouper"
        assert step.job == "irmsd"
        assert step.extra_option_tokens == ("-T", "0.3", "-N", "10")
        assert step.suffix_option_tokens == ("--inversion", "on")

    def test_missing_colon_errors(self):
        with pytest.raises(ValueError, match="expected PROGRAM:JOB"):
            parse_chain_step("gaussian")

    def test_empty_parts_error(self):
        with pytest.raises(ValueError, match="expected PROGRAM:JOB"):
            parse_chain_step("gaussian:")


class TestChainStepRegistry:
    def test_chain_program_registry_matches_project_settings(self):
        assert tuple(CHAIN_PROGRAM_SETTINGS) == ChainProjectSettings.PROGRAMS
        assert tuple(_program_commands()) == tuple(CHAIN_PROGRAM_SETTINGS)

    def test_pipeline_jobs_are_program_commands(self):
        for program, jobs in _EXPECTED_SPECS.items():
            names = _pipeline_job_names(program)
            assert set(jobs).issubset(names)
            assert CHAIN_NESTED_ONLY_JOBS.isdisjoint(names)
        assert "irc" in _pipeline_job_names("gaussian")
        validate_chain_step("gaussian", "irc")
        assert "com" not in _pipeline_job_names("gaussian")
        assert "inp" not in _pipeline_job_names("orca")

    def test_file_and_trajectory_commands_are_rejected(self):
        for program, job in (
            ("gaussian", "com"),
            ("gaussian", "traj"),
            ("gaussian", "dias"),
            ("orca", "inp"),
        ):
            with pytest.raises(
                ValueError,
                match=f"do not support {program} '{job}'",
            ):
                validate_chain_step(program, job)

    def test_grouper_strategy_comes_from_the_grouper_command(self):
        validate_chain_step("grouper", "irmsd")
        with pytest.raises(ValueError, match="job 'nope'"):
            validate_chain_step("grouper", "nope")

    def test_nested_only_jobs_are_explicit_workflow_names(self):
        assert CHAIN_NESTED_ONLY_JOBS == frozenset(
            {"pka", "fukui", "redox", "reaction"}
        )
        for program in CHAIN_PROGRAM_SETTINGS:
            assert CHAIN_NESTED_ONLY_JOBS.isdisjoint(
                _pipeline_job_names(program)
            )

    def test_nested_only_pka_errors(self):
        with pytest.raises(
            ValueError,
            match="do not support gaussian 'pka'",
        ):
            validate_chain_step("gaussian", "pka")
        assert "pka" in CHAIN_NESTED_ONLY_JOBS
        assert "redox" in CHAIN_NESTED_ONLY_JOBS

    def test_nested_only_redox_errors(self):
        with pytest.raises(
            ValueError,
            match="do not support gaussian 'redox'",
        ):
            validate_chain_step("gaussian", "redox")

    def test_removed_yaml_job_errors_as_nested(self):
        with pytest.raises(
            ValueError,
            match="do not support gaussian 'fukui'",
        ):
            validate_chain_step("gaussian", "fukui")

    def test_nested_qmmm_is_not_a_gaussian_step(self):
        with pytest.raises(
            ValueError,
            match="Unsupported chain step job 'qmmm' for program 'gaussian'",
        ):
            validate_chain_step("gaussian", "qmmm")

    def test_unknown_job_lists_supported_jobs(self):
        with pytest.raises(
            ValueError,
            match="Unsupported chain step job 'foo' for program 'gaussian'",
        ):
            validate_chain_step("gaussian", "foo")

    def test_unsupported_program_errors(self):
        with pytest.raises(
            ValueError, match="Unsupported chain step program 'foo'"
        ):
            validate_chain_step("foo", "opt")


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
            if program == "gaussian" and job == "crest":
                continue
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
            assert child.label.startswith("mol_00_")
            assert child.settings.charge == 0
            assert child.settings.multiplicity == 1

    def test_any_program_command_can_be_a_step(
        self, isolated_config_dir, temporary_working_dir
    ):
        settings = _chain_settings(gaussian="gas_solv")
        chain = build_chain_job(
            settings,
            steps=(ChainStep(program="gaussian", job="irc"),),
            molecule=_h2(),
            label="mol",
            charge=0,
            multiplicity=1,
        )
        child = chain.phases[0].resolve_jobs()[0]
        assert isinstance(child, GaussianIRCJob)
        assert child.label == "mol_00_gaussian_irc"
        assert child.settings.charge == 0

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
        best.energy = -2.0
        dummy = DummyCRESTJob(_h2(), best, complete=True)
        best.write_xyz(dummy.best_file, mode="w")
        chain.phases[0]._resolved_jobs = [dummy]
        children = chain.phases[1].resolve_jobs()
        assert len(children) == 1
        assert isinstance(children[0], XTBOptJob)
        assert np.allclose(children[0].molecule.positions, best.positions)


def _write_xyz(path, molecules):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as handle:
        for molecule in molecules:
            handle.write(f"{molecule.num_atoms}\ncomment\n")
            for symbol, position in zip(
                molecule.chemical_symbols, molecule.positions
            ):
                handle.write(
                    f"{symbol} {position[0]} {position[1]} {position[2]}\n"
                )


def _crest_with_conformers(molecules, charge=0, multiplicity=1):
    source = _h2(charge=charge, multiplicity=multiplicity)
    job = DummyCRESTJob(source, molecules[0], complete=True)
    _write_xyz(job.conformers_file, molecules)
    return job


class StubGaussianCrest(GaussianCrestJob):
    def __init__(self, molecules):
        self._outs = list(molecules)

    def is_complete(self):
        return True

    @property
    def num_confs_to_opt(self):
        return len(self._outs)

    @property
    def all_conformers_jobs(self):
        return [
            DummyHandoffJob(f"c{index}", molecule, complete=True)
            for index, molecule in enumerate(self._outs, start=1)
        ]


class TestCrestEnsembleHandoff:
    def test_gaussian_crest_reads_conformers_file(
        self, isolated_config_dir, temporary_working_dir
    ):
        settings = _chain_settings(crest="test", gaussian="gas_solv")
        crest_step = parse_chain_step("crest:conformers -nc 2")
        gaussian_step = parse_chain_step(
            "gaussian: -l structure_from_lowest -c 0 -m 1 crest -j opt -nc 10"
        )
        chain = build_chain_job(
            settings,
            steps=(crest_step, gaussian_step),
            molecule=_h2(),
            label="mol",
            charge=0,
            multiplicity=1,
        )
        frames = [_h2(z=0.70), _h2(z=0.90)]
        chain.phases[0]._resolved_jobs = [_crest_with_conformers(frames)]
        children = chain.phases[1].resolve_jobs()
        assert len(children) == 1
        job = children[0]
        assert isinstance(job, GaussianCrestJob)
        assert job.label == "structure_from_lowest_opt"
        assert job.num_confs_to_opt == 10
        assert len(job.all_conformers) == 2

    def test_omitted_num_confs_uses_best_conformer(
        self, isolated_config_dir, temporary_working_dir
    ):
        settings = _chain_settings(crest="test", gaussian="gas_solv")
        chain = build_chain_job(
            settings,
            steps=(
                ChainStep(program="crest", job="conformers"),
                parse_chain_step("gaussian: crest -j opt"),
            ),
            molecule=_h2(),
            label="mol",
            charge=0,
            multiplicity=1,
        )
        best = _h2(z=0.82)
        crest = _crest_with_conformers([_h2(z=0.70), _h2(z=0.90)])
        best.write_xyz(crest.best_file, mode="w")
        chain.phases[0]._resolved_jobs = [crest]
        job = chain.phases[1].resolve_jobs()[0]
        assert job.num_confs_to_opt == 1
        assert job.label == "mol_opt"
        assert np.allclose(job.all_conformers[0].positions, best.positions)

    def test_grouping_and_num_confs_are_rejected(
        self, isolated_config_dir, temporary_working_dir
    ):
        settings = _chain_settings(crest="test", gaussian="gas_solv")
        chain = build_chain_job(
            settings,
            steps=(
                ChainStep(program="crest", job="conformers"),
                parse_chain_step("gaussian: crest -j opt -nc 2 -g rmsd"),
            ),
            molecule=_h2(),
            label="mol",
            charge=0,
            multiplicity=1,
        )
        crest = _crest_with_conformers([_h2(z=0.70), _h2(z=0.90)])
        _h2(z=0.70).write_xyz(crest.best_file, mode="w")
        chain.phases[0]._resolved_jobs = [crest]
        with pytest.raises(ValueError, match="mutually exclusive"):
            chain.phases[1].resolve_jobs()

    def test_gaussian_crest_must_follow_crest_conformers(
        self, isolated_config_dir
    ):
        settings = _chain_settings(gaussian="gas_solv")
        with pytest.raises(ValueError, match="must follow crest:conformers"):
            build_chain_job(
                settings,
                steps=(ChainStep(program="gaussian", job="crest"),),
                molecule=_h2(),
                label="mol",
            )

    def test_next_step_uses_lowest_energy_conformer(
        self, isolated_config_dir, temporary_working_dir
    ):
        settings = _chain_settings(crest="test", gaussian="gas_solv")
        chain = build_chain_job(
            settings,
            steps=(
                ChainStep(program="crest", job="conformers"),
                ChainStep(program="gaussian", job="crest"),
                ChainStep(program="gaussian", job="opt"),
            ),
            molecule=_h2(),
            label="mol",
            charge=0,
            multiplicity=1,
        )
        higher = _h2(z=0.80)
        higher.energy = -1.0
        lowest = _h2(z=1.10)
        lowest.energy = -2.5
        chain.phases[1]._resolved_jobs = [StubGaussianCrest([higher, lowest])]
        children = chain.phases[2].resolve_jobs()
        assert len(children) == 1
        assert children[0].label == "mol_02_gaussian_opt"
        assert np.allclose(children[0].molecule.positions, lowest.positions)

    def test_next_step_waits_without_conformer_energies(
        self, isolated_config_dir, temporary_working_dir
    ):
        settings = _chain_settings(crest="test", gaussian="gas_solv")
        chain = build_chain_job(
            settings,
            steps=(
                ChainStep(program="crest", job="conformers"),
                ChainStep(program="gaussian", job="crest"),
                ChainStep(program="gaussian", job="opt"),
            ),
            molecule=_h2(),
            label="mol",
            charge=0,
            multiplicity=1,
        )
        chain.phases[1]._resolved_jobs = [StubGaussianCrest([_h2(z=0.80)])]
        assert chain.phases[2].resolve_jobs() == []

    def test_grouper_skips_a_single_conformer(
        self, isolated_config_dir, temporary_working_dir
    ):
        settings = _chain_settings(crest="test", xtb="test")
        chain = build_chain_job(
            settings,
            steps=(
                ChainStep(program="crest", job="conformers"),
                ChainStep(program="grouper", job="formula"),
                ChainStep(program="xtb", job="opt"),
            ),
            molecule=_h2(),
            label="mol",
        )
        frame = _h2(z=0.82, charge=3, multiplicity=2)
        crest = _crest_with_conformers([frame], charge=3, multiplicity=2)
        frame.write_xyz(crest.best_file, mode="w")
        chain.phases[0]._resolved_jobs = [crest]
        grouped = chain.phases[1].resolve_jobs()
        assert len(grouped) == 1
        assert grouped[0].is_complete()
        children = chain.phases[2].resolve_jobs()
        assert len(children) == 1
        assert children[0].label == "mol_02_xtb_opt"
        assert children[0].molecule.charge == 3
        assert children[0].molecule.multiplicity == 2
        assert np.allclose(children[0].molecule.positions, frame.positions)

    def test_grouper_representatives_are_one_job_each(
        self, isolated_config_dir, temporary_working_dir
    ):
        settings = _chain_settings(crest="test", gaussian="gas_solv")
        chain = build_chain_job(
            settings,
            steps=(
                ChainStep(program="crest", job="conformers"),
                ChainStep(program="grouper", job="irmsd"),
                ChainStep(program="gaussian", job="opt"),
            ),
            molecule=_h2(),
            label="mol",
            charge=0,
            multiplicity=1,
        )
        first = _h2(z=0.70, charge=1, multiplicity=2)
        second = _h2(z=1.20, charge=1, multiplicity=2)
        grouper = GrouperJob(
            molecules=[first, second],
            grouping_strategy="irmsd",
            label="grouped",
        )
        _write_xyz(
            os.path.join(grouper.output_dir, "grouped_group_1.xyz"),
            [first, second],
        )
        _write_xyz(
            os.path.join(grouper.output_dir, "grouped_group_2.xyz"),
            [second],
        )
        assert grouper.is_complete()
        chain.phases[1]._resolved_jobs = [grouper]
        children = chain.phases[2].resolve_jobs()
        assert [child.label for child in children] == [
            "mol_02_gaussian_opt_g1",
            "mol_02_gaussian_opt_g2",
        ]
        assert np.allclose(children[0].molecule.positions, first.positions)
        assert np.allclose(children[1].molecule.positions, second.positions)
        assert children[0].molecule.charge == 1
        assert children[1].molecule.multiplicity == 2

    def test_grouper_does_not_require_a_chain_alias(
        self, isolated_config_dir, temporary_working_dir
    ):
        settings = _chain_settings(crest="test", xtb="test")
        chain = build_chain_job(
            settings,
            steps=(
                parse_chain_step("crest:conformers -nc 2"),
                parse_chain_step("grouper: formula"),
                ChainStep(program="xtb", job="opt"),
            ),
            molecule=_h2(),
            label="mol",
        )
        frames = [_h2(z=0.70), _h2(z=0.90)]
        chain.phases[0]._resolved_jobs = [_crest_with_conformers(frames)]
        children = chain.phases[1].resolve_jobs()
        assert len(children) == 1
        assert isinstance(children[0], GrouperJob)
        assert children[0].grouping_strategy == "formula"
        assert len(children[0].molecules) == 2


class TestCrestNumConformers:
    def test_omitted_num_confs_uses_best_file(
        self, isolated_config_dir, temporary_working_dir
    ):
        higher = _h2(z=0.70)
        higher.energy = -1.0
        lowest = _h2(z=1.20)
        lowest.energy = -3.0
        best = _h2(z=0.95)
        settings = _chain_settings(crest="test", gaussian="gas_solv")
        chain = build_chain_job(
            settings,
            steps=(
                parse_chain_step("crest:conformers"),
                ChainStep(program="gaussian", job="opt"),
            ),
            molecule=_h2(),
            label="mol",
            charge=0,
            multiplicity=1,
        )
        crest = _crest_with_conformers([higher, lowest])
        best.write_xyz(crest.best_file, mode="w")
        chain.phases[0]._resolved_jobs = [crest]
        children = chain.phases[1].resolve_jobs()
        assert len(children) == 1
        assert children[0].label == "mol_01_gaussian_opt"
        assert np.allclose(children[0].molecule.positions, best.positions)

    def test_one_selects_the_lowest_conformer(
        self, isolated_config_dir, temporary_working_dir
    ):
        higher = _h2(z=0.70)
        higher.energy = -1.0
        lowest = _h2(z=1.20)
        lowest.energy = -3.0
        settings = _chain_settings(crest="test", gaussian="gas_solv")
        chain = build_chain_job(
            settings,
            steps=(
                parse_chain_step("crest:conformers -nc 1"),
                ChainStep(program="gaussian", job="opt"),
            ),
            molecule=_h2(),
            label="mol",
            charge=0,
            multiplicity=1,
        )
        frames = [higher, lowest]
        crest = _crest_with_conformers(frames)
        for index, molecule in enumerate(frames):
            molecule.write_xyz(
                crest.conformers_file, mode="w" if index == 0 else "a"
            )
        chain.phases[0]._resolved_jobs = [crest]
        children = chain.phases[1].resolve_jobs()
        assert len(children) == 1
        assert children[0].label == "mol_01_gaussian_opt"
        assert np.allclose(children[0].molecule.positions, lowest.positions)

    def test_selects_lowest_conformers(
        self, isolated_config_dir, temporary_working_dir
    ):
        settings = _chain_settings(crest="test", gaussian="gas_solv")
        chain = build_chain_job(
            settings,
            steps=(
                parse_chain_step("crest:conformers -nc 2"),
                ChainStep(program="gaussian", job="opt"),
            ),
            molecule=_h2(),
            label="mol",
            charge=0,
            multiplicity=1,
        )
        higher = _h2(z=0.70, charge=1, multiplicity=2)
        higher.energy = -1.0
        lowest = _h2(z=1.20, charge=1, multiplicity=2)
        lowest.energy = -3.0
        middle = _h2(z=0.90, charge=1, multiplicity=2)
        middle.energy = -2.0
        frames = [higher, lowest, middle]
        crest = _crest_with_conformers(frames, charge=1, multiplicity=2)
        for index, molecule in enumerate(frames):
            molecule.write_xyz(
                crest.conformers_file, mode="w" if index == 0 else "a"
            )
        chain.phases[0]._resolved_jobs = [crest]
        children = chain.phases[1].resolve_jobs()
        assert [child.label for child in children] == [
            "mol_01_gaussian_opt_c1",
            "mol_01_gaussian_opt_c2",
        ]
        assert np.allclose(children[0].molecule.positions, lowest.positions)
        assert np.allclose(children[1].molecule.positions, middle.positions)
        assert children[0].molecule.charge == 1
        assert children[1].molecule.multiplicity == 2

    def test_uses_available_conformers_when_fewer_than_requested(
        self, isolated_config_dir, temporary_working_dir
    ):
        settings = _chain_settings(crest="test", xtb="test")
        chain = build_chain_job(
            settings,
            steps=(
                parse_chain_step("crest:conformers -nc 10"),
                ChainStep(program="xtb", job="opt"),
            ),
            molecule=_h2(),
            label="mol",
        )
        frames = [_h2(z=0.70), _h2(z=0.90)]
        chain.phases[0]._resolved_jobs = [_crest_with_conformers(frames)]
        children = chain.phases[1].resolve_jobs()
        assert len(children) == 2
        assert children[0].label == "mol_01_xtb_opt_c1"

    def test_waits_until_conformers_file_exists(
        self, isolated_config_dir, temporary_working_dir
    ):
        settings = _chain_settings(crest="test", gaussian="gas_solv")
        chain = build_chain_job(
            settings,
            steps=(
                parse_chain_step("crest:conformers -nc 2"),
                ChainStep(program="gaussian", job="opt"),
            ),
            molecule=_h2(),
            label="mol",
            charge=0,
            multiplicity=1,
        )
        chain.phases[0]._resolved_jobs = [
            DummyCRESTJob(_h2(), _h2(z=0.80), complete=True)
        ]
        assert chain.phases[1].resolve_jobs() == []

    def test_gaussian_crest_reads_selected_conformers(
        self, isolated_config_dir, temporary_working_dir
    ):
        settings = _chain_settings(crest="test", gaussian="gas_solv")
        chain = build_chain_job(
            settings,
            steps=(
                parse_chain_step("crest:conformers -nc 2"),
                parse_chain_step("gaussian: crest -j opt"),
            ),
            molecule=_h2(),
            label="mol",
            charge=0,
            multiplicity=1,
        )
        higher = _h2(z=0.70)
        higher.energy = -1.0
        lowest = _h2(z=1.20)
        lowest.energy = -3.0
        middle = _h2(z=0.90)
        middle.energy = -2.0
        frames = [higher, lowest, middle]
        crest = _crest_with_conformers(frames)
        for index, molecule in enumerate(frames):
            molecule.write_xyz(
                crest.conformers_file, mode="w" if index == 0 else "a"
            )
        chain.phases[0]._resolved_jobs = [crest]
        job = chain.phases[1].resolve_jobs()[0]
        assert isinstance(job, GaussianCrestJob)
        assert job.num_confs_to_opt == 2
        assert len(job.all_conformers) == 2
        assert np.allclose(job.all_conformers[0].positions, lowest.positions)
        assert np.allclose(job.all_conformers[1].positions, middle.positions)


class TestIRCHandoff:
    @staticmethod
    def _irc_with_endpoints(monkeypatch, endpoints):
        job = GaussianIRCJob(
            molecule=_h2(),
            settings=GaussianIRCJobSettings(),
            label="mol_00_gaussian_irc",
        )
        monkeypatch.setattr(job, "endpoint_jobs", lambda: list(endpoints))
        return job

    def test_endpoint_jobs_are_the_directional_logs(
        self, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        job = GaussianIRCJob(
            molecule=_h2(),
            settings=GaussianIRCJobSettings(),
            label="mol_02_gaussian_irc",
        )
        logs = [
            os.path.basename(item.outputfile) for item in job.endpoint_jobs()
        ]
        assert logs == [
            "mol_02_gaussian_ircf.log",
            "mol_02_gaussian_ircr.log",
        ]

        forward = GaussianIRCJob(
            molecule=_h2(),
            settings=GaussianIRCJobSettings(direction="forward"),
            label="mol_02_gaussian_irc",
        )
        assert len(forward.endpoint_jobs()) == 1
        assert forward.endpoint_jobs()[0].outputfile.endswith(
            "mol_02_gaussian_irc.log"
        )

    def test_both_endpoints_become_following_jobs(
        self, isolated_config_dir, temporary_working_dir, monkeypatch
    ):
        settings = _chain_settings(gaussian="gas_solv")
        chain = build_chain_job(
            settings,
            steps=(
                ChainStep(program="gaussian", job="irc"),
                ChainStep(program="gaussian", job="opt"),
            ),
            molecule=_h2(),
            label="mol",
            charge=0,
            multiplicity=1,
        )
        forward = _h2(z=0.90)
        reverse = _h2(z=1.20)
        chain.phases[0]._resolved_jobs = [
            self._irc_with_endpoints(
                monkeypatch,
                [
                    DummyHandoffJob("ircf", forward, complete=True),
                    DummyHandoffJob("ircr", reverse, complete=True),
                ],
            )
        ]
        children = chain.phases[1].resolve_jobs()
        assert [child.label for child in children] == [
            "mol_01_gaussian_opt_c1",
            "mol_01_gaussian_opt_c2",
        ]
        assert np.allclose(children[0].molecule.positions, forward.positions)
        assert np.allclose(children[1].molecule.positions, reverse.positions)

    def test_missing_endpoint_geometry_waits(
        self, isolated_config_dir, temporary_working_dir, monkeypatch
    ):
        settings = _chain_settings(gaussian="gas_solv")
        chain = build_chain_job(
            settings,
            steps=(
                ChainStep(program="gaussian", job="irc"),
                ChainStep(program="gaussian", job="opt"),
            ),
            molecule=_h2(),
            label="mol",
            charge=0,
            multiplicity=1,
        )
        chain.phases[0]._resolved_jobs = [
            self._irc_with_endpoints(
                monkeypatch,
                [
                    DummyHandoffJob("ircf", _h2(z=0.90), complete=True),
                    DummySPJob("ircr", _h2(), None, complete=True),
                ],
            )
        ]
        assert chain.phases[1].resolve_jobs() == []


class TestCompositeHandoff:
    def _chain(self, **aliases):
        return build_chain_job(
            _chain_settings(**aliases),
            steps=(
                ChainStep(program="gaussian", job="opt"),
                ChainStep(program="gaussian", job="opt"),
            ),
            molecule=_h2(),
            label="mol",
            charge=0,
            multiplicity=1,
        )

    def test_qrc_directions_become_following_jobs(
        self,
        isolated_config_dir,
        temporary_working_dir,
        gaussian_outputs_test_directory,
        monkeypatch,
    ):
        chain = self._chain(gaussian="gas_solv")
        molecule = Molecule.from_filepath(
            os.path.join(gaussian_outputs_test_directory, "KOH.log")
        )
        forward = _h2(z=0.90)
        reverse = _h2(z=1.20)
        qrc = GaussianQRCJob(
            molecule=molecule,
            settings=GaussianJobSettings.default(),
            label="mol_00_gaussian_qrc",
        )
        monkeypatch.setattr(
            qrc,
            "_prepare_both_qrc_jobs",
            lambda: [
                DummyHandoffJob("qrcf", forward, complete=True),
                DummyHandoffJob("qrcr", reverse, complete=True),
            ],
        )
        chain.phases[0]._resolved_jobs = [qrc]
        children = chain.phases[1].resolve_jobs()
        assert [child.label for child in children] == [
            "mol_01_gaussian_opt_c1",
            "mol_01_gaussian_opt_c2",
        ]
        assert np.allclose(children[0].molecule.positions, forward.positions)
        assert np.allclose(children[1].molecule.positions, reverse.positions)

    def test_orca_qrc_directions_become_following_jobs(
        self,
        isolated_config_dir,
        temporary_working_dir,
        gaussian_outputs_test_directory,
        monkeypatch,
    ):
        chain = build_chain_job(
            _chain_settings(orca="gas_solv"),
            steps=(
                ChainStep(program="orca", job="qrc"),
                ChainStep(program="orca", job="opt"),
            ),
            molecule=_h2(),
            label="mol",
            charge=0,
            multiplicity=1,
        )
        molecule = Molecule.from_filepath(
            os.path.join(gaussian_outputs_test_directory, "KOH.log")
        )
        forward = _h2(z=0.80)
        reverse = _h2(z=1.10)
        qrc = ORCAQRCJob(
            molecule=molecule,
            settings=ORCAJobSettings.default(),
            label="mol_00_orca_qrc",
        )
        monkeypatch.setattr(
            qrc,
            "_prepare_both_qrc_jobs",
            lambda: [
                DummyHandoffJob("qrcf", forward, complete=True),
                DummyHandoffJob("qrcr", reverse, complete=True),
            ],
        )
        chain.phases[0]._resolved_jobs = [qrc]
        children = chain.phases[1].resolve_jobs()
        assert [child.label for child in children] == [
            "mol_01_orca_opt_c1",
            "mol_01_orca_opt_c2",
        ]
        assert np.allclose(children[0].molecule.positions, forward.positions)
        assert np.allclose(children[1].molecule.positions, reverse.positions)

    def test_link_irc_directions_become_following_jobs(
        self, isolated_config_dir, temporary_working_dir, monkeypatch
    ):
        chain = self._chain(gaussian="gas_solv")
        forward = _h2(z=0.70)
        reverse = _h2(z=1.30)
        link = GaussianLinkJob(
            molecule=_h2(),
            settings=GaussianLinkJobSettings(jobtype="irc"),
            label="mol_00_gaussian_link",
        )
        monkeypatch.setattr(link, "is_complete", lambda: True)
        monkeypatch.setattr(
            link,
            "_get_irc_jobs",
            lambda: [
                DummyHandoffJob("ircf", forward, complete=True),
                DummyHandoffJob("ircr", reverse, complete=True),
            ],
        )
        chain.phases[0]._resolved_jobs = [link]
        children = chain.phases[1].resolve_jobs()
        assert [child.label for child in children] == [
            "mol_01_gaussian_opt_c1",
            "mol_01_gaussian_opt_c2",
        ]
        assert np.allclose(children[0].molecule.positions, forward.positions)
        assert np.allclose(children[1].molecule.positions, reverse.positions)


class TestOutputFileHandoff:
    def test_gaussian_log_keeps_vibrational_modes(
        self,
        isolated_config_dir,
        temporary_working_dir,
        gaussian_outputs_test_directory,
    ):
        log = os.path.join(gaussian_outputs_test_directory, "KOH.log")
        chain = build_chain_job(
            _chain_settings(gaussian="gas_solv"),
            steps=(
                ChainStep(program="gaussian", job="ts"),
                ChainStep(program="gaussian", job="qrc"),
            ),
            molecule=Molecule.from_filepath(log),
            label="mol",
            charge=0,
            multiplicity=1,
        )
        ts_job = chain.phases[0].resolve_jobs()[0]
        shutil.copy(log, ts_job.outputfile)
        qrc_job = chain.phases[1].resolve_jobs()[0]
        assert isinstance(qrc_job, GaussianQRCJob)
        assert len(qrc_job.molecule.vibrational_modes) == 3
        assert not os.path.exists("mol_01_gaussian_qrc.xyz")

    def test_orca_output_feeds_gaussian(
        self,
        isolated_config_dir,
        temporary_working_dir,
        orca_outputs_directory,
    ):
        out = os.path.join(orca_outputs_directory, "phenol_pka_B_sp.out")
        expected = Molecule.from_filepath(out)
        chain = build_chain_job(
            _chain_settings(orca="gas_solv", gaussian="gas_solv"),
            steps=(
                ChainStep(program="orca", job="sp"),
                ChainStep(program="gaussian", job="opt"),
            ),
            molecule=expected,
            label="mol",
            charge=-1,
            multiplicity=1,
        )
        orca_job = chain.phases[0].resolve_jobs()[0]
        shutil.copy(out, orca_job.outputfile)
        child = chain.phases[1].resolve_jobs()[0]
        assert isinstance(child, GaussianOptJob)
        assert child.settings.charge == -1
        assert np.allclose(child.molecule.positions, expected.positions)
        assert not os.path.exists("mol_01_gaussian_opt.xyz")

    def test_xtb_output_feeds_orca(
        self, isolated_config_dir, temporary_working_dir, xtb_water_outfolder
    ):
        chain = build_chain_job(
            _chain_settings(xtb="test", orca="gas_solv"),
            steps=(
                ChainStep(program="xtb", job="opt"),
                ChainStep(program="orca", job="sp"),
            ),
            molecule=Molecule.from_filepath(
                os.path.join(xtb_water_outfolder, "water.xyz")
            ),
            label="mol",
            charge=0,
            multiplicity=1,
        )
        xtb_job = chain.phases[0].resolve_jobs()[0]
        shutil.copytree(
            xtb_water_outfolder, xtb_job.folder, dirs_exist_ok=True
        )
        os.rename(
            os.path.join(xtb_job.folder, "water_ohess.out"), xtb_job.outputfile
        )
        expected = xtb_job.optimized_structure()
        child = chain.phases[1].resolve_jobs()[0]
        assert isinstance(child, ORCASinglePointJob)
        assert np.allclose(
            child.molecule.positions, expected.positions, atol=1e-4
        )
        assert child.molecule.energy == expected.energy

    def test_later_step_reads_charge_from_the_log(
        self,
        isolated_config_dir,
        temporary_working_dir,
        gaussian_outputs_test_directory,
    ):
        log = os.path.join(
            gaussian_outputs_test_directory, "collidine-H_opt.log"
        )
        chain = build_chain_job(
            _chain_settings(gaussian="gas_solv"),
            steps=(
                ChainStep(program="gaussian", job="opt"),
                ChainStep(program="gaussian", job="sp"),
            ),
            molecule=Molecule.from_filepath(log),
            label="mol",
            charge=0,
            multiplicity=1,
        )
        opt_job = chain.phases[0].resolve_jobs()[0]
        assert opt_job.settings.charge == 0
        shutil.copy(log, opt_job.outputfile)
        sp_job = chain.phases[1].resolve_jobs()[0]
        assert sp_job.settings.charge == 1
        assert sp_job.settings.multiplicity == 1

    def test_step_charge_overrides_the_handed_off_log(
        self,
        isolated_config_dir,
        temporary_working_dir,
        gaussian_outputs_test_directory,
    ):
        log = os.path.join(
            gaussian_outputs_test_directory, "collidine-H_opt.log"
        )
        chain = build_chain_job(
            _chain_settings(gaussian="gas_solv"),
            steps=(
                ChainStep(program="gaussian", job="opt"),
                parse_chain_step("gaussian: -c 3 -m 2 sp"),
            ),
            molecule=Molecule.from_filepath(log),
            label="mol",
            charge=0,
            multiplicity=1,
        )
        opt_job = chain.phases[0].resolve_jobs()[0]
        shutil.copy(log, opt_job.outputfile)
        sp_job = chain.phases[1].resolve_jobs()[0]
        assert sp_job.settings.charge == 3
        assert sp_job.settings.multiplicity == 2
        assert not os.path.exists("mol_01_orca_sp.xyz")
