from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from click.testing import CliRunner

from chemsmart.cli.chain.chain import chain
from chemsmart.cli.gaussian.gaussian import gaussian
from chemsmart.cli.run import run
from chemsmart.cli.sub import sub
from chemsmart.jobs.chain import ChainJob
from chemsmart.jobs.crest.conformers import CRESTConformerSearchJob
from chemsmart.jobs.gaussian.fukui import GaussianFukuiJob
from chemsmart.jobs.gaussian.pka import GaussianpKaJob
from chemsmart.jobs.gaussian.reaction import GaussianReactionJob
from chemsmart.jobs.gaussian.redox import GaussianRedoxJob
from chemsmart.jobs.orca.fukui import ORCAFukuiJob
from chemsmart.jobs.orca.reaction import ORCAReactionJob
from chemsmart.jobs.orca.redox import ORCARedoxJob
from chemsmart.settings.gaussian import GaussianProjectSettings
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


def _write_chain_yaml(isolated_config_dir, name, content):
    chain_dir = isolated_config_dir / "chain"
    chain_dir.mkdir(exist_ok=True)
    (chain_dir / f"{name}.yaml").write_text(content)
    return name


def _invoke_chain(args, obj=None, **extra):
    runner = CliRunner()
    if obj is None:
        obj = {"jobrunner": MagicMock()}
    return runner.invoke(chain, args, obj=obj, catch_exceptions=False, **extra)


class TestChainHelp:
    def test_run_help_lists_chain(self):
        result = CliRunner().invoke(run, ["--help"])
        assert result.exit_code == 0, result.output
        assert "chain" in result.output

    def test_sub_help_lists_chain(self):
        result = CliRunner().invoke(sub, ["--help"])
        assert result.exit_code == 0, result.output
        assert "chain" in result.output

    def test_chain_help_lists_custom_and_workflow_subcommands(self):
        result = CliRunner().invoke(run, ["chain", "--help"])
        assert result.exit_code == 0, result.output
        assert "\n  custom" in result.output
        assert "\n  pka" in result.output
        assert "\n  fukui" in result.output
        assert "\n  redox" in result.output
        assert "\n  reaction" in result.output
        assert "\n  gaussian" not in result.output
        assert "\n  orca" not in result.output
        assert "\n  xtb" not in result.output
        assert "\n  crest" not in result.output

    def test_run_help_lists_pka_fukui_and_redox_analysis(self):
        result = CliRunner().invoke(run, ["--help"])
        assert result.exit_code == 0, result.output
        assert "\n  pka" in result.output
        assert "\n  fukui" in result.output
        assert "\n  redox" in result.output
        assert "\n  reaction" not in result.output

    def test_sub_help_hides_chain_workflow_subcommands(self):
        result = CliRunner().invoke(sub, ["--help"])
        assert result.exit_code == 0, result.output
        assert "\n  pka" not in result.output
        assert "\n  fukui" not in result.output
        assert "\n  redox" not in result.output
        assert "\n  reaction" not in result.output


class TestChainPipeline:
    def test_pipeline_builds_chain_job(
        self, isolated_config_dir, single_molecule_xyz_file
    ):
        _write_chain_yaml(
            isolated_config_dir,
            "combined",
            """
crest: test
xtb: test
gaussian: gas_solv
orca: gas_solv
""",
        )
        result = _invoke_chain(
            [
                "-p",
                "combined",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "-s",
                "crest:conformers",
                "-s",
                "xtb:opt",
                "-s",
                "gaussian:opt",
                "-s",
                "orca:sp",
            ],
            standalone_mode=False,
        )
        assert result.exit_code == 0, result.output
        job = result.return_value
        assert isinstance(job, ChainJob)
        assert job.PROGRAM == "chain"
        assert job.programs == ("crest", "xtb", "gaussian", "orca")
        assert job.label == "crest_best"
        assert [phase.name for phase in job.phases] == [
            "00_crest_conformers",
            "01_xtb_opt",
            "02_gaussian_opt",
            "03_orca_sp",
        ]
        children = job.phases[0].resolve_jobs()
        assert len(children) == 1
        assert isinstance(children[0], CRESTConformerSearchJob)
        assert children[0].label == "crest_best_00_crest_conformers"

    def test_sub_test_writes_chain_pipeline_script(
        self,
        isolated_config_dir,
        tmp_path,
        server_yaml_file,
        single_molecule_xyz_file,
        monkeypatch,
    ):
        monkeypatch.chdir(tmp_path)
        _write_chain_yaml(
            isolated_config_dir,
            "combined",
            "crest: test\nxtb: test\ngaussian: gas_solv\n",
        )
        result = CliRunner().invoke(
            sub,
            [
                "-s",
                server_yaml_file,
                "--test",
                "chain",
                "-p",
                "combined",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "-s",
                "crest:conformers",
                "-s",
                "xtb:opt",
                "-s",
                "gaussian:opt",
            ],
            catch_exceptions=False,
        )
        assert result.exit_code == 0, result.output
        submit_script = tmp_path / "chemsmart_sub_crest_best.sh"
        run_script = tmp_path / "chemsmart_run_crest_best.py"
        assert submit_script.exists()
        assert run_script.exists()
        submit_text = submit_script.read_text()
        assert submit_text.count("conda activate") == 1
        assert "module load craype-x86-rome" in submit_text
        assert "g16.login" in submit_text
        assert "GAUSS_EXEDIR" in submit_text
        assert "'chain'" in run_script.read_text()

    def test_pipeline_without_steps_is_usage_error(
        self, isolated_config_dir, single_molecule_xyz_file
    ):
        _write_chain_yaml(
            isolated_config_dir,
            "combined",
            "gaussian: gas_solv\n",
        )
        result = CliRunner().invoke(
            chain,
            [
                "-p",
                "combined",
                "-f",
                single_molecule_xyz_file,
            ],
            obj={"jobrunner": MagicMock()},
        )
        assert result.exit_code == 2, result.output
        assert "requires at least one -s/--steps" in result.output

    def test_pipeline_without_filename_is_usage_error(
        self, isolated_config_dir
    ):
        _write_chain_yaml(
            isolated_config_dir,
            "combined",
            "gaussian: gas_solv\n",
        )
        result = CliRunner().invoke(
            chain,
            ["-p", "combined", "-s", "gaussian:opt"],
            obj={"jobrunner": MagicMock()},
        )
        assert result.exit_code == 2, result.output
        assert "requires a structure file" in result.output

    def test_unknown_nested_command_is_usage_error(
        self, isolated_config_dir, single_molecule_xyz_file
    ):
        _write_chain_yaml(
            isolated_config_dir,
            "combined",
            "gaussian: gas_solv\n",
        )
        result = CliRunner().invoke(
            chain,
            [
                "-p",
                "combined",
                "-f",
                single_molecule_xyz_file,
                "opt",
            ],
            obj={"jobrunner": MagicMock()},
        )
        assert result.exit_code == 2, result.output
        assert "No such command 'opt'" in result.output

    def test_invalid_step_token_is_usage_error(
        self, isolated_config_dir, single_molecule_xyz_file
    ):
        _write_chain_yaml(
            isolated_config_dir,
            "combined",
            "gaussian: gas_solv\n",
        )
        result = CliRunner().invoke(
            chain,
            [
                "-p",
                "combined",
                "-f",
                single_molecule_xyz_file,
                "-s",
                "gaussian",
            ],
            obj={"jobrunner": MagicMock()},
        )
        assert result.exit_code == 2, result.output
        assert "expected PROGRAM:JOB" in result.output

    def test_nested_only_step_is_usage_error(
        self, isolated_config_dir, single_molecule_xyz_file
    ):
        _write_chain_yaml(
            isolated_config_dir,
            "combined",
            "gaussian: gas_solv\n",
        )
        result = CliRunner().invoke(
            chain,
            [
                "-p",
                "combined",
                "-f",
                single_molecule_xyz_file,
                "-s",
                "gaussian:pka",
            ],
            obj={"jobrunner": MagicMock()},
        )
        assert result.exit_code == 2, result.output
        assert "do not support gaussian 'pka'" in result.output


class TestChainPipelineExtras:
    def test_pipeline_quoted_step_applies_opt_options(
        self, isolated_config_dir, single_molecule_xyz_file
    ):
        _write_chain_yaml(
            isolated_config_dir,
            "combined",
            "gaussian: gas_solv\n",
        )
        result = _invoke_chain(
            [
                "-p",
                "combined",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "-s",
                "gaussian: -o maxstep=8,maxsize=12 opt",
            ],
            standalone_mode=False,
        )
        assert result.exit_code == 0, result.output
        job = result.return_value
        assert isinstance(job, ChainJob)
        child = job.phases[0].resolve_jobs()[0]
        assert (
            child.settings.additional_opt_options_in_route
            == "maxstep=8,maxsize=12"
        )

    def test_custom_subcommand_builds_pipeline(
        self, isolated_config_dir, single_molecule_xyz_file
    ):
        _write_chain_yaml(
            isolated_config_dir,
            "combined",
            "gaussian: gas_solv\n",
        )
        result = _invoke_chain(
            [
                "-p",
                "combined",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "-s",
                "gaussian:opt",
                "custom",
            ],
            standalone_mode=False,
        )
        assert result.exit_code == 0, result.output
        assert isinstance(result.return_value, ChainJob)
        assert [phase.name for phase in result.return_value.phases] == [
            "00_gaussian_opt"
        ]

    def test_missing_chain_project_is_usage_error(self, isolated_config_dir):
        result = CliRunner().invoke(
            chain,
            ["-p", "missing", "custom"],
            obj={"jobrunner": MagicMock()},
        )
        assert result.exit_code == 2, result.output

    def test_run_gaussian_without_chain_uses_own_project(
        self, single_molecule_xyz_file
    ):
        args = [
            "-p",
            "test",
            "-f",
            single_molecule_xyz_file,
            "-c",
            "0",
            "-m",
            "1",
            "opt",
        ]
        with (
            patch(
                "chemsmart.settings.gaussian.GaussianProjectSettings."
                "from_project",
                wraps=GaussianProjectSettings.from_project,
            ) as gaussian_from_project,
            patch(
                "chemsmart.jobs.gaussian.opt.GaussianOptJob",
                return_value=MagicMock(),
            ),
        ):
            result = CliRunner().invoke(
                gaussian,
                args,
                obj={"jobrunner": MagicMock()},
                catch_exceptions=False,
            )
        assert result.exit_code == 0, result.output
        gaussian_from_project.assert_called_once_with("test")


def _water_xyz(isolated_config_dir):
    path = isolated_config_dir / "water.xyz"
    path.write_text(
        "3\nwater\n"
        "O  0.000  0.000  0.000\n"
        "H  0.957  0.000  0.000\n"
        "H -0.240  0.927  0.000\n"
    )
    return str(path)


def _combined_yaml(isolated_config_dir):
    return _write_chain_yaml(
        isolated_config_dir,
        "combined",
        "gaussian: gas_solv\norca: gas_solv\n",
    )


def _ref_xyz_pair(directory):
    ox = directory / "ref_ox.xyz"
    red = directory / "ref_red.xyz"
    ox.write_text("2\nref ox\nH 0.0 0.0 0.0\nH 0.0 0.0 0.74\n")
    red.write_text("2\nref red\nH 0.0 0.0 0.0\nH 0.0 0.0 0.74\n")
    return str(ox), str(red)


class TestChainWorkflows:
    def test_pka_program_gaussian_uses_chain_alias_project(
        self, isolated_config_dir
    ):
        _combined_yaml(isolated_config_dir)
        water = _water_xyz(isolated_config_dir)
        with patch(
            "chemsmart.settings.gaussian.GaussianProjectSettings.from_project",
            wraps=GaussianProjectSettings.from_project,
        ) as gaussian_from_project:
            result = _invoke_chain(
                [
                    "-p",
                    "combined",
                    "-f",
                    water,
                    "-c",
                    "0",
                    "-m",
                    "1",
                    "pka",
                    "--program",
                    "gaussian",
                    "-s",
                    "direct",
                    "-pi",
                    "2",
                ],
                standalone_mode=False,
            )
        assert result.exit_code == 0, result.output
        assert isinstance(result.return_value, GaussianpKaJob)
        gaussian_from_project.assert_called_with("gas_solv")
        assert result.return_value.skip_completed is True

    def test_pka_chain_level_rerun_sets_skip_completed_false(
        self, isolated_config_dir
    ):
        _combined_yaml(isolated_config_dir)
        water = _water_xyz(isolated_config_dir)
        result = _invoke_chain(
            [
                "-R",
                "-p",
                "combined",
                "-f",
                water,
                "-c",
                "0",
                "-m",
                "1",
                "pka",
                "--program",
                "gaussian",
                "-s",
                "direct",
                "-pi",
                "2",
            ],
            standalone_mode=False,
        )
        assert result.exit_code == 0, result.output
        assert isinstance(result.return_value, GaussianpKaJob)
        assert result.return_value.skip_completed is False

    def test_pka_workflow_level_rerun_sets_skip_completed_false(
        self, isolated_config_dir
    ):
        _combined_yaml(isolated_config_dir)
        water = _water_xyz(isolated_config_dir)
        result = _invoke_chain(
            [
                "-p",
                "combined",
                "-f",
                water,
                "-c",
                "0",
                "-m",
                "1",
                "pka",
                "-R",
                "--program",
                "gaussian",
                "-s",
                "direct",
                "-pi",
                "2",
            ],
            standalone_mode=False,
        )
        assert result.exit_code == 0, result.output
        assert isinstance(result.return_value, GaussianpKaJob)
        assert result.return_value.skip_completed is False

    def test_pka_submit_requires_program(
        self, isolated_config_dir, single_molecule_xyz_file
    ):
        _combined_yaml(isolated_config_dir)
        result = CliRunner().invoke(
            chain,
            [
                "-p",
                "combined",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "pka",
                "-s",
                "direct",
                "-pi",
                "1",
            ],
            obj={"jobrunner": MagicMock()},
        )
        assert result.exit_code == 2, result.output
        assert "--program is required for submit" in result.output

    def test_pka_analyze_works_without_project(self):
        result = CliRunner().invoke(
            chain,
            ["pka", "analyze"],
            obj={"jobrunner": MagicMock()},
        )
        assert result.exit_code == 2, result.output
        assert "all 8 output files are required" in result.output
        assert "-p/--project" not in result.output

    def test_steps_cannot_combine_with_pka(
        self, isolated_config_dir, single_molecule_xyz_file
    ):
        _combined_yaml(isolated_config_dir)
        result = CliRunner().invoke(
            chain,
            [
                "-p",
                "combined",
                "-f",
                single_molecule_xyz_file,
                "-s",
                "gaussian:opt",
                "pka",
                "--program",
                "gaussian",
            ],
            obj={"jobrunner": MagicMock()},
        )
        assert result.exit_code == 2, result.output
        assert "Cannot combine -s/--steps" in result.output

    def test_fukui_program_orca_builds_orca_job(
        self, isolated_config_dir, single_molecule_xyz_file
    ):
        _combined_yaml(isolated_config_dir)
        result = _invoke_chain(
            [
                "-p",
                "combined",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "fukui",
                "--program",
                "orca",
            ],
            standalone_mode=False,
        )
        assert result.exit_code == 0, result.output
        assert isinstance(result.return_value, ORCAFukuiJob)

    def test_fukui_program_gaussian_builds_gaussian_job(
        self, isolated_config_dir, single_molecule_xyz_file
    ):
        _combined_yaml(isolated_config_dir)
        result = _invoke_chain(
            [
                "-p",
                "combined",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "fukui",
                "--program",
                "gaussian",
            ],
            standalone_mode=False,
        )
        assert result.exit_code == 0, result.output
        assert isinstance(result.return_value, GaussianFukuiJob)

    def test_fukui_analyze_works_without_project(self):
        result = CliRunner().invoke(
            chain,
            ["fukui", "analyze"],
            obj={"jobrunner": MagicMock()},
        )
        assert result.exit_code == 2, result.output
        assert "Missing option" in result.output
        assert "--neutral-filename" in result.output or "-n" in result.output
        assert "-p/--project" not in result.output

    def test_reaction_program_gaussian_builds_gaussian_job(
        self, isolated_config_dir, single_molecule_xyz_file
    ):
        _combined_yaml(isolated_config_dir)
        with patch(
            "chemsmart.settings.gaussian.GaussianProjectSettings.from_project",
            wraps=GaussianProjectSettings.from_project,
        ) as gaussian_from_project:
            result = _invoke_chain(
                [
                    "-p",
                    "combined",
                    "-f",
                    single_molecule_xyz_file,
                    "-c",
                    "0",
                    "-m",
                    "1",
                    "reaction",
                    "--program",
                    "gaussian",
                ],
                standalone_mode=False,
            )
        assert result.exit_code == 0, result.output
        assert isinstance(result.return_value, GaussianReactionJob)
        gaussian_from_project.assert_called_with("gas_solv")

    def test_reaction_program_orca_builds_orca_job(
        self, isolated_config_dir, single_molecule_xyz_file
    ):
        _combined_yaml(isolated_config_dir)
        result = _invoke_chain(
            [
                "-p",
                "combined",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "reaction",
                "--program",
                "orca",
            ],
            standalone_mode=False,
        )
        assert result.exit_code == 0, result.output
        assert isinstance(result.return_value, ORCAReactionJob)

    def test_reaction_has_no_analyze_command(self):
        result = CliRunner().invoke(
            chain,
            ["reaction", "--help"],
            obj={"jobrunner": MagicMock()},
        )
        assert result.exit_code == 0, result.output
        assert "\n  submit" in result.output
        assert "\n  batch" in result.output
        assert "\n  analyze" not in result.output

    def test_sub_test_chain_reaction_product_uses_cli_flag(
        self, isolated_config_dir, tmp_path, server_yaml_file, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        _combined_yaml(isolated_config_dir)
        reactant = tmp_path / "r.xyz"
        product = tmp_path / "p.xyz"
        reactant.write_text("2\nh2\nH 0 0 0\nH 0 0 0.74\n")
        product.write_text("2\nh2\nH 0 0 0\nH 0 0 0.90\n")
        result = CliRunner().invoke(
            sub,
            [
                "-s",
                server_yaml_file,
                "--test",
                "chain",
                "-p",
                "combined",
                "-f",
                str(reactant),
                "-c",
                "0",
                "-m",
                "1",
                "reaction",
                "--program",
                "gaussian",
                "--product",
                str(product),
            ],
            catch_exceptions=False,
        )
        assert result.exit_code == 0, result.output
        run_script = tmp_path / "chemsmart_run_r_reaction.py"
        assert run_script.exists(), sorted(tmp_path.glob("chemsmart_run_*.py"))
        text = run_script.read_text()
        assert "'--product'" in text
        assert "'--products'" not in text

    def test_redox_submit_requires_program(
        self, isolated_config_dir, single_molecule_xyz_file
    ):
        _combined_yaml(isolated_config_dir)
        result = CliRunner().invoke(
            chain,
            [
                "-p",
                "combined",
                "-f",
                single_molecule_xyz_file,
                "redox",
            ],
            obj={"jobrunner": MagicMock()},
        )
        assert result.exit_code == 2, result.output
        assert "--program is required for submit" in result.output

    def test_redox_program_gaussian_builds_gaussian_job(
        self, isolated_config_dir, single_molecule_xyz_file
    ):
        _combined_yaml(isolated_config_dir)
        ref_ox, ref_red = _ref_xyz_pair(isolated_config_dir)
        with patch(
            "chemsmart.settings.gaussian.GaussianProjectSettings.from_project",
            wraps=GaussianProjectSettings.from_project,
        ) as gaussian_from_project:
            result = _invoke_chain(
                [
                    "-p",
                    "combined",
                    "-f",
                    single_molecule_xyz_file,
                    "-c",
                    "1",
                    "-m",
                    "2",
                    "redox",
                    "--program",
                    "gaussian",
                    "--ref-ox",
                    ref_ox,
                    "--ref-red",
                    ref_red,
                ],
                standalone_mode=False,
            )
        assert result.exit_code == 0, result.output
        assert isinstance(result.return_value, GaussianRedoxJob)
        gaussian_from_project.assert_called_with("gas_solv")

    def test_redox_program_orca_builds_orca_job(
        self, isolated_config_dir, single_molecule_xyz_file
    ):
        _combined_yaml(isolated_config_dir)
        ref_ox, ref_red = _ref_xyz_pair(isolated_config_dir)
        result = _invoke_chain(
            [
                "-p",
                "combined",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "1",
                "-m",
                "2",
                "redox",
                "--program",
                "orca",
                "--ref-ox",
                ref_ox,
                "--ref-red",
                ref_red,
            ],
            standalone_mode=False,
        )
        assert result.exit_code == 0, result.output
        assert isinstance(result.return_value, ORCARedoxJob)

    def test_redox_analyze_is_registered_without_project(self):
        result = CliRunner().invoke(
            chain,
            ["redox", "--help"],
            obj={"jobrunner": MagicMock()},
        )
        assert result.exit_code == 0, result.output
        assert "\n  analyze" in result.output

        analyze_help = CliRunner().invoke(
            chain,
            ["redox", "analyze", "--help"],
            obj={"jobrunner": MagicMock()},
        )
        assert analyze_help.exit_code == 0, analyze_help.output
        assert "--ox-gas" in analyze_help.output
        assert "--ref-red-solv" in analyze_help.output
        assert "--e-ref" in analyze_help.output

    def test_redox_analyze_works_without_project(self):
        result = CliRunner().invoke(
            chain,
            ["redox", "analyze"],
            obj={"jobrunner": MagicMock()},
        )
        assert result.exit_code == 2, result.output
        assert "Missing option" in result.output
        assert "--e-ref" in result.output or "--ox-gas" in result.output
        assert "-p/--project" not in result.output

    @pytest.mark.parametrize(
        "step,snippet",
        [
            ("gaussian:fukui", "do not support gaussian 'fukui'"),
            ("gaussian:redox", "do not support gaussian 'redox'"),
            ("gaussian:reaction", "do not support gaussian 'reaction'"),
            ("orca:pka", "do not support orca 'pka'"),
            ("orca:fukui", "do not support orca 'fukui'"),
            ("orca:redox", "do not support orca 'redox'"),
            ("orca:reaction", "do not support orca 'reaction'"),
        ],
    )
    def test_nested_only_workflow_steps_are_usage_errors(
        self, isolated_config_dir, single_molecule_xyz_file, step, snippet
    ):
        _combined_yaml(isolated_config_dir)
        result = CliRunner().invoke(
            chain,
            [
                "-p",
                "combined",
                "-f",
                single_molecule_xyz_file,
                "-s",
                step,
            ],
            obj={"jobrunner": MagicMock()},
        )
        assert result.exit_code == 2, result.output
        assert snippet in result.output

    def test_pka_analyze_computes_without_project(self, tmp_path, monkeypatch):
        files = {
            name: str(tmp_path / name)
            for name in (
                "ha.log",
                "a.log",
                "hb.log",
                "b.log",
                "has.log",
                "as.log",
                "hbs.log",
                "bs.log",
            )
        }
        for path in files.values():
            Path(path).write_text("Gaussian, Inc.\n")
        called = {}

        def _fake_print(*args, **kwargs):
            called["kwargs"] = kwargs

        monkeypatch.setattr("chemsmart.cli.pka.print_pka_summary", _fake_print)
        result = CliRunner().invoke(
            chain,
            [
                "pka",
                "analyze",
                "-ha",
                files["ha.log"],
                "-a",
                files["a.log"],
                "-hr",
                files["hb.log"],
                "-r",
                files["b.log"],
                "-has",
                files["has.log"],
                "-as",
                files["as.log"],
                "--href-solv",
                files["hbs.log"],
                "--ref-solv",
                files["bs.log"],
                "-rp",
                "6.75",
            ],
            obj={"jobrunner": MagicMock()},
        )
        assert result.exit_code == 0, result.output
        assert called["kwargs"]["ha_gas_file"] == files["ha.log"]
        assert called["kwargs"]["pka_reference"] == 6.75
        assert "-p/--project" not in result.output

    def test_fukui_analyze_computes_without_project(
        self, tmp_path, monkeypatch
    ):
        files = {
            "n": tmp_path / "mol_n.log",
            "rc": tmp_path / "mol_rc.log",
            "ra": tmp_path / "mol_ra.log",
        }
        for path in files.values():
            path.write_text("n")
        called = {}

        def _fake_analyze(**kwargs):
            called.update(kwargs)

        monkeypatch.setattr("chemsmart.cli.fukui.analyze_fukui", _fake_analyze)
        result = CliRunner().invoke(
            chain,
            [
                "fukui",
                "analyze",
                "-n",
                str(files["n"]),
                "-c",
                str(files["rc"]),
                "-a",
                str(files["ra"]),
            ],
            obj={"jobrunner": MagicMock()},
        )
        assert result.exit_code == 0, result.output
        assert called["neutral_filename"] == str(files["n"])
        assert called["radical_cation_filename"] == str(files["rc"])
        assert "-p/--project" not in result.output


class TestExistingProgramCLI:
    def test_gaussian_pka_still_uses_gaussian_project(self, tmp_path):
        water = tmp_path / "water.xyz"
        water.write_text(
            "3\nwater\n"
            "O  0.000  0.000  0.000\n"
            "H  0.957  0.000  0.000\n"
            "H -0.240  0.927  0.000\n"
        )
        with patch(
            "chemsmart.settings.gaussian.GaussianProjectSettings.from_project",
            wraps=GaussianProjectSettings.from_project,
        ) as gaussian_from_project:
            result = CliRunner().invoke(
                gaussian,
                [
                    "-p",
                    "test",
                    "-f",
                    str(water),
                    "-c",
                    "0",
                    "-m",
                    "1",
                    "pka",
                    "-s",
                    "direct",
                    "-pi",
                    "2",
                ],
                obj={"jobrunner": MagicMock()},
                catch_exceptions=False,
                standalone_mode=False,
            )
        assert result.exit_code == 0, result.output
        assert isinstance(result.return_value, GaussianpKaJob)
        gaussian_from_project.assert_called_with("test")

    def test_run_pka_analyze_still_works(self, tmp_path, monkeypatch):
        files = {
            name: tmp_path / name
            for name in (
                "ha.log",
                "a.log",
                "hb.log",
                "b.log",
                "has.log",
                "as.log",
                "hbs.log",
                "bs.log",
            )
        }
        for path in files.values():
            path.write_text("Gaussian, Inc.\n")
        called = {}

        def _fake_print(*args, **kwargs):
            called["kwargs"] = kwargs

        monkeypatch.setattr("chemsmart.cli.pka.print_pka_summary", _fake_print)
        result = CliRunner().invoke(
            run,
            [
                "pka",
                "analyze",
                "-ha",
                str(files["ha.log"]),
                "-a",
                str(files["a.log"]),
                "-hr",
                str(files["hb.log"]),
                "-r",
                str(files["b.log"]),
                "-has",
                str(files["has.log"]),
                "-as",
                str(files["as.log"]),
                "--href-solv",
                str(files["hbs.log"]),
                "--ref-solv",
                str(files["bs.log"]),
                "-rp",
                "6.75",
            ],
        )
        assert result.exit_code == 0, result.output
        assert called["kwargs"]["ha_gas_file"] == str(files["ha.log"])
        assert called["kwargs"]["pka_reference"] == 6.75

    def test_run_fukui_still_works(self, tmp_path, monkeypatch):
        files = {
            "n": tmp_path / "mol_n.log",
            "rc": tmp_path / "mol_rc.log",
            "ra": tmp_path / "mol_ra.log",
        }
        for path in files.values():
            path.write_text("n")
        called = {}

        def _fake_analyze(**kwargs):
            called.update(kwargs)

        monkeypatch.setattr("chemsmart.cli.fukui.analyze_fukui", _fake_analyze)
        result = CliRunner().invoke(
            run,
            [
                "fukui",
                "-n",
                str(files["n"]),
                "-c",
                str(files["rc"]),
                "-a",
                str(files["ra"]),
            ],
        )
        assert result.exit_code == 0, result.output
        assert called["neutral_filename"] == str(files["n"])
