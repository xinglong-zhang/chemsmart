import os
import shutil

import click
import numpy as np
import pytest
from click.testing import CliRunner

from chemsmart.cli.main import entry_point as main_entry_point
from chemsmart.cli.run import run as run_entry_point
from chemsmart.cli.sub import sub as sub_entry_point
from chemsmart.io.gaussian.input import Gaussian16Input
from chemsmart.jobs.gaussian.settings import GaussianJobSettings
from chemsmart.utils.cli import CtxObjArguments, MyCommand, MyGroup


@pytest.fixture
def cli_runner():
    return CliRunner()


def test_cli_entry_point():
    @click.group(cls=MyGroup, chain=True)
    @click.pass_context
    def test_group(ctx):
        pass

    @test_group.result_callback()
    @click.pass_context
    def group_callback(ctx, *args, **kwargs):
        commands = ctx.obj["subcommand"]
        args = CtxObjArguments(commands)
        ctx.obj["cmd_line"] = args.reconstruct_command_line()

    @click.command("command2", cls=MyCommand)
    def test_command1():
        pass

    @click.command("command1", cls=MyCommand)
    def test_command2():
        pass

    test_group.add_command(test_command1, name="addition")
    test_group.add_command(test_command2, name="subtraction")
    return test_group


class TestCLI:
    """Parent class with useful fixtures and methods for testing the CLI."""

    @pytest.fixture(autouse=True)
    def _prevent_click_io_error(self, tmpdir, monkeypatch, caplog):
        """Prevents click from raising an IOError when pytest logging is enabled.

        See https://github.com/pallets/click/issues/824
        """
        caplog.set_level(10000000)

    @pytest.fixture(autouse=True)
    def _cd_tmpdir(self, tmpdir, monkeypatch):
        """Changes the working directory to the temporary directory."""
        monkeypatch.chdir(tmpdir)
        yield
        monkeypatch.undo()

    @staticmethod
    def run_and_check(
        command, runner=CliRunner(), entry_point=None, obj=None, input=None
    ):
        """Runs a command and raise an error if the command fails.

        Returns:
            results (click.testing.Result): The results of the command.
        """
        if obj is None:
            obj = {}

        if entry_point is None:
            entry_point = main_entry_point

        results = runner.invoke(entry_point, command, obj=obj, input=input)

        if results.exception:
            raise results.exception
        return results


class TestGaussianCLI(TestCLI):
    def test_output_command_string_with_space_in_args(
        self, tmpdir, gaussian_singlet_opt_outfile, monkeypatch
    ):
        shutil.copy(src=gaussian_singlet_opt_outfile, dst=tmpdir)

        # Redirect the server config used by the CLI to a temporary copy
        from pathlib import Path

        from chemsmart.settings import server as server_module

        src_server = (
            Path(__file__).resolve().parent.parent
            / "chemsmart"
            / "settings"
            / "templates"
            / "server.yaml"
        )
        server_dir = Path(tmpdir) / "chemsmart_templates" / "server"
        server_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy(src=str(src_server), dst=str(server_dir / "PBS.yaml"))

        class _StubUserSettings:
            def __init__(self, server_path):
                self._server_path = str(server_path)

            @property
            def user_server_dir(self):
                return self._server_path

        monkeypatch.setattr(
            server_module, "user_settings", _StubUserSettings(server_dir)
        )

        command = (
            "-s PBS --test --print-command --no-scratch  gaussian -p test2 -f nhc_neutral_singlet.log "
            "-l nhc userjob -r 'opt freq mn15 def2svp' "
        )

        results = self.run_and_check(
            command, entry_point=sub_entry_point, obj={"verbose": False}
        )
        assert (
            "opt freq mn15 def2svp" in results.output
        )  # not broken up into indidual strings

    def test_gaussian_cli_sub(self, tmpdir, gaussian_singlet_opt_outfile):
        shutil.copy(src=gaussian_singlet_opt_outfile, dst=tmpdir)

        command = "-s SLURM --test --print-command --no-scratch gaussian -p test -f nhc_neutral_singlet.log opt"
        results = self.run_and_check(command, entry_point=sub_entry_point)
        assert "m062x" in results.output
        assert "def2svp" in results.output
        assert "smd" in results.output
        assert "dichloroethane" in results.output

    def test_gaussian_cli_run(self, tmpdir, gaussian_singlet_opt_outfile):
        shutil.copy(src=gaussian_singlet_opt_outfile, dst=tmpdir)

        command = "--no-scratch --fake gaussian -p test -f nhc_neutral_singlet.log opt"
        self.run_and_check(command, entry_point=run_entry_point)

        # input file created
        inputfile = os.path.join(tmpdir, "nhc_neutral_singlet_opt_fake.com")
        assert os.path.exists(inputfile)
        gaussian_input = Gaussian16Input(inputfile)
        assert (
            gaussian_input.route_string
            == "# opt freq m062x def2svp scrf=(smd,solvent=dichloroethane)"
        )

        settings = GaussianJobSettings.from_filepath(
            "nhc_neutral_singlet_opt_fake.com"
        )
        assert settings.charge == 0
        assert settings.multiplicity == 1

    def test_gaussian_cli_dieze_tag(
        self, tmpdir, gaussian_singlet_opt_outfile
    ):
        shutil.copy(src=gaussian_singlet_opt_outfile, dst=tmpdir)

        command = "--fake --no-scratch gaussian -p test -f nhc_neutral_singlet.log -d n opt"
        self.run_and_check(command, entry_point=run_entry_point)

        # input file created
        inputfile = os.path.join(tmpdir, "nhc_neutral_singlet_opt_fake.com")
        assert os.path.exists(inputfile)
        gaussian_input = Gaussian16Input(inputfile)
        assert (
            gaussian_input.route_string
            == "#n opt freq m062x def2svp scrf=(smd,solvent=dichloroethane)"
        )

    def test_gaussian_cli_charge_multiplicity_append_name(
        self, tmpdir, gaussian_singlet_opt_outfile
    ):
        shutil.copy(src=gaussian_singlet_opt_outfile, dst=tmpdir)

        command = (
            "--fake --no-scratch gaussian -p test -f nhc_neutral_singlet.log -d t -c 1 -m 2 -a charge_multiplicity "
            "scan -c [[1,2]] -n 10 -s 0.2"
        )
        self.run_and_check(command, entry_point=run_entry_point)

        # input file created
        inputfile = os.path.join(
            tmpdir, "nhc_neutral_singlet_charge_multiplicity_fake.com"
        )
        assert os.path.exists(inputfile)
        gaussian_input = Gaussian16Input(inputfile)
        assert (
            gaussian_input.route_string
            == "#t opt=modredundant m062x def2svp scrf=(smd,solvent=dichloroethane)"
        )
        assert gaussian_input.content_groups[-1] == ["B 1 2 S 10 0.2"]

        settings = GaussianJobSettings.from_filepath(inputfile)
        assert settings.charge == 1
        assert settings.multiplicity == 2

    def test_gaussian_cli_modify_settings(
        self, tmpdir, gaussian_singlet_opt_outfile
    ):
        shutil.copy(src=gaussian_singlet_opt_outfile, dst=tmpdir)

        command = (
            '--no-scratch --fake gaussian -p test -f nhc_neutral_singlet.log -x b3lyp -b "6-31G*" -m 3 '
            "-l completely_new_filename -o maxstep=8,maxsize=12 -r nosymm modred -c [[2,3,6],[6,7]] "
        )
        self.run_and_check(command, entry_point=run_entry_point)

        # input file created
        inputfile = os.path.join(tmpdir, "completely_new_filename_fake.com")
        assert os.path.exists(inputfile)
        gaussian_input = Gaussian16Input(inputfile)
        assert (
            gaussian_input.route_string
            == "# opt=(modredundant,maxstep=8,maxsize=12) freq b3lyp 6-31g* scrf=(smd,solvent=dichloroethane) nosymm"
        )  # modred job
        assert gaussian_input.content_groups[-1] == ["A 2 3 6 F", "B 6 7 F"]

        settings = GaussianJobSettings.from_filepath(inputfile)
        assert settings.charge == 0
        assert settings.multiplicity == 3

    def test_gaussian_sp_job_solvent_settings(
        self, tmpdir, gaussian_singlet_opt_outfile
    ):
        shutil.copy(src=gaussian_singlet_opt_outfile, dst=tmpdir)

        command = "--no-scratch gaussian -p test -f nhc_neutral_singlet.log sp"
        self.run_and_check(command, entry_point=run_entry_point)

        inputfile = os.path.join(
            tmpdir, "nhc_neutral_singlet_sp_smd_dichloroethane.com"
        )
        assert os.path.exists(inputfile)
        gaussian_input = Gaussian16Input(inputfile)
        assert (
            gaussian_input.route_string
            == "# m062x def2tzvp scrf=(smd,solvent=dichloroethane)"
        )
        assert len(gaussian_input.append_info_after_coords_as_string) == 0

        settings = GaussianJobSettings.from_filepath(inputfile)
        assert settings.charge == 0
        assert settings.multiplicity == 1
        assert settings.freq is False
        assert settings.solvent_model == "smd"
        assert settings.solvent_id == "dichloroethane"

    def test_gaussian_sp_job_custom_solvent_settings(
        self, tmpdir, gaussian_singlet_opt_outfile
    ):
        shutil.copy(src=gaussian_singlet_opt_outfile, dst=tmpdir)

        command = (
            "--no-scratch gaussian -p test3 -f nhc_neutral_singlet.log sp"
        )
        self.run_and_check(command, entry_point=run_entry_point)

        # input file created
        inputfile = os.path.join(
            tmpdir, "nhc_neutral_singlet_sp_smd_generic.com"
        )
        assert os.path.exists(inputfile)
        gaussian_input = Gaussian16Input(inputfile)
        assert (
            gaussian_input.route_string
            == "# mn15 def2qzvp scrf=(smd,solvent=generic,read)"
        )

        settings = GaussianJobSettings.from_filepath(inputfile)
        assert settings.charge == 0
        assert settings.multiplicity == 1
        assert settings.freq is False
        assert settings.solvent_model == "smd"
        assert settings.solvent_id == "generic,read"

    def test_gaussian_sp_job_modify_solvent_from_project_solvent(
        self, tmpdir, gaussian_singlet_opt_outfile
    ):
        shutil.copy(src=gaussian_singlet_opt_outfile, dst=tmpdir)

        command = "--no-scratch gaussian -p test -f nhc_neutral_singlet.log sp -sm cpcm -si chloroform"
        self.run_and_check(command, entry_point=run_entry_point)

        # input file created
        inputfile = os.path.join(
            tmpdir, "nhc_neutral_singlet_sp_cpcm_chloroform.com"
        )  # note that file name is also changed from test above
        assert os.path.exists(inputfile)
        gaussian_input = Gaussian16Input(inputfile)
        assert (
            gaussian_input.route_string
            == "# m062x def2tzvp scrf=(cpcm,solvent=chloroform)"
        )

        settings = GaussianJobSettings.from_filepath(inputfile)
        assert settings.charge == 0
        assert settings.multiplicity == 1
        assert settings.freq is False
        assert settings.solvent_model == "cpcm"
        assert settings.solvent_id == "chloroform"

    def test_gaussian_sp_job_modify_one_param_from_project_solvent(
        self, tmpdir, gaussian_singlet_opt_outfile
    ):
        shutil.copy(src=gaussian_singlet_opt_outfile, dst=tmpdir)

        command = "--no-scratch gaussian -p test -f nhc_neutral_singlet.log sp -sm cpcm"
        self.run_and_check(command, entry_point=run_entry_point)

        # input file created
        inputfile = os.path.join(
            tmpdir, "nhc_neutral_singlet_sp_cpcm_dichloroethane.com"
        )
        assert os.path.exists(inputfile)
        gaussian_input = Gaussian16Input(inputfile)
        assert (
            gaussian_input.route_string
            == "# m062x def2tzvp scrf=(cpcm,solvent=dichloroethane)"
        )
        settings = GaussianJobSettings.from_filepath(inputfile)
        assert settings.charge == 0
        assert settings.multiplicity == 1
        assert settings.freq is False
        assert settings.solvent_model == "cpcm"
        assert settings.solvent_id == "dichloroethane"

    def test_gaussian_sp_job_add_solvent_from_project_no_solvent(
        self, tmpdir, gaussian_singlet_opt_outfile
    ):
        shutil.copy(src=gaussian_singlet_opt_outfile, dst=tmpdir)

        command = " --no-scratch gaussian -p test4 -f nhc_neutral_singlet.log sp -sm dipole -si water"
        self.run_and_check(command, entry_point=run_entry_point)

        # input file created
        inputfile = os.path.join(
            tmpdir, "nhc_neutral_singlet_sp_dipole_water.com"
        )  # note that file name is also changed from test above
        assert os.path.exists(inputfile)
        gaussian_input = Gaussian16Input(inputfile)
        assert (
            gaussian_input.route_string
            == "# mn15 def2qzvp scrf=(dipole,solvent=water)"
        )

        settings = GaussianJobSettings.from_filepath(inputfile)
        assert settings.charge == 0
        assert settings.multiplicity == 1
        assert settings.freq is False
        assert settings.solvent_model == "dipole"
        assert settings.solvent_id == "water"

    @pytest.mark.parametrize(
        "remove_solvent_option", ["--remove-solvent", "-r"]
    )
    def test_gaussian_sp_job_turn_off_solvent(
        self, tmpdir, gaussian_singlet_opt_outfile, remove_solvent_option
    ):
        shutil.copy(src=gaussian_singlet_opt_outfile, dst=tmpdir)

        command = f"--no-scratch gaussian -p test -f nhc_neutral_singlet.log sp {remove_solvent_option}"
        self.run_and_check(command, entry_point=run_entry_point)

        inputfile = os.path.join(
            tmpdir, "nhc_neutral_singlet_sp_gas_phase.com"
        )  # Note: file name is also changed from test above
        assert os.path.exists(inputfile)
        gaussian_input = Gaussian16Input(inputfile)
        assert gaussian_input.route_string == "# m062x def2tzvp"

    def test_gaussian_job_from_genecp(
        self, tmpdir, gaussian_ts_genecp_outfile
    ):
        shutil.copy(src=gaussian_ts_genecp_outfile, dst=tmpdir)

        command = " --no-scratch gaussian -p test4 -f pd_genecp_ts.log opt"
        self.run_and_check(command, entry_point=run_entry_point)

        # input file created
        inputfile = os.path.join(tmpdir, "pd_genecp_ts_opt.com")
        assert os.path.exists(inputfile)
        gaussian_input = Gaussian16Input(inputfile)
        assert gaussian_input.route_string == "# opt freq mn15 def2svp"

    def test_gaussian_job_from_premature_termination_with_opt_points(
        self, tmpdir, gaussian_failed_modred_outfile
    ):
        shutil.copy(src=gaussian_failed_modred_outfile, dst=tmpdir)

        command = "--no-scratch gaussian -p test4 -f cage_free_failed_modred.log -x m062x -l cage_free_from_failed_modred_opt opt"
        self.run_and_check(command, entry_point=run_entry_point)

        # input file created
        inputfile = os.path.join(
            tmpdir, "cage_free_from_failed_modred_opt.com"
        )
        assert os.path.exists(inputfile)
        gaussian_input = Gaussian16Input(inputfile)
        assert gaussian_input.route_string == "# opt freq m062x def2svp"
        mol = gaussian_input.molecule

        # last structure from prematurely terminated outputfile (last point is not excluded herein)
        assert np.allclose(
            mol.positions[0],
            np.array([-2.04970900, -0.92164400, 0.37502900]),
            atol=1e-4,
        )

        command2 = "--no-scratch gaussian -p test4 -f cage_free_failed_modred.log -x m062x -i -2 -l cage_free_from_failed_modred2_opt opt"
        self.run_and_check(command2, entry_point=run_entry_point)

        # input file created
        inputfile = os.path.join(
            tmpdir, "cage_free_from_failed_modred2_opt.com"
        )
        assert os.path.exists(inputfile)
        gaussian_input = Gaussian16Input(inputfile)
        assert gaussian_input.route_string == "# opt freq m062x def2svp"
        mol = gaussian_input.molecule

        # use second last structure from prematurely terminated outputfile (via -i -2 in cli)
        assert np.allclose(
            mol.positions[0],
            np.array([-2.06360000, -0.93860900, 0.37364500]),
            atol=1e-4,
        )

    def test_gaussian_freeze_atoms_job(
        self, tmpdir, gaussian_singlet_opt_outfile
    ):
        shutil.copy(src=gaussian_singlet_opt_outfile, dst=tmpdir)

        command = "--no-scratch gaussian -p test4 -f nhc_neutral_singlet.log opt -f 1-10"  # 1-indexed labels as on GaussView
        self.run_and_check(command, entry_point=run_entry_point)

        # input file created
        inputfile = os.path.join(tmpdir, "nhc_neutral_singlet_opt.com")
        assert os.path.exists(inputfile)
        gaussian_input = Gaussian16Input(inputfile)
        assert gaussian_input.route_string == "# opt freq mn15 def2svp"
        assert gaussian_input.has_frozen_coordinates
        assert gaussian_input.frozen_coordinate_indices == list(range(1, 11))
        inputlines = gaussian_input.contents
        assert (
            inputlines[8]
            == "C         -1    0.9200700000    0.5470510000   -0.9928680000"
        )  # first line of coords
        assert (
            inputlines[-2]
            == "F          0   -5.8786510000   -1.1600890000    1.1374510000"
        )  # last line of coords

        assert gaussian_input.molecule.chemical_formula == "C19H12F3I2N3O"
        settings = GaussianJobSettings.from_filepath(
            "nhc_neutral_singlet_opt.com"
        )
        assert settings.charge == 0
        assert settings.multiplicity == 1

    def test_qmmm_job(self, tmpdir, gaussian_singlet_opt_outfile):
        shutil.copy(src=gaussian_singlet_opt_outfile, dst=tmpdir)

        command = "--no-scratch --fake gaussian -p test4 -f nhc_neutral_singlet.log qmmm"  # 1-indexed labels as on GaussView
        self.run_and_check(command, entry_point=run_entry_point)
        # input file created
        inputfile = os.path.join(tmpdir, "nhc_neutral_singlet_qmmm_fake.com")
        print(inputfile)
        assert os.path.exists(inputfile)
