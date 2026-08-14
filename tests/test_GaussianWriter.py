import os
from filecmp import cmp
from shutil import copy

import pytest

from chemsmart.io.gaussian.output import Gaussian16Output
from chemsmart.io.molecules.structure import Molecule, QMMMMolecule
from chemsmart.jobs.gaussian import (
    GaussianModredJob,
    GaussianOptJob,
    GaussianQMMMJob,
    GaussianScanJob,
    GaussianSinglePointJob,
    GaussianTSJob,
)
from chemsmart.jobs.gaussian.settings import (
    GaussianJobSettings,
    GaussianQMMMJobSettings,
)
from chemsmart.jobs.gaussian.writer import GaussianInputWriter
from chemsmart.settings.gaussian import GaussianProjectSettings
from chemsmart.utils.utils import cmp_with_ignore


class TestGaussianInputWriter:
    def test_write_opt_job(
        self,
        tmpdir,
        single_molecule_xyz_file,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_jobrunner_no_scratch,
        gaussian_written_opt_file,
    ):
        # get project settings
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        settings = project_settings.opt_settings()
        settings.charge = 0
        settings.multiplicity = 1
        job = GaussianOptJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="gaussian_opt",
            jobrunner=gaussian_jobrunner_no_scratch,
        )

        assert isinstance(job, GaussianOptJob)
        g16_writer = GaussianInputWriter(job=job)

        # write input file
        g16_writer.write(target_directory=tmpdir)
        g16_file = os.path.join(tmpdir, "gaussian_opt.com")
        assert os.path.isfile(g16_file)
        assert cmp(
            g16_file, gaussian_written_opt_file, shallow=False
        )  # writes input file as expected

        # job run will result in the job being run and
        # the output file copied back to run folder
        # job.run()
        # assert job.is_complete()

    def test_write_semiempirical_opt_job(
        self,
        tmpdir,
        single_molecule_xyz_file,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_jobrunner_no_scratch,
        gaussian_written_pm6_opt_file,
    ):
        # get project settings
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        settings = project_settings.opt_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.semiempirical = "PM6"
        job = GaussianOptJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="gaussian_pm6_opt",
            jobrunner=gaussian_jobrunner_no_scratch,
        )

        assert isinstance(job, GaussianOptJob)
        g16_writer = GaussianInputWriter(job=job)

        # write input file
        g16_writer.write(target_directory=tmpdir)
        g16_file = os.path.join(tmpdir, "gaussian_pm6_opt.com")
        assert os.path.isfile(g16_file)
        assert cmp(
            g16_file, gaussian_written_pm6_opt_file, shallow=False
        )  # writes input file as expected

        # job run will result in the job being run and
        # the output file copied back to run folder
        # job.run()
        # assert job.is_complete()

    def test_write_opt_job_with_route(
        self,
        tmpdir,
        single_molecule_xyz_file,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_jobrunner_no_scratch,
        gaussian_written_opt_file_with_route,
    ):
        # get project settings
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        settings = project_settings.opt_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.route_to_be_written = "#p pbepbe/6-31g(d) opt"
        settings.title = "Optimisation job with supplied route"
        job = GaussianOptJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="gaussian_opt",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert isinstance(job, GaussianOptJob)
        g16_writer = GaussianInputWriter(job=job)

        # write input file
        g16_writer.write(target_directory=tmpdir)
        g16_file = os.path.join(tmpdir, "gaussian_opt.com")
        assert os.path.isfile(g16_file)
        assert cmp(
            g16_file, gaussian_written_opt_file_with_route, shallow=False
        )

    def test_write_modred_job(
        self,
        tmpdir,
        single_molecule_xyz_file,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_jobrunner_no_scratch,
        gaussian_written_modred_file,
    ):
        # get project settings
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        settings = project_settings.modred_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.modred = [[1, 2], [3, 4, 5]]
        job = GaussianModredJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="gaussian_modred",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert isinstance(job, GaussianModredJob)
        g16_writer = GaussianInputWriter(job=job)

        # write input file
        g16_writer.write(target_directory=tmpdir)
        g16_file = os.path.join(tmpdir, "gaussian_modred.com")
        assert os.path.isfile(g16_file)
        assert cmp(g16_file, gaussian_written_modred_file, shallow=False)

    def test_write_scan_job_multiple_degrees_of_freedom(
        self,
        tmpdir,
        single_molecule_xyz_file,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_jobrunner_no_scratch,
        gaussian_written_scan_multiple_degrees_of_freedom_file,
    ):
        # get project settings
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        settings = project_settings.scan_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.modred = {
            "coords": [[1, 2], [3, 4, 5]],
            "num_steps": [10, 18],
            "step_size": [0.1, 5.0],
        }
        job = GaussianScanJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="gaussian_scan_multiple_degrees_of_freedom",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert isinstance(job, GaussianScanJob)
        g16_writer = GaussianInputWriter(job=job)

        # write input file
        g16_writer.write(target_directory=tmpdir)
        g16_file = os.path.join(
            tmpdir, "gaussian_scan_multiple_degrees_of_freedom.com"
        )
        assert os.path.isfile(g16_file)
        assert cmp(
            g16_file,
            gaussian_written_scan_multiple_degrees_of_freedom_file,
            shallow=False,
        )

    def test_write_scan_job_single_degree_of_freedom(
        self,
        tmpdir,
        single_molecule_xyz_file,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_jobrunner_no_scratch,
        gaussian_written_scan_single_degree_of_freedom_file,
    ):
        # get project settings
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        settings = project_settings.scan_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.modred = {
            "coords": [1, 2],
            "num_steps": [10],
            "step_size": [0.1],
        }
        job = GaussianScanJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="gaussian_scan_single_degree_of_freedom",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert isinstance(job, GaussianScanJob)
        g16_writer = GaussianInputWriter(job=job)

        # write input file
        g16_writer.write(target_directory=tmpdir)
        g16_file = os.path.join(
            tmpdir, "gaussian_scan_single_degree_of_freedom.com"
        )
        assert os.path.isfile(g16_file)
        assert cmp(
            g16_file,
            gaussian_written_scan_single_degree_of_freedom_file,
            shallow=False,
        )

    def test_write_scan_job_multiple_degrees_of_freedom_with_constraints(
        self,
        tmpdir,
        single_molecule_xyz_file,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_jobrunner_no_scratch,
        gaussian_written_scan_multiple_degrees_of_freedom_with_constraints_file,
    ):
        # get project settings
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        settings = project_settings.scan_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.modred = {
            "coords": [
                [1, 2],
                [3, 4],
                [5, 6, 7],
                [8, 9, 10],
                [11, 12, 13, 14],
            ],
            "num_steps": [10, 12, 18, 15, 36],
            "step_size": [0.1, 0.2, 3.0, 4.0, 10.0],
        }
        job = GaussianScanJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="gaussian_scan_multiple_degrees_of_freedom_with_constraints",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert isinstance(job, GaussianScanJob)
        g16_writer = GaussianInputWriter(job=job)

        # write input file
        g16_writer.write(target_directory=tmpdir)
        g16_file = os.path.join(
            tmpdir,
            "gaussian_scan_multiple_degrees_of_freedom_with_constraints.com",
        )
        assert os.path.isfile(g16_file)
        assert cmp(
            g16_file,
            gaussian_written_scan_multiple_degrees_of_freedom_with_constraints_file,
            shallow=False,
        )

    def test_write_ts_job(
        self,
        tmpdir,
        single_molecule_xyz_file,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_jobrunner_no_scratch,
        gaussian_written_ts_file,
    ):
        # get project settings
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        settings = project_settings.ts_settings()
        settings.charge = 0
        settings.multiplicity = 1
        job = GaussianTSJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="gaussian_ts",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert isinstance(job, GaussianTSJob)
        g16_writer = GaussianInputWriter(job=job)

        # write input file
        g16_writer.write(target_directory=tmpdir)
        g16_file = os.path.join(tmpdir, "gaussian_ts.com")
        assert os.path.isfile(g16_file)
        assert cmp(g16_file, gaussian_written_ts_file, shallow=False)

    def test_write_qmmm_job(
        self,
        tmpdir,
        single_molecule_xyz_file,
        gaussian_yaml_settings_qmmm_project_name,
        gaussian_written_qmmm_file,
        gaussian_jobrunner_no_scratch,
    ):
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_qmmm_project_name
        )
        qmmm_settings = project_settings.qmmm_settings()
        qmmm_settings.charge = 0
        qmmm_settings.multiplicity = 1
        qmmm_settings.charge_total = 0
        qmmm_settings.mult_total = 1
        qmmm_settings.high_level_atoms = [1, 2, 3]
        job = GaussianQMMMJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=qmmm_settings,
            label="gaussian_qmmm",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert isinstance(job, GaussianQMMMJob)
        g16_writer = GaussianInputWriter(job=job)
        g16_writer.write(target_directory=tmpdir)
        g16_file = os.path.join(tmpdir, "gaussian_qmmm.com")
        assert qmmm_settings.high_level_atoms == [1, 2, 3]
        with open(g16_file, "r") as f:
            g16_content = f.read()
            # check that the qmmm section is present in the written file
        with open(gaussian_written_qmmm_file, "r") as f:
            expected_content = f.read()
        assert g16_content == expected_content
        assert os.path.isfile(g16_file)
        assert cmp(g16_file, gaussian_written_qmmm_file, shallow=False)

    def test_write_qmmm_job_with_scale_factors(
        self,
        tmpdir,
        single_molecule_xyz_file,
        gaussian_yaml_settings_qmmm_project_name,
        gaussian_jobrunner_no_scratch,
    ):
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_qmmm_project_name
        )
        qmmm_settings = project_settings.qmmm_settings()
        qmmm_settings.charge = 0
        qmmm_settings.multiplicity = 1
        qmmm_settings.charge_total = 0
        qmmm_settings.mult_total = 1
        qmmm_settings.high_level_atoms = [1, 2, 3]
        qmmm_settings.bonded_atoms = [(3, 4)]
        qmmm_settings.scale_factors = {(3, 4): [0.709]}
        job = GaussianQMMMJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=qmmm_settings,
            label="gaussian_qmmm_scale",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        g16_writer = GaussianInputWriter(job=job)
        g16_writer.write(target_directory=tmpdir)
        g16_file = os.path.join(tmpdir, "gaussian_qmmm_scale.com")
        with open(g16_file, "r") as f:
            content = f.read()
        assert "H 3" in content
        assert "0.709" in content
        assert "geom=connectivity" in content.lower()

    def test_write_qmmm_link_atoms_share_high_boundary(
        self,
        tmpdir,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
    ):
        """Both MM boundary atoms may link to the same high-layer atom."""
        from chemsmart.jobs.gaussian.settings import GaussianQMMMJobSettings

        qmmm_settings = GaussianQMMMJobSettings(
            parent_jobtype="opt",
            jobtype="opt",
            freq=False,
            high_level_functional="mn15",
            high_level_basis="def2svp",
            low_level_force_field="UFF",
            charge_total=0,
            mult_total=1,
            charge_high=0,
            mult_high=1,
            high_level_atoms=[1, 2, 3],
            bonded_atoms=[(3, 4), (3, 5)],
        )
        job = GaussianQMMMJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=qmmm_settings,
            label="gaussian_qmmm_shared_link",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        GaussianInputWriter(job=job).write(target_directory=tmpdir)
        with open(
            os.path.join(tmpdir, "gaussian_qmmm_shared_link.com")
        ) as handle:
            content = handle.read()
        assert content.count(" L H 3") == 2

    def test_write_qmmm_auto_assigns_link_atoms_for_cut_bonds(
        self,
        tmpdir,
        gaussian_jobrunner_no_scratch,
    ):
        """Cut covalent bonds are assigned as link atoms when bonded_atoms is omitted."""
        import numpy as np

        from chemsmart.io.molecules.structure import Molecule
        from chemsmart.jobs.gaussian.settings import GaussianQMMMJobSettings

        molecule = Molecule(
            symbols=["C", "H", "H", "H", "C", "H", "H", "H"],
            positions=np.array(
                [
                    [-0.48611108, -0.34722222, 0.00000000],
                    [-0.12945666, -1.35603222, 0.00000000],
                    [-0.12943824, 0.15717597, -0.87365150],
                    [-1.55611108, -0.34720903, 0.00000000],
                    [0.02723114, 0.37873406, 1.25740497],
                    [-0.32782521, 1.38810715, 1.25642745],
                    [-0.33103797, -0.12453438, 2.13105486],
                    [1.09722933, 0.37702737, 1.25838372],
                ],
                dtype=float,
            ),
        )
        qmmm_settings = GaussianQMMMJobSettings(
            parent_jobtype="opt",
            jobtype="opt",
            freq=False,
            high_level_functional="mn15",
            high_level_basis="def2svp",
            low_level_force_field="UFF",
            charge_total=0,
            mult_total=1,
            charge_high=0,
            mult_high=1,
            high_level_atoms=[1, 2, 3, 4],
        )
        job = GaussianQMMMJob(
            molecule=molecule,
            settings=qmmm_settings,
            label="gaussian_qmmm_auto_link",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        GaussianInputWriter(job=job).write(target_directory=tmpdir)
        with open(
            os.path.join(tmpdir, "gaussian_qmmm_auto_link.com")
        ) as handle:
            content = handle.read()
        assert "L H 1" in content
        assert job.molecule.bonded_atoms == [(1, 5)]

    def test_write_qmmm_rejects_multiple_link_atoms_on_same_atom(
        self,
        tmpdir,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
    ):
        """Gaussian permits only one link-atom specification per atom."""
        from chemsmart.jobs.gaussian.settings import GaussianQMMMJobSettings

        qmmm_settings = GaussianQMMMJobSettings(
            parent_jobtype="opt",
            jobtype="opt",
            freq=False,
            high_level_functional="mn15",
            high_level_basis="def2svp",
            low_level_force_field="UFF",
            charge_total=0,
            mult_total=1,
            charge_high=0,
            mult_high=1,
            high_level_atoms=[1, 2, 3],
            bonded_atoms=[(2, 4), (3, 4)],
        )
        job = GaussianQMMMJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=qmmm_settings,
            label="gaussian_qmmm_double_link",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        with pytest.raises(ValueError, match="only one link-atom"):
            GaussianInputWriter(job=job).write(target_directory=tmpdir)

    def test_write_qmmm_genecp_follows_conventional_basis_path(
        self,
        tmpdir,
        gaussian_jobrunner_no_scratch,
    ):
        """ONIOM gen/genecp uses the same shared basis path as non-QMMM jobs."""
        import numpy as np

        from chemsmart.io.molecules.structure import Molecule
        from chemsmart.jobs.gaussian.settings import GaussianQMMMJobSettings

        molecule = Molecule(
            symbols=["C", "H", "H", "H", "Cu", "H"],
            positions=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [-1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [2.0, 0.0, 0.0],
                    [3.0, 0.0, 0.0],
                ],
                dtype=float,
            ),
        )
        # Organic high layer: choose a light basis (no heavy_elements), same as
        # a conventional job that does not request GenECP.
        qmmm_light = GaussianQMMMJobSettings(
            parent_jobtype="opt",
            jobtype="opt",
            freq=False,
            high_level_functional="mn15",
            high_level_basis="def2svp",
            low_level_force_field="UFF",
            charge_total=2,
            mult_total=2,
            charge_high=2,
            mult_high=1,
            high_level_atoms=[1, 2, 3, 4],
            bonded_atoms=[(4, 5)],
        )
        job = GaussianQMMMJob(
            molecule=molecule,
            settings=qmmm_light,
            label="gaussian_qmmm_light_basis",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        GaussianInputWriter(job=job).write(target_directory=tmpdir)
        with open(
            os.path.join(tmpdir, "gaussian_qmmm_light_basis.com")
        ) as handle:
            content = handle.read()
        route = next(
            line for line in content.splitlines() if line.startswith("#")
        )
        assert "oniom(mn15/def2svp:UFF)" in route
        assert "genecp" not in route.lower()
        assert "2 2 2 1 2 1" in content
        assert "L H 4" in content
        assert "****" not in content

        # Cu in the high layer with genecp: shared path demotes to gen (Z<=36)
        # and writes the Cu basis block, same as conventional GenECP jobs.
        qmmm_cu = GaussianQMMMJobSettings(
            parent_jobtype="opt",
            jobtype="opt",
            freq=False,
            high_level_functional="mn15",
            high_level_basis="genecp",
            low_level_force_field="UFF",
            heavy_elements=["Cu"],
            heavy_elements_basis="def2-SVPD",
            light_elements_basis="def2SVP",
            charge_total=0,
            mult_total=1,
            charge_high=0,
            mult_high=1,
            high_level_atoms=[1, 2, 3, 4, 5],
            bonded_atoms=[(5, 6)],
        )
        assert qmmm_cu.basis == "genecp"
        job = GaussianQMMMJob(
            molecule=molecule,
            settings=qmmm_cu,
            label="gaussian_qmmm_cu_high",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        GaussianInputWriter(job=job).write(target_directory=tmpdir)
        with open(os.path.join(tmpdir, "gaussian_qmmm_cu_high.com")) as handle:
            content = handle.read()
        route = next(
            line for line in content.splitlines() if line.startswith("#")
        )
        assert "oniom(mn15/gen:UFF)" in route
        assert "Cu     0" in content or "Cu 0" in content

    def test_write_qmmm_amber_requires_mm_atom_info(
        self,
        tmpdir,
        single_molecule_xyz_file,
        gaussian_yaml_settings_qmmm_project_name,
        gaussian_jobrunner_no_scratch,
    ):
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_qmmm_project_name
        )
        qmmm_settings = project_settings.qmmm_settings()
        qmmm_settings.low_level_force_field = "AMBER=HardFirst"
        qmmm_settings.charge_total = 0
        qmmm_settings.mult_total = 1
        qmmm_settings.high_level_atoms = [1, 2, 3]
        job = GaussianQMMMJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=qmmm_settings,
            label="gaussian_qmmm_amber_missing",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        with pytest.raises(
            ValueError, match="mm_atom_info_file|Element-Type-Charge"
        ):
            GaussianInputWriter(job=job).write(target_directory=tmpdir)

    def test_write_qmmm_amber_incomplete_mm_atom_info_raises(
        self,
        tmpdir,
        gaussian_jobrunner_no_scratch,
    ):
        molecule = QMMMMolecule(
            symbols=["C", "O", "H"],
            positions=[
                [0.0, 0.0, 0.0],
                [1.4, 0.0, 0.0],
                [1.9, 0.9, 0.0],
            ],
            high_level_atoms=[1],
            low_level_atoms=[2, 3],
            bonded_atoms=[(1, 2)],
            mm_atom_info=[
                ("CT", 0.03, None, None),
                None,
                ("HC", 0.09, None, None),
            ],
        )
        settings = GaussianQMMMJobSettings(
            high_level_functional="b3lyp",
            high_level_basis="sto-3g",
            low_level_force_field="AMBER=HardFirst",
            charge_total=0,
            mult_total=1,
            high_level_atoms=[1],
            bonded_atoms=[(1, 2)],
        )
        job = GaussianQMMMJob(
            molecule=molecule,
            settings=settings,
            label="gaussian_qmmm_incomplete_mm",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        with pytest.raises(
            ValueError, match="mm_atom_info_file|Element-Type-Charge"
        ):
            GaussianInputWriter(job=job).write(target_directory=tmpdir)

    def test_write_qmmm_inherits_layers_from_molecule(
        self,
        tmpdir,
        gaussian_jobrunner_no_scratch,
    ):
        molecule = QMMMMolecule(
            symbols=["C", "O", "H"],
            positions=[
                [0.0, 0.0, 0.0],
                [1.4, 0.0, 0.0],
                [1.9, 0.9, 0.0],
            ],
            high_level_atoms=[1],
            medium_level_atoms=[2],
            low_level_atoms=[3],
            bonded_atoms=[(1, 2)],
            scale_factors={(1, 2): [0.7]},
            mm_parameters="NonBon 3 1 0 0 0.0 0.0 0.5 0.0 0.0 0.0\n",
        )
        settings = GaussianQMMMJobSettings(
            high_level_functional="b3lyp",
            high_level_basis="sto-3g",
            medium_level_functional="hf",
            medium_level_basis="sto-3g",
            low_level_force_field="UFF",
            charge_total=0,
            mult_total=1,
            high_level_atoms=None,
            medium_level_atoms=None,
            low_level_atoms=None,
            bonded_atoms=None,
            scale_factors=None,
        )
        job = GaussianQMMMJob(
            molecule=molecule,
            settings=settings,
            label="gaussian_qmmm_inherit",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        GaussianInputWriter(job=job).write(target_directory=tmpdir)
        out = os.path.join(tmpdir, "gaussian_qmmm_inherit.com")
        with open(out) as handle:
            content = handle.read()
        assert " H" in content
        assert " M" in content or " L" in content
        assert "NonBon 3 1 0 0 0.0 0.0 0.5 0.0 0.0 0.0" in content

    def test_append_mm_parameters_skips_plain_molecule(self):
        """False branch: molecule is not a QMMMMolecule and has no params file."""
        from io import StringIO

        molecule = Molecule(symbols=["C"], positions=[[0.0, 0.0, 0.0]])
        settings = GaussianQMMMJobSettings(
            high_level_functional="b3lyp",
            high_level_basis="sto-3g",
            low_level_force_field="UFF",
            charge_total=0,
            mult_total=1,
            high_level_atoms=[1],
        )
        job = GaussianQMMMJob(
            molecule=molecule,
            settings=settings,
            label="gaussian_qmmm_plain_mm",
            jobrunner=None,
        )
        buf = StringIO()
        GaussianInputWriter(job=job)._append_mm_parameters(buf)
        assert buf.getvalue() == ""

    def test_write_qmmm_pm6_skips_connectivity_section(
        self,
        tmpdir,
        single_molecule_xyz_file,
        gaussian_yaml_settings_qmmm_project_name,
        gaussian_jobrunner_no_scratch,
    ):
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_qmmm_project_name
        )
        qmmm_settings = project_settings.qmmm_settings()
        qmmm_settings.low_level_force_field = "PM6"
        qmmm_settings.charge_total = 0
        qmmm_settings.mult_total = 1
        qmmm_settings.high_level_atoms = [1, 2, 3]
        job = GaussianQMMMJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=qmmm_settings,
            label="gaussian_qmmm_pm6",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        GaussianInputWriter(job=job).write(target_directory=tmpdir)
        out = os.path.join(tmpdir, "gaussian_qmmm_pm6.com")
        with open(out) as handle:
            content = handle.read()
        assert "geom=connectivity" not in content.lower()
        assert "\n1 2 1.0" not in content

    def test_write_qmmm_mm_parameters_file_overrides_molecule_params(
        self,
        tmpdir,
        gaussian_yaml_settings_qmmm_project_name,
        gaussian_jobrunner_no_scratch,
    ):
        com_path = os.path.join(tmpdir, "prepared_oniom.com")
        with open(com_path, "w") as handle:
            handle.write(
                "%chk=prepared_oniom.chk\n"
                "%nprocshared=2\n"
                "%mem=4GB\n"
                "# opt oniom(b3lyp/sto-3g:AMBER=HardFirst) geom=connectivity\n"
                "\n"
                "prepared ONIOM\n"
                "\n"
                "0 1 0 1 0 1\n"
                "C-CT-0.03      0.0000000000    0.0000000000    0.0000000000 H\n"
                "O-OH--0.65     1.4000000000    0.0000000000    0.0000000000 L H-HC-0.09 1\n"
                "H-HC-0.09      1.9000000000    0.9000000000    0.0000000000 L\n"
                "\n"
                "1 2 1.0\n"
                "2 3 1.0\n"
                "3\n"
                "\n"
                "HrmBnd1 CT OH HC 50.0 109.5\n"
                "\n"
            )
        override = os.path.join(tmpdir, "override_params.dat")
        with open(override, "w") as handle:
            handle.write("NonBon 3 1 0 0 0.0 0.0 0.5 0.0 0.0 0.0")

        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_qmmm_project_name
        )
        qmmm_settings = project_settings.qmmm_settings()
        qmmm_settings.low_level_force_field = "AMBER=HardFirst"
        qmmm_settings.charge_total = 0
        qmmm_settings.mult_total = 1
        qmmm_settings.high_level_atoms = [1]
        qmmm_settings.bonded_atoms = [(1, 2)]
        qmmm_settings.mm_parameters_file = override
        job = GaussianQMMMJob.from_filename(
            filename=com_path,
            settings=qmmm_settings,
            label="gaussian_qmmm_params_override",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        GaussianInputWriter(job=job).write(target_directory=tmpdir)
        out = os.path.join(tmpdir, "gaussian_qmmm_params_override.com")
        with open(out) as handle:
            content = handle.read()
        assert "NonBon 3 1 0 0 0.0 0.0 0.5 0.0 0.0 0.0" in content
        assert "HrmBnd1 CT OH HC 50.0 109.5" not in content

    def test_write_qmmm_amber_from_oniom_com_without_sidecar_files(
        self,
        tmpdir,
        gaussian_yaml_settings_qmmm_project_name,
        gaussian_jobrunner_no_scratch,
    ):
        com_path = os.path.join(tmpdir, "prepared_oniom.com")
        with open(com_path, "w") as handle:
            handle.write(
                "%chk=prepared_oniom.chk\n"
                "%nprocshared=2\n"
                "%mem=4GB\n"
                "# opt oniom(b3lyp/sto-3g:AMBER=HardFirst) geom=connectivity\n"
                "\n"
                "prepared ONIOM\n"
                "\n"
                "0 1 0 1 0 1\n"
                "C-CT-0.03      0.0000000000    0.0000000000    0.0000000000 H\n"
                "O-OH--0.65     1.4000000000    0.0000000000    0.0000000000 L H-HC-0.09 1\n"
                "H-HC-0.09      1.9000000000    0.9000000000    0.0000000000 L\n"
                "\n"
                "1 2 1.0\n"
                "2 3 1.0\n"
                "3\n"
                "\n"
                "HrmBnd1 CT OH HC 50.0 109.5\n"
                "NonBon 3 1 0 0 0.0 0.0 0.5 0.0 0.0 0.0\n"
                "\n"
            )

        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_qmmm_project_name
        )
        qmmm_settings = project_settings.qmmm_settings()
        qmmm_settings.low_level_force_field = "AMBER=HardFirst"
        qmmm_settings.charge_total = 0
        qmmm_settings.mult_total = 1
        qmmm_settings.high_level_atoms = [1]
        qmmm_settings.bonded_atoms = [(1, 2)]
        job = GaussianQMMMJob.from_filename(
            filename=com_path,
            settings=qmmm_settings,
            label="gaussian_qmmm_from_com",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert job.molecule.mm_atom_info is not None
        assert job.molecule.mm_parameters is not None
        GaussianInputWriter(job=job).write(target_directory=tmpdir)
        out = os.path.join(tmpdir, "gaussian_qmmm_from_com.com")
        with open(out) as handle:
            content = handle.read()
        assert "C-CT-0.03" in content
        assert "O-OH--0.65" in content
        assert "H-HC-0.09 1" in content
        assert "HrmBnd1 CT OH HC 50.0 109.5" in content
        assert "NonBon 3 1 0 0 0.0 0.0 0.5 0.0 0.0 0.0" in content

    def test_write_qmmm_uses_layers_from_oniom_com(
        self,
        tmpdir,
        gaussian_qmmm_inputfile_2layer,
        gaussian_yaml_settings_qmmm_project_name,
        gaussian_jobrunner_no_scratch,
    ):
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_qmmm_project_name
        )
        qmmm_settings = project_settings.qmmm_settings()
        qmmm_settings.charge_total = 0
        qmmm_settings.mult_total = 1
        qmmm_settings.high_level_atoms = None
        job = GaussianQMMMJob.from_filename(
            filename=gaussian_qmmm_inputfile_2layer,
            settings=qmmm_settings,
            label="gaussian_qmmm_from_layers",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        GaussianInputWriter(job=job).write(target_directory=tmpdir)
        assert job.molecule.high_level_atoms == [1, 2, 3, 4]
        out = os.path.join(tmpdir, "gaussian_qmmm_from_layers.com")
        with open(out) as handle:
            content = handle.read()
        assert content.count(" H\n") >= 4
        assert " L" in content

    def test_write_qmmm_amber_with_mm_atom_info(
        self,
        tmpdir,
        single_molecule_xyz_file,
        gaussian_yaml_settings_qmmm_project_name,
        gaussian_jobrunner_no_scratch,
    ):
        molecule = Molecule.from_filepath(single_molecule_xyz_file)
        mm_info = os.path.join(tmpdir, "mm_atoms.dat")
        with open(mm_info, "w") as handle:
            for i, symbol in enumerate(molecule.chemical_symbols, start=1):
                if i == 4:
                    handle.write(f"{i} {symbol} 0.03 HC 0.09\n")
                else:
                    handle.write(f"{i} {symbol} 0.0\n")
        params = os.path.join(tmpdir, "mm_params.dat")
        with open(params, "w") as handle:
            handle.write("NonBon 3 1 0 0 0.0 0.0 0.5 0.0 0.0 0.0\n")

        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_qmmm_project_name
        )
        qmmm_settings = project_settings.qmmm_settings()
        qmmm_settings.low_level_force_field = "AMBER=HardFirst"
        qmmm_settings.charge = 0
        qmmm_settings.multiplicity = 1
        qmmm_settings.charge_total = 0
        qmmm_settings.mult_total = 1
        qmmm_settings.high_level_atoms = [1, 2, 3]
        qmmm_settings.bonded_atoms = [(3, 4)]
        qmmm_settings.mm_atom_info_file = mm_info
        qmmm_settings.mm_parameters_file = params
        job = GaussianQMMMJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=qmmm_settings,
            label="gaussian_qmmm_amber",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        GaussianInputWriter(job=job).write(target_directory=tmpdir)
        g16_file = os.path.join(tmpdir, "gaussian_qmmm_amber.com")
        with open(g16_file, "r") as handle:
            content = handle.read()
        assert "AMBER=HardFirst" in content
        assert "geom=connectivity" in content.lower()
        assert "C-C-0.0" in content
        assert "H-HC-0.09 3" in content
        assert "NonBon 3 1 0 0 0.0 0.0 0.5 0.0 0.0 0.0" in content
        assert "\n1 2 " in content

    def test_write_qmmm_missing_mm_parameters_file_raises(
        self,
        tmpdir,
        single_molecule_xyz_file,
        gaussian_yaml_settings_qmmm_project_name,
        gaussian_jobrunner_no_scratch,
    ):
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_qmmm_project_name
        )
        qmmm_settings = project_settings.qmmm_settings()
        qmmm_settings.charge_total = 0
        qmmm_settings.mult_total = 1
        qmmm_settings.high_level_atoms = [1, 2, 3]
        qmmm_settings.mm_parameters_file = os.path.join(
            tmpdir, "missing_mm_params.dat"
        )
        job = GaussianQMMMJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=qmmm_settings,
            label="gaussian_qmmm_missing_params",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        with pytest.raises(FileNotFoundError, match="MM parameters file"):
            GaussianInputWriter(job=job).write(target_directory=tmpdir)

    def test_write_qmmm_input_from_logfile(
        self,
        tmpdir,
        gaussian_yaml_settings_qmmm_project_name,
        gaussian_singlet_opt_outfile,
        gaussian_jobrunner_no_scratch,
        gaussian_written_qmmm_log_file,
    ):
        """Taking the Gaussian nhc_neutral_singlet.log
        output and write qmmm .com"""
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_qmmm_project_name
        )
        qmmm_settings = project_settings.qmmm_settings()
        qmmm_settings.charge = 0
        qmmm_settings.multiplicity = 1
        qmmm_settings.charge_total = 0
        qmmm_settings.real_multiplicity = 1
        qmmm_settings.high_level_atoms = [3, 12, 14, 7, 9]
        qmmm_settings.medium_level_atoms = [8, 17, 19, 20, 21, 22, 23, 24, 25]
        qmmm_settings.bonded_atoms = [(1, 3)]
        job_settings = GaussianJobSettings.from_logfile(
            gaussian_singlet_opt_outfile
        )
        keywords = ("charge", "multiplicity", "title")
        qmmm_settings = qmmm_settings.merge(job_settings, keywords=keywords)
        job = GaussianQMMMJob.from_filename(
            filename=gaussian_singlet_opt_outfile,
            settings=qmmm_settings,
            label="gaussian_qmmm_from_log",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert isinstance(job, GaussianQMMMJob)
        g16_writer = GaussianInputWriter(job=job)
        # write input file
        g16_writer.write(target_directory=tmpdir)
        g16_file = os.path.join(tmpdir, "gaussian_qmmm_from_log.com")
        assert os.path.isfile(g16_file)
        assert cmp(
            g16_file,
            gaussian_written_qmmm_log_file,
            shallow=False,
        )

    def test_write_opt_input_from_logfile(
        self,
        tmpdir,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_singlet_opt_outfile,
        gaussian_jobrunner_no_scratch,
        gaussian_written_ts_from_nhc_singlet_log_file,
    ):
        """Taking the Gaussian nhc_neutral_singlet.log output
        and write aldehyde_opt.com using the settings from the .log file."""
        # get project settings
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        ts_settings = project_settings.ts_settings()
        job_settings = GaussianJobSettings.from_logfile(
            gaussian_singlet_opt_outfile
        )
        # also merge the title keywords
        keywords = ("charge", "multiplicity", "title")
        ts_settings = ts_settings.merge(job_settings, keywords=keywords)
        job = GaussianTSJob.from_filename(
            filename=gaussian_singlet_opt_outfile,
            settings=ts_settings,
            label="gaussian_ts_from_log",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert isinstance(job, GaussianTSJob)
        g16_writer = GaussianInputWriter(job=job)
        # write input file
        g16_writer.write(target_directory=tmpdir)
        g16_file = os.path.join(tmpdir, "gaussian_ts_from_log.com")
        assert os.path.isfile(g16_file)
        assert cmp(
            g16_file,
            gaussian_written_ts_from_nhc_singlet_log_file,
            shallow=False,
        )

    def test_write_sp_input_with_solvation_from_logfile(
        self,
        tmpdir,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_singlet_opt_outfile,
        gaussian_jobrunner_no_scratch,
        gaussian_written_sp_from_nhc_singlet_log_with_solvent_file,
    ):
        """Test writing simple .com input file using settings from .log file,
        including solvation."""

        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        sp_settings = project_settings.sp_settings()
        job_settings = GaussianJobSettings.from_logfile(
            gaussian_singlet_opt_outfile
        )
        sp_settings = sp_settings.merge(job_settings)
        job = GaussianSinglePointJob.from_filename(
            filename=gaussian_singlet_opt_outfile,
            settings=sp_settings,
            label="gaussian_sp_from_log_with_solvent",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert isinstance(job, GaussianSinglePointJob)
        g16_writer = GaussianInputWriter(job=job)
        # write input file
        g16_writer.write(target_directory=tmpdir)
        g16_file = os.path.join(
            tmpdir, "gaussian_sp_from_log_with_solvent.com"
        )
        assert os.path.isfile(g16_file)
        assert cmp(
            g16_file,
            gaussian_written_sp_from_nhc_singlet_log_with_solvent_file,
            shallow=False,
        )

    def test_write_sp_with_custom_solvation_from_logfile(
        self,
        tmpdir,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_singlet_opt_outfile,
        gaussian_jobrunner_no_scratch,
        smd_TBME_solvent_parameters_text_file,
        gaussian_written_sp_from_nhc_singlet_log_with_custom_solvent_file,
    ):
        """Test writing input file from log file.
        Simply taking the Gaussian nhc_neutral_singlet.log output and write
        gaussian_sp_custom_solv.com using the settings from the .log
        file and including custom solvation parameters from file smd_TBME.
        """
        smd_TBME_tmp_path = os.path.join(tmpdir, "smd_TBME")
        copy(smd_TBME_solvent_parameters_text_file, smd_TBME_tmp_path)

        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        sp_settings = project_settings.sp_settings()
        sp_settings.solvent_model = "smd"
        sp_settings.custom_solvent = smd_TBME_tmp_path
        job_settings = GaussianJobSettings.from_logfile(
            gaussian_singlet_opt_outfile
        )
        sp_settings = sp_settings.merge(job_settings)

        job = GaussianSinglePointJob.from_filename(
            filename=gaussian_singlet_opt_outfile,
            settings=sp_settings,
            label="gaussian_sp_from_log_with_custom_solvent",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert isinstance(job, GaussianSinglePointJob)
        g16_writer = GaussianInputWriter(job=job)
        # write input file
        g16_writer.write(target_directory=tmpdir)
        g16_file = os.path.join(
            tmpdir, "gaussian_sp_from_log_with_custom_solvent.com"
        )
        assert os.path.isfile(g16_file)
        assert cmp(
            g16_file,
            gaussian_written_sp_from_nhc_singlet_log_with_custom_solvent_file,
            shallow=False,
        )

    def test_write_ts_with_custom_basis_from_logfile(
        self,
        tmpdir,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_singlet_opt_outfile,
        gaussian_jobrunner_no_scratch,
        Ni_def2tzvp_PCHOSi_svp_text_file,
        gaussian_written_sp_from_nhc_singlet_log_with_custom_basis_file,
    ):
        custom_basis_tmp_path = os.path.join(tmpdir, "custom_basis.txt")
        copy(Ni_def2tzvp_PCHOSi_svp_text_file, custom_basis_tmp_path)

        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        ts_settings = project_settings.ts_settings()
        ts_settings.gen_genecp_file = custom_basis_tmp_path
        job_settings = GaussianJobSettings.from_logfile(
            gaussian_singlet_opt_outfile
        )
        ts_settings = ts_settings.merge(job_settings)

        job = GaussianTSJob.from_filename(
            filename=gaussian_singlet_opt_outfile,
            settings=ts_settings,
            label="gaussian_sp_from_log_with_custom_basis",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert isinstance(job, GaussianTSJob)
        g16_writer = GaussianInputWriter(job=job)
        # write input file
        g16_writer.write(target_directory=tmpdir)
        g16_file = os.path.join(
            tmpdir, "gaussian_sp_from_log_with_custom_basis.com"
        )
        assert os.path.isfile(g16_file)
        assert cmp(
            g16_file,
            gaussian_written_sp_from_nhc_singlet_log_with_custom_basis_file,
            shallow=False,
        )

    def test_write_ts_with_custom_basis_using_api(
        self,
        tmpdir,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_ts_genecp_outfile,
        gaussian_jobrunner_no_scratch,
        gaussian_written_sp_from_nhc_singlet_log_with_custom_basis_from_api_file,
    ):
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        ts_settings = project_settings.ts_settings()
        job_settings = GaussianJobSettings.from_logfile(
            gaussian_ts_genecp_outfile
        )
        ts_settings = ts_settings.merge(job_settings)

        # update settings for heavy elements
        ts_settings.heavy_elements = ["Pd"]
        ts_settings.heavy_elements_basis = "def2-TZVPPD"
        ts_settings.light_elements_basis = "def2-SVP"

        molecule = Molecule.from_filepath(gaussian_ts_genecp_outfile)

        job = GaussianTSJob(
            molecule=molecule,
            settings=ts_settings,
            label="gaussian_sp_from_log_with_custom_basis_from_api",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert isinstance(job, GaussianTSJob)
        g16_writer = GaussianInputWriter(job=job)
        # write input file
        g16_writer.write(target_directory=tmpdir)
        g16_file = os.path.join(
            tmpdir, "gaussian_sp_from_log_with_custom_basis_from_api.com"
        )

        # compare the written input file with the expected input file, except
        # the line containing the version number
        # (basis set exchange api may be different)
        assert cmp_with_ignore(
            g16_file,
            gaussian_written_sp_from_nhc_singlet_log_with_custom_basis_from_api_file,
            ignore_string="Version",
        )

    def test_write_modred_with_custom_basis_for_all_elements_in_structure_using_api(
        self,
        tmpdir,
        gaussian_yaml_settings_gas_solv_project_name,
        modred_genecp_inputfile,
        gaussian_jobrunner_no_scratch,
        gaussian_modred_with_custom_basis_for_all_atoms_from_api,
    ):
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        modred_settings = project_settings.modred_settings()
        job_settings = GaussianJobSettings.from_filepath(
            modred_genecp_inputfile
        )
        modred_settings = modred_settings.merge(job_settings)

        modred_settings.modred = [[1, 2], [3, 4, 5]]
        modred_settings.basis = "genecp"
        modred_settings.heavy_elements = [
            "C",
            "H",
            "O",
            "N",
            "Pd",
            "P",
            "S",
            "F",
            "Br",
            "I",
        ]
        # more than all elements in the system but will be filtered
        # to only those in the system for input preparation
        modred_settings.heavy_elements_basis = "def2-TZVPPD"
        modred_settings.light_elements_basis = None  # light element basis not specified as all use custom basis from heavy_elements_basis

        molecule = Molecule.from_filepath(modred_genecp_inputfile)
        job = GaussianModredJob(
            molecule=molecule,
            settings=modred_settings,
            label="gaussian_modred_with_custom_basis_for_all_atoms_from_api",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert isinstance(job, GaussianModredJob)
        g16_writer = GaussianInputWriter(job=job)

        # write input file
        g16_writer.write(target_directory=tmpdir)
        g16_file = os.path.join(
            tmpdir,
            "gaussian_modred_with_custom_basis_for_all_atoms_from_api.com",
        )
        assert os.path.isfile(g16_file)

        # compare the written input file with the expected input file, except
        # the line containing the version number
        # (basis set exchange api may be different)
        assert cmp_with_ignore(
            g16_file,
            gaussian_modred_with_custom_basis_for_all_atoms_from_api,
            ignore_string="Version",
        )

    def test_write_gaussian_input_from_pbc_logfile(
        self,
        tmpdir,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_pbc_2d_outputfile,
        gaussian_jobrunner_no_scratch,
        gaussian_written_opt_from_graphite_2d_pbc_log,
    ):
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        opt_settings = project_settings.opt_settings()
        file_settings = GaussianJobSettings.from_filepath(
            filepath=gaussian_pbc_2d_outputfile
        )
        opt_settings = opt_settings.merge(file_settings)

        molecule = Molecule.from_filepath(gaussian_pbc_2d_outputfile)
        g16_outputfile = Gaussian16Output(gaussian_pbc_2d_outputfile)
        assert g16_outputfile.list_of_pbc_conditions == [1, 1, 0]
        assert g16_outputfile.input_translation_vectors == [
            [2.47532, 0.0, 0.0],
            [-1.21995, 2.13345, 0.0],
        ]

        job = GaussianSinglePointJob(
            molecule=molecule,
            settings=opt_settings,
            label="graphite_2d_opt_from_log",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert isinstance(job, GaussianSinglePointJob)

        g16_writer = GaussianInputWriter(job=job)
        g16_writer.write(target_directory=tmpdir)
        g16_file = os.path.join(tmpdir, "graphite_2d_opt_from_log.com")
        assert os.path.isfile(g16_file)
        assert cmp(
            g16_file,
            gaussian_written_opt_from_graphite_2d_pbc_log,
            shallow=False,
        )
