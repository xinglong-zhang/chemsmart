import os
from shutil import copy

from chemsmart.jobs.nciplot import NCIPLOTJob
from chemsmart.jobs.nciplot.settings import NCIPLOTJobSettings
from chemsmart.jobs.nciplot.writer import NCIPLOTInputWriter


class TestNCIPLOTInputWriter:
    def test_write_nci_for_one_file(
        self,
        tmpdir,
        single_molecule_xyz_file,
        nciplot_jobrunner_no_scratch,
    ):
        job_settings = NCIPLOTJobSettings(
            rthres=10.0,
            ligand_file_number=1,
            ligand_radius=1.0,
            radius_positions="(1.0, 1.1, 1.2)",
            radius_r=1.5,
            intercut1=0.5,
            intercut2=0.6,
            increments="0.1,0.1,0.1",
            fragments={1: [1, 2, 3], 2: [4, 5, 6]},
            cutoff_density_dat=0.01,
            cutoff_rdg_dat=0.02,
            cutoff_density_cube=0.03,
            cutoff_rdg_cube=0.04,
            dgrid=True,
            integrate=True,
            ranges=[[-0.1, -0.02], [-0.02, 0.02], [0.02, 0.1]],
            grid_quality="coarse",
        )

        # copy file to tmpdir
        tmpdir_xyz_file = os.path.join(tmpdir, "single_molecule.xyz")
        copy(single_molecule_xyz_file, tmpdir_xyz_file)

        # change to tmpdir and run test
        os.chdir(tmpdir)

        # create nciplot job
        job = NCIPLOTJob(
            filenames=("single_molecule.xyz",),
            settings=job_settings,
            jobrunner=nciplot_jobrunner_no_scratch,
        )
        nciplot_writer = NCIPLOTInputWriter(job=job)

        # write input file
        nciplot_writer.write(target_directory=tmpdir)

        # if no job label, then base filename is used
        nci_file = os.path.join(tmpdir, "single_molecule.nci")
        assert os.path.exists(nci_file)
        lines = open(nci_file, "r").readlines()
        assert lines[0] == "1\n"
        assert lines[1] == "single_molecule.xyz\n"
        assert lines[2] == "RTHRES 10.0\n"
        assert lines[3] == "LIGAND 1 1.0\n"
        assert lines[4] == "RADIUS 1.0  1.1  1.2 1.5\n"
        assert lines[5] == "INTERMOLECULAR\n"
        assert lines[6] == "INTERCUT 0.5 0.6\n"
        assert lines[7] == "INCREMENTS 0.1 0.1 0.1\n"
        assert lines[8] == "FRAGMENTS\n"
        assert lines[9] == "  1 1 2 3\n"
        assert lines[10] == "  2 4 5 6\n"
        assert lines[11] == "END\n"
        assert lines[12] == "CUTOFFS 0.01 0.02\n"
        assert lines[13] == "CUTPLOT 0.03 0.04\n"
        assert lines[14] == "DGRID\n"
        assert lines[15] == "INTEGRATE\n"
        assert lines[16] == "RANGE 3\n"
        assert lines[17] == "-0.1 -0.02\n"
        assert lines[18] == "-0.02 0.02\n"
        assert lines[19] == "0.02 0.1\n"
        assert lines[20] == "COARSE\n"

    def test_write_nci_for_two_files(
        self,
        tmpdir,
        single_molecule_xyz_file,
        nciplot_jobrunner_no_scratch,
    ):
        job_settings = NCIPLOTJobSettings(
            rthres=10.0,
            ligand_file_number=1,
            ligand_radius=1.0,
            radius_positions="(1.0, 1.1, 1.2)",
            radius_r=1.5,
            intercut1=0.5,
            intercut2=0.6,
            increments="0.1,0.1,0.1",
            fragments={1: [1, 2, 3], 2: [4, 5, 6]},
            cutoff_density_dat=0.01,
            cutoff_rdg_dat=0.02,
            cutoff_density_cube=0.03,
            cutoff_rdg_cube=0.04,
            dgrid=True,
            integrate=True,
            ranges=[[-0.1, -0.02], [-0.02, 0.02], [0.02, 0.1]],
            grid_quality="ultrafine",
        )

        # copy file to tmpdir
        tmpdir_xyz_file = os.path.join(tmpdir, "single_molecule.xyz")
        copy(single_molecule_xyz_file, tmpdir_xyz_file)

        tmpdir_xyz_file2 = os.path.join(tmpdir, "single_molecule2.xyz")
        copy(single_molecule_xyz_file, tmpdir_xyz_file2)

        # change to tmpdir and run test
        os.chdir(tmpdir)

        # create nciplot job
        job = NCIPLOTJob(
            filenames=("single_molecule.xyz", "single_molecule2.xyz"),
            settings=job_settings,
            label="nci_two_files",
            jobrunner=nciplot_jobrunner_no_scratch,
        )
        nciplot_writer = NCIPLOTInputWriter(job=job)

        # write input file
        nciplot_writer.write(target_directory=tmpdir)

        # if job label is given, then filename takes job label
        # NCIJobSettings.label is not used
        nci_file = os.path.join(tmpdir, "nci_two_files.nci")
        assert os.path.exists(nci_file)
        lines = open(nci_file, "r").readlines()
        assert lines[0] == "2\n"
        assert lines[1] == "single_molecule.xyz\n"
        assert lines[2] == "single_molecule2.xyz\n"
        assert lines[3] == "RTHRES 10.0\n"
        assert lines[4] == "LIGAND 1 1.0\n"
        assert lines[5] == "RADIUS 1.0  1.1  1.2 1.5\n"
        assert lines[6] == "INTERMOLECULAR\n"
        assert lines[7] == "INTERCUT 0.5 0.6\n"
        assert lines[8] == "INCREMENTS 0.1 0.1 0.1\n"
        assert lines[9] == "FRAGMENTS\n"
        assert lines[10] == "  1 1 2 3\n"
        assert lines[11] == "  2 4 5 6\n"
        assert lines[12] == "END\n"
        assert lines[13] == "CUTOFFS 0.01 0.02\n"
        assert lines[14] == "CUTPLOT 0.03 0.04\n"
        assert lines[15] == "DGRID\n"
        assert lines[16] == "INTEGRATE\n"
        assert lines[17] == "RANGE 3\n"
        assert lines[18] == "-0.1 -0.02\n"
        assert lines[19] == "-0.02 0.02\n"
        assert lines[20] == "0.02 0.1\n"
        assert lines[21] == "ULTRAFINE\n"

        # job run will result in the job being run and the output file copied back to run folder
        # job.run()
        # assert job.is_complete()
