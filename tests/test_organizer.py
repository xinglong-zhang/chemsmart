from chemsmart.io.organizer import FileOrganizer


class TestOrganizer:

    def test_file_organizer(self, tmpdir, excel_file):
        file_organizer = FileOrganizer(
            filename=excel_file, sheetname="co2", cols="B:D", skip=2, row=100
        )
        df = file_organizer.load_data()
        assert df.shape == (98, 3)
        folders, new_filenames, old_filenames = file_organizer.clean_data()
        assert len(folders) == len(new_filenames) == len(old_filenames) == 31
