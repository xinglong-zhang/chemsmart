
from chemsmart.io.organizer import FileOrganizer


class TestOrganizer:

    def test_file_organizer(self, tmpdir, excel_file):
        file_organizer = FileOrganizer(filename=excel_file, sheetname="co2")
        df = file_organizer.load_data()
        print(df)
