from chemsmart.utils.mixins import FileMixin

class CubeFile(FileMixin):
    def __init__(self, filename):
        super().__init__(filename=filename)

