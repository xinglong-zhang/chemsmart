import os
from functools import cached_property
class GaussianFileMixin:
    """Mixin class for files that can be supplied to prepare Gaussian input files.

    Args:
        file to be given to be parsed.
    """

    def __init__(self, filename):
        self.filename = filename

    @property
    def filepath(self):
        return os.path.abspath(self.filename)

    @property
    def base_filename_with_extension(self):
        return self.filepath.split('/')[-1]

    @property
    def basename(self):
        return self.base_filename_with_extension.split('.')[0]

    @cached_property
    def contents(self):
        with open(self.filepath) as f:
            return [line.strip() for line in f.readlines()]

    @cached_property
    def content_lines_string(self):
        with open(self.filepath) as f:
            return f.read()