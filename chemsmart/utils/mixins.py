import os
from functools import cached_property

class FileMixin:
    """Mixin class for files that can be opened and read"""

    @property
    def filepath(self):
        return os.path.abspath(self.filename)

    @property
    def filepath_directory(self):
        return os.path.split(self.filepath)[0]

    @property
    def base_filename_with_extension(self):
        return os.path.split(self.filepath)[1]

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
