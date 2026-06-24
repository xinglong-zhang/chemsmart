"""
GROMACS input writer.

At this stage, prepared workflows use user-provided MDP, structure and topology
files directly. This writer is kept as a small extension point for future MDP
generation.
"""

from __future__ import annotations

from pathlib import Path


class GromacsInputWriter:
    """
    Minimal writer interface for future GROMACS input generation.
    """

    def __init__(self, job):
        self.job = job

    def write(self):
        """
        Write generated GROMACS input files.

        The current implementation does not generate files automatically because
        prepared workflows rely on user-provided files.
        """
        return []

    def write_text_file(self, filename, content):
        """
        Utility method for writing a text file into the job folder.
        """
        path = Path(self.job.folder) / filename
        path.write_text(content, encoding="utf-8")
        return path
