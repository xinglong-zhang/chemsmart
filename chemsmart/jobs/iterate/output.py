"""Incremental XYZ output delivery for CHEMSMART Iterate jobs."""

from __future__ import annotations

import os
import stat
import tempfile
import uuid
from typing import TYPE_CHECKING, BinaryIO, Optional

from chemsmart.jobs.iterate.report import (
    STAGE_WRITE,
    STATUS_SUCCESS,
    STATUS_WRITE_FAILED,
    CombinationResult,
)

if TYPE_CHECKING:
    from chemsmart.io.molecules.structure import Molecule
    from chemsmart.jobs.iterate.job import IterateJob


def _xyz_bytes(label: str, molecule: "Molecule") -> bytes:
    """Render one labelled molecule using CHEMSMART's XYZ layout."""
    lines = [f"{molecule.num_atoms}\n", f"       {label}\n"]
    lines.extend(
        f"{symbol:2s}  {position[0]:15.10f}  "
        f"{position[1]:15.10f}  {position[2]:15.10f}\n"
        for symbol, position in zip(
            molecule.chemical_symbols, molecule.positions
        )
    )
    return "".join(lines).encode("utf-8")


def _write_bytes_atomically(path: str, content: bytes) -> None:
    """Atomically replace ``path`` with ``content``."""
    descriptor, temporary_path = _open_output_temporary(path)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_path, path)
    except Exception:
        try:
            os.close(descriptor)
        except OSError:
            pass
        try:
            os.remove(temporary_path)
        except OSError:
            pass
        raise


def _open_output_temporary(path: str) -> tuple[int, str]:
    """Create a same-directory output temporary with normal file modes."""
    directory = os.path.dirname(os.path.abspath(path))
    os.makedirs(directory, exist_ok=True)
    try:
        existing_mode = stat.S_IMODE(os.stat(path).st_mode)
    except FileNotFoundError:
        existing_mode = None

    while True:
        temporary_path = os.path.join(
            directory,
            f".{os.path.basename(path)}.{uuid.uuid4().hex}.tmp",
        )
        try:
            descriptor = os.open(
                temporary_path,
                os.O_WRONLY | os.O_CREAT | os.O_EXCL,
                0o666,
            )
        except FileExistsError:
            continue
        if existing_mode is not None:
            os.fchmod(descriptor, existing_mode)
        return descriptor, temporary_path


class IterateOutputWriter:
    """Write successful results incrementally in the parent process.

    Separate outputs are committed atomically as results arrive. Merged
    outputs are spooled incrementally, then assembled in original combination
    order and atomically committed when execution finishes. In both modes the
    generated :class:`Molecule` is released immediately after serialization.
    """

    def __init__(self, job: "IterateJob") -> None:
        self.job = job
        self._spool: Optional[BinaryIO] = None
        self._spool_path: Optional[str] = None
        self._blocks: dict[int, tuple[int, int]] = {}
        self._claimed_paths: set[str] = set()

    def consume(self, result: CombinationResult) -> None:
        """Serialize one successful result and release its molecule."""
        if (
            result.execution_status != STATUS_SUCCESS
            or result.molecule is None
        ):
            return

        molecule = result.molecule
        try:
            content = _xyz_bytes(result.label, molecule)
            if self.job.separate_outputs:
                self._write_separate(result, content)
            else:
                self._write_spool(result, content)
        except Exception as error:
            if self.job.separate_outputs:
                output_path = result.output_path
            else:
                output_path = self.job.outputfile
            self._mark_write_failed(result, error, output_path)
        finally:
            result.molecule = None

    def finalize(self, results: list[CombinationResult]) -> list[str]:
        """Commit output and return paths in original combination order."""
        if self.job.separate_outputs:
            written = [
                result
                for result in results
                if result.execution_status == STATUS_SUCCESS
                and result.output_path
            ]
            for index, result in enumerate(written, start=1):
                result.structure_index = index
            return [result.output_path for result in written]

        if not self._blocks:
            self.abort()
            return []

        output_path = self.job.outputfile
        temporary_path: Optional[str] = None
        descriptor: Optional[int] = None
        written_results: list[CombinationResult] = []

        try:
            assert self._spool is not None
            self._spool.flush()
            os.fsync(self._spool.fileno())
            descriptor, temporary_path = _open_output_temporary(output_path)
            with os.fdopen(descriptor, "wb") as output:
                descriptor = None
                for result in results:
                    if result.execution_status != STATUS_SUCCESS:
                        continue
                    offset, length = self._blocks[result.combination_number]
                    self._spool.seek(offset)
                    content = self._spool.read(length)
                    if len(content) != length:
                        raise OSError(
                            "Incomplete temporary XYZ block for combination "
                            f"{result.combination_number}."
                        )
                    output.write(content)
                    written_results.append(result)
                output.flush()
                os.fsync(output.fileno())
            os.replace(temporary_path, output_path)
            temporary_path = None
        except Exception as error:
            for result in results:
                if result.execution_status == STATUS_SUCCESS:
                    self._mark_write_failed(result, error, output_path)
            return []
        finally:
            if descriptor is not None:
                try:
                    os.close(descriptor)
                except OSError:
                    pass
            if temporary_path is not None:
                try:
                    os.remove(temporary_path)
                except OSError:
                    pass
            self.abort()

        for index, result in enumerate(written_results, start=1):
            result.output_path = output_path
            result.structure_index = index
        return [output_path]

    def abort(self) -> None:
        """Discard only uncommitted merged-output temporary data."""
        if self._spool is not None:
            try:
                self._spool.close()
            except OSError:
                pass
            self._spool = None
        if self._spool_path is not None:
            try:
                os.remove(self._spool_path)
            except OSError:
                pass
            self._spool_path = None
        self._blocks.clear()

    def _write_separate(
        self, result: CombinationResult, content: bytes
    ) -> None:
        output_directory = self.job.output_directory or "."
        output_path = os.path.join(output_directory, f"{result.label}.xyz")
        if output_path in self._claimed_paths:
            raise ValueError(
                f"Duplicate output filename '{output_path}' within this run."
            )
        self._claimed_paths.add(output_path)
        result.output_path = output_path
        _write_bytes_atomically(output_path, content)

    def _write_spool(self, result: CombinationResult, content: bytes) -> None:
        if result.combination_number in self._blocks:
            raise ValueError(
                "Duplicate combination number "
                f"{result.combination_number} in merged output."
            )
        spool = self._ensure_spool()
        offset = spool.tell()
        try:
            spool.write(content)
        except Exception:
            try:
                spool.seek(offset)
                spool.truncate()
            except OSError:
                pass
            raise
        self._blocks[result.combination_number] = (offset, len(content))

    def _ensure_spool(self) -> BinaryIO:
        if self._spool is not None:
            return self._spool
        output_path = self.job.outputfile
        directory = os.path.dirname(os.path.abspath(output_path))
        os.makedirs(directory, exist_ok=True)
        descriptor, self._spool_path = tempfile.mkstemp(
            prefix=f".{os.path.basename(output_path)}.",
            suffix=".spool",
            dir=directory,
        )
        self._spool = os.fdopen(descriptor, "w+b")
        return self._spool

    @staticmethod
    def _mark_write_failed(
        result: CombinationResult,
        error: Exception,
        output_path: Optional[str] = None,
    ) -> None:
        result.execution_status = STATUS_WRITE_FAILED
        result.failure_stage = STAGE_WRITE
        result.error_type = type(error).__name__
        result.error_message = str(error)
        if output_path is not None:
            result.output_path = output_path
