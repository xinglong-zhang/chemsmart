"""DGL QM9 and revQM9 geometry input."""

import os
import pickle
from importlib import import_module
from unittest.mock import MagicMock, patch
from zipfile import ZipFile

import numpy as np
import pytest
from click.testing import CliRunner
from numpy.lib.format import write_array_header_1_0

from chemsmart.io.molecules.structure import Molecule
from chemsmart.io.npz import NPZFile


def _write_revqm9(path):
    coords = np.empty(2, dtype=object)
    numbers = np.empty(2, dtype=object)
    elements = np.empty(2, dtype=object)
    coords[:] = [np.zeros((1, 3)), np.array([[0, 0, 0], [0, 0, 0.74]])]
    numbers[:] = [np.array([2]), np.array([1, 1])]
    elements[:] = [np.array(["He"]), np.array(["H", "H"])]
    np.savez_compressed(
        path, coords=coords, charges=numbers, elements=elements
    )
    return coords


def test_revqm9_ragged_arrays_and_selection(tmp_path):
    path = tmp_path / "revQM9.npz"
    coords = _write_revqm9(path)
    molecules = Molecule.from_filepath(path, index=":", return_list=True)
    assert [list(m.symbols) for m in molecules] == [["He"], ["H", "H"]]
    assert [m.structure_index_in_file for m in molecules] == [1, 2]
    assert molecules[0].charge is None and molecules[0].energy is None
    np.testing.assert_array_equal(molecules[1].positions, coords[1])
    assert list(Molecule.from_filepath(path, index="-1").symbols) == ["H", "H"]
    assert [
        m.structure_index_in_file
        for m in Molecule.from_filepath(path, index="2,1")
    ] == [2, 1]
    assert len(Molecule.from_filepath(path, index="1:2")) == 1
    assert len(Molecule.from_filepath(path, index="1", return_list=True)) == 1


def test_dgl_concatenated_molecules_and_selection(tmp_path):
    path = tmp_path / "qm9_eV.npz"
    coords = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]])
    np.savez(path, R=coords, Z=[2, 1, 1, 8], N=[1, 2, 1])
    molecules = Molecule.from_filepath(path, index=":", return_list=True)
    assert [list(m.symbols) for m in molecules] == [["He"], ["H", "H"], ["O"]]
    assert [m.structure_index_in_file for m in molecules] == [1, 2, 3]
    np.testing.assert_array_equal(molecules[1].positions, coords[1:3])
    assert list(NPZFile(path).molecule.symbols) == ["O"]
    assert [
        m.structure_index_in_file
        for m in Molecule.from_filepath(path, index="3,1-2")
    ] == [3, 1, 2]
    with pytest.raises(ValueError, match="1-based"):
        Molecule.from_filepath(path, index="0")
    with pytest.raises(IndexError):
        Molecule.from_filepath(path, index="4")


@pytest.mark.parametrize(
    "counts, numbers, message",
    [
        ([1, 1], [1, 1, 1], "sum to the number"),
        ([1, 0, 2], [1, 1, 1], "positive integer"),
        ([1, -1, 3], [1, 1, 1], "positive integer"),
        ([1.0, 2.0], [1, 1, 1], "positive integer"),
        ([[1, 2]], [1, 1, 1], "positive integer"),
        ([1, 2], [1, 0, 1], "zero-padded"),
        ([1, 2], [1, 1], "matching atom counts"),
    ],
)
def test_dgl_rejects_invalid_data(tmp_path, counts, numbers, message):
    path = tmp_path / "invalid.npz"
    np.savez(path, R=np.zeros((3, 3)), Z=numbers, N=counts)
    with pytest.raises(ValueError, match=message):
        Molecule.from_filepath(path, index=":")


def test_npz_rejects_unsupported_layout_and_nonfinite_coordinates(tmp_path):
    path = tmp_path / "bad.npz"
    np.savez(path, positions=[[0, 0, 0]], atomic_numbers=[2])
    with pytest.raises(ValueError, match="DGL QM9 fields"):
        Molecule.from_filepath(path)
    np.savez(path, R=[[float("nan"), 0, 0]], Z=[2], N=[1])
    with pytest.raises(ValueError, match="finite"):
        Molecule.from_filepath(path)


@pytest.mark.parametrize(
    "extra", [{"elements": ["C"]}, {"atomic_numbers": [6]}]
)
def test_npz_rejects_conflicting_atom_identity(tmp_path, extra):
    path = tmp_path / "conflict.npz"
    np.savez(path, coords=[[0, 0, 0]], charges=[8], **extra)
    with pytest.raises(ValueError, match="conflicting atomic identities"):
        Molecule.from_filepath(path)


def test_npz_rejects_mismatched_atom_numbers_and_elements(tmp_path):
    path = tmp_path / "conflict.npz"
    np.savez(
        path,
        coords=[[0, 0, 0]],
        charges=[8],
        atomic_numbers=[8],
        elements=["C"],
    )
    with pytest.raises(ValueError, match="conflicting atomic identities"):
        Molecule.from_filepath(path)


def test_npz_rejects_padded_atom_arrays(tmp_path):
    path = tmp_path / "padded.npz"
    np.savez(path, R=np.zeros((2, 3, 3)), Z=[[6, 1, 0], [8, 0, 0]], N=[2, 1])
    with pytest.raises(ValueError, match="zero-padded"):
        Molecule.from_filepath(path, index=":")


def test_npz_rejects_arbitrary_pickle_globals(tmp_path):
    marker = tmp_path / "unpickled"

    class Payload:
        def __reduce__(self):
            return os.mkdir, (str(marker),)

    path = tmp_path / "untrusted.npz"
    np.savez(path, coords=np.array([Payload()], dtype=object), charges=[[1]])
    with pytest.raises(ValueError, match="Disallowed pickle global"):
        Molecule.from_filepath(path)
    assert not marker.exists()


def test_npz_rejects_pickled_array_header_mismatch(tmp_path):
    path = tmp_path / "mismatch.npz"
    with ZipFile(path, "w") as archive:
        with archive.open("coords.npy", "w") as stream:
            write_array_header_1_0(
                stream, {"descr": "|O", "fortran_order": False, "shape": (2,)}
            )
            pickle.dump(np.array([np.zeros((1, 3))], dtype=object), stream)
        with archive.open("charges.npy", "w") as stream:
            np.lib.format.write_array(stream, np.array([[2]]))
    with pytest.raises(ValueError, match="does not match its NPY header"):
        Molecule.from_filepath(path)


@pytest.mark.parametrize(
    "backend, project, job_path, subcommand",
    [
        (
            "gaussian",
            "gas_solv",
            "chemsmart.jobs.gaussian.opt.GaussianOptJob",
            "opt",
        ),
        ("orca", "gas_solv", "chemsmart.jobs.orca.opt.ORCAOptJob", "opt"),
        ("xtb", "test", "chemsmart.jobs.xtb.opt.XTBOptJob", "opt"),
        (
            "crest",
            "test",
            "chemsmart.jobs.crest.conformers.CRESTConformerSearchJob",
            "conformers",
        ),
    ],
)
def test_cli_selects_npz_molecule(
    tmp_path,
    backend,
    project,
    job_path,
    subcommand,
    make_cli_ctx_obj,
    gaussian_jobrunner_no_scratch,
):
    path = tmp_path / "revQM9.npz"
    _write_revqm9(path)
    command = getattr(
        import_module(f"chemsmart.cli.{backend}.{backend}"), backend
    )
    original_reader = NPZFile.get_molecules
    with (
        patch(job_path, return_value=MagicMock()) as job,
        patch.object(
            NPZFile,
            "get_molecules",
            autospec=True,
            side_effect=original_reader,
        ),
    ):
        result = CliRunner().invoke(
            command,
            [
                "-p",
                project,
                "-f",
                str(path),
                "-i",
                "2",
                "-c",
                "0",
                "-m",
                "1",
                subcommand,
            ],
            obj=make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
    assert result.exit_code == 0, (result.output, result.exception)
    job.assert_called_once()
    assert list(job.call_args.kwargs["molecule"].symbols) == ["H", "H"]
