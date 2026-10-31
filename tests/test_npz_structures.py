"""DGL QM9 and revQM9 geometry input."""

import os
import pickle
import sqlite3
from importlib import import_module
from pathlib import Path
from unittest.mock import MagicMock, patch
from zipfile import ZipFile

import numpy as np
import pytest
from ase.units import Hartree
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


class TestNPZStructures:
    @pytest.mark.parametrize(
        "filename", ["revQM9_subset.npz", "qm9_eV_subset.npz"]
    )
    def test_npz_database_import(self, tmp_path, filename):
        from chemsmart.cli.database.assemble import assemble
        from chemsmart.database.database import Database

        source = Path(__file__).parent / "data" / "IOTests" / filename
        output = tmp_path / "dataset.db"
        args = ["-f", str(source), "-i", "3,1", "-o", str(output)]
        result = CliRunner().invoke(assemble, args, obj={})
        assert result.exit_code == 0, (result.output, result.exception)
        with sqlite3.connect(output) as conn:
            rows = conn.execute(
                "SELECT r.program, r.total_energy, r.homo_energy, "
                "r.normal_termination, rs.structure_index_in_file, rs.energy "
                "FROM records r JOIN record_structures rs USING (record_id) "
                "ORDER BY r.record_index"
            ).fetchall()
        assert len(rows) == 2
        assert [row[4] for row in rows] == [3, 1]
        assert all(row[3] is None for row in rows)
        record = Database(str(output)).get_record(record_index=2)
        properties = record["molecules"][0]["dataset_properties"]
        restored = Molecule.from_filepath(str(output), record_index=2)
        assert (restored.charge, restored.multiplicity) == (0, 1)
        assert restored.structure_index_in_file == 1
        if filename.startswith("rev"):
            assert rows[1][1] == pytest.approx(-40.47574250955533)
            assert rows[1][5] == pytest.approx(rows[1][1])
            assert rows[1][2] == pytest.approx(-0.4093204819908599 * Hartree)
            assert properties["homo"]["value"] == pytest.approx(-0.40932048)
            assert properties["homo"]["unit"] is None
            assert restored.homo == pytest.approx(-0.40932048)
            assert restored.atomization == pytest.approx(-0.66538218)
            assert restored.dataset_source == "revQM9"
            assert properties["homo"]["unit_status"] == "inferred"
            assert properties["dipole"]["suggested_unit"] == "Debye"
            assert record["molecules"][0][
                "dipole_moment_magnitude"
            ] == pytest.approx(1.2762092032023672e-05)
            assert record["results"]["lumo_energy"] == pytest.approx(
                0.06115874826417237 * Hartree
            )
            assert properties["atomization"]["value"] == pytest.approx(
                -0.66538218
            )
        else:
            assert rows[1][1] is None
            assert rows[1][2] == pytest.approx(-10.54985436)
            assert properties["U0"]["value"] == pytest.approx(-17.17218053)
            assert properties["U0"]["unit"] == "eV"
            assert restored.U0 == pytest.approx(-17.17218053)
            assert restored.dataset_source == "qm9_eV"

    def test_revqm9_subset_fixture(self):
        # First five entries from the official revQM9 dataset.
        path = Path(__file__).parent / "data" / "IOTests" / "revQM9_subset.npz"
        molecules = Molecule.from_filepath(path, index=":", return_list=True)
        assert len(molecules) == 5
        assert [m.empirical_formula for m in molecules[:3]] == [
            "CH4",
            "H3N",
            "H2O",
        ]
        assert molecules[0].energy == pytest.approx(-40.47574250955533)
        assert molecules[0].homo == pytest.approx(-0.40932048)
        assert molecules[0].atomization == pytest.approx(-0.66538218)
        assert molecules[2].structure_index_in_file == 3
        selected = Molecule.from_filepath(path, index="3,1")
        assert [m.empirical_formula for m in selected] == ["H2O", "CH4"]
        assert selected[1].energy == pytest.approx(molecules[0].energy)
        assert selected[1].homo == molecules[0].homo

    def test_dgl_qm9_ev_subset_fixture(self):
        # First five molecules from DGL's qm9_eV.npz, with original arrays.
        path = Path(__file__).parent / "data" / "IOTests" / "qm9_eV_subset.npz"
        molecules = Molecule.from_filepath(path, index=":", return_list=True)
        assert len(molecules) == 5
        assert [m.num_atoms for m in molecules] == [5, 4, 3, 4, 3]
        assert [m.empirical_formula for m in molecules[:3]] == [
            "CH4",
            "H3N",
            "H2O",
        ]
        np.testing.assert_allclose(
            molecules[0].positions[0],
            [-0.0126981359, 1.0858041578, 0.0080009958],
        )
        assert molecules[0].structure_index_in_file == 1
        assert molecules[0].energy is None
        assert molecules[0].U0 == pytest.approx(-17.17218053)
        selected = Molecule.from_filepath(path, index="3,1")
        assert [m.empirical_formula for m in selected] == ["H2O", "CH4"]

    def test_rejects_misaligned_properties(self, tmp_path):
        path = tmp_path / "invalid_properties.npz"
        np.savez(
            path, coords=np.zeros((2, 1, 3)), charges=[[1], [8]], homo=[-0.2]
        )
        with pytest.raises(ValueError, match="homo.*one value per molecule"):
            Molecule.from_filepath(path, index=":")

    def test_revqm9_ragged_arrays_and_selection(self, tmp_path):
        path = tmp_path / "revQM9.npz"
        coords = _write_revqm9(path)
        molecules = Molecule.from_filepath(path, index=":", return_list=True)
        assert [list(m.symbols) for m in molecules] == [["He"], ["H", "H"]]
        assert [m.structure_index_in_file for m in molecules] == [1, 2]
        assert molecules[0].charge == 0 and molecules[0].multiplicity == 1
        assert molecules[0].energy is None
        np.testing.assert_array_equal(molecules[1].positions, coords[1])
        assert list(Molecule.from_filepath(path, index="-1").symbols) == [
            "H",
            "H",
        ]
        assert [
            m.structure_index_in_file
            for m in Molecule.from_filepath(path, index="2,1")
        ] == [2, 1]
        assert len(Molecule.from_filepath(path, index="1:2")) == 1
        assert (
            len(Molecule.from_filepath(path, index="1", return_list=True)) == 1
        )

    def test_dgl_concatenated_molecules_and_selection(self, tmp_path):
        path = tmp_path / "qm9_eV.npz"
        coords = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]])
        np.savez(path, R=coords, Z=[2, 1, 1, 8], N=[1, 2, 1])
        molecules = Molecule.from_filepath(path, index=":", return_list=True)
        assert [list(m.symbols) for m in molecules] == [
            ["He"],
            ["H", "H"],
            ["O"],
        ]
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

    def test_revqm9_reads_total_energy_for_selected_molecules(self, tmp_path):
        path = tmp_path / "revqm9.npz"
        np.savez(
            path,
            coords=np.zeros((2, 1, 3)),
            charges=[[1], [8]],
            energies=[-1.2, -3.4],
        )
        molecules = Molecule.from_filepath(path, index="2,1")
        assert [m.energy for m in molecules] == pytest.approx([-3.4, -1.2])

    def test_dgl_u0_is_not_assigned_as_total_energy(self, tmp_path):
        path = tmp_path / "dgl.npz"
        np.savez(path, R=np.zeros((1, 3)), Z=[1], N=[1], U0=[-1.2])
        assert Molecule.from_filepath(path).energy is None

    @pytest.mark.parametrize(
        "energies", [[1.0], [1.0, float("nan")], [[1.0], [2.0]]]
    )
    def test_npz_rejects_invalid_energy_values(self, tmp_path, energies):
        path = tmp_path / "invalid_energy.npz"
        np.savez(
            path,
            coords=np.zeros((2, 1, 3)),
            charges=[[1], [8]],
            energies=energies,
        )
        with pytest.raises(ValueError, match="energies"):
            Molecule.from_filepath(path, index=":")

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
    def test_dgl_rejects_invalid_data(
        self, tmp_path, counts, numbers, message
    ):
        path = tmp_path / "invalid.npz"
        np.savez(path, R=np.zeros((3, 3)), Z=numbers, N=counts)
        with pytest.raises(ValueError, match=message):
            Molecule.from_filepath(path, index=":")

    def test_npz_rejects_unsupported_layout_and_nonfinite_coordinates(
        self, tmp_path
    ):
        path = tmp_path / "bad.npz"
        np.savez(path, positions=[[0, 0, 0]], atomic_numbers=[2])
        with pytest.raises(ValueError, match="DGL QM9 fields"):
            Molecule.from_filepath(path)
        np.savez(path, R=[[float("nan"), 0, 0]], Z=[2], N=[1])
        with pytest.raises(ValueError, match="finite"):
            Molecule.from_filepath(path)

    @pytest.mark.parametrize(
        "coords, charges, message",
        [
            (np.zeros((2, 1, 3)), [[1]], "matching molecule counts"),
            (np.empty((0, 1, 3)), np.empty((0, 1)), "no molecular structures"),
            (np.zeros((1, 1, 2)), [[1]], r"coordinates must have shape"),
            (np.zeros((1, 2, 3)), [[1]], "atom and coordinate counts differ"),
            (np.zeros((1, 1, 3)), [[-1]], "invalid atomic numbers"),
            (np.zeros((1, 1, 3)), [[1.5]], "invalid atomic numbers"),
            (np.zeros((1, 1, 3)), [[119]], "invalid atomic numbers"),
            (np.zeros((1, 1, 3)), [[float("nan")]], "invalid atomic numbers"),
        ],
    )
    def test_revqm9_rejects_invalid_structures(
        self, tmp_path, coords, charges, message
    ):
        path = tmp_path / "invalid.npz"
        np.savez(path, coords=coords, charges=charges)
        with pytest.raises(ValueError, match=message):
            Molecule.from_filepath(path, index=":")

    def test_revqm9_rejects_invalid_element_symbol(self, tmp_path):
        path = tmp_path / "invalid_element.npz"
        np.savez(path, coords=[[0, 0, 0]], charges=[1], elements=["Invalid"])
        with pytest.raises(ValueError, match="invalid elements"):
            Molecule.from_filepath(path)

    @pytest.mark.parametrize(
        "extra", [{"elements": ["C"]}, {"atomic_numbers": [6]}]
    )
    def test_npz_rejects_conflicting_atom_identity(self, tmp_path, extra):
        path = tmp_path / "conflict.npz"
        np.savez(path, coords=[[0, 0, 0]], charges=[8], **extra)
        with pytest.raises(ValueError, match="conflicting atomic identities"):
            Molecule.from_filepath(path)

    def test_npz_rejects_mismatched_atom_numbers_and_elements(self, tmp_path):
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

    def test_npz_rejects_padded_atom_arrays(self, tmp_path):
        path = tmp_path / "padded.npz"
        np.savez(
            path, R=np.zeros((2, 3, 3)), Z=[[6, 1, 0], [8, 0, 0]], N=[2, 1]
        )
        with pytest.raises(ValueError, match="zero-padded"):
            Molecule.from_filepath(path, index=":")

    def test_npz_rejects_arbitrary_pickle_globals(self, tmp_path):
        marker = tmp_path / "unpickled"

        class Payload:
            def __reduce__(self):
                return os.mkdir, (str(marker),)

        path = tmp_path / "untrusted.npz"
        np.savez(
            path, coords=np.array([Payload()], dtype=object), charges=[[1]]
        )
        with pytest.raises(ValueError, match="Disallowed pickle global"):
            Molecule.from_filepath(path)
        assert not marker.exists()

    def test_npz_rejects_pickled_array_header_mismatch(self, tmp_path):
        path = tmp_path / "mismatch.npz"
        with ZipFile(path, "w") as archive:
            with archive.open("coords.npy", "w") as stream:
                write_array_header_1_0(
                    stream,
                    {"descr": "|O", "fortran_order": False, "shape": (2,)},
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
        self,
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
