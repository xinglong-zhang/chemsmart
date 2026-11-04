"""DGL QM9 and revQM9 geometry input."""

import logging
import os
import pickle
import shutil
import sqlite3
from importlib import import_module
from unittest.mock import MagicMock, patch
from zipfile import ZipFile

import numpy as np
import pytest
from ase.units import Hartree
from click.testing import CliRunner
from numpy.lib.format import write_array_header_1_0

from chemsmart.cli.database.assemble import assemble
from chemsmart.database.assemble import NPZAssembler
from chemsmart.database.database import Database
from chemsmart.database.inspect import DatabaseInspector
from chemsmart.io.molecules.structure import Molecule
from chemsmart.io.npz import NPZFile
from tests.helpers import make_cli_ctx_obj, write_revqm9_npz


class TestNPZStructures:
    def test_failed_npz_insert_rolls_back_record(
        self, tmp_path, io_test_directory
    ):
        source = os.path.join(io_test_directory, "revQM9_subset.npz")
        db = Database(str(tmp_path / "records.db"))
        db.create()
        with db.get_connection() as conn:
            conn.execute("DROP TABLE dataset_properties")
        records = NPZAssembler(
            str(source), index="1", dataset="revqm9"
        ).assemble()
        assert db.insert_records(records) == 0
        assert db.count_records() == 0

    def test_failed_npz_update_preserves_existing_record(
        self, tmp_path, io_test_directory
    ):
        source = os.path.join(io_test_directory, "revQM9_subset.npz")
        db = Database(str(tmp_path / "records.db"))
        db.create()
        assert (
            db.insert_records(
                NPZAssembler(
                    str(source), index="1", dataset="revqm9"
                ).assemble()
            )
            == 1
        )
        with db.get_connection() as conn:
            conn.execute("DROP TABLE dataset_properties")
        records = NPZAssembler(
            str(source), index="1", dataset="revqm9", method="PBE0"
        ).assemble()
        assert db.insert_records(records) == 0
        record = db.get_record(record_index=1)
        assert record["meta"]["method"] == "aPBE0"
        assert len(record["molecules"]) == 1

    def test_npz_reimport_preserves_record_index_and_declared_data(
        self, tmp_path, io_test_directory
    ):
        source = os.path.join(io_test_directory, "revQM9_subset.npz")
        output = tmp_path / "records.db"
        base = ["-f", str(source), "-i", "1", "-o", str(output)]

        def run(extra):
            result = CliRunner().invoke(assemble, base + extra, obj={})
            assert result.exit_code == 0, (result.output, result.exception)

        run([])
        db = Database(str(output))
        record = db.get_record(record_index=1)
        assert record["meta"]["method"] is None
        assert record["results"]["total_energy"] == pytest.approx(
            -40.47574250955533
        )
        restored = Molecule.from_filepath(str(output), record_index=1)
        assert restored.energy == pytest.approx(-40.47574250955533)
        run(["--dataset", "revqm9"])
        record = db.get_record(record_index=1)
        assert record["meta"]["method"] == "aPBE0"
        homo = record["results"]["homo_energy"]
        assert homo is not None
        run([])
        record = db.get_record(record_index=1)
        assert db.count_records() == 1
        assert record["meta"]["method"] == "aPBE0"
        assert record["results"]["homo_energy"] == homo
        run(["--dataset", "revqm9", "--method", "PBE0"])
        assert db.get_record(record_index=1)["meta"]["method"] == "PBE0"

    def test_npz_database_cli_metadata_override(
        self, tmp_path, io_test_directory
    ):
        source = os.path.join(io_test_directory, "revQM9_subset.npz")
        output = tmp_path / "custom.db"
        result = CliRunner().invoke(
            assemble,
            [
                "-f",
                str(source),
                "--dataset",
                "revqm9",
                "--method",
                "PBE0",
                "--basis",
                "def2-SVP",
                "-i",
                "1",
                "-o",
                str(output),
            ],
            obj={},
        )
        assert result.exit_code == 0, (result.output, result.exception)
        record = Database(str(output)).get_record(record_index=1)
        assert record["meta"]["method"] == "PBE0"
        assert record["meta"]["basis"] == "def2-SVP"

    @pytest.mark.parametrize(
        "dataset,method,basis,expected_source,expected_method,expected_basis",
        [
            (None, None, None, None, None, None),
            ("revqm9", None, None, "revQM9", "aPBE0", "cc-pVTZ"),
            ("qm9", None, None, "qm9_eV", "B3LYP", "6-31G(2df,p)"),
            ("qm9", "PBE0", "def2-SVP", "qm9_eV", "PBE0", "def2-SVP"),
            (
                "revqm9",
                "custom",
                "custom-basis",
                "revQM9",
                "custom",
                "custom-basis",
            ),
        ],
    )
    def test_npz_database_metadata(
        self,
        tmp_path,
        dataset,
        method,
        basis,
        expected_source,
        expected_method,
        expected_basis,
        io_test_directory,
    ):
        filename = (
            "qm9_eV_subset.npz" if dataset == "qm9" else "revQM9_subset.npz"
        )
        source = os.path.join(io_test_directory, filename)
        record = next(
            NPZAssembler(
                source, index="1", dataset=dataset, method=method, basis=basis
            ).assemble()
        )
        assert record.provenance["program"] == expected_source
        assert record.meta["method"] == expected_method
        assert record.meta["basis"] == expected_basis
        if dataset is None:
            properties = record.molecules[0]["dataset_properties"]
            assert properties["homo"]["unit"] is None

    def test_npz_rejects_mismatched_dataset(self, io_test_directory):
        source = os.path.join(io_test_directory, "revQM9_subset.npz")
        with pytest.raises(ValueError, match="layout does not match"):
            list(NPZAssembler(source, index="1", dataset="qm9").assemble())

    @pytest.mark.parametrize(
        "filename", ["revQM9_subset.npz", "qm9_eV_subset.npz"]
    )
    def test_npz_database_import(self, tmp_path, filename, io_test_directory):
        source = os.path.join(io_test_directory, filename)
        output = tmp_path / "dataset.db"
        dataset = "revqm9" if filename.startswith("rev") else "qm9"
        args = [
            "-f",
            str(source),
            "--dataset",
            dataset,
            "-i",
            "3,1",
            "-o",
            str(output),
        ]
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
            assert restored.energy == pytest.approx(-40.47574250955533)
            assert rows[1][2] == pytest.approx(-0.4093204819908599 * Hartree)
            assert properties["homo"]["value"] == pytest.approx(-0.40932048)
            assert properties["homo"]["unit"] == "Hartree"
            assert restored.homo == pytest.approx(-0.40932048)
            assert restored.atomization == pytest.approx(-0.66538218)
            assert restored.dataset_source == "revQM9"
            assert properties["dipole"]["unit"] == "Debye"
            assert properties["atomization"]["unit"] == "Hartree"
            assert properties["lumo"]["unit"] == "Hartree"
            assert record["molecules"][0][
                "dipole_moment_magnitude"
            ] == pytest.approx(1.2762092032023672e-05)
            assert record["results"]["lumo_energy"] == pytest.approx(
                0.06115874826417237 * Hartree
            )
            detail = DatabaseInspector(
                str(output), index=2
            ).format_record_detail()
            assert "HOMO Energy (eV)" in detail
            assert "LUMO Energy (eV)" in detail
            assert "-11.138" in detail
            assert "1.664" in detail
            assert "α-HOMO Energy" not in detail
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

    def test_revqm9_subset_fixture(self, io_test_directory):
        # First five entries from the official revQM9 dataset.
        path = os.path.join(io_test_directory, "revQM9_subset.npz")
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

    def test_dgl_qm9_ev_subset_fixture(self, io_test_directory):
        # First five molecules from DGL's qm9_eV.npz, with original arrays.
        path = os.path.join(io_test_directory, "qm9_eV_subset.npz")
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
        coords = write_revqm9_npz(path)
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

    def test_npz_rejects_padded_qm9_arrays(self, tmp_path):
        path = tmp_path / "padded.npz"
        np.savez(
            path, R=np.zeros((2, 3, 3)), Z=[[6, 1, 0], [8, 0, 0]], N=[2, 1]
        )
        with pytest.raises(ValueError, match="must concatenate all atoms"):
            Molecule.from_filepath(path, index=":")

    def test_npz_rejects_padded_qm9_arrays_that_contradict_their_counts(
        self, tmp_path
    ):
        """No zeros reveal the padding, and N disagrees with the atom count."""
        path = tmp_path / "padded.npz"
        np.savez(path, R=np.zeros((2, 2, 3)), Z=[[1, 1], [8, 8]], N=[1, 2])
        with pytest.raises(ValueError, match="must concatenate all atoms"):
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

    @pytest.mark.parametrize("selection", ["2", "-1"])
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
        selection,
        gaussian_jobrunner_no_scratch,
    ):
        path = tmp_path / "revQM9.npz"
        write_revqm9_npz(path)
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
            ) as reader,
        ):
            ctx_obj = make_cli_ctx_obj(gaussian_jobrunner_no_scratch)
            result = CliRunner().invoke(
                command,
                [
                    "-p",
                    project,
                    "-f",
                    str(path),
                    "-i",
                    selection,
                    "-c",
                    "0",
                    "-m",
                    "1",
                    subcommand,
                ],
                obj=ctx_obj,
            )
        assert result.exit_code == 0, (result.output, result.exception)
        job.assert_called_once()
        assert list(job.call_args.kwargs["molecule"].symbols) == ["H", "H"]
        if backend != "crest":
            assert reader.call_args.kwargs["index"] == selection
            assert ctx_obj["molecule_indices"] == [2]

    @pytest.mark.parametrize(
        "backend,project,job_path",
        [
            (
                "gaussian",
                "gas_solv",
                "chemsmart.jobs.gaussian.opt.GaussianOptJob",
            ),
            ("orca", "gas_solv", "chemsmart.jobs.orca.opt.ORCAOptJob"),
            ("xtb", "test", "chemsmart.jobs.xtb.opt.XTBOptJob"),
        ],
    )
    def test_cli_npz_preserves_selected_order(
        self,
        tmp_path,
        backend,
        project,
        job_path,
        gaussian_jobrunner_no_scratch,
    ):
        path = tmp_path / "revQM9.npz"
        write_revqm9_npz(path)
        command = getattr(
            import_module(f"chemsmart.cli.{backend}.{backend}"), backend
        )
        ctx_obj = make_cli_ctx_obj(gaussian_jobrunner_no_scratch)
        with patch(job_path, return_value=MagicMock()) as job:
            result = CliRunner().invoke(
                command,
                [
                    "-p",
                    project,
                    "-f",
                    str(path),
                    "-i",
                    "2,1",
                    "-c",
                    "0",
                    "-m",
                    "1",
                    "opt",
                ],
                obj=ctx_obj,
            )
        assert result.exit_code == 0, (result.output, result.exception)
        assert ctx_obj["molecule_indices"] == [2, 1]
        assert job.call_count == 2

    @pytest.mark.parametrize(
        "extra",
        [["-d", "."], ["-t", "log"], ["-p", "orca"], ["--include-failed"]],
    )
    def test_npz_import_rejects_directory_options(self, tmp_path, extra):
        source = tmp_path / "input.npz"
        write_revqm9_npz(source)
        output = tmp_path / "output.db"
        result = CliRunner().invoke(
            assemble, ["-f", str(source), "-o", str(output), *extra], obj={}
        )
        assert result.exit_code != 0
        assert "cannot be combined" in result.output
        assert not output.exists()

    def test_dipole_magnitude_does_not_need_a_declared_dataset(
        self, io_test_directory
    ):
        source = os.path.join(io_test_directory, "qm9_eV_subset.npz")
        for dataset in ("qm9", None):
            record = next(
                NPZAssembler(source, index="2", dataset=dataset).assemble()
            )
            molecule = record.molecules[0]
            assert molecule["dipole_moment_magnitude"] == pytest.approx(1.6256)


class TestNPZImportRobustness:
    """Long imports must keep their progress and report what they skipped."""

    def _flaky(self, fails_on, error=RuntimeError("perception failed")):
        original = NPZAssembler.get_molecule_info
        calls = {"count": 0}

        def get_molecule_info(assembler, molecule):
            calls["count"] += 1
            if calls["count"] == fails_on:
                raise error
            return original(assembler, molecule)

        return patch.object(
            NPZAssembler, "get_molecule_info", get_molecule_info
        )

    def test_failed_structure_is_skipped_and_reported(
        self, tmp_path, io_test_directory
    ):
        source = os.path.join(io_test_directory, "qm9_eV_subset.npz")
        db = Database(str(tmp_path / "records.db"))
        db.create()
        assembler = NPZAssembler(source, dataset="qm9")
        with self._flaky(fails_on=4):
            assert db.insert_records(assembler.assemble()) == 4
        assert db.count_records() == 4
        assert assembler.failed_structures == {4: "perception failed"}

    def test_missing_dependency_stops_the_import(
        self, tmp_path, io_test_directory
    ):
        """It would fail every structure, so skipping each hides the cause."""
        source = os.path.join(io_test_directory, "qm9_eV_subset.npz")
        db = Database(str(tmp_path / "records.db"))
        db.create()
        assembler = NPZAssembler(source, dataset="qm9")
        with self._flaky(fails_on=2, error=ImportError("needs Open Babel")):
            with pytest.raises(ImportError, match="Open Babel"):
                db.insert_records(assembler.assemble())
        assert db.count_records() == 1

    def test_failing_record_source_keeps_earlier_records(
        self, tmp_path, io_test_directory
    ):
        source = os.path.join(io_test_directory, "qm9_eV_subset.npz")
        db = Database(str(tmp_path / "records.db"))
        db.create()
        records = NPZAssembler(source, dataset="qm9").assemble()

        def stops_early():
            yield next(records)
            yield next(records)
            raise RuntimeError("source stopped")

        with pytest.raises(RuntimeError, match="source stopped"):
            db.insert_records(stops_early())
        assert db.count_records() == 2

    def test_progress_is_committed_in_batches(
        self, tmp_path, io_test_directory
    ):
        source = os.path.join(io_test_directory, "qm9_eV_subset.npz")
        output = tmp_path / "records.db"
        db = Database(str(output))
        db.create()
        records = NPZAssembler(source, dataset="qm9").assemble()
        seen_by_another_connection = []

        def observing():
            for record in records:
                with sqlite3.connect(output) as other:
                    seen_by_another_connection.append(
                        other.execute(
                            "SELECT COUNT(*) FROM records"
                        ).fetchone()[0]
                    )
                yield record

        with patch("chemsmart.database.database.COMMIT_BATCH_SIZE", 2):
            db.insert_records(observing())
        assert seen_by_another_connection == [0, 0, 2, 2, 4]

    def test_undeclared_import_warns_about_assumed_units(
        self, tmp_path, io_test_directory, caplog
    ):
        source = os.path.join(io_test_directory, "revQM9_subset.npz")
        with caplog.at_level(logging.WARNING):
            list(NPZAssembler(source, index="1").assemble())
        assert "--dataset revqm9" in caplog.text

    def test_declared_import_does_not_warn(self, io_test_directory, caplog):
        source = os.path.join(io_test_directory, "revQM9_subset.npz")
        with caplog.at_level(logging.WARNING):
            list(NPZAssembler(source, index="1", dataset="revqm9").assemble())
        assert "--dataset" not in caplog.text

    def test_cli_stores_an_absolute_source_path(
        self, tmp_path, io_test_directory, monkeypatch
    ):
        shutil.copy(
            os.path.join(io_test_directory, "qm9_eV_subset.npz"),
            tmp_path / "input.npz",
        )
        monkeypatch.chdir(tmp_path)
        result = CliRunner().invoke(
            assemble, ["-f", "input.npz", "-i", "1", "-o", "x.db"], obj={}
        )
        assert result.exit_code == 0, (result.output, result.exception)
        record = Database(str(tmp_path / "x.db")).get_record(record_index=1)
        assert os.path.isabs(record["provenance"]["source"])

    def test_cli_reports_skipped_structures(
        self, tmp_path, io_test_directory, caplog
    ):
        source = os.path.join(io_test_directory, "qm9_eV_subset.npz")
        with self._flaky(fails_on=2), caplog.at_level(logging.WARNING):
            result = CliRunner().invoke(
                assemble, ["-f", source, "-o", str(tmp_path / "x.db")], obj={}
            )
        assert result.exit_code == 0, (result.output, result.exception)
        assert "Skipped 1 structure(s)" in caplog.text


class TestNPZDatasetSourceRoundTrip:
    @pytest.mark.parametrize(
        ("filename", "dataset", "source"),
        [
            ("revQM9_subset.npz", None, "revQM9"),
            ("revQM9_subset.npz", "revqm9", "revQM9"),
            ("qm9_eV_subset.npz", None, "qm9_eV"),
            ("qm9_eV_subset.npz", "qm9", "qm9_eV"),
        ],
    )
    def test_source_survives_the_database(
        self, tmp_path, io_test_directory, filename, dataset, source
    ):
        path = os.path.join(io_test_directory, filename)
        output = tmp_path / "records.db"
        db = Database(str(output))
        db.create()
        db.insert_records(
            NPZAssembler(path, index="1", dataset=dataset).assemble()
        )
        direct = Molecule.from_filepath(path, index="1", dataset=dataset)
        restored = Molecule.from_filepath(str(output), record_index=1)
        assert direct.dataset_source == restored.dataset_source == source

    def test_coordinates_only_archive_keeps_its_source(self, tmp_path):
        path = tmp_path / "input.npz"
        np.savez(path, coords=np.zeros((1, 1, 3)), charges=[[2]])
        output = tmp_path / "records.db"
        db = Database(str(output))
        db.create()
        assert db.insert_records(NPZAssembler(str(path)).assemble()) == 1
        restored = Molecule.from_filepath(str(output), record_index=1)
        assert restored.dataset_source == "revQM9"


class TestNPZSelection:
    @pytest.mark.parametrize("index", ["10:12", "99", "1,99"])
    def test_selection_outside_the_file_is_explained(
        self, io_test_directory, index
    ):
        path = os.path.join(io_test_directory, "revQM9_subset.npz")
        with pytest.raises(IndexError, match="selects no structures"):
            Molecule.from_filepath(path, index=index)

    def test_selection_error_names_the_file_size(self, io_test_directory):
        path = os.path.join(io_test_directory, "revQM9_subset.npz")
        with pytest.raises(IndexError, match="which has 5"):
            Molecule.from_filepath(path, index="99")

    def test_indices_are_the_positions_in_the_file(self, io_test_directory):
        path = os.path.join(io_test_directory, "qm9_eV_subset.npz")
        molecules, indices = Molecule.from_filepath_with_indices(path, "3,1")
        assert indices == [3, 1]
        assert [m.chemical_formula for m in molecules] == ["H2O", "CH4"]

    def test_no_index_reports_no_positions(self, io_test_directory):
        path = os.path.join(io_test_directory, "qm9_eV_subset.npz")
        molecules, indices = Molecule.from_filepath_with_indices(path)
        assert len(molecules) == 5
        assert indices is None

    def test_other_formats_report_positions_too(
        self, multiple_molecules_xyz_file
    ):
        molecules, indices = Molecule.from_filepath_with_indices(
            multiple_molecules_xyz_file, "2"
        )
        assert len(molecules) == 1
        assert indices == [2]


class TestNPZEnergiesWithoutADeclaredDataset:
    """Energies are recorded either way, in each column's fixed unit."""

    def _record(self, io_test_directory, filename, dataset=None):
        source = os.path.join(io_test_directory, filename)
        return next(
            NPZAssembler(source, index="1", dataset=dataset).assemble()
        )

    def test_revqm9_energies_match_a_declared_import(self, io_test_directory):
        undeclared = self._record(io_test_directory, "revQM9_subset.npz")
        declared = self._record(
            io_test_directory, "revQM9_subset.npz", "revqm9"
        )
        assert undeclared.results == declared.results
        assert undeclared.molecules[0]["energy"] == pytest.approx(
            -40.47574250955533
        )

    def test_revqm9_total_energy_stays_in_hartree_and_orbitals_in_ev(
        self, io_test_directory
    ):
        results = self._record(io_test_directory, "revQM9_subset.npz").results
        assert results["total_energy"] == pytest.approx(-40.47574250955533)
        assert results["homo_energy"] == pytest.approx(
            -0.4093204819908599 * Hartree
        )
        assert results["fmo_gap"] == pytest.approx(
            (0.06115874826417237 + 0.4093204819908599) * Hartree
        )

    def test_qm9_orbital_energies_are_recorded(self, io_test_directory):
        results = self._record(io_test_directory, "qm9_eV_subset.npz").results
        assert results["homo_energy"] == pytest.approx(-10.54985436)
        # U0 is atom-reference-adjusted, so it is not a total energy.
        assert results["total_energy"] is None

    def test_declaring_still_adds_what_it_alone_can_state(
        self, io_test_directory
    ):
        undeclared = self._record(io_test_directory, "revQM9_subset.npz")
        declared = self._record(
            io_test_directory, "revQM9_subset.npz", "revqm9"
        )
        assert undeclared.provenance["program"] is None
        assert undeclared.meta["method"] is None
        assert (
            undeclared.molecules[0]["dataset_properties"]["homo"]["unit"]
            is None
        )
        assert declared.provenance["program"] == "revQM9"
        assert declared.meta["method"] == "aPBE0"
        assert (
            declared.molecules[0]["dataset_properties"]["homo"]["unit"]
            == "Hartree"
        )

    def test_values_reach_the_database_and_export(
        self, tmp_path, io_test_directory
    ):
        source = os.path.join(io_test_directory, "revQM9_subset.npz")
        output = tmp_path / "records.db"
        result = CliRunner().invoke(
            assemble, ["-f", source, "-o", str(output)], obj={}
        )
        assert result.exit_code == 0, (result.output, result.exception)
        record = Database(str(output)).get_record(record_index=1)
        assert record["results"]["total_energy"] == pytest.approx(
            -40.47574250955533
        )
        assert record["results"]["lumo_energy"] == pytest.approx(
            0.06115874826417237 * Hartree
        )


class TestNPZValuesInDatabaseColumns:
    """Calculated values go in the column with the same meaning and unit."""

    def _stored(self, tmp_path, io_test_directory, filename, *options):
        source = os.path.join(io_test_directory, filename)
        output = tmp_path / "records.db"
        result = CliRunner().invoke(
            assemble,
            ["-f", source, "-i", "1", "-o", str(output), *options],
            obj={},
        )
        assert result.exit_code == 0, (result.output, result.exception)
        return Database(str(output)).get_record(record_index=1)

    @pytest.mark.parametrize("options", [[], ["--dataset", "qm9"]])
    def test_qm9_zero_point_energy_is_in_hartree(
        self, tmp_path, io_test_directory, options
    ):
        record = self._stored(
            tmp_path, io_test_directory, "qm9_eV_subset.npz", *options
        )
        # Methane's zero-point energy is about 28 kcal/mol.
        assert record["results"]["zero_point_energy"] == pytest.approx(
            0.04475, abs=1e-5
        )

    @pytest.mark.parametrize("options", [[], ["--dataset", "qm9"]])
    def test_qm9_rotational_constants_are_in_hz(
        self, tmp_path, io_test_directory, options
    ):
        record = self._stored(
            tmp_path, io_test_directory, "qm9_eV_subset.npz", *options
        )
        expected = [157.7118e9, 157.70997e9, 157.70699e9]
        assert record["results"][
            "rotational_constants_in_Hz"
        ] == pytest.approx(expected, rel=1e-5)
        assert record["molecules"][0]["rotational_constants"] == pytest.approx(
            expected, rel=1e-5
        )

    def test_atom_referenced_energies_stay_out_of_absolute_columns(
        self, tmp_path, io_test_directory
    ):
        record = self._stored(
            tmp_path,
            io_test_directory,
            "qm9_eV_subset.npz",
            "--dataset",
            "qm9",
        )
        for column in ("internal_energy", "enthalpy", "gibbs_free_energy"):
            assert record["results"].get(column) is None
        properties = record["molecules"][0]["dataset_properties"]
        assert properties["U0"]["value"] == pytest.approx(-17.17218053)
        assert properties["G"]["value"] == pytest.approx(-16.15192, abs=1e-4)

    @pytest.mark.parametrize(
        ("filename", "fields"),
        [
            (
                "qm9_eV_subset.npz",
                {
                    "id",
                    "A",
                    "B",
                    "C",
                    "mu",
                    "alpha",
                    "homo",
                    "lumo",
                    "gap",
                    "r2",
                    "zpve",
                    "U0",
                    "U",
                    "H",
                    "G",
                    "Cv",
                },
            ),
            (
                "revQM9_subset.npz",
                {"inds", "atomization", "homo", "lumo", "dipole"},
            ),
        ],
    )
    def test_every_per_molecule_value_is_kept(
        self, tmp_path, io_test_directory, filename, fields
    ):
        record = self._stored(tmp_path, io_test_directory, filename)
        assert set(record["molecules"][0]["dataset_properties"]) == fields

    def test_revqm9_has_no_values_for_qm9_only_columns(
        self, tmp_path, io_test_directory
    ):
        record = self._stored(
            tmp_path,
            io_test_directory,
            "revQM9_subset.npz",
            "--dataset",
            "revqm9",
        )
        assert record["results"].get("zero_point_energy") is None
        assert record["results"].get("rotational_constants_in_Hz") is None
        assert record["molecules"][0][
            "dipole_moment_magnitude"
        ] == pytest.approx(1.2762092032023672e-05)
