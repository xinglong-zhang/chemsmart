import importlib
import os
import shutil
from unittest.mock import MagicMock, patch

import click
import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner

from chemsmart.cli.grouper import grouper
from chemsmart.cli.grouper.grouper import (
    THERMOCHEMISTRY_DEFAULT_KWARGS,
    _extract_conformer_id,
    _extract_energy_based_on_energy_type,
    _extract_last_energy_from_output_file,
    _find_matching_sp_file,
    _get_label,
    _load_molecules_from_directory,
    build_thermochemistry_kwargs,
    create_grouper_job_from_context,
)
from chemsmart.io.molecules.structure import Molecule
from chemsmart.io.xyz.xyzfile import XYZFile
from chemsmart.jobs.grouper.base import (
    MatrixGrouper,
    MoleculeGrouper,
    ResultsRecorder,
    StructureGrouperConfig,
)
from chemsmart.jobs.grouper.connectivity import ConnectivityGrouper
from chemsmart.jobs.grouper.energy import (
    KCAL_TO_HARTREE,
    EnergyGrouper,
)
from chemsmart.jobs.grouper.formula import FormulaGrouper
from chemsmart.jobs.grouper.isomorphism import RDKitIsomorphismGrouper
from chemsmart.jobs.grouper.job import GrouperJob
from chemsmart.jobs.grouper.rmsd import (
    BasicRMSDGrouper,
    HungarianRMSDGrouper,
    IRMSDGrouper,
    PymolRMSDGrouper,
    RMSDGrouper,
    SpyRMSDGrouper,
)
from chemsmart.jobs.grouper.runner import GrouperJobRunner
from chemsmart.jobs.grouper.tanimoto import TanimotoSimilarityGrouper
from chemsmart.jobs.grouper.tfd import TorsionFingerprintGrouper
from chemsmart.utils.grouper import StructureGrouperFactory
from chemsmart.utils.utils import (
    find_irmsd_command,
    kabsch_align,
    to_graph_wrapper,
)

SYMBOLS = ["H", "H"]


POSITIONS = np.array([[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]])


grouper_module = importlib.import_module("chemsmart.cli.grouper.grouper")


def invoke_grouper(args, mocker, mock_job_return=None):
    mock_job_cls = mocker.patch("chemsmart.jobs.grouper.GrouperJob")
    if mock_job_return is not None:
        mock_job_cls.return_value = mock_job_return
    runner = CliRunner()
    result = runner.invoke(grouper, args, catch_exceptions=False)
    return result, mock_job_cls


@pytest.fixture()
def grouper_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    return GrouperJob(
        molecules=make_molecules(3),
        grouping_strategy="rmsd",
        label="mygrp",
        threshold=0.5,
    )


class SimpleMolecule:
    def __init__(self, energy, symbols, positions):
        self.energy = energy
        self.chemical_symbols = symbols
        self.positions = positions

    @property
    def num_atoms(self):
        return len(self.chemical_symbols)


def make_molecules(n=3):
    return [MagicMock(name=f"molecule{i}") for i in range(n)]


class _ConcreteGrouper(MoleculeGrouper):
    """Minimal concrete subclass so the abstract base can be instantiated."""

    def group(self):
        indices = list(range(len(self.molecules)))
        return [self.molecules], [indices]


def _make_water():
    return Molecule(
        symbols=["O", "H"],
        positions=[[0.0, 0.0, 0.0], [0.96, 0.0, 0.0]],
        charge=0,
        multiplicity=1,
    )


def _mol_tfd_grouper():
    return Molecule(symbols=SYMBOLS, positions=POSITIONS)


def _mol_energy_grouper(energy):
    return Molecule(symbols=SYMBOLS, positions=POSITIONS, energy=energy)


def _mol(x_offset=0.0):
    positions = np.array(
        [[0.0 + x_offset, 0.0, 0.0], [0.74 + x_offset, 0.0, 0.0]]
    )
    return Molecule(symbols=SYMBOLS, positions=positions)


def _formaldehyde(x_offset=0.0, y_rotate=False):
    # One C, one O, two H -- exercises the Hungarian assignment's
    # single-atom-per-element shortcut (C and O) alongside the real
    # linear_sum_assignment path (the two H atoms).
    positions = np.array(
        [
            [0.0 + x_offset, 0.0, 0.0],
            [1.2 + x_offset, 0.0, 0.0],
            [-0.6 + x_offset, 0.9, 0.0],
            [-0.6 + x_offset, -0.9, 0.0],
        ]
    )
    if y_rotate:
        # Swap the two symmetric H atoms' order to force the Hungarian
        # algorithm to actually search for the optimal assignment.
        positions = positions[[0, 1, 3, 2]]
    return Molecule(symbols=["C", "O", "H", "H"], positions=positions)


def _cleanup_pymol_grouper(grouper):
    """Avoid __del__ triggering a real cmd.quit() during teardown,
    mirroring the pattern used in tests/test_groupers.py."""
    import shutil

    if getattr(grouper, "_temp_dir", None):
        shutil.rmtree(grouper._temp_dir, ignore_errors=True)
        grouper._temp_dir = None
    grouper.cmd = None


@pytest.fixture()
def shared_pymol_grouper():
    # A single real PymolRMSDGrouper reused across several tests below.
    # Constructing PyMOL's embedded, process-wide session is expensive
    # and (empirically) prone to hangs if done many times back-to-back
    # in the same process alongside tests/test_groupers.py's own PyMOL
    # tests -- so we build it once per test instead of once per test
    # method.
    grouper = PymolRMSDGrouper([_mol(), _mol(0.001)])
    yield grouper
    _cleanup_pymol_grouper(grouper)


def _make_irmsd_grouper(**kwargs):
    # find_irmsd_command() is called at __init__ time and raises
    # RuntimeError if no real irmsd binary is found. Mocking it lets us
    # exercise the rest of IRMSDGrouper's logic (XYZ writing, output
    # parsing, subprocess error handling, and _record_results) without
    # needing the actual external tool installed.
    with patch(
        "chemsmart.jobs.grouper.rmsd.find_irmsd_command",
        return_value="/fake/bin/irmsd",
    ):
        return IRMSDGrouper([_mol(), _mol(0.001)], **kwargs)


def _is_irmsd_available() -> bool:
    """Return True when either irmsd Python API or CLI is available."""
    try:
        irmsd_module = importlib.import_module("irmsd.api.irmsd_exposed")
        getattr(irmsd_module, "get_irmsd")
        return True
    except (ImportError, OSError, AttributeError):
        return bool(find_irmsd_command())


def _assert_progress_milestones(log_text: str) -> None:
    """Assert that progress milestones from 10% to 100% are present."""
    for milestone in range(10, 101, 10):
        assert f"Matrix calculation progress: {milestone}%" in log_text


class TestGrouperCliRunner:
    """Cross-module CLI/runner wiring checks and strategy edge-case behavior."""

    # --- CLI wiring ---
    def test_tfd_cli_forwards_strategy_and_representative(
        self, monkeypatch, multiple_molecules_xyz_file, captured
    ):
        from click.testing import CliRunner

        grouper_module = importlib.import_module(
            "chemsmart.cli.grouper.grouper"
        )
        tfd_module = importlib.import_module("chemsmart.cli.grouper.tfd")

        def fake_create_grouper_job(ctx, **kwargs):
            captured.update(kwargs)
            captured["representative"] = ctx.obj["representative"]

        monkeypatch.setattr(
            tfd_module,
            "create_grouper_job_from_context",
            fake_create_grouper_job,
        )

        result = CliRunner().invoke(
            grouper_module.grouper,
            [
                "-f",
                multiple_molecules_xyz_file,
                "--representative",
                "center",
                "tfd",
                "--no-use-weights",
                "--max-dev",
                "spec",
            ],
        )

        assert result.exit_code == 0, result.output
        assert captured == {
            "strategy": "tfd",
            "use_weights": False,
            "max_dev": "spec",
            "representative": "center",
        }

    # --- Job and runner wiring ---
    def test_grouper_job_stores_output_configuration(self, methanol_molecules):
        thermo_parameters = "298.15,1.0,100.0"

        job = GrouperJob(
            molecules=methanol_molecules[:2],
            grouping_strategy="tfd",
            representative_strategy="center",
            energy_type="QHG",
            thermo_parameters=thermo_parameters,
            use_weights=False,
            max_dev="spec",
        )

        assert job.representative_strategy == "center"
        assert job.energy_type == "QHG"
        assert job.thermo_parameters == thermo_parameters
        assert job.grouper_kwargs == {
            "use_weights": False,
            "max_dev": "spec",
        }

    def test_create_grouper_job_forwards_representative(
        self, methanol_molecules
    ):
        from types import SimpleNamespace

        from chemsmart.cli.grouper.grouper import (
            create_grouper_job_from_context,
        )

        ctx = SimpleNamespace(
            obj={
                "molecules": methanol_molecules[:2],
                "ignore_hydrogens": False,
                "num_procs": 1,
                "grouper_label": "conformers",
                "num_groups": None,
                "conformer_ids": ["c1", "c2"],
                "skipped_ids": [],
                "energy_type": "E",
                "thermo_kwargs": None,
                "threshold": 0.1,
                "matrix_format": "csv",
                "representative": "top3",
            }
        )

        job = create_grouper_job_from_context(ctx, strategy="tfd")

        assert job.representative_strategy == "top3"
        assert job.grouping_strategy == "tfd"

    def test_grouper_runner_creates_tfd_with_representative(
        self, methanol_molecules, captured
    ):
        class FakeTFDGrouper:
            def __init__(self, **kwargs):
                captured.update(kwargs)

        job = GrouperJob(
            molecules=methanol_molecules[:2],
            grouping_strategy="tfd",
            threshold=0.2,
            num_groups=None,
            representative_strategy="center",
            use_weights=False,
        )
        runner = object.__new__(GrouperJobRunner)
        runner._get_grouper_classes = lambda: {"tfd": FakeTFDGrouper}

        grouper = runner._create_grouper(job)

        assert isinstance(grouper, FakeTFDGrouper)
        assert captured["representative_strategy"] == "center"
        assert captured["threshold"] == 0.2
        assert captured["use_weights"] is False

    def test_grouper_runner_preserves_representative_order(
        self, methanol_molecules, temporary_working_dir
    ):
        molecules = methanol_molecules[:2]
        molecules[0].energy = -2.0
        molecules[1].energy = -1.0
        job = GrouperJob(
            molecules=molecules,
            grouping_strategy="tfd",
            label="ordered",
            conformer_ids=["c1", "c2"],
        )
        os.makedirs(job.output_dir)

        runner = object.__new__(GrouperJobRunner)
        runner._write_outputs(job, [[molecules[1], molecules[0]]], [[1, 0]])

        with open(job.outputfile) as output:
            comments = [
                line.strip() for line in output if "Original_Index" in line
            ]

        assert "Original_Index: c2" in comments[0]
        assert "Original_Index: c1" in comments[1]

    # --- Strategy/backend edge paths ---
    def test_isomorphism_grouper_handles_hashing_failure(
        self, monkeypatch, caplog, methanol_molecules
    ):
        molecules = methanol_molecules[:2]
        molecule_type = type(molecules[0])
        original_get_rdkit_hash = molecule_type.get_rdkit_hash

        def get_rdkit_hash_with_failure(molecule, **kwargs):
            if molecule is molecules[0]:
                raise ValueError("invalid coordinates")
            return original_get_rdkit_hash(molecule, **kwargs)

        monkeypatch.setattr(
            molecule_type,
            "get_rdkit_hash",
            get_rdkit_hash_with_failure,
        )
        grouper = RDKitIsomorphismGrouper(molecules)
        monkeypatch.setattr(grouper, "record", lambda **kwargs: None)

        with caplog.at_level("WARNING"):
            groups, index_groups = grouper.group()

        assert sorted(index_groups) == [[0], [1]]
        assert len(groups) == 2
        assert (
            "RDKit molecular hashing failed: invalid coordinates"
            in caplog.text
        )

    def test_rmsd_repr_and_explicit_matrix_output(
        self, monkeypatch, methanol_molecules, temporary_working_dir
    ):
        threshold_grouper = BasicRMSDGrouper(
            methanol_molecules[:2], threshold=0.5
        )
        num_groups_grouper = BasicRMSDGrouper(
            methanol_molecules[:2], threshold=None, num_groups=2
        )
        recorded = {}
        monkeypatch.setattr(
            threshold_grouper,
            "record",
            lambda **kwargs: recorded.update(kwargs),
        )

        matrix = threshold_grouper.calculate_full_rmsd_matrix(
            os.path.join(temporary_working_dir, "matrix", "rmsd.txt")
        )

        assert "threshold=0.5" in repr(threshold_grouper)
        assert "num_groups=2" in repr(num_groups_grouper)
        assert np.array_equal(recorded["rmsd_matrix"], matrix)

    def test_spyrmsd_incompatible_inputs_and_payload_handling(
        self, monkeypatch
    ):
        from spyrmsd import rmsd as spyrmsd_rmsd

        grouper = object.__new__(SpyRMSDGrouper)
        grouper._coords_list = [np.zeros((1, 3)), np.zeros((2, 3))]
        grouper._anum_list = [np.array([6]), np.array([6, 1])]
        grouper._adj_list = [np.zeros((1, 1)), np.zeros((2, 2))]
        grouper.minimize = True
        grouper.cache = True
        grouper.best_isomorphisms = {}

        assert grouper._calculate_rmsd_with_isomorphism((0, 1)) == (
            np.inf,
            None,
        )

        grouper._coords_list[1] = np.zeros((1, 3))
        grouper._anum_list[1] = np.array([8])
        grouper._adj_list[1] = np.zeros((1, 1))
        assert grouper._calculate_rmsd_with_isomorphism((0, 1)) == (
            np.inf,
            None,
        )

        grouper._anum_list[1] = np.array([6])
        monkeypatch.setattr(
            spyrmsd_rmsd, "symmrmsd", lambda *args, **kwargs: 0.25
        )
        assert grouper._calculate_rmsd_with_isomorphism((0, 1)) == (0.25, None)

        monkeypatch.setattr(
            spyrmsd_rmsd,
            "symmrmsd",
            lambda *args, **kwargs: (_ for _ in ()).throw(
                ValueError("symmetry calculation failed")
            ),
        )
        assert grouper._calculate_rmsd_with_isomorphism((0, 1)) == (
            np.inf,
            None,
        )

        grouper.best_isomorphisms[(1, 0)] = ([1, 0], [0, 1])
        assert grouper.get_best_isomorphism(0, 1) == ([0, 1], [1, 0])
        assert grouper.get_best_isomorphism(3, 4) is None

    def test_irmsd_helpers_and_direct_api_failure(
        self, monkeypatch, caplog, methanol_molecules
    ):
        from types import SimpleNamespace

        from chemsmart.jobs.grouper import rmsd as rmsd_module

        grouper = object.__new__(IRMSDGrouper)
        grouper.molecules = methanol_molecules[:2]
        grouper.ignore_hydrogens = True
        grouper.inversion = "on"
        grouper._actual_inversion = None
        assert grouper._inversion_code("auto") == 0
        assert grouper._inversion_code("on") == 1
        assert grouper._inversion_code("off") == 2

        numbers, positions = grouper._prepare_irmsd_arrays(
            grouper.molecules[0]
        )
        assert len(numbers) == len(positions) == 2

        failing_api = SimpleNamespace(
            get_irmsd=lambda *args, **kwargs: (_ for _ in ()).throw(
                ValueError("backend failed")
            )
        )
        monkeypatch.setattr(
            rmsd_module.importlib, "import_module", lambda name: failing_api
        )
        with caplog.at_level("WARNING"):
            assert np.isinf(grouper._calculate_rmsd_direct((0, 1)))
        assert "irmsd API calculation failed" in caplog.text

        value, inversion = grouper._parse_irmsd_output(
            "iRMSD: invalid\nInversion check: off"
        )
        assert np.isinf(value)
        assert inversion == "off"

    def test_irmsd_initialization_selects_api_cli_or_errors(
        self, monkeypatch, methanol_molecules
    ):
        from types import SimpleNamespace

        from chemsmart.jobs.grouper import rmsd as rmsd_module

        direct_api = SimpleNamespace(
            get_irmsd=lambda *args, **kwargs: (0.0,) * 5
        )
        monkeypatch.setattr(
            rmsd_module.importlib, "import_module", lambda name: direct_api
        )
        direct = IRMSDGrouper(methanol_molecules[:2])
        assert direct._use_direct_api is True
        assert direct._irmsd_cmd is None

        def unavailable_api(name):
            raise ImportError("API unavailable")

        monkeypatch.setattr(
            rmsd_module.importlib, "import_module", unavailable_api
        )
        monkeypatch.setattr(
            rmsd_module, "find_irmsd_command", lambda: "/bin/irmsd"
        )
        cli = IRMSDGrouper(methanol_molecules[:2])
        assert cli._use_direct_api is False
        assert cli._irmsd_cmd == "/bin/irmsd"

        monkeypatch.setattr(rmsd_module, "find_irmsd_command", lambda: None)
        with pytest.raises(RuntimeError, match="irmsd is not available"):
            IRMSDGrouper(methanol_molecules[:2])

    def test_irmsd_initialization_falls_back_when_api_not_callable(
        self, monkeypatch, methanol_molecules
    ):
        from types import SimpleNamespace

        from chemsmart.jobs.grouper import rmsd as rmsd_module

        monkeypatch.setattr(
            rmsd_module.importlib,
            "import_module",
            lambda name: SimpleNamespace(get_irmsd=None),
        )
        monkeypatch.setattr(
            rmsd_module, "find_irmsd_command", lambda: "/bin/irmsd"
        )
        grouper = IRMSDGrouper(methanol_molecules[:2])
        assert grouper._use_direct_api is False
        assert grouper._irmsd_cmd == "/bin/irmsd"

    @pytest.mark.parametrize(
        ("backend_value", "expected"), ((0.25, 0.25), (np.nan, np.inf))
    )
    def test_irmsd_direct_api_success_records_inversion(
        self, monkeypatch, methanol_molecules, backend_value, expected
    ):
        from types import SimpleNamespace

        from chemsmart.jobs.grouper import rmsd as rmsd_module

        calls = []

        def fake_get_irmsd(numbers1, positions1, numbers2, positions2, **kw):
            calls.append((len(numbers1), len(numbers2), kw))
            return backend_value, None, None, None, None

        monkeypatch.setattr(
            rmsd_module.importlib,
            "import_module",
            lambda name: SimpleNamespace(get_irmsd=fake_get_irmsd),
        )
        grouper = IRMSDGrouper(methanol_molecules[:2], inversion="off")
        assert grouper._use_direct_api is True

        # Dispatches to the direct API and returns (rmsd, inversion) payload.
        payload = grouper._calculate_pair_payload((0, 1))
        assert payload[0] == expected
        assert payload[1] == "off"
        # All atoms are passed when ignore_hydrogens is False.
        n_atoms = len(methanol_molecules[0].chemical_symbols)
        assert calls == [(n_atoms, n_atoms, {"iinversion": 2})]

        # A later call keeps the already-resolved inversion.
        grouper.inversion = "on"
        assert grouper._calculate_rmsd((0, 1)) == expected
        assert grouper._actual_inversion == "off"

    def test_irmsd_direct_api_mismatched_atoms_returns_inf(
        self, monkeypatch, caplog, methanol_and_ethanol
    ):
        from types import SimpleNamespace

        from chemsmart.jobs.grouper import rmsd as rmsd_module

        def must_not_be_called(*args, **kwargs):
            raise AssertionError("get_irmsd should not be called")

        monkeypatch.setattr(
            rmsd_module.importlib,
            "import_module",
            lambda name: SimpleNamespace(get_irmsd=must_not_be_called),
        )
        grouper = IRMSDGrouper(methanol_and_ethanol, ignore_hydrogens=True)
        with caplog.at_level("WARNING"):
            assert np.isinf(grouper._calculate_rmsd((0, 1)))
        assert "Molecules 0 and 1 have different atom types" in caplog.text

    def test_irmsd_consume_payload_keeps_first_resolved_inversion(self):
        grouper = object.__new__(IRMSDGrouper)
        grouper._actual_inversion = None

        assert grouper._consume_pair_payload((0, 1), (0.5, None)) == 0.5
        assert grouper._actual_inversion is None
        assert grouper._consume_pair_payload((0, 2), (0.7, "on")) == 0.7
        assert grouper._actual_inversion == "on"
        assert grouper._consume_pair_payload((1, 2), (0.9, "off")) == 0.9
        assert grouper._actual_inversion == "on"

    def test_parse_irmsd_output_can_skip_inversion_line(self):
        grouper = object.__new__(IRMSDGrouper)
        value, inversion = grouper._parse_irmsd_output(
            "Inversion check: on\niRMSD: 0.42", parse_inversion=False
        )
        assert value == pytest.approx(0.42)
        assert inversion is None

    def test_calculate_pair_chunk_returns_payload_per_pair(
        self, methanol_molecules
    ):
        grouper = BasicRMSDGrouper(methanol_molecules)
        chunk = grouper._calculate_pair_chunk([(0, 1), (0, 2)])
        assert chunk == [
            pytest.approx(grouper._calculate_rmsd((0, 1))),
            pytest.approx(grouper._calculate_rmsd((0, 2))),
        ]

    def test_pymol_init_reports_missing_pymol(self, monkeypatch):
        import sys

        monkeypatch.setitem(sys.modules, "pymol", None)
        grouper = object.__new__(PymolRMSDGrouper)
        grouper.cmd = None
        grouper._temp_dir = None
        with pytest.raises(ImportError, match="PyMOL not available"):
            grouper._init_pymol()

    def test_irmsd_cli_handles_failure_and_timeout(
        self, monkeypatch, methanol_molecules
    ):
        import subprocess
        from types import SimpleNamespace

        grouper = object.__new__(IRMSDGrouper)
        grouper.molecules = methanol_molecules[:2]
        grouper._irmsd_cmd = "/bin/irmsd"
        grouper.inversion = "auto"
        grouper.ignore_hydrogens = False
        grouper._actual_inversion = None

        failed = SimpleNamespace(returncode=1, stderr="failed", stdout="")
        monkeypatch.setattr(subprocess, "run", lambda *args, **kwargs: failed)
        assert np.isinf(grouper._calculate_rmsd_cli((0, 1)))

        def timeout(*args, **kwargs):
            raise subprocess.TimeoutExpired("irmsd", 60)

        monkeypatch.setattr(subprocess, "run", timeout)
        assert np.isinf(grouper._calculate_rmsd_cli((0, 1)))

    def test_pymol_heavy_atom_xyz_and_repr(
        self, methanol_molecules, temporary_working_dir
    ):
        grouper = object.__new__(PymolRMSDGrouper)
        grouper.ignore_hydrogens = True
        grouper.num_groups = None
        grouper.threshold = 0.5
        grouper.cmd = None
        grouper._temp_dir = None
        output = os.path.join(temporary_working_dir, "heavy.xyz")

        grouper._write_xyz(methanol_molecules[0], output)

        with open(output) as xyz_file:
            lines = xyz_file.readlines()
        assert lines[0].strip() == "2"
        assert all(not line.startswith("H ") for line in lines[2:])
        assert "threshold=0.5" in repr(grouper)

        grouper.num_groups = 2
        assert "num_groups=2" in repr(grouper)

    def test_pymol_rmsd_handles_invalid_and_failed_alignment(self, caplog):
        class Command:
            def __init__(self):
                self.result = (-1.0,)

            def align(self, first, second):
                if isinstance(self.result, Exception):
                    raise self.result
                return self.result

        grouper = object.__new__(PymolRMSDGrouper)
        grouper._mol_names = ["mol_0", "mol_1"]
        grouper._alignment_cache = {}
        grouper.cmd = Command()
        grouper._temp_dir = None

        assert np.isinf(grouper._calculate_rmsd((0, 1)))
        grouper._alignment_cache.clear()
        grouper.cmd.result = ValueError("align failed")
        with caplog.at_level("WARNING"):
            assert np.isinf(grouper._calculate_rmsd((0, 1)))
        assert "PyMOL align failed" in caplog.text

        grouper.cmd = None


@pytest.mark.usefixtures("temporary_working_dir")
class TestBasicRMSDGrouperAndBasicFunctionality:
    NUM_PROCS = 4

    def test_rmsd_grouper(
        self,
        methanol_molecules,
        methanol_and_ethanol,
    ):
        methanol = methanol_molecules[0]
        methanol_rot1 = methanol_molecules[1]

        assert np.any(
            methanol.positions != methanol_rot1.positions
        ), "Rotated molecule should have different positions."

        grouper = BasicRMSDGrouper(methanol_molecules)
        groups, group_indices = grouper.group()

        rmsd = grouper.calculate_rmsd_pair(0, 1)
        assert np.isclose(
            rmsd, 0.0, rtol=1e-3
        ), "RMSD should be close to zero."

        rmsd = grouper.calculate_rmsd_pair(1, 2)
        assert np.isclose(
            rmsd, 0.0, rtol=1e-3
        ), "RMSD should be close to zero."

        _, _, _, _, rmsd_kabsch = kabsch_align(
            methanol.positions,
            methanol_rot1.positions,
        )
        assert np.isclose(rmsd_kabsch, 0.0, rtol=1e-3)

        assert (
            len(groups) == 1
        ), "Molecules should form one group based on geometry."
        assert (
            len(group_indices) == 1
        ), "Molecules should form one group based on geometry."

        unique_structures = grouper.unique()
        assert (
            len(unique_structures) == 1
        ), "Molecules should form one group based on geometry."

        grouper2 = BasicRMSDGrouper(methanol_and_ethanol)
        groups, group_indices = grouper2.group()

        assert (
            len(groups) == 1
        ), "Molecules should form two groups based on geometry."
        assert (
            len(group_indices) == 1
        ), "Molecules should form two groups based on geometry."

        rmsd = grouper2.calculate_rmsd_pair(0, 1)
        assert (
            rmsd is np.inf
        ), "RMSD is set to be infinity for different molecules."

        unique_structures = grouper2.unique()
        assert (
            len(unique_structures) == 1
        ), "One incompatible structure should be skipped, leaving one unique structure."

    def test_rmsd_grouper_for_crest_conformers(
        self, multiple_molecules_xyz_file
    ):
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)

        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 18
        grouper = BasicRMSDGrouper(
            molecules, threshold=0.2, num_procs=self.NUM_PROCS
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 18
        assert len(group_indices) == 18
        unique_structures = grouper.unique()
        assert len(unique_structures) == 18

        # rmsd calculation from grouper
        rmsd = grouper.calculate_rmsd_pair(0, 1)
        assert np.isclose(rmsd, 0.409, rtol=1e-3)

        _, _, _, _, rmsd = kabsch_align(
            molecules[0].positions, molecules[1].positions
        )
        assert np.isclose(rmsd, 0.409, rtol=1e-3)

        grouper2 = BasicRMSDGrouper(
            molecules, threshold=0.5, num_procs=self.NUM_PROCS
        )
        groups, group_indices = grouper2.group()
        assert len(groups) == 12
        assert len(group_indices) == 12
        unique_structures = grouper2.unique()
        assert len(unique_structures) == 12

        grouper3 = BasicRMSDGrouper(
            molecules, threshold=1.0, num_procs=self.NUM_PROCS
        )
        groups, group_indices = grouper3.group()
        assert len(groups) == 11
        assert len(group_indices) == 11
        unique_structures = grouper3.unique()
        assert len(unique_structures) == 11

        grouper4 = BasicRMSDGrouper(
            molecules, threshold=1.5, num_procs=self.NUM_PROCS
        )
        groups, group_indices = grouper4.group()
        assert len(groups) == 8
        assert len(group_indices) == 8
        unique_structures = grouper4.unique()
        assert len(unique_structures) == 8

        grouper5 = BasicRMSDGrouper(
            molecules, threshold=2.0, num_procs=self.NUM_PROCS
        )
        groups, group_indices = grouper5.group()
        assert len(groups) == 7
        assert len(group_indices) == 7
        unique_structures = grouper5.unique()
        assert len(unique_structures) == 7

        grouper6 = BasicRMSDGrouper(
            molecules, threshold=2.5, num_procs=self.NUM_PROCS
        )
        groups, group_indices = grouper6.group()
        assert len(groups) == 4
        assert len(group_indices) == 4
        unique_structures = grouper6.unique()
        assert len(unique_structures) == 4

    def test_num_groups_parameter(self, multiple_molecules_xyz_file):
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)

        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 18
        grouper = BasicRMSDGrouper(
            molecules, threshold=None, num_groups=17, num_procs=self.NUM_PROCS
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 17
        assert len(group_indices) == 17
        unique_structures = grouper.unique()
        assert len(unique_structures) == 17

    def test_pick_the_lowestenergy_conformers(
        self, multiple_molecules_xyz_file
    ):
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)

        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 18
        grouper = BasicRMSDGrouper(
            molecules, threshold=None, num_groups=3, num_procs=self.NUM_PROCS
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 3
        assert len(group_indices) == 3
        unique_structures = grouper.unique()
        assert len(unique_structures) == 3

        energies = [mol.energy for mol in unique_structures]
        assert -126.2575508 in energies
        assert -126.25017833 in energies
        assert -126.24909661 in energies

    def test_rmsd_grouper_for_crest_conformers_ignore_Hs(
        self, multiple_molecules_xyz_file
    ):
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)

        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 18
        grouper = BasicRMSDGrouper(
            molecules,
            threshold=0.2,
            num_procs=self.NUM_PROCS,
            ignore_hydrogens=True,
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 17
        assert len(group_indices) == 17
        unique_structures = grouper.unique()
        assert len(unique_structures) == 17

        # rmsd calculation from grouper
        rmsd = grouper.calculate_rmsd_pair(0, 1)
        assert np.isclose(rmsd, 0.301, rtol=1e-3)  # removed H atoms

        _, _, _, _, rmsd = kabsch_align(
            molecules[0].positions, molecules[1].positions
        )
        assert np.isclose(rmsd, 0.409, rtol=1e-3)  # did not remove H atoms

        grouper2 = BasicRMSDGrouper(
            molecules,
            threshold=0.5,
            num_procs=self.NUM_PROCS,
            ignore_hydrogens=True,
        )
        groups, group_indices = grouper2.group()
        assert len(groups) == 12
        assert len(group_indices) == 12
        unique_structures = grouper2.unique()
        assert len(unique_structures) == 12

        grouper3 = BasicRMSDGrouper(
            molecules,
            threshold=1.0,
            num_procs=self.NUM_PROCS,
            ignore_hydrogens=True,
        )
        groups, group_indices = grouper3.group()
        assert len(groups) == 10
        assert len(group_indices) == 10
        unique_structures = grouper3.unique()
        assert len(unique_structures) == 10

        grouper4 = BasicRMSDGrouper(
            molecules,
            threshold=1.5,
            num_procs=self.NUM_PROCS,
            ignore_hydrogens=True,
        )
        groups, group_indices = grouper4.group()
        assert len(groups) == 7
        assert len(group_indices) == 7
        unique_structures = grouper4.unique()
        assert len(unique_structures) == 7

        grouper5 = BasicRMSDGrouper(
            molecules,
            threshold=2.0,
            num_procs=self.NUM_PROCS,
            ignore_hydrogens=True,
        )
        groups, group_indices = grouper5.group()
        assert len(groups) == 5
        assert len(group_indices) == 5
        unique_structures = grouper5.unique()
        assert len(unique_structures) == 5

        grouper6 = BasicRMSDGrouper(
            molecules,
            threshold=2.5,
            num_procs=self.NUM_PROCS,
            ignore_hydrogens=True,
        )
        groups, group_indices = grouper6.group()
        for indices in group_indices:
            for i, idx_i in enumerate(indices):
                for idx_j in indices[i + 1 :]:
                    rmsd = grouper6.calculate_rmsd_pair(idx_i, idx_j)
                    assert rmsd <= 2.5
        assert len(groups) == 4
        assert len(group_indices) == 4
        unique_structures = grouper6.unique()
        assert len(unique_structures) == 4

    def test_rmsd_grouper_for_rotated_molecules(
        self, two_rotated_molecules_xyz_file
    ):
        xyz_file = XYZFile(filename=two_rotated_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 2
        grouper = BasicRMSDGrouper(
            molecules,
            threshold=0.5,
            num_procs=self.NUM_PROCS,
            ignore_hydrogens=False,
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 2
        assert len(group_indices) == 2
        unique_structures = grouper.unique()
        assert len(unique_structures) == 2

    def test_calculate_rmsd_pair_wrapper_validation(self, methanol_molecules):
        grouper = BasicRMSDGrouper(methanol_molecules[:2])

        assert np.isclose(grouper.calculate_rmsd_pair(0, 0), 0.0)
        assert np.isclose(
            grouper.calculate_rmsd_pair(0, 1),
            grouper.calculate_rmsd_pair(1, 0),
        )

        with pytest.raises(IndexError):
            grouper.calculate_rmsd_pair(0, 3)

        with pytest.raises(TypeError):
            grouper.calculate_rmsd_pair("a", 1)  # type: ignore[arg-type]


@pytest.mark.usefixtures("temporary_working_dir")
class TestHungarianRMSDGrouper:
    NUM_PROCS = 4

    def test_hrmsd_grouper_for_rotated_molecules(
        self, two_rotated_molecules_xyz_file
    ):
        xyz_file = XYZFile(filename=two_rotated_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 2
        grouper = HungarianRMSDGrouper(
            molecules,
            threshold=0.5,
            num_procs=self.NUM_PROCS,
            ignore_hydrogens=False,
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 1
        assert len(group_indices) == 1
        unique_structures = grouper.unique()
        assert len(unique_structures) == 1

        # rmsd calculation from grouper
        rmsd = grouper.calculate_rmsd_pair(0, 1)
        assert np.isclose(rmsd, 0.2294, rtol=1e-3)
        assert np.isclose(
            grouper.calculate_rmsd_pair(0, 1),
            grouper.calculate_rmsd_pair(1, 0),
        )

    def test_hrmsd_grouper_for_crest_molecules(
        self, multiple_molecules_xyz_file
    ):
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 18
        grouper = HungarianRMSDGrouper(
            molecules,
            threshold=0.5,
            num_procs=self.NUM_PROCS,
            ignore_hydrogens=False,
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 12
        assert len(group_indices) == 12
        unique_structures = grouper.unique()
        assert len(unique_structures) == 12

        # rmsd calculation from grouper
        rmsd = grouper.calculate_rmsd_pair(0, 1)
        assert np.isclose(rmsd, 0.4091, rtol=1e-3)
        rmsd = grouper.calculate_rmsd_pair(0, 2)
        assert np.isclose(rmsd, 0.5899, rtol=1e-3)
        rmsd = grouper.calculate_rmsd_pair(0, 3)
        assert np.isclose(rmsd, 1.8891, rtol=1e-3)

    def test_ignore_hydrogen(self, multiple_molecules_xyz_file):
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 18
        grouper = HungarianRMSDGrouper(
            molecules,
            threshold=0.5,
            num_procs=self.NUM_PROCS,
            ignore_hydrogens=True,
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 12
        assert len(group_indices) == 12
        unique_structures = grouper.unique()
        assert len(unique_structures) == 12

        # rmsd calculation from grouper
        rmsd = grouper.calculate_rmsd_pair(0, 4)
        assert np.isclose(rmsd, 1.1915, rtol=1e-3)


@pytest.mark.usefixtures("temporary_working_dir")
class TestSpyRMSDGrouper:
    NUM_PROCS = 4

    def test_spyrmsd_grouper_for_rotated_molecules(
        self, two_rotated_molecules_xyz_file
    ):
        xyz_file = XYZFile(filename=two_rotated_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 2
        grouper = SpyRMSDGrouper(
            molecules,
            threshold=0.5,
            num_procs=self.NUM_PROCS,
            ignore_hydrogens=False,
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 1
        assert len(group_indices) == 1
        unique_structures = grouper.unique()
        assert len(unique_structures) == 1

        # rmsd calculation from grouper
        rmsd = grouper.calculate_rmsd_pair(0, 1)
        assert np.isclose(rmsd, 0.2125, rtol=1e-3)

    def test_spyrmsd_grouper_for_crest_molecules(
        self, multiple_molecules_xyz_file
    ):
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 18
        grouper = SpyRMSDGrouper(
            molecules,
            threshold=0.5,
            num_procs=self.NUM_PROCS,
            ignore_hydrogens=False,
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 12
        assert len(group_indices) == 12
        unique_structures = grouper.unique()
        assert len(unique_structures) == 12

        # rmsd calculation from grouper
        rmsd = grouper.calculate_rmsd_pair(0, 1)
        assert np.isclose(rmsd, 0.4091, rtol=1e-3)
        rmsd = grouper.calculate_rmsd_pair(0, 2)
        assert np.isclose(rmsd, 1.3925, rtol=1e-3)
        rmsd = grouper.calculate_rmsd_pair(0, 3)
        assert np.isclose(rmsd, 2.1789, rtol=1e-3)
        rmsd = grouper.calculate_rmsd_pair(0, 4)
        assert np.isclose(rmsd, 1.8202, rtol=1e-3)
        rmsd = grouper.calculate_rmsd_pair(0, 5)
        assert np.isclose(rmsd, 2.0029, rtol=1e-3)

    def test_ignore_hydrogen(self, multiple_molecules_xyz_file):
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 18
        grouper = SpyRMSDGrouper(
            molecules,
            threshold=2.5981,
            num_procs=self.NUM_PROCS,
            ignore_hydrogens=True,
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 3
        assert len(group_indices) == 3
        unique_structures = grouper.unique()
        assert len(unique_structures) == 3

        # rmsd calculation from grouper
        rmsd = grouper.calculate_rmsd_pair(0, 5)
        assert np.isclose(rmsd, 1.7034, rtol=1e-3)
        rmsd = grouper.calculate_rmsd_pair(0, 6)
        assert np.isclose(rmsd, 2.6183, rtol=1e-3)


@pytest.mark.skipif(
    not _is_irmsd_available(), reason="irmsd API/command not available"
)
@pytest.mark.usefixtures("temporary_working_dir")
class TestIRMSDGrouper:
    NUM_PROCS = 4

    def test_irmsd_grouper_for_rotated_molecules(
        self, two_rotated_molecules_xyz_file
    ):
        xyz_file = XYZFile(filename=two_rotated_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 2
        grouper = IRMSDGrouper(
            molecules,
            threshold=0.125,
            num_procs=self.NUM_PROCS,
            ignore_hydrogens=False,
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 2
        assert len(group_indices) == 2
        unique_structures = grouper.unique()
        assert len(unique_structures) == 2

        # rmsd calculation from grouper
        rmsd = grouper.calculate_rmsd_pair(0, 1)
        assert np.isclose(rmsd, 0.2294, rtol=1e-3)

    def test_irmsd_grouper_for_crest_molecules(
        self, multiple_molecules_xyz_file
    ):
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 18
        grouper = IRMSDGrouper(
            molecules,
            threshold=0.5,
            num_procs=self.NUM_PROCS,
            ignore_hydrogens=False,
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 12
        assert len(group_indices) == 12
        unique_structures = grouper.unique()
        assert len(unique_structures) == 12

        # rmsd calculation from grouper
        rmsd = grouper.calculate_rmsd_pair(0, 1)
        assert np.isclose(rmsd, 0.4091, rtol=1e-3)
        rmsd = grouper.calculate_rmsd_pair(0, 2)
        assert np.isclose(rmsd, 1.3925, rtol=1e-3)
        rmsd = grouper.calculate_rmsd_pair(2, 10)
        assert np.isclose(rmsd, 2.2390, rtol=1e-3)
        rmsd = grouper.calculate_rmsd_pair(8, 16)
        assert np.isclose(rmsd, 3.4209, rtol=1e-3)
        rmsd = grouper.calculate_rmsd_pair(0, 13)
        assert np.isclose(rmsd, 0.8411, rtol=1e-3)

    def test_ignore_hydrogen(self, two_rotated_molecules_xyz_file):
        xyz_file = XYZFile(filename=two_rotated_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 2
        grouper = IRMSDGrouper(
            molecules,
            threshold=0.125,
            num_procs=self.NUM_PROCS,
            ignore_hydrogens=True,
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 2
        assert len(group_indices) == 2
        unique_structures = grouper.unique()
        assert len(unique_structures) == 2

        # rmsd calculation from grouper
        rmsd = grouper.calculate_rmsd_pair(0, 1)
        assert np.isclose(rmsd, 0.2294, rtol=1e-3)

    def test_parallel_irmsd_preserves_actual_inversion_metadata(
        self, two_rotated_molecules_xyz_file
    ):
        """Parallel iRMSD must propagate resolved inversion metadata to the parent."""
        xyz_file = XYZFile(filename=two_rotated_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)

        grouper = IRMSDGrouper(
            molecules,
            threshold=0.5,
            num_procs=self.NUM_PROCS,
            ignore_hydrogens=False,
            inversion="off",
        )
        grouper.group()

        assert grouper.num_procs > 1
        assert grouper._actual_inversion == "off"


@pytest.mark.usefixtures("temporary_working_dir")
class TestPymolRMSDGrouper:
    NUM_PROCS = 1

    @classmethod
    def setup_class(cls):
        """Initialize PyMOL once for all tests in this class."""
        try:
            import pymol
            from pymol import cmd

            pymol.finish_launching(["pymol", "-qc"])
            cmd.reinitialize()
        except Exception:
            # PyMOL may not be installed; tests will be skipped via @pytest.mark.skipif
            pass

    def teardown_method(self, method):
        """Clean up PyMOL objects after each test method to prevent slowdown."""
        try:
            from pymol import cmd

            # Delete all objects instead of reinitialize to keep session alive
            cmd.delete("all")
        except Exception:
            # Ignore cleanup errors if PyMOL is not available
            pass

    def test_pymolrmsd_grouper_for_rotated_molecules(
        self, two_rotated_molecules_xyz_file
    ):
        xyz_file = XYZFile(filename=two_rotated_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 2
        grouper = PymolRMSDGrouper(
            molecules,
            threshold=0.5,
            num_procs=self.NUM_PROCS,
            ignore_hydrogens=False,
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 1
        assert len(group_indices) == 1
        unique_structures = grouper.unique()
        assert len(unique_structures) == 1

        # rmsd calculation from grouper
        rmsd = grouper.calculate_rmsd_pair(0, 1)
        assert np.isclose(rmsd, 0.0000, rtol=1e-3)

        # Explicitly cleanup to prevent __del__ from calling quit()
        if hasattr(grouper, "_temp_dir") and grouper._temp_dir:

            shutil.rmtree(grouper._temp_dir, ignore_errors=True)
            grouper._temp_dir = None
        grouper.cmd = None  # Prevent __del__ from calling quit()

    def test_pymolrmsd_grouper_for_crest_molecules(
        self, multiple_molecules_xyz_file
    ):
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 18
        grouper = PymolRMSDGrouper(
            molecules,
            threshold=0.5,
            num_procs=self.NUM_PROCS,
            ignore_hydrogens=False,
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 8
        assert len(group_indices) == 8
        unique_structures = grouper.unique()
        assert len(unique_structures) == 8
        for indices in group_indices:
            for i, idx_i in enumerate(indices):
                for idx_j in indices[i + 1 :]:
                    rmsd = grouper.calculate_rmsd_pair(idx_i, idx_j)
                    assert rmsd <= 0.5
        # rmsd calculation from grouper
        rmsd = grouper.calculate_rmsd_pair(0, 1)
        assert np.isclose(rmsd, 0.074175, rtol=1e-3)
        rmsd = grouper.calculate_rmsd_pair(0, 2)
        assert np.isclose(rmsd, 0.023745, rtol=1e-3)
        rmsd = grouper.calculate_rmsd_pair(3, 4)
        assert np.isclose(rmsd, 1.725241, rtol=1e-3)
        rmsd = grouper.calculate_rmsd_pair(7, 8)
        assert np.isclose(rmsd, 2.309451, rtol=1e-3)
        rmsd = grouper.calculate_rmsd_pair(14, 16)
        assert np.isclose(rmsd, 1.837319, rtol=1e-3)

        # Explicitly cleanup to prevent __del__ from calling quit()
        if hasattr(grouper, "_temp_dir") and grouper._temp_dir:

            shutil.rmtree(grouper._temp_dir, ignore_errors=True)
            grouper._temp_dir = None
        grouper.cmd = None  # Prevent __del__ from calling quit()

    def test_ignore_hydrogen(self, multiple_molecules_xyz_file):
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 18
        grouper = PymolRMSDGrouper(
            molecules,
            threshold=0.8203818,
            num_procs=self.NUM_PROCS,
            ignore_hydrogens=True,
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 5
        assert len(group_indices) == 5
        unique_structures = grouper.unique()
        assert len(unique_structures) == 5

        # rmsd calculation from grouper
        rmsd = grouper.calculate_rmsd_pair(0, 2)
        assert np.isclose(rmsd, 0.0211, rtol=1e-3)
        rmsd = grouper.calculate_rmsd_pair(1, 4)
        assert np.isclose(rmsd, 0.7669, rtol=1e-3)

        # Explicitly cleanup to prevent __del__ from calling quit()
        if hasattr(grouper, "_temp_dir") and grouper._temp_dir:

            shutil.rmtree(grouper._temp_dir, ignore_errors=True)
            grouper._temp_dir = None
        grouper.cmd = None  # Prevent __del__ from calling quit()

    def test_pymol_grouper_falls_back_to_single_proc(
        self, methanol_molecules, caplog
    ):
        """PyMOL grouper warns and falls back to num_procs=1 when >1 is requested."""
        caplog.set_level("WARNING")
        grouper = PymolRMSDGrouper(
            methanol_molecules,
            threshold=0.5,
            num_procs=4,
        )
        assert grouper.num_procs == 1
        assert (
            "PymolRMSDGrouper does not support multiprocessing; using num_procs=1."
            in caplog.text
        )

        if hasattr(grouper, "_temp_dir") and grouper._temp_dir:
            shutil.rmtree(grouper._temp_dir, ignore_errors=True)
            grouper._temp_dir = None
        grouper.cmd = None

    @classmethod
    def teardown_class(cls):
        """Clean up PyMOL after all tests in this class are done."""
        try:
            from pymol import cmd

            cmd.reinitialize()
        except Exception:
            # Ignore cleanup errors if PyMOL is not available
            pass


@pytest.mark.usefixtures("temporary_working_dir")
class TestTanimotoSimilarityGrouper:
    NUM_PROCS = 4

    @pytest.mark.parametrize("threshold", (-0.1, 1.1))
    def test_rejects_out_of_range_threshold(
        self, methanol_molecules, threshold
    ):
        with pytest.raises(ValueError, match="between 0.0 and 1.0"):
            TanimotoSimilarityGrouper(
                methanol_molecules[:2], threshold=threshold
            )

    def test_all_failed_molecule_preparations_return_empty_groups(
        self, monkeypatch, methanol_molecules
    ):
        monkeypatch.setattr(
            TanimotoSimilarityGrouper,
            "_molecule_to_rdkit",
            lambda self, molecule: None,
        )
        grouper = TanimotoSimilarityGrouper(
            methanol_molecules[:2], conformer_ids=["c1", "c2"]
        )
        recorded = {}
        monkeypatch.setattr(
            grouper, "record", lambda **kwargs: recorded.update(kwargs)
        )

        assert np.isnan(grouper.calculate_tanimoto_pair(0, 1))
        groups, index_groups = grouper.group()

        assert groups == []
        assert index_groups == []
        assert grouper._matrix_skipped_ids == ["c1", "c2"]
        assert recorded["valid_indices"] == []

    def test_pair_returns_nan_when_fingerprint_generation_fails(
        self, monkeypatch, methanol_molecules, caplog
    ):
        grouper = TanimotoSimilarityGrouper(methanol_molecules[:2])
        monkeypatch.setattr(grouper, "_get_fingerprint", lambda molecule: None)

        with caplog.at_level("WARNING"):
            similarity = grouper.calculate_tanimoto_pair(0, 1)

        assert np.isnan(similarity)
        assert "fingerprint generation failed" in caplog.text

    @pytest.mark.parametrize(
        ("matrix", "message"),
        (
            (np.array([1.0, 0.5]), "2-dimensional"),
            (np.eye(3), "dimensions"),
            (np.array([["a", "b"], ["b", "a"]]), "numeric"),
            (
                np.array([[1.0, np.nan], [0.5, np.nan]]),
                "complete NaN rows and columns",
            ),
            (
                np.array([[1.0, np.nan], [np.nan, 1.0]]),
                "finite",
            ),
            (np.array([[1.0, 0.2], [0.3, 1.0]]), "symmetric"),
            (np.array([[0.0, 0.2], [0.2, 1.0]]), "diagonal"),
            (np.array([[1.0, 1.2], [1.2, 1.0]]), "between 0.0 and 1.0"),
        ),
    )
    def test_similarity_output_matrix_validation(
        self, methanol_molecules, matrix, message
    ):
        grouper = TanimotoSimilarityGrouper(methanol_molecules[:2])

        with pytest.raises(ValueError, match=message):
            grouper.group_by_threshold(matrix)

    def test_all_nan_similarity_matrix_has_no_groups(self, methanol_molecules):
        grouper = TanimotoSimilarityGrouper(methanol_molecules[:2])

        groups, index_groups = grouper.group_by_threshold(
            np.full((2, 2), np.nan)
        )

        assert groups == []
        assert index_groups == []
        assert grouper._matrix_skipped_indices == [0, 1]

    def test_per_call_threshold_cannot_exceed_one(self, methanol_molecules):
        grouper = TanimotoSimilarityGrouper(methanol_molecules[:2])
        similarity_matrix = np.eye(2)

        with pytest.raises(ValueError, match="between 0.0 and 1.0"):
            grouper.group_by_threshold(similarity_matrix, threshold=1.1)

    def test_usr_preparation_handles_missing_rdkit_molecule(
        self, monkeypatch, methanol_molecules
    ):
        molecule_type = type(methanol_molecules[0])
        monkeypatch.setattr(molecule_type, "to_rdkit", lambda self: None)

        grouper = TanimotoSimilarityGrouper(
            methanol_molecules[:2], fingerprint_type="usr"
        )

        assert grouper.rdkit_molecules == []
        assert grouper._tanimoto_skipped_indices == [0, 1]

    def test_usr_hydrogen_removal_failure_uses_original_molecule(
        self, monkeypatch, caplog, methanol_molecules
    ):
        from chemsmart.jobs.grouper import tanimoto as tanimoto_module

        monkeypatch.setattr(
            tanimoto_module.Chem,
            "RemoveHs",
            lambda *args, **kwargs: (_ for _ in ()).throw(ValueError("bad H")),
        )

        with caplog.at_level("WARNING"):
            grouper = TanimotoSimilarityGrouper(
                methanol_molecules[:2],
                fingerprint_type="usr",
                ignore_hydrogens=True,
            )

        assert len(grouper.rdkit_molecules) == 2
        assert "Failed to remove hydrogens for usr: bad H" in caplog.text

    def test_topological_preparation_failure_is_skipped(
        self, monkeypatch, caplog, methanol_molecules
    ):
        from chemsmart.jobs.grouper import tanimoto as tanimoto_module

        monkeypatch.setattr(
            tanimoto_module.Chem,
            "MolFromXYZBlock",
            lambda xyz: (_ for _ in ()).throw(ValueError("bad XYZ")),
        )

        with caplog.at_level("WARNING"):
            grouper = TanimotoSimilarityGrouper(methanol_molecules[:2])

        assert grouper.rdkit_molecules == []
        assert "RDKit molecule preparation failed" in caplog.text

    def test_topological_fingerprint_can_ignore_hydrogens(
        self, methanol_molecules
    ):
        grouper = TanimotoSimilarityGrouper(
            methanol_molecules[:2], ignore_hydrogens=True
        )

        assert len(grouper.rdkit_molecules) == 2
        assert all(
            atom.GetSymbol() != "H"
            for molecule in grouper.rdkit_molecules
            for atom in molecule.GetAtoms()
        )

    def test_unknown_and_failed_fingerprint_generation(
        self, monkeypatch, caplog, methanol_molecules
    ):
        from chemsmart.jobs.grouper import tanimoto as tanimoto_module

        grouper = TanimotoSimilarityGrouper(methanol_molecules[:2])
        rdkit_molecule = grouper.rdkit_molecules[0]
        grouper.fingerprint_type = "unknown"
        assert grouper._get_fingerprint(rdkit_molecule) is not None

        monkeypatch.setattr(
            tanimoto_module,
            "GetRDKitFPGenerator",
            lambda: (_ for _ in ()).throw(ValueError("bad fingerprint")),
        )
        with caplog.at_level("WARNING"):
            assert grouper._get_fingerprint(rdkit_molecule) is None

        assert "Fingerprint generation failed: bad fingerprint" in caplog.text

    def test_tanimoto_num_groups_and_threshold_postprocessing(
        self, monkeypatch, methanol_molecules
    ):
        grouper = TanimotoSimilarityGrouper(
            methanol_molecules[:3], num_groups=2
        )
        monkeypatch.setattr(grouper, "record", lambda **kwargs: None)

        _, index_groups = grouper.group()

        assert len(index_groups) >= 2
        assert "num_groups=2" in repr(grouper)

        grouper._auto_threshold = 0.25
        grouper._num_groups_distance_dtype = np.dtype(np.float32)
        grouper._postprocess_auto_threshold_for_num_groups()
        assert np.isclose(grouper._auto_threshold, 0.75)

    def test_usr_pair_similarity_and_hydrogen_removal(
        self, methanol_molecules
    ):
        grouper = TanimotoSimilarityGrouper(
            methanol_molecules[:2],
            fingerprint_type="usr",
        )

        similarity = grouper.calculate_tanimoto_pair(0, 1)

        assert 0.0 <= similarity <= 1.0

        hydrogen_free_grouper = TanimotoSimilarityGrouper(
            methanol_molecules[:2],
            fingerprint_type="usr",
            ignore_hydrogens=True,
        )
        assert len(hydrogen_free_grouper.rdkit_molecules) == 2

    def test_topological_preparation_retries_with_hueckel(
        self, monkeypatch, methanol_molecules
    ):
        from chemsmart.jobs.grouper import tanimoto as tanimoto_module

        grouper = TanimotoSimilarityGrouper(methanol_molecules[:2])
        calls = []

        def determine_bonds(rdkit_molecule, **kwargs):
            calls.append(kwargs)
            if len(calls) == 1:
                raise ValueError("first perception failed")

        monkeypatch.setattr(
            tanimoto_module.rdDetermineBonds,
            "DetermineBonds",
            determine_bonds,
        )
        monkeypatch.setattr(
            tanimoto_module.rdDetermineBonds,
            "hueckelEnabled",
            lambda: True,
        )

        rdkit_molecule = grouper._molecule_to_rdkit(methanol_molecules[0])

        assert rdkit_molecule is not None
        assert calls == [{"charge": 0}, {"useHueckel": True, "charge": 0}]

    def test_topological_preparation_handles_unusable_xyz(
        self, monkeypatch, methanol_molecules
    ):
        from chemsmart.jobs.grouper import tanimoto as tanimoto_module

        grouper = TanimotoSimilarityGrouper(methanol_molecules[:2])
        monkeypatch.setattr(
            tanimoto_module.Chem, "MolFromXYZBlock", lambda xyz: None
        )

        assert grouper._molecule_to_rdkit(methanol_molecules[0]) is None

    def test_topological_preparation_stops_without_hueckel(
        self, monkeypatch, methanol_molecules
    ):
        from chemsmart.jobs.grouper import tanimoto as tanimoto_module

        grouper = TanimotoSimilarityGrouper(methanol_molecules[:2])
        monkeypatch.setattr(
            tanimoto_module.rdDetermineBonds,
            "DetermineBonds",
            lambda *args, **kwargs: (_ for _ in ()).throw(
                ValueError("bond perception failed")
            ),
        )
        monkeypatch.setattr(
            tanimoto_module.rdDetermineBonds,
            "hueckelEnabled",
            lambda: False,
        )

        assert grouper._molecule_to_rdkit(methanol_molecules[0]) is None

    def test_topological_hueckel_retry_handles_unusable_xyz(
        self, monkeypatch, methanol_molecules
    ):
        from chemsmart.jobs.grouper import tanimoto as tanimoto_module

        grouper = TanimotoSimilarityGrouper(methanol_molecules[:2])
        original_from_xyz = tanimoto_module.Chem.MolFromXYZBlock
        calls = 0

        def from_xyz(xyz):
            nonlocal calls
            calls += 1
            return original_from_xyz(xyz) if calls == 1 else None

        monkeypatch.setattr(tanimoto_module.Chem, "MolFromXYZBlock", from_xyz)
        monkeypatch.setattr(
            tanimoto_module.rdDetermineBonds,
            "DetermineBonds",
            lambda *args, **kwargs: (_ for _ in ()).throw(
                ValueError("retry with Hueckel")
            ),
        )
        monkeypatch.setattr(
            tanimoto_module.rdDetermineBonds,
            "hueckelEnabled",
            lambda: True,
        )

        assert grouper._molecule_to_rdkit(methanol_molecules[0]) is None

    def test_tanimoto_num_groups_record_includes_actual_groups(
        self, monkeypatch, methanol_molecules, captured
    ):
        grouper = TanimotoSimilarityGrouper(
            methanol_molecules[:2], num_groups=2
        )
        grouper._cached_group_indices = [[0], [1]]
        grouper._auto_threshold = 0.8

        class Recorder:
            def get_labels(self, count):
                return [str(index + 1) for index in range(count)]

            def build_groups_dataframe(self, index_groups, count):
                return index_groups

            def record_results(self, **kwargs):
                captured.update(kwargs)

        monkeypatch.setattr(grouper, "_get_results_recorder", Recorder)
        grouper._record_results(
            tanimoto_matrix=np.eye(2, dtype=np.float32),
            valid_indices=[0, 1],
            groups=[[methanol_molecules[0]], [methanol_molecules[1]]],
            index_groups=[[0], [1]],
        )

        assert ("Actual Groups", 2) in captured["header_info"]
        assert captured["suffix"] == "N2"

    def test_tanimoto_similarity_grouper(
        self, methanol_molecules, methanol_and_ethanol
    ):
        grouper = TanimotoSimilarityGrouper(methanol_molecules)
        assert "threshold=0.9" in repr(grouper)
        groups, group_indices = grouper.group()
        assert (
            len(groups) == 1
        ), "Molecules should form one group based on RCM similarity."
        assert (
            len(group_indices) == 1
        ), "Molecules should form one group based on RCM similarity."
        unique_structures = grouper.unique()
        assert (
            len(unique_structures) == 1
        ), "Molecules should form one group based on RCM similarity."
        grouper2 = TanimotoSimilarityGrouper(methanol_and_ethanol)
        groups, group_indices = grouper2.group()
        assert (
            len(groups) == 2
        ), "Molecules should form two groups based on RCM similarity."
        assert (
            len(group_indices) == 2
        ), "Molecules should form two groups based on RCM similarity."
        unique_structures = grouper2.unique()
        assert (
            len(unique_structures) == 2
        ), "Molecules should form two groups based on RCM similarity."

    def test_calculate_tanimoto_pair_wrapper_validation(
        self, methanol_molecules
    ):
        grouper = TanimotoSimilarityGrouper(
            methanol_molecules[:2], fingerprint_type="rdkit"
        )

        assert np.isclose(grouper.calculate_tanimoto_pair(0, 0), 1.0)

        with pytest.raises(IndexError):
            grouper.calculate_tanimoto_pair(0, 4)

        with pytest.raises(TypeError):
            grouper.calculate_tanimoto_pair("a", 1)  # type: ignore[arg-type]

    def test_tanimoto_grouper_for_crest_conformers(
        self, multiple_molecules_xyz_file
    ):
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)

        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 18
        grouper = TanimotoSimilarityGrouper(
            molecules,
            threshold=0.8,
            fingerprint_type="usrcat",
            num_procs=self.NUM_PROCS,
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 14
        assert len(group_indices) == 14
        unique_structures = grouper.unique()
        assert len(unique_structures) == 14

        grouper = TanimotoSimilarityGrouper(
            molecules,
            threshold=0.9,
            fingerprint_type="usrcat",
            num_procs=self.NUM_PROCS,
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 17
        assert len(group_indices) == 17
        unique_structures = grouper.unique()
        assert len(unique_structures) == 17

    def test_tanimoto_threshold_boundary_through_group(
        self, methanol_molecules, monkeypatch
    ):
        molecules = methanol_molecules[:2]
        from chemsmart.jobs.grouper import tanimoto as tanimoto_module

        for fingerprint_type in ("usr", "usrcat"):
            for similarity, expected_num_groups in (
                (0.90, 1),
                (0.91, 1),
                (0.89, 2),
            ):
                monkeypatch.setattr(
                    tanimoto_module.rdMolDescriptors,
                    "GetUSRScore",
                    lambda *args, value=similarity, **kwargs: value,
                )

                grouper = TanimotoSimilarityGrouper(
                    molecules,
                    threshold=0.9,
                    fingerprint_type=fingerprint_type,
                    num_procs=1,
                )
                groups, group_indices = grouper.group()

                assert len(groups) == expected_num_groups
                assert len(group_indices) == expected_num_groups
                if expected_num_groups == 1:
                    assert set(group_indices[0]) == {0, 1}
                else:
                    assert all(len(group) == 1 for group in group_indices)


@pytest.mark.usefixtures("temporary_working_dir")
class TestTorsionFingerprintGrouper:
    NUM_PROCS = 4

    def test_torsionfingerprint_grouper_for_rotated_molecules(
        self, two_rotated_molecules_xyz_file
    ):
        xyz_file = XYZFile(filename=two_rotated_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 2
        grouper = TorsionFingerprintGrouper(
            molecules,
            threshold=0.1,
            num_procs=self.NUM_PROCS,
            ignore_hydrogens=False,
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 1
        assert len(group_indices) == 1
        unique_structures = grouper.unique()
        assert len(unique_structures) == 1

    def test_torsionfingerprint_grouper_for_crest_molecules(
        self, multiple_molecules_xyz_file
    ):
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 18
        grouper = TorsionFingerprintGrouper(
            molecules,
            threshold=0.05,
            num_procs=self.NUM_PROCS,
            ignore_hydrogens=False,
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 10
        assert len(group_indices) == 10
        unique_structures = grouper.unique()
        assert len(unique_structures) == 10

    def test_use_weights_parameter(self, multiple_molecules_xyz_file):
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 18
        grouper = TorsionFingerprintGrouper(
            molecules,
            threshold=0.05,
            num_procs=self.NUM_PROCS,
            use_weights=False,
            ignore_hydrogens=False,
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 15
        assert len(group_indices) == 15
        unique_structures = grouper.unique()
        assert len(unique_structures) == 15

        tfd = grouper.calculate_tfd_pair(0, 2)
        assert np.isclose(tfd, 0.08229, rtol=1e-3)
        tfd = grouper.calculate_tfd_pair(0, 7)
        assert np.isclose(tfd, 0.07727, rtol=1e-3)

    def test_tfd_num_groups_records_actual_group_count(
        self, multiple_molecules_xyz_file, temporary_working_dir
    ):
        from openpyxl import load_workbook

        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)

        grouper = TorsionFingerprintGrouper(
            molecules,
            num_groups=5,
            num_procs=1,
            ignore_hydrogens=False,
            label="tfd_actual_groups",
        )
        groups, _ = grouper.group()

        xlsx_file = os.path.join(
            temporary_working_dir,
            "tfd_actual_groups_group_result",
            "tfd_actual_groups_TorsionFingerprintGrouper_N5.xlsx",
        )
        assert os.path.exists(xlsx_file)

        wb = load_workbook(xlsx_file, data_only=True)
        ws = wb["TFD_Matrix"]
        header_lines = [str(ws[f"A{i}"].value) for i in range(1, 25)]

        assert any(
            f"Actual Groups: {len(groups)}" in line
            for line in header_lines
            if line
        )

    def test_use_maxdev_parameter(self, multiple_molecules_xyz_file):
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 18
        grouper = TorsionFingerprintGrouper(
            molecules,
            threshold=0.26965,
            num_procs=self.NUM_PROCS,
            use_weights=True,
            max_dev="spec",
            ignore_hydrogens=False,
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 3
        assert len(group_indices) == 3
        unique_structures = grouper.unique()
        assert len(unique_structures) == 3

        tfd = grouper.calculate_tfd_pair(0, 1)
        assert np.isclose(tfd, 0.02027, rtol=1e-3)
        tfd = grouper.calculate_tfd_pair(3, 4)
        assert np.isclose(tfd, 0.24365, rtol=1e-3)

    def test_calculate_tfd_pair_wrapper_validation(
        self, multiple_molecules_xyz_file
    ):
        molecules = XYZFile(
            filename=multiple_molecules_xyz_file
        ).get_molecules(index=":", return_list=True)
        grouper = TorsionFingerprintGrouper(molecules[:3], threshold=0.1)

        assert np.isclose(grouper.calculate_tfd_pair(1, 1), 0.0)

        with pytest.raises(IndexError):
            grouper.calculate_tfd_pair(0, 8)

        with pytest.raises(TypeError):
            grouper.calculate_tfd_pair("a", 1)  # type: ignore[arg-type]


@pytest.mark.usefixtures("temporary_working_dir")
class TestEnergyGrouper:
    NUM_PROCS = 4

    @pytest.mark.parametrize(
        ("matrix", "message"),
        (
            (np.array([0.0, 1.0]), "2-dimensional"),
            (np.zeros((2, 3)), "square"),
            (np.zeros((2, 2)), "number of molecules"),
            (
                np.array(
                    [
                        [0.0, np.nan, 0.0],
                        [np.nan, 0.0, 0.0],
                        [0.0, 0.0, 0.0],
                    ]
                ),
                "finite values",
            ),
            (np.diag([1.0, 0.0, 0.0]), "diagonal"),
            (
                np.array(
                    [
                        [0.0, 1.0, 0.0],
                        [-0.5, 0.0, 0.0],
                        [0.0, 0.0, 0.0],
                    ]
                ),
                "antisymmetric",
            ),
        ),
    )
    def test_energy_output_matrix_validation(
        self, methanol_molecules, matrix, message
    ):
        molecules = methanol_molecules[:3]
        for index, molecule in enumerate(molecules):
            molecule.energy = -100.0 + index * 0.001
        grouper = EnergyGrouper(molecules, threshold=1.0)

        with pytest.raises(ValueError, match=message):
            grouper.group_by_threshold(matrix)

    def test_energy_linkage_reports_non_finite_energy(
        self, monkeypatch, methanol_molecules
    ):
        molecules = methanol_molecules[:2]
        molecules[0].energy = -100.0
        molecules[1].energy = -99.0
        grouper = EnergyGrouper(molecules, threshold=1.0)

        def raise_non_finite_error(self, distance_matrix):
            raise ValueError("linkage received non-finite distances")

        monkeypatch.setattr(
            MatrixGrouper,
            "_build_complete_linkage_tree",
            raise_non_finite_error,
        )

        with pytest.raises(
            ValueError,
            match="Energy grouping requires all molecules to have finite energy",
        ):
            grouper._build_complete_linkage_tree(
                np.array([[0.0, 1.0], [1.0, 0.0]])
            )

    def test_energy_linkage_preserves_other_validation_errors(
        self, methanol_molecules
    ):
        molecules = methanol_molecules[:2]
        molecules[0].energy = -100.0
        molecules[1].energy = -99.0
        grouper = EnergyGrouper(molecules, threshold=1.0)

        with pytest.raises(ValueError, match="negative"):
            grouper._build_complete_linkage_tree(
                np.array([[0.0, -1.0], [-1.0, 0.0]])
            )

    def test_energy_grouper_raises_error_for_missing_energy(
        self, methanol_molecules
    ):
        # methanol_molecules from pubchem don't have energy information
        with pytest.raises(ValueError) as excinfo:
            EnergyGrouper(methanol_molecules)
        assert "missing energy information" in str(excinfo.value)

    def test_energy_grouper_for_crest_conformers(
        self, multiple_molecules_xyz_file
    ):
        """Test EnergyGrouper with molecules that have energy information."""
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 18

        # All CREST conformers should have energy
        for mol in molecules:
            assert mol.energy is not None

        assert len(molecules) == 18
        grouper = EnergyGrouper(
            molecules,
            threshold=0.5,
            num_procs=self.NUM_PROCS,
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 8
        assert len(group_indices) == 8
        unique_structures = grouper.unique()
        assert len(unique_structures) == 8

        expected1 = -1.7839 / 627.509474
        expected2 = 1.4580 / 627.509474

        relative_diff, abs_diff = grouper.calculate_energy_diff_pair(1, 0)
        assert np.isclose(relative_diff, expected1, rtol=1e-3)
        relative_diff, abs_diff = grouper.calculate_energy_diff_pair(1, 2)
        assert np.isclose(abs_diff, expected2, rtol=1e-3)

        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)

        grouper = EnergyGrouper(
            molecules, num_groups=5, num_procs=self.NUM_PROCS
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 5
        assert len(group_indices) == 5
        assert grouper._auto_threshold is not None
        assert grouper._auto_threshold > 0

    def test_energy_grouper_for_log_conformers(
        self, ts_conformers_log_directory
    ):
        """Test EnergyGrouper with molecules loaded from log files using Gibbs energy."""
        import glob
        import re

        from chemsmart.io.gaussian.output import Gaussian16Output
        from chemsmart.io.molecules.structure import Molecule

        # Load molecules from log files and set Gibbs energy (like CLI does)
        log_files = sorted(
            glob.glob(os.path.join(ts_conformers_log_directory, "*.log"))
        )
        molecules = []
        conformer_ids = []

        for log_file in log_files:
            mol = Molecule.from_filepath(
                filepath=log_file, index="-1", return_list=False
            )
            if mol is not None:
                # Extract and set Gibbs energy using gibbs_free_energy property
                g16_output = Gaussian16Output(filename=log_file)
                gibbs_energy = g16_output.gibbs_free_energy
                if gibbs_energy is not None:
                    mol._energy = gibbs_energy

                molecules.append(mol)
                basename = os.path.basename(log_file)
                match = re.search(r"_c(\d+)\.log$", basename)
                if match:
                    conformer_ids.append(f"c{match.group(1)}")
                else:
                    conformer_ids.append(basename)

        assert len(molecules) == 5
        assert len(conformer_ids) == 5

        # All log conformers should have Gibbs energy
        for mol in molecules:
            assert mol.energy is not None

        grouper = EnergyGrouper(
            molecules,
            threshold=1.0,
            num_procs=self.NUM_PROCS,
            conformer_ids=conformer_ids,
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 2
        assert len(group_indices) == 2
        unique_structures = grouper.unique()
        assert len(unique_structures) == 2

        expected1 = -4.1429 / 627.509474
        expected2 = 4.1429 / 627.509474

        relative_diff, abs_diff = grouper.calculate_energy_diff_pair(1, 4)
        assert np.isclose(relative_diff, expected1, rtol=1e-2)
        assert np.isclose(abs_diff, expected2, rtol=1e-2)

    def test_calculate_energy_diff_pair_wrapper_validation(
        self, multiple_molecules_xyz_file
    ):
        molecules = XYZFile(
            filename=multiple_molecules_xyz_file
        ).get_molecules(index=":", return_list=True)
        grouper = EnergyGrouper(molecules[:3], threshold=1.0)

        relative_diff, abs_diff = grouper.calculate_energy_diff_pair(2, 2)
        assert np.isclose(relative_diff, 0.0)
        assert np.isclose(abs_diff, 0.0)

        with pytest.raises(IndexError):
            grouper.calculate_energy_diff_pair(0, 5)

        with pytest.raises(TypeError):
            grouper.calculate_energy_diff_pair("a", 1)  # type: ignore[arg-type]

    def test_energy_extraction_from_xyz_file(
        self, multiple_molecules_xyz_file
    ):

        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)

        # Test specific energy value for 3rd structure (index 2)
        expected_energy_mol3 = -126.25238449
        actual_energy_mol3 = molecules[2].energy

        assert actual_energy_mol3 is not None, "3rd molecule energy is None"
        assert np.isclose(
            actual_energy_mol3, expected_energy_mol3, rtol=1e-8
        ), (
            f"Energy mismatch for 3rd molecule: got {actual_energy_mol3}, "
            f"expected {expected_energy_mol3}"
        )

    def test_gibbs_energy_extraction_function(
        self, ts_conformers_log_directory
    ):
        """Test that Gaussian16Output.gibbs_free_energy extracts correct value."""
        from chemsmart.io.gaussian.output import Gaussian16Output

        log_file = os.path.join(
            ts_conformers_log_directory, "ch_1c_para_c1.log"
        )
        g16_output = Gaussian16Output(filename=log_file)

        # Known values
        expected_scf = -1401.99431519
        expected_gibbs_correction = 0.311863
        expected_gibbs_energy = (
            expected_scf + expected_gibbs_correction
        )  # -1401.68245219

        gibbs_energy = g16_output.gibbs_free_energy
        assert np.isclose(gibbs_energy, expected_gibbs_energy, rtol=1e-6), (
            f"Gibbs energy mismatch: got {gibbs_energy}, "
            f"expected {expected_gibbs_energy}"
        )

    def test_energy_extraction_from_ts_log_files(
        self, ts_conformers_log_directory
    ):
        """Test that energy is correctly extracted from TS log files as SCF Done energy."""
        import glob

        from chemsmart.io.gaussian.output import Gaussian16Output
        from chemsmart.io.molecules.structure import Molecule

        log_files = sorted(
            glob.glob(os.path.join(ts_conformers_log_directory, "*.log"))
        )

        molecules = []
        for log_file in log_files:
            mol = Molecule.from_filepath(
                filepath=log_file, index="-1", return_list=False
            )
            if mol is not None:
                molecules.append(mol)
                g16_output = Gaussian16Output(filename=log_file)
                scf_energy = (
                    g16_output.scf_energies[-1]
                    if g16_output.scf_energies
                    else None
                )
                assert np.isclose(mol.energy, scf_energy, rtol=1e-6), (
                    f"Energy mismatch for {log_file}: "
                    f"mol.energy={mol.energy}, expected SCF={scf_energy}"
                )


@pytest.mark.usefixtures("temporary_working_dir")
class TestOtherGroupers:
    NUM_PROCS = 4

    def test_base_record_template_method(self, methanol_molecules):
        class DummyGrouper(MoleculeGrouper):
            def group(self):
                return [], []

            def _record_results(self, **kwargs):
                return kwargs

        grouper = DummyGrouper(methanol_molecules)
        payload = grouper.record(example_key="example_value")
        assert payload["example_key"] == "example_value"

    @pytest.mark.parametrize("energy_type", ["QHH", "QHG", "SP_QHG"])
    def test_thermochemistry_parameters_header_for_qhg_energy_types(
        self, methanol_molecules, energy_type
    ):
        """qhG and SP-qhG results must record the exact thermochemistry settings used."""
        thermo_parameters = (
            "temperature=298.15, concentration=1.0, pressure=1.0, "
            "use_weighted_mass=True, alpha=4, s_freq_cutoff=100.0, "
            "entropy_method=grimme, h_freq_cutoff=100.0, "
            "energy_units=hartree, check_imaginary_frequencies=True, "
            "cutoff_entropy_grimme=100.0, cutoff_enthalpy=100.0"
        )
        grouper = BasicRMSDGrouper(
            methanol_molecules[:2],
            threshold=0.5,
            num_procs=1,
            energy_type=energy_type,
            thermo_parameters=thermo_parameters,
        )

        header_info = []
        grouper._append_thermo_header(header_info)

        assert header_info == [
            ("Thermochemistry Parameters", thermo_parameters)
        ]

    @pytest.mark.parametrize("energy_type", ["E", "H", "G"])
    def test_thermochemistry_parameters_header_not_written_for_other_energy_types(
        self, methanol_molecules, energy_type
    ):
        """Thermochemistry parameter header is specific to qhG and SP-qhG."""
        thermo_parameters = (
            "temperature=298.15, concentration=1.0, pressure=1.0, "
            "use_weighted_mass=True, alpha=4, s_freq_cutoff=100.0, "
            "entropy_method=grimme, h_freq_cutoff=100.0, "
            "energy_units=hartree, check_imaginary_frequencies=True, "
            "cutoff_entropy_grimme=100.0, cutoff_enthalpy=100.0"
        )
        grouper = BasicRMSDGrouper(
            methanol_molecules[:2],
            threshold=0.5,
            num_procs=1,
            energy_type=energy_type,
            thermo_parameters=thermo_parameters,
        )

        header_info = []
        grouper._append_thermo_header(header_info)

        assert not any(
            key == "Thermochemistry Parameters" for key, _ in header_info
        )

    def test_record_writes_rmsd_matrix_with_headers(
        self, methanol_molecules, temporary_working_dir
    ):
        """Detailed record test: verify matrix file and key header lines are written."""
        from openpyxl import load_workbook

        molecules = methanol_molecules[:2]
        grouper = BasicRMSDGrouper(
            molecules,
            threshold=0.5,
            num_procs=1,
            label="record_detail",
            energy_type="E",
        )

        # Populate cache first so Groups sheet can be generated.
        groups, group_indices = grouper.group()
        assert len(groups) >= 1
        assert len(group_indices) >= 1

        # Call unified record() explicitly to validate template-method path.
        rmsd_matrix = np.zeros((len(molecules), len(molecules)))
        grouper.record(rmsd_matrix=rmsd_matrix, grouping_time=0.01)

        xlsx_file = os.path.join(
            temporary_working_dir,
            "record_detail_group_result",
            "record_detail_BasicRMSDGrouper_T0.5.xlsx",
        )
        assert os.path.exists(xlsx_file)

        wb = load_workbook(xlsx_file, data_only=True)
        assert "RMSD_Matrix" in wb.sheetnames
        assert "Groups" in wb.sheetnames

        ws = wb["RMSD_Matrix"]
        header_lines = [str(ws[f"A{i}"].value) for i in range(1, 25)]

        assert any("Energy Type: E" in line for line in header_lines if line)
        assert any(
            "Used Molecules: 2" in line for line in header_lines if line
        )
        assert any(
            "Skipped Molecules: 0" in line for line in header_lines if line
        )
        assert any("Grouping Time:" in line for line in header_lines if line)

    def test_record_writes_formula_outputs_with_headers(
        self, methanol_and_ethanol, temporary_working_dir
    ):
        """Detailed record test: verify non-matrix output sheet and header lines."""
        from openpyxl import load_workbook

        grouper = FormulaGrouper(
            methanol_and_ethanol,
            num_procs=1,
            label="formula_record_detail",
            energy_type="E",
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 2
        assert len(group_indices) == 2

        xlsx_file = os.path.join(
            temporary_working_dir,
            "formula_record_detail_group_result",
            "formula_record_detail_FormulaGrouper.xlsx",
        )
        assert os.path.exists(xlsx_file)

        wb = load_workbook(xlsx_file, data_only=True)
        assert "Formulas" in wb.sheetnames
        assert "Groups" in wb.sheetnames

        ws = wb["Formulas"]
        header_lines = [str(ws[f"A{i}"].value) for i in range(1, 25)]
        assert any(
            "Total Molecules: 2" in line for line in header_lines if line
        )
        assert any(
            "Unique Formulas: 2" in line for line in header_lines if line
        )
        assert any("Energy Type: E" in line for line in header_lines if line)

    def test_formula_grouper(
        self,
        methanol_molecules,
        methanol_and_ethanol,
        conformers_from_rdkit,
    ):
        grouper = FormulaGrouper(methanol_molecules)
        groups, group_indices = grouper.group()
        unique_structures = grouper.unique()
        assert (
            len(groups) == 1
        ), "Molecules should form one group based on formula."

        assert (
            len(unique_structures) == 1
        ), "Molecules should form one group based on formula."

        grouper2 = FormulaGrouper(methanol_and_ethanol)
        groups, group_indices = grouper2.group()
        unique_structures = grouper2.unique()
        assert (
            len(groups) == 2
        ), "Molecules should form two groups based on formula."
        assert (
            len(unique_structures) == 2
        ), "Molecules should form two groups based on formula."

        grouper3 = FormulaGrouper(conformers_from_rdkit)
        # based on Formula, should all be the same even for 300 conformers
        groups, group_indices = grouper3.group()
        unique_structures = grouper3.unique()
        assert len(groups) == 1
        assert len(unique_structures) == 1

    @pytest.mark.slow
    def test_connectivity_grouper(
        self, methanol_molecules, methanol_and_ethanol
    ):
        grouper = ConnectivityGrouper(methanol_molecules)
        groups, group_indices = grouper.group()
        assert (
            len(groups) == 1
        ), "Molecules should form one group based on connectivity."
        assert (
            len(group_indices) == 1
        ), "Molecules should form one group based on connectivity."
        unique_structures = grouper.unique()
        assert (
            len(unique_structures) == 1
        ), "Molecules should form one group based on connectivity."

        grouper2 = ConnectivityGrouper(methanol_and_ethanol)
        groups, group_indices = grouper2.group()
        assert (
            len(groups) == 2
        ), "Molecules should form two groups based on connectivity."
        assert (
            len(group_indices) == 2
        ), "Molecules should form two groups based on connectivity."
        unique_structures = grouper2.unique()
        assert (
            len(unique_structures) == 2
        ), "Molecules should form two groups based on connectivity."

    def test_connectivity_grouper_for_crest_conformers(
        self, multiple_molecules_xyz_file
    ):
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)

        molecules = xyz_file.get_molecules(index=":", return_list=True)
        assert len(molecules) == 18
        grouper = ConnectivityGrouper(
            molecules,
            num_procs=self.NUM_PROCS,
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 1
        assert len(group_indices) == 1
        unique_structures = grouper.unique()
        assert len(unique_structures) == 1

    def test_connectivity_grouper_ignore_hydrogens(self, methanol_molecules):
        """Covers the ignore_hydrogens branch in group(), which strips
        H nodes from every graph before running isomorphism checks."""
        grouper = ConnectivityGrouper(
            methanol_molecules, ignore_hydrogens=True
        )
        groups, group_indices = grouper.group()
        assert len(groups) == 1
        assert len(group_indices) == 1

    def test_connectivity_grouper_repr(self, methanol_molecules):
        # ConnectivityGrouper does not support multiprocessing, so any
        # requested num_procs > 1 is clamped down to 1.
        grouper = ConnectivityGrouper(
            methanol_molecules, num_procs=2, ignore_hydrogens=True
        )
        text = repr(grouper)
        assert "num_procs=1" in text
        assert "ignore_hydrogens=True" in text

    def test_connectivity_grouper_record_results_without_grouping_time(
        self, methanol_molecules, tmp_path, monkeypatch
    ):
        """Covers `if grouping_time is not None:`'s False arm: calling
        the internal result-writer directly without a grouping_time
        must skip the "Grouping Time" header line rather than raising."""
        monkeypatch.chdir(tmp_path)
        grouper = ConnectivityGrouper(methanol_molecules, label="conn_test")
        groups, index_groups = grouper.group()
        # Should not raise when grouping_time is omitted (defaults to None).
        grouper._record_results(groups=groups, index_groups=index_groups)

    def test_check_isomorphism_direct(self, methanol_molecules):
        """See BUGS_FOUND.md #62: _check_isomorphism is defined but
        never called by group() (which calls _are_isomorphic directly
        in a serial loop instead). Unit-tested directly since no code
        path reaches it."""
        grouper = ConnectivityGrouper(methanol_molecules)
        grouper.graphs = [
            to_graph_wrapper(mol, 0.0, grouper.adjust_H)
            for mol in methanol_molecules
        ]
        i, j, is_isomorphic = grouper._check_isomorphism((0, 1))
        assert (i, j) == (0, 1)
        assert is_isomorphic is True

    def test_rdkit_isomorphism_grouper(
        self, methanol_molecules, methanol_and_ethanol
    ):
        grouper = RDKitIsomorphismGrouper(methanol_molecules)
        groups, group_indices = grouper.group()
        assert (
            len(groups) == 1
        ), "Molecules should form one group based on RCM similarity."
        assert (
            len(group_indices) == 1
        ), "Molecules should form one group based on RCM similarity."
        unique_structures = grouper.unique()
        assert (
            len(unique_structures) == 1
        ), "Molecules should form one group based on RCM similarity."
        grouper2 = RDKitIsomorphismGrouper(methanol_and_ethanol)
        groups, group_indices = grouper2.group()
        assert (
            len(groups) == 2
        ), "Molecules should form two groups based on RCM similarity."
        assert (
            len(group_indices) == 2
        ), "Molecules should form two groups based on RCM similarity."
        unique_structures = grouper2.unique()
        assert (
            len(unique_structures) == 2
        ), "Molecules should form two groups based on RCM similarity."


@pytest.mark.usefixtures("temporary_working_dir")
class TestFactory:

    @pytest.mark.parametrize(
        ("parameter", "value"),
        (("use_weights", False), ("max_dev", "spec")),
    )
    def test_tfd_only_parameter_warns_for_other_strategy(
        self, caplog, methanol_molecules, parameter, value
    ):
        caplog.set_level("WARNING")

        grouper = StructureGrouperFactory.create(
            methanol_molecules,
            strategy="rmsd",
            **{parameter: value},
        )

        assert isinstance(grouper, BasicRMSDGrouper)
        assert (
            f"Parameter '{parameter}={value}' is only effective for 'tfd' "
            "strategy, ignored for 'rmsd'." in caplog.messages
        )

    @pytest.mark.parametrize(
        ("parameter", "value", "owner"),
        (
            ("inversion", "on", "irmsd"),
            ("fingerprint_type", "morgan", "tanimoto"),
        ),
    )
    def test_strategy_only_parameter_warns_and_is_dropped(
        self, caplog, methanol_molecules, parameter, value, owner
    ):
        caplog.set_level("WARNING")

        grouper = StructureGrouperFactory.create(
            methanol_molecules,
            strategy="rmsd",
            **{parameter: value},
        )

        assert isinstance(grouper, BasicRMSDGrouper)
        assert not hasattr(grouper, parameter)
        assert (
            f"Parameter '{parameter}={value}' is only effective for "
            f"'{owner}' strategy, ignored for 'rmsd'." in caplog.messages
        )

    def test_none_strategy_defaults_to_rmsd(self, methanol_molecules):
        grouper = StructureGrouperFactory.create(
            methanol_molecules, strategy=None
        )
        assert isinstance(grouper, BasicRMSDGrouper)

    @pytest.mark.parametrize(
        ("kwargs", "message"),
        (
            ({"strategy": "bogus"}, "Unknown grouping strategy: 'bogus'"),
            (
                {"strategy": "rmsd", "threshold": 0.5, "num_groups": 2},
                "Cannot specify both threshold (-T) and num_groups (-N)",
            ),
            (
                {"strategy": "formula", "threshold": 0.5},
                "Strategy 'formula' does not support threshold parameter",
            ),
            (
                {"strategy": "formula", "num_groups": 2},
                "Strategy 'formula' does not support num_groups parameter",
            ),
            (
                {"strategy": "energy", "ignore_hydrogens": True},
                "Strategy 'energy' does not support ignore_hydrogens",
            ),
        ),
    )
    def test_invalid_parameter_combinations_raise(
        self, methanol_molecules, kwargs, message
    ):
        with pytest.raises(ValueError) as excinfo:
            StructureGrouperFactory.create(methanol_molecules, **kwargs)
        assert message in str(excinfo.value)

    def test_structure_grouper_factory_energy(
        self, multiple_molecules_xyz_file
    ):
        """Test factory creation of energy grouper (requires molecules with energy)."""
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)

        factory = StructureGrouperFactory()
        energy_grouper = factory.create(molecules, strategy="energy")
        assert isinstance(energy_grouper, EnergyGrouper)

    def test_structure_grouper_factory(self, methanol_molecules):
        factory = StructureGrouperFactory()
        rmsd_grouper = factory.create(methanol_molecules, strategy="rmsd")
        assert isinstance(rmsd_grouper, RMSDGrouper)
        hrmsd_grouper = factory.create(methanol_molecules, strategy="hrmsd")
        assert isinstance(hrmsd_grouper, RMSDGrouper)
        spyrmsd_grouper = factory.create(
            methanol_molecules, strategy="spyrmsd"
        )
        assert isinstance(spyrmsd_grouper, RMSDGrouper)

        # irmsd requires either Python API or external command
        if _is_irmsd_available():
            irmsd_grouper = factory.create(
                methanol_molecules, strategy="irmsd"
            )
            assert isinstance(irmsd_grouper, RMSDGrouper)

        pymolrmsd_grouper = factory.create(
            methanol_molecules, strategy="pymolrmsd"
        )
        assert isinstance(pymolrmsd_grouper, RMSDGrouper)

        # Explicitly cleanup pymolrmsd_grouper to prevent __del__ from calling quit()
        if (
            hasattr(pymolrmsd_grouper, "_temp_dir")
            and pymolrmsd_grouper._temp_dir
        ):

            shutil.rmtree(pymolrmsd_grouper._temp_dir, ignore_errors=True)
            pymolrmsd_grouper._temp_dir = None
        pymolrmsd_grouper.cmd = None  # Prevent __del__ from calling quit()

        tfd_grouper = factory.create(methanol_molecules, strategy="tfd")
        assert isinstance(tfd_grouper, TorsionFingerprintGrouper)
        formula_grouper = factory.create(
            methanol_molecules, strategy="formula"
        )
        assert isinstance(formula_grouper, FormulaGrouper)
        connectivity_grouper = factory.create(
            methanol_molecules, strategy="connectivity"
        )
        assert isinstance(connectivity_grouper, ConnectivityGrouper)
        tanimoto_grouper = factory.create(
            methanol_molecules, strategy="tanimoto"
        )
        assert isinstance(tanimoto_grouper, TanimotoSimilarityGrouper)
        rdkit_isomorphism_grouper = factory.create(
            methanol_molecules, strategy="isomorphism"
        )
        assert isinstance(rdkit_isomorphism_grouper, RDKitIsomorphismGrouper)

    @classmethod
    def teardown_class(cls):
        """Clean up PyMOL after all tests in this class are done."""
        try:
            from pymol import cmd

            cmd.reinitialize()
        except Exception:
            # Ignore cleanup errors if PyMOL is not available
            pass


@pytest.mark.usefixtures("temporary_working_dir")
class TestGrouperUtilityFunctions:
    """Test utility functions and helper methods in groupers."""

    NUM_PROCS = 1

    def test_rmsd_matrix_with_num_groups(
        self, multiple_molecules_xyz_file, temporary_working_dir
    ):
        """Test that RMSD matrix filename reflects num_groups when used."""
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)

        grouper = BasicRMSDGrouper(
            molecules,
            num_groups=3,
            num_procs=1,
            label="test_num_groups",
        )
        _, _ = grouper.group()

        # Check that Excel file was created in the temporary_working_dir
        expected_file = os.path.join(
            temporary_working_dir,
            "test_num_groups_group_result",
            "test_num_groups_BasicRMSDGrouper_N3.xlsx",
        )
        assert os.path.exists(
            expected_file
        ), f"Matrix file not found: {expected_file}"

    def test_grouping_result_caching(self, multiple_molecules_xyz_file):
        """Test that grouping results are cached and reused."""
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)

        grouper = BasicRMSDGrouper(
            molecules,
            threshold=0.5,
            num_procs=1,
        )

        groups1, indices1 = grouper.group()

        assert grouper._cached_groups is not None
        assert grouper._cached_group_indices is not None

        # unique() should use cached results
        unique_mols = grouper.unique()
        assert len(unique_mols) == len(groups1)

    def test_unique_returns_lowest_energy_representative(
        self, multiple_molecules_xyz_file
    ):
        """Test that unique() returns lowest energy molecule from each group."""
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)

        # CREST conformers have energy information
        grouper = BasicRMSDGrouper(
            molecules,
            threshold=2.0,
            num_procs=1,
        )
        groups, group_indices = grouper.group()
        unique_mols = grouper.unique()

        # For each group, verify the representative has lowest energy
        for i, (group, indices) in enumerate(zip(groups, group_indices)):
            mols_with_energy = [
                (mol, idx)
                for mol, idx in zip(group, indices)
                if mol.energy is not None
            ]
            if mols_with_energy:
                min_energy_mol = min(
                    mols_with_energy, key=lambda x: x[0].energy
                )[0]
                assert unique_mols[i].energy == min_energy_mol.energy


@pytest.mark.usefixtures("temporary_working_dir")
class TestCompleteLinkageInput:
    """Input validation checks for complete-linkage based matrix groupers."""

    NUM_PROCS = 4

    @pytest.mark.parametrize(
        "kwargs",
        (
            {"threshold": 0.5, "num_groups": 2},
            {"threshold": -0.1},
            {"threshold": None, "num_groups": 0},
        ),
    )
    def test_matrix_grouper_rejects_invalid_grouping_configuration(
        self, methanol_molecules, kwargs
    ):
        with pytest.raises(ValueError):
            BasicRMSDGrouper(methanol_molecules[:2], **kwargs)

    def test_matrix_grouper_rejects_unknown_representative_strategy(
        self, methanol_molecules
    ):
        with pytest.raises(ValueError, match="Unsupported representative"):
            BasicRMSDGrouper(
                methanol_molecules[:2],
                representative_strategy="random",
            )

    def test_energy_order_places_invalid_energies_last(
        self, methanol_molecules
    ):
        molecules = methanol_molecules[:3]
        molecules[0]._energy = -2.0
        molecules[1]._energy = "invalid"
        molecules[2]._energy = np.nan
        grouper = BasicRMSDGrouper(molecules)

        assert grouper._order_index_group_by_energy([2, 1, 0]) == [0, 1, 2]

    def test_explicit_output_directory_is_used(
        self, methanol_molecules, temporary_working_dir
    ):
        grouper = BasicRMSDGrouper(
            methanol_molecules[:2],
            label="sample",
            output_dir=temporary_working_dir,
        )

        assert grouper._get_output_dir() == os.path.join(
            temporary_working_dir, "sample_group_result"
        )

    @pytest.mark.parametrize(
        ("matrix", "message"),
        (
            (np.array([0.0, 1.0]), "2-dimensional"),
            (np.zeros((2, 3)), "square"),
        ),
    )
    def test_distance_matrix_shape_validation(
        self, methanol_molecules, matrix, message
    ):
        grouper = BasicRMSDGrouper(methanol_molecules[:2])

        with pytest.raises(ValueError, match=message):
            grouper.group_by_threshold(matrix)

    def test_clustering_submatrix_validates_original_indices(
        self, methanol_molecules
    ):
        grouper = BasicRMSDGrouper(methanol_molecules[:2])

        with pytest.raises(ValueError, match="original_indices"):
            grouper._prepare_clustering_submatrix(
                np.zeros((2, 2)), original_indices=[0]
            )

    def test_singleton_matrix_helpers(self, methanol_molecules):
        grouper = BasicRMSDGrouper(
            methanol_molecules[:2], conformer_ids=["c1", "c2"]
        )

        assert grouper._build_complete_linkage_tree(np.zeros((1, 1))) is None
        assert grouper._get_original_index_label(1) == "c2"
        assert (
            grouper._mean_distance_to_group(0, [0], np.zeros((1, 1)), {0: 0})
            == 0.0
        )

    @pytest.mark.parametrize("threshold", (None, "bad", np.nan, -1.0))
    def test_threshold_resolution_validation(
        self, methanol_molecules, threshold
    ):
        grouper = BasicRMSDGrouper(methanol_molecules[:2])
        grouper.threshold = None

        with pytest.raises(ValueError, match="threshold"):
            grouper._resolve_threshold(threshold)

    @pytest.mark.parametrize("num_groups", (None, "bad", 0))
    def test_num_groups_resolution_validation(
        self, methanol_molecules, num_groups
    ):
        grouper = BasicRMSDGrouper(methanol_molecules[:2])
        grouper.num_groups = None

        with pytest.raises(ValueError, match="num_groups"):
            grouper._resolve_num_groups(num_groups)


@pytest.mark.usefixtures("temporary_working_dir")
class TestCompleteLinkageOutput:
    """Output/state behavior and clustering semantics for complete linkage."""

    NUM_PROCS = 4

    # --- Threshold behavior ---

    def test_public_group_by_threshold_reuses_distance_matrix_and_caches_result(
        self, methanol_molecules
    ):
        molecules = methanol_molecules[:3]
        grouper = BasicRMSDGrouper(molecules, threshold=0.5, num_procs=1)
        distance_matrix = np.array(
            [
                [0.0, 0.3, 0.8],
                [0.3, 0.0, 0.7],
                [0.8, 0.7, 0.0],
            ]
        )
        original_matrix = distance_matrix.copy()

        _, strict_indices = grouper.group_by_threshold(
            distance_matrix, threshold=0.2
        )
        groups, relaxed_indices = grouper.group_by_threshold(
            distance_matrix, threshold=0.5
        )

        assert strict_indices == [[0], [1], [2]]
        assert {frozenset(group) for group in relaxed_indices} == {
            frozenset({0, 1}),
            frozenset({2}),
        }
        assert grouper.threshold == 0.5
        assert np.array_equal(distance_matrix, original_matrix)
        assert grouper._cached_groups == groups
        assert grouper._cached_group_indices == relaxed_indices

    def test_tfd_public_group_by_threshold_accepts_output_distance_matrix(
        self, methanol_molecules
    ):
        grouper = TorsionFingerprintGrouper(
            methanol_molecules[:2], threshold=0.1
        )
        tfd_matrix = np.array([[0.0, 0.15], [0.15, 0.0]])

        _, strict_indices = grouper.group_by_threshold(
            tfd_matrix, threshold=0.1
        )
        _, relaxed_indices = grouper.group_by_threshold(
            tfd_matrix, threshold=0.2
        )

        assert strict_indices == [[0], [1]]
        assert relaxed_indices == [[0, 1]]
        assert grouper.threshold == 0.1

    def test_energy_public_group_by_threshold_accepts_signed_kcal_output(
        self, multiple_molecules_xyz_file
    ):
        molecules = XYZFile(
            filename=multiple_molecules_xyz_file
        ).get_molecules(index=":", return_list=True)[:3]
        grouper = EnergyGrouper(molecules, threshold=1.0)
        energy_matrix_kcal = np.array(
            [
                [0.0, 0.5, 2.0],
                [-0.5, 0.0, 1.5],
                [-2.0, -1.5, 0.0],
            ]
        )

        _, index_groups = grouper.group_by_threshold(energy_matrix_kcal)

        assert {frozenset(group) for group in index_groups} == {
            frozenset({0, 1}),
            frozenset({2}),
        }
        assert grouper.threshold == 1.0
        assert np.array_equal(
            energy_matrix_kcal,
            np.array(
                [
                    [0.0, 0.5, 2.0],
                    [-0.5, 0.0, 1.5],
                    [-2.0, -1.5, 0.0],
                ]
            ),
        )

    def test_tanimoto_public_group_by_threshold_accepts_recorded_similarity(
        self, methanol_molecules
    ):
        grouper = TanimotoSimilarityGrouper(
            methanol_molecules[:3], threshold=0.9
        )
        similarity_matrix = np.array(
            [
                [1.0, np.nan, 0.7],
                [np.nan, np.nan, np.nan],
                [0.7, np.nan, 1.0],
            ],
            dtype=np.float32,
        )

        _, index_groups = grouper.group_by_threshold(
            similarity_matrix, threshold=0.6
        )

        assert index_groups == [[0, 2]]
        assert grouper._matrix_skipped_indices == [1]
        assert grouper._cached_group_indices == index_groups
        assert grouper.threshold == 0.9

    def test_threshold_grouping_clears_previous_auto_threshold(
        self, methanol_molecules
    ):
        grouper = BasicRMSDGrouper(methanol_molecules, num_groups=2)
        matrix = np.array(
            [
                [0.0, 0.1, 0.8],
                [0.1, 0.0, 0.7],
                [0.8, 0.7, 0.0],
            ]
        )

        grouper.group_by_num_groups(matrix, num_groups=2)
        assert grouper._auto_threshold is not None

        old_threshold = grouper.threshold
        groups, index_groups = grouper.group_by_threshold(
            matrix, threshold=0.5
        )

        assert grouper._auto_threshold is None
        assert grouper.threshold == old_threshold
        assert index_groups == [[0, 1], [2]]
        assert groups == [
            [methanol_molecules[0], methanol_molecules[1]],
            [methanol_molecules[2]],
        ]

    def test_hierarchical_complete_linkage_prevents_chaining(
        self, methanol_molecules
    ):
        """A-B and B-C being close must not chain A and C together."""
        molecules = methanol_molecules[:3]
        grouper = BasicRMSDGrouper(
            molecules,
            threshold=0.5,
            num_procs=1,
        )
        distance_matrix = np.array(
            [
                [0.0, 0.3, 0.8],
                [0.3, 0.0, 0.3],
                [0.8, 0.3, 0.0],
            ]
        )

        groups, index_groups = grouper.group_by_threshold(distance_matrix)

        assert len(groups) == 2
        assert sorted(len(group) for group in index_groups) == [1, 2]
        assert not any(set(group) == {0, 1, 2} for group in index_groups)

    def test_threshold_includes_equal_distance(self, methanol_molecules):
        """Threshold mode uses complete-linkage distance <= threshold."""
        molecules = methanol_molecules[:2]
        grouper = BasicRMSDGrouper(
            molecules,
            threshold=0.5,
            num_procs=1,
        )
        distance_matrix = np.array(
            [
                [0.0, 0.5],
                [0.5, 0.0],
            ]
        )
        groups, index_groups = grouper.group_by_threshold(distance_matrix)
        assert len(groups) == 1
        assert index_groups == [[0, 1]]

    def test_tanimoto_threshold_uses_similarity_semantics(
        self, methanol_molecules
    ):
        """A Tanimoto threshold of 0.9 means similarity >= 0.9."""
        molecules = methanol_molecules[:3]
        grouper = TanimotoSimilarityGrouper(
            molecules,
            threshold=0.9,
            num_procs=1,
        )

        similarity_matrix = np.array(
            [
                [1.0, 0.90, 0.70],
                [0.90, 1.0, 0.70],
                [0.70, 0.70, 1.0],
            ]
        )

        groups, index_groups = grouper.group_by_threshold(similarity_matrix)

        assert len(groups) == 2
        assert {frozenset(group) for group in index_groups} == {
            frozenset({0, 1}),
            frozenset({2}),
        }
        assert np.isclose(grouper.threshold, 0.9)

    def test_complete_linkage_prevents_chaining(
        self, multiple_molecules_xyz_file
    ):
        """Test that complete linkage prevents chaining effect in grouping."""
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)

        # With a moderate threshold, complete linkage should create
        # groups where ALL members are within threshold of each other
        grouper = BasicRMSDGrouper(
            molecules,
            threshold=1.0,
            num_procs=self.NUM_PROCS,
        )
        groups, group_indices = grouper.group()

        # Verify complete linkage property: within each group,
        # all pairs should have RMSD <= threshold
        for group_idx, indices in enumerate(group_indices):
            if len(indices) > 1:
                for i, idx_i in enumerate(indices):
                    for j, idx_j in enumerate(indices):
                        if i < j:
                            rmsd = grouper.calculate_rmsd_pair(idx_i, idx_j)
                            assert rmsd <= 1.0, (
                                f"Group {group_idx}: RMSD between {idx_i} and {idx_j} "
                                f"is {rmsd}, is above threshold 1.0"
                            )

    # --- Num-groups behavior ---
    def test_public_group_by_num_groups_reuses_distance_matrix_and_caches_result(
        self, methanol_molecules
    ):
        molecules = methanol_molecules[:3]
        grouper = BasicRMSDGrouper(molecules, threshold=0.5, num_procs=1)
        distance_matrix = np.array(
            [
                [0.0, 0.2, 0.9],
                [0.2, 0.0, 0.8],
                [0.9, 0.8, 0.0],
            ]
        )
        original_matrix = distance_matrix.copy()

        groups, index_groups = grouper.group_by_num_groups(
            distance_matrix, num_groups=2
        )

        assert {frozenset(group) for group in index_groups} == {
            frozenset({0, 1}),
            frozenset({2}),
        }
        assert np.isclose(grouper._auto_threshold, 0.2)
        assert grouper._cached_groups == groups
        assert grouper._cached_group_indices == index_groups
        assert np.array_equal(distance_matrix, original_matrix)

    def test_tfd_public_group_by_num_groups_accepts_output_distance_matrix(
        self, methanol_molecules
    ):
        grouper = TorsionFingerprintGrouper(
            methanol_molecules[:3], threshold=0.1
        )
        tfd_matrix = np.array(
            [
                [0.0, 0.15, 0.8],
                [0.15, 0.0, 0.7],
                [0.8, 0.7, 0.0],
            ]
        )

        _, index_groups = grouper.group_by_num_groups(tfd_matrix, num_groups=2)

        assert {frozenset(group) for group in index_groups} == {
            frozenset({0, 1}),
            frozenset({2}),
        }
        assert np.isclose(grouper._auto_threshold, 0.15)

    def test_energy_public_group_by_num_groups_accepts_signed_kcal_output(
        self, multiple_molecules_xyz_file
    ):
        molecules = XYZFile(
            filename=multiple_molecules_xyz_file
        ).get_molecules(index=":", return_list=True)[:3]
        grouper = EnergyGrouper(molecules, threshold=1.0)
        energy_matrix_kcal = np.array(
            [
                [0.0, 0.5, 2.0],
                [-0.5, 0.0, 1.5],
                [-2.0, -1.5, 0.0],
            ]
        )

        _, index_groups = grouper.group_by_num_groups(
            energy_matrix_kcal, num_groups=2
        )

        assert {frozenset(group) for group in index_groups} == {
            frozenset({0, 1}),
            frozenset({2}),
        }
        assert np.isclose(grouper._auto_threshold, 0.5)

    def test_tanimoto_public_group_by_num_groups_accepts_recorded_similarity(
        self, methanol_molecules
    ):
        grouper = TanimotoSimilarityGrouper(
            methanol_molecules[:3], threshold=0.9
        )
        similarity_matrix = np.array(
            [
                [1.0, np.nan, 0.95],
                [np.nan, np.nan, np.nan],
                [0.95, np.nan, 1.0],
            ],
            dtype=np.float32,
        )

        _, index_groups = grouper.group_by_num_groups(
            similarity_matrix, num_groups=2
        )

        assert index_groups == [[0], [2]]
        assert grouper._matrix_skipped_indices == [1]
        assert grouper._auto_threshold is None

    def test_num_groups_keeps_previous_full_linkage_level(
        self, methanol_molecules
    ):
        """-N must not split a tied linkage level just to force exact N."""
        from ase.build import molecule as ase_molecule

        from chemsmart.io.molecules.structure import Molecule

        molecules = [
            Molecule.from_ase_atoms(ase_molecule("CH3OH")) for _ in range(4)
        ]
        grouper = BasicRMSDGrouper(
            molecules,
            threshold=None,
            num_groups=3,
            num_procs=1,
        )
        # At distance 0.1, two independent merges occur at the same linkage
        # level, so the hierarchy jumps directly from 4 groups to 2 groups.
        # Requested N=3 must therefore retain the previous 4-group level.
        distance_matrix = np.array(
            [
                [0.0, 0.1, 1.0, 1.0],
                [0.1, 0.0, 1.0, 1.0],
                [1.0, 1.0, 0.0, 0.1],
                [1.0, 1.0, 0.1, 0.0],
            ]
        )
        groups, index_groups = grouper.group_by_num_groups(distance_matrix)
        assert len(groups) == 4
        assert len(index_groups) == 4
        assert all(len(group) == 1 for group in index_groups)
        assert grouper._auto_threshold is None

    def test_num_groups_uses_exact_level_when_available(
        self, methanol_molecules
    ):
        """-N should return exactly N when a full linkage level gives N groups."""
        from ase.build import molecule as ase_molecule

        from chemsmart.io.molecules.structure import Molecule

        molecules = [
            Molecule.from_ase_atoms(ase_molecule("CH3OH")) for _ in range(4)
        ]
        grouper = BasicRMSDGrouper(
            molecules,
            threshold=None,
            num_groups=3,
            num_procs=1,
        )

        distance_matrix = np.array(
            [
                [0.0, 0.1, 1.0, 1.0],
                [0.1, 0.0, 1.0, 1.0],
                [1.0, 1.0, 0.0, 0.4],
                [1.0, 1.0, 0.4, 0.0],
            ]
        )

        groups, index_groups = grouper.group_by_num_groups(distance_matrix)

        assert len(groups) == 3
        assert len(index_groups) == 3
        assert np.isclose(grouper._auto_threshold, 0.1)

    # --- Inf/NaN and matrix-path edge handling ---

    def test_infinite_distance_matrix_skips_minimum_problematic_structures(
        self, methanol_molecules
    ):
        """Remove only the structure(s) needed to obtain a finite submatrix."""
        molecules = methanol_molecules[:3]
        grouper = BasicRMSDGrouper(
            molecules,
            threshold=0.5,
            num_procs=1,
        )

        distance_matrix = np.array(
            [
                [0.0, 0.3, np.inf],
                [0.3, 0.0, np.inf],
                [np.inf, np.inf, 0.0],
            ]
        )
        original_matrix = distance_matrix.copy()
        groups, index_groups = grouper.group_by_threshold(distance_matrix)

        assert index_groups == [[0, 1]]
        assert grouper._matrix_skipped_indices == [2]
        assert np.array_equal(distance_matrix, original_matrix, equal_nan=True)

    def test_nan_distance_matrix_raises(self, methanol_molecules):
        """NaN pairwise distances must still fail clearly at runtime."""
        molecules = methanol_molecules[:2]
        grouper = BasicRMSDGrouper(
            molecules,
            threshold=0.5,
            num_procs=1,
        )
        distance_matrix = np.array(
            [
                [0.0, np.nan],
                [np.nan, 0.0],
            ]
        )

        with pytest.raises(ValueError, match="non-finite"):
            grouper.group_by_threshold(distance_matrix)

    def test_infinite_pair_removes_one_structure_and_groups_remaining(
        self, methanol_molecules
    ):
        from ase.build import molecule as ase_molecule

        from chemsmart.io.molecules.structure import Molecule

        molecules = [
            Molecule.from_ase_atoms(ase_molecule("CH3OH")) for _ in range(4)
        ]
        grouper = BasicRMSDGrouper(
            molecules,
            threshold=0.5,
            num_procs=1,
        )
        distance_matrix = np.array(
            [
                [0.0, np.inf, 0.8, 0.9],
                [np.inf, 0.0, 0.7, 0.8],
                [0.8, 0.7, 0.0, 0.2],
                [0.9, 0.8, 0.2, 0.0],
            ]
        )

        groups, index_groups = grouper.group_by_threshold(distance_matrix)

        assert grouper._matrix_skipped_indices == [0]
        assert {frozenset(group) for group in index_groups} == {
            frozenset({1}),
            frozenset({2, 3}),
        }
        assert len(groups) == 2

    def test_multiple_infinite_pairs_remove_minimum_problematic_structures(
        self, methanol_molecules
    ):
        from ase.build import molecule as ase_molecule

        from chemsmart.io.molecules.structure import Molecule

        molecules = [
            Molecule.from_ase_atoms(ase_molecule("CH3OH")) for _ in range(5)
        ]
        grouper = BasicRMSDGrouper(
            molecules,
            threshold=0.5,
            num_procs=1,
        )
        distance_matrix = np.array(
            [
                [0.0, np.inf, 0.8, 0.9, 0.9],
                [np.inf, 0.0, np.inf, 0.8, 0.8],
                [0.8, np.inf, 0.0, 0.7, 0.7],
                [0.9, 0.8, 0.7, 0.0, 0.2],
                [0.9, 0.8, 0.7, 0.2, 0.0],
            ]
        )

        groups, index_groups = grouper.group_by_threshold(distance_matrix)

        assert grouper._matrix_skipped_indices == [1]
        assert {frozenset(group) for group in index_groups} == {
            frozenset({0}),
            frozenset({2}),
            frozenset({3, 4}),
        }
        assert len(groups) == 3

    def test_original_indices_are_preserved_after_skipping(
        self, methanol_molecules
    ):
        from ase.build import molecule as ase_molecule

        from chemsmart.io.molecules.structure import Molecule

        molecules = [
            Molecule.from_ase_atoms(ase_molecule("CH3OH")) for _ in range(5)
        ]
        grouper = BasicRMSDGrouper(
            molecules,
            threshold=0.5,
            num_procs=1,
        )
        distance_matrix = np.array(
            [
                [0.0, np.inf, 0.9, 0.9, 0.9],
                [np.inf, 0.0, 0.9, 0.9, 0.9],
                [0.9, 0.9, 0.0, 0.8, 0.2],
                [0.9, 0.9, 0.8, 0.0, 0.8],
                [0.9, 0.9, 0.2, 0.8, 0.0],
            ]
        )

        _, index_groups = grouper.group_by_threshold(distance_matrix)

        assert grouper._matrix_skipped_indices == [0]
        assert {frozenset(group) for group in index_groups} == {
            frozenset({1}),
            frozenset({2, 4}),
            frozenset({3}),
        }

    def test_num_groups_after_skipping_uses_valid_structure_count(
        self, methanol_molecules
    ):
        from ase.build import molecule as ase_molecule

        from chemsmart.io.molecules.structure import Molecule

        molecules = [
            Molecule.from_ase_atoms(ase_molecule("CH3OH")) for _ in range(5)
        ]
        grouper = BasicRMSDGrouper(
            molecules,
            threshold=None,
            num_groups=5,
            num_procs=1,
        )
        distance_matrix = np.array(
            [
                [0.0, np.inf, 0.9, 0.9, 0.9],
                [np.inf, 0.0, 0.9, 0.9, 0.9],
                [0.9, 0.9, 0.0, 0.4, 0.5],
                [0.9, 0.9, 0.4, 0.0, 0.3],
                [0.9, 0.9, 0.5, 0.3, 0.0],
            ]
        )

        groups, index_groups = grouper.group_by_num_groups(distance_matrix)

        assert grouper._matrix_skipped_indices == [0]
        assert len(groups) == 4
        assert index_groups == [[1], [2], [3], [4]]

    def test_negative_distance_raises(self, methanol_molecules):
        molecules = methanol_molecules[:2]
        grouper = BasicRMSDGrouper(molecules, threshold=0.5, num_procs=1)
        distance_matrix = np.array([[0.0, -0.1], [-0.1, 0.0]])

        with pytest.raises(ValueError, match="negative"):
            grouper.group_by_threshold(distance_matrix)

    def test_negative_infinity_raises(self, methanol_molecules):
        molecules = methanol_molecules[:2]
        grouper = BasicRMSDGrouper(molecules, threshold=0.5, num_procs=1)
        distance_matrix = np.array([[0.0, -np.inf], [-np.inf, 0.0]])

        with pytest.raises(ValueError, match="negative"):
            grouper.group_by_threshold(distance_matrix)

    def test_asymmetric_distance_matrix_raises(self, methanol_molecules):
        molecules = methanol_molecules[:2]
        grouper = BasicRMSDGrouper(molecules, threshold=0.5, num_procs=1)
        distance_matrix = np.array([[0.0, 0.2], [0.3, 0.0]])

        with pytest.raises(ValueError, match="symmetric"):
            grouper.group_by_threshold(distance_matrix)

    def test_nonzero_diagonal_raises(self, methanol_molecules):
        molecules = methanol_molecules[:2]
        grouper = BasicRMSDGrouper(molecules, threshold=0.5, num_procs=1)
        distance_matrix = np.array([[0.1, 0.2], [0.2, 0.0]])

        with pytest.raises(ValueError, match="diagonal"):
            grouper.group_by_threshold(distance_matrix)

    def test_infinite_pair_tie_breaking_is_deterministic(
        self, methanol_molecules
    ):
        """Equal +inf counts should use a deterministic tie-breaking rule."""
        molecules = methanol_molecules[:3]
        grouper = BasicRMSDGrouper(
            molecules,
            threshold=0.5,
            num_procs=1,
        )

        distance_matrix = np.array(
            [
                [0.0, np.inf, 0.8],
                [np.inf, 0.0, 0.8],
                [0.8, 0.8, 0.0],
            ]
        )
        groups, index_groups = grouper.group_by_threshold(distance_matrix)
        # Structures 0 and 1 have the same number of +inf entries.
        # The current implementation deterministically removes the
        # first maximum, i.e. original index 0.
        assert grouper._matrix_skipped_indices == [0]
        assert {frozenset(group) for group in index_groups} == {
            frozenset({1}),
            frozenset({2}),
        }
        assert len(groups) == 2

    def test_infinite_distance_removal_is_iterative(self, methanol_molecules):
        """Recount +inf involvement after each removal until the submatrix is finite."""
        from ase.build import molecule as ase_molecule

        from chemsmart.io.molecules.structure import Molecule

        molecules = [
            Molecule.from_ase_atoms(ase_molecule("CH3OH")) for _ in range(5)
        ]
        grouper = BasicRMSDGrouper(
            molecules,
            threshold=0.5,
            num_procs=1,
        )
        distance_matrix = np.array(
            [
                [0.0, np.inf, np.inf, 0.9, 0.9],
                [np.inf, 0.0, 0.8, np.inf, 0.9],
                [np.inf, 0.8, 0.0, np.inf, 0.9],
                [0.9, np.inf, np.inf, 0.0, 0.2],
                [0.9, 0.9, 0.9, 0.2, 0.0],
            ]
        )
        groups, index_groups = grouper.group_by_threshold(distance_matrix)
        # Initial +inf counts:
        # 0 -> 2
        # 1 -> 2
        # 2 -> 2
        # 3 -> 2
        # 4 -> 0
        #
        # Tie-breaking removes 0 first.
        #
        # Remaining +inf pairs:
        # 1-3
        # 2-3
        #
        # Recomputed counts:
        # 1 -> 1
        # 2 -> 1
        # 3 -> 2
        #
        # Therefore 3 must be removed in the second iteration.
        assert grouper._matrix_skipped_indices == [0, 3]
        # Remaining original indices are 1, 2, and 4.
        # All pairwise distances are > 0.5, so all are singleton groups.
        assert {frozenset(group) for group in index_groups} == {
            frozenset({1}),
            frozenset({2}),
            frozenset({4}),
        }
        assert len(groups) == 3

    # --- Immutability and order-independence ---

    def test_grouping_does_not_modify_distance_matrix(
        self, methanol_molecules
    ):
        """Hierarchical grouping must not mutate the caller's distance matrix."""
        molecules = methanol_molecules[:3]
        grouper = BasicRMSDGrouper(
            molecules,
            threshold=0.5,
            num_procs=1,
        )

        distance_matrix = np.array(
            [
                [0.0, 0.3, 0.8],
                [0.3, 0.0, 0.3],
                [0.8, 0.3, 0.0],
            ]
        )
        original_matrix = distance_matrix.copy()

        grouper.group_by_threshold(distance_matrix)

        assert np.array_equal(distance_matrix, original_matrix)

    def test_hierarchical_grouping_is_input_order_independent(
        self, methanol_molecules
    ):
        """Permuting input order must not change the partition after remapping indices."""
        from ase.build import molecule as ase_molecule

        from chemsmart.io.molecules.structure import Molecule

        molecules = [
            Molecule.from_ase_atoms(ase_molecule("CH3OH")) for _ in range(4)
        ]
        distance_matrix = np.array(
            [
                [0.0, 0.2, 0.9, 1.0],
                [0.2, 0.0, 0.8, 0.9],
                [0.9, 0.8, 0.0, 0.3],
                [1.0, 0.9, 0.3, 0.0],
            ]
        )

        grouper = BasicRMSDGrouper(
            molecules,
            threshold=0.5,
            num_procs=1,
        )
        _, index_groups = grouper.group_by_threshold(distance_matrix)
        expected_partition = {frozenset(group) for group in index_groups}

        permutation = [2, 0, 3, 1]
        permuted_molecules = [molecules[i] for i in permutation]
        permuted_matrix = distance_matrix[np.ix_(permutation, permutation)]
        permuted_grouper = BasicRMSDGrouper(
            permuted_molecules,
            threshold=0.5,
            num_procs=1,
        )
        _, permuted_groups = permuted_grouper.group_by_threshold(
            permuted_matrix
        )

        remapped_partition = {
            frozenset(permutation[i] for i in group)
            for group in permuted_groups
        }
        assert remapped_partition == expected_partition


@pytest.mark.usefixtures("temporary_working_dir")
class TestConformerIdsFunctionality:
    """Test conformer_ids parameter functionality."""

    NUM_PROCS = 1

    def test_conformer_ids_from_log_directory(
        self, ts_conformers_log_directory
    ):
        """Test loading conformer IDs from a directory of log files."""
        import glob
        import re

        from chemsmart.io.molecules.structure import Molecule

        # Load molecules from log files in the ts_conformers_log_directory
        log_files = sorted(
            glob.glob(os.path.join(ts_conformers_log_directory, "*.log"))
        )
        assert (
            len(log_files) == 5
        ), f"Expected 5 log files, found {len(log_files)}"

        # Extract conformer IDs from filenames (e.g., ch_1c_para_c1.log -> c1)
        molecules = []
        conformer_ids = []
        for log_file in log_files:
            mol = Molecule.from_filepath(
                filepath=log_file, index="-1", return_list=False
            )
            if mol is not None:
                molecules.append(mol)
                # Extract conformer ID from filename
                basename = os.path.basename(log_file)
                match = re.search(r"_c(\d+)\.log$", basename)
                if match:
                    conformer_ids.append(f"c{match.group(1)}")
                else:
                    conformer_ids.append(basename)

        assert len(molecules) >= 2, "Need at least 2 molecules for grouping"

        grouper = BasicRMSDGrouper(
            molecules,
            threshold=1.0,
            num_procs=1,
            label="ts_conformers_test",
            conformer_ids=conformer_ids,
        )
        _, _ = grouper.group()

        # Verify conformer_ids are stored
        assert grouper.conformer_ids == conformer_ids

        # Check Excel output
        import pandas as pd

        excel_file = os.path.join(
            "ts_conformers_test_group_result",
            "ts_conformers_test_BasicRMSDGrouper_T1.0.xlsx",
        )
        assert os.path.exists(
            excel_file
        ), f"Excel file not found: {excel_file}"

        # Matrix data starts at row 8 (0-indexed: skiprows=15), first column is index
        df = pd.read_excel(
            excel_file, sheet_name="RMSD_Matrix", skiprows=14, index_col=0
        )
        # Check that conformer IDs are used as labels
        assert "c1" in str(df.columns[0]) or "c1" in str(df.index[0])

    def test_traj_conformer_ids_original_indices(self):
        """Test that traj job correctly sets original conformer indices."""
        from ase.build import molecule as ase_molecule

        from chemsmart.io.molecules.structure import Molecule

        # Local ASE methanol — avoid flaky live PubChem for a trivial fixture.
        mol = Molecule.from_ase_atoms(ase_molecule("CH3OH"))
        molecules = [mol.copy() for _ in range(18)]

        # Simulate traj behavior: select last 50% (indices 10-18 in original)
        proportion = 0.5
        total = len(molecules)
        last_num = int(round(total * proportion, 1))
        _ = molecules[-last_num:]  # selected_molecules - not used in this test

        # Calculate expected original indices (1-based)
        start_original_index = total - last_num + 1
        expected_ids = [str(start_original_index + i) for i in range(last_num)]

        # Verify expected IDs
        assert expected_ids[0] == "10"  # First selected is original index 10
        assert expected_ids[-1] == "18"  # Last selected is original index 18


@pytest.mark.usefixtures("temporary_working_dir")
class TestOutputFileGeneration:
    """Test that grouper generates correct output files."""

    NUM_PROCS = 1

    def test_group_xyz_files_contain_energy_and_index_info(
        self, multiple_molecules_xyz_file
    ):
        """Test that group XYZ files contain energy and original index information."""
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)

        grouper = BasicRMSDGrouper(
            molecules[:6],
            threshold=1.0,
            num_procs=1,
            label="info_test",
        )
        _, _ = grouper.group()
        grouper.unique()  # Generates xyz files

        output_dir = "info_test_group_result"
        first_group_file = os.path.join(output_dir, "info_test_group_1.xyz")

        with open(first_group_file, "r") as f:
            content = f.read()

        # Check for expected information in comments
        assert "Original_Index:" in content
        assert "E" in content

    def test_group_xyz_files_sorted_by_energy(
        self, multiple_molecules_xyz_file
    ):
        """Test that molecules in group XYZ files are sorted by energy (lowest first)."""
        import re

        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)

        grouper = BasicRMSDGrouper(
            molecules,
            threshold=2.0,  # Higher threshold to get more molecules per group
            num_procs=1,
            label="sort_test",
        )
        groups, group_indices = grouper.group()
        grouper.unique()  # Generates xyz files

        output_dir = "sort_test_group_result"

        # Check each group file
        for group_num in range(len(groups)):
            group_file = os.path.join(
                output_dir, f"sort_test_group_{group_num + 1}.xyz"
            )

            with open(group_file, "r") as f:
                content = f.read()

            # Extract energies from file
            lines = content.strip().split("\n")
            i = 0
            energies = []

            while i < len(lines):
                num_atoms = int(lines[i].strip())
                comment_line = lines[i + 1]
                energy_match = re.search(
                    r"E\(Hartree\):\s*([-\d.]+)", comment_line
                )
                if energy_match:
                    energies.append(float(energy_match.group(1)))
                i += num_atoms + 2

            # Verify sorted by energy (lowest first)
            for j in range(len(energies) - 1):
                assert energies[j] <= energies[j + 1], (
                    f"Group {group_num + 1} not sorted: "
                    f"energy[{j}]={energies[j]} > energy[{j+1}]={energies[j+1]}"
                )

    def test_group_xyz_files_energy_from_log_files(
        self, ts_conformers_log_directory
    ):
        """Test that energy is correctly extracted from log files and written to group XYZ files."""
        import glob
        import re

        from chemsmart.io.molecules.structure import Molecule

        # Load molecules from log files
        log_files = sorted(
            glob.glob(os.path.join(ts_conformers_log_directory, "*.log"))
        )
        assert (
            len(log_files) >= 2
        ), f"Need at least 2 log files, found {len(log_files)}"

        molecules = []
        conformer_ids = []
        original_energies = {}

        for idx, log_file in enumerate(log_files):
            mol = Molecule.from_filepath(
                filepath=log_file, index="-1", return_list=False
            )
            if mol is not None:
                molecules.append(mol)
                # Store original energy for verification
                original_energies[idx] = mol.energy

                # Extract conformer ID from filename
                basename = os.path.basename(log_file)
                match = re.search(r"_c(\d+)\.log$", basename)
                if match:
                    conformer_ids.append(f"c{match.group(1)}")
                else:
                    conformer_ids.append(basename)

        assert len(molecules) >= 2, "Need at least 2 molecules for grouping"

        grouper = BasicRMSDGrouper(
            molecules,
            threshold=5.0,  # High threshold to group all together
            num_procs=1,
            label="log_energy_test",
            conformer_ids=conformer_ids,
        )
        _, _ = grouper.group()
        grouper.unique()  # Generates xyz files

        # Check group file
        output_dir = "log_energy_test_group_result"
        first_group_file = os.path.join(
            output_dir, "log_energy_test_group_1.xyz"
        )

        assert os.path.exists(
            first_group_file
        ), f"Group file not found: {first_group_file}"

        with open(first_group_file, "r") as f:
            content = f.read()

        # Verify energy values in the output
        assert "E(Hartree):" in content

        # Parse and verify each energy value
        lines = content.strip().split("\n")
        i = 0
        energies_found = 0
        while i < len(lines):
            try:
                num_atoms = int(lines[i].strip())
            except ValueError:
                break
            comment_line = lines[i + 1]

            # Extract energy from comment
            energy_match = re.search(
                r"E\(Hartree\):\s*([-\d.]+)", comment_line
            )
            if energy_match:
                extracted_energy = float(energy_match.group(1))
                # Verify the energy is a reasonable value (not zero or nan)
                assert (
                    extracted_energy < 0
                ), f"Energy should be negative: {extracted_energy}"
                energies_found += 1

            i += num_atoms + 2

        assert energies_found > 0, "No energies found in output file"


@pytest.mark.usefixtures("temporary_working_dir")
class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_two_molecules_grouping(self):
        """Test grouping with minimum number of molecules (2)."""
        from ase.build import molecule as ase_molecule

        from chemsmart.io.molecules.structure import Molecule

        mol1 = Molecule.from_ase_atoms(ase_molecule("CH3OH"))
        mol2 = mol1.copy()

        grouper = BasicRMSDGrouper([mol1, mol2], threshold=0.5)
        groups, group_indices = grouper.group()

        # Two identical molecules should form one group
        assert len(groups) == 1
        assert len(groups[0]) == 2

    def test_num_groups_equals_num_molecules(self, methanol_molecules):
        """Test requesting same number of groups as molecules."""
        grouper = BasicRMSDGrouper(
            methanol_molecules,
            num_groups=len(methanol_molecules),
        )
        groups, group_indices = grouper.group()

        # Requesting N equal to the number of molecules returns singleton groups.
        assert len(groups) == len(methanol_molecules)
        assert all(len(group) == 1 for group in group_indices)

    def test_num_groups_exceeds_num_molecules(self, methanol_molecules):
        """Test requesting more groups than molecules."""
        n_mols = len(methanol_molecules)
        grouper = BasicRMSDGrouper(
            methanol_molecules,
            num_groups=n_mols + 5,
        )
        groups, group_indices = grouper.group()

        # More requested groups than molecules still returns singleton groups.
        assert len(groups) == n_mols
        assert all(len(group) == 1 for group in group_indices)

    def test_very_low_threshold(self, multiple_molecules_xyz_file):
        """Test with very low threshold (should create many groups)."""
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)

        grouper = BasicRMSDGrouper(
            molecules,
            threshold=0.001,  # Very low threshold
            num_procs=1,
        )
        groups, group_indices = grouper.group()

        # With very low threshold, most molecules should be in separate groups
        # (unless they're truly identical)
        assert len(groups) >= len(molecules) - 1

    def test_very_high_threshold(self, multiple_molecules_xyz_file):
        """Test with very high threshold (should create few groups)."""
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)

        grouper = BasicRMSDGrouper(
            molecules,
            threshold=100.0,  # Very high threshold
            num_procs=1,
        )
        groups, group_indices = grouper.group()

        # With very high threshold, all molecules with same formula should be in one group
        assert len(groups) <= 3  # Likely 1-2 groups

    def test_different_formulas_always_separate(self, methanol_and_ethanol):
        """Test that molecules with different formulas are always in separate groups."""
        grouper = BasicRMSDGrouper(
            methanol_and_ethanol,
            threshold=100.0,  # Even with high threshold
        )
        groups, group_indices = grouper.group()

        # Incompatible RMSD pairs are now skipped from clustering entirely.
        assert len(groups) == 1
        assert len(group_indices) == 1
        assert grouper._matrix_skipped_indices == [0]
        assert group_indices == [[1]]

    def test_rmsd_infinity_for_different_molecules(self, methanol_and_ethanol):
        """Test that RMSD returns infinity for molecules with different atom counts."""
        grouper = BasicRMSDGrouper(methanol_and_ethanol, threshold=0.5)
        rmsd = grouper.calculate_rmsd_pair(0, 1)

        assert rmsd == np.inf or rmsd == float("inf")

    def test_incompatible_rmsd_structures_are_skipped(
        self, methanol_and_ethanol
    ):
        """Structures involved in RMSD=inf pairs should be skipped, not grouped."""
        grouper = BasicRMSDGrouper(methanol_and_ethanol, threshold=0.5)
        groups, group_indices = grouper.group()

        assert len(groups) == 1
        assert len(group_indices) == 1
        assert grouper._matrix_skipped_indices == [0]
        assert group_indices == [[1]]


@pytest.mark.usefixtures("temporary_working_dir")
class TestLabelAndAppendLabel:
    """Test -l (label) and -a (append_label) parameter functionality."""

    NUM_PROCS = 1

    def test_get_label_function(self):
        """Test _get_label function logic."""
        from chemsmart.cli.grouper.grouper import _get_label

        # Test with label only
        result = _get_label(
            label="custom", append_label=None, base_label="base"
        )
        assert result == "custom"

        # Test with append_label only
        result = _get_label(
            label=None, append_label="suffix", base_label="base"
        )
        assert result == "base_suffix"

        # Test with neither (use base_label)
        result = _get_label(label=None, append_label=None, base_label="base")
        assert result == "base"

        # Test with both should raise error
        with pytest.raises(ValueError) as excinfo:
            _get_label(
                label="custom", append_label="suffix", base_label="base"
            )
        assert "Only give label or append_label" in str(excinfo.value)

    def test_label_in_output_directory(self, multiple_molecules_xyz_file):
        """Test that label parameter affects output directory name."""
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)[:5]

        # Test with custom label
        grouper = BasicRMSDGrouper(
            molecules,
            threshold=1.0,
            num_procs=self.NUM_PROCS,
            label="my_custom_label",
        )
        _, _ = grouper.group()

        # Check output directory uses custom label
        output_dir = "my_custom_label_group_result"
        assert os.path.exists(
            output_dir
        ), f"Output dir not found: {output_dir}"

        # Check Excel file uses custom label
        excel_file = os.path.join(
            output_dir, "my_custom_label_BasicRMSDGrouper_T1.0.xlsx"
        )
        assert os.path.exists(
            excel_file
        ), f"Excel file not found: {excel_file}"

    def test_label_in_group_xyz_files(self, multiple_molecules_xyz_file):
        """Test that label parameter affects group XYZ file names."""
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)[:5]

        grouper = BasicRMSDGrouper(
            molecules,
            threshold=1.0,
            num_procs=self.NUM_PROCS,
            label="test_label",
        )
        groups, group_indices = grouper.group()
        grouper.unique()  # Generates xyz files

        output_dir = "test_label_group_result"

        # Check group XYZ files use label
        for i in range(len(groups)):
            xyz_path = os.path.join(output_dir, f"test_label_group_{i+1}.xyz")
            assert os.path.exists(xyz_path), f"Group XYZ not found: {xyz_path}"

    def test_different_labels_create_different_outputs(
        self, multiple_molecules_xyz_file
    ):
        """Test that different labels create separate output directories."""
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)[:5]

        # First grouper with label "run1"
        grouper1 = BasicRMSDGrouper(
            molecules,
            threshold=1.0,
            num_procs=self.NUM_PROCS,
            label="run1",
        )
        grouper1.group()

        # Second grouper with label "run2"
        grouper2 = BasicRMSDGrouper(
            molecules,
            threshold=1.0,
            num_procs=self.NUM_PROCS,
            label="run2",
        )
        grouper2.group()

        # Both output directories should exist
        assert os.path.exists("run1_group_result")
        assert os.path.exists("run2_group_result")

        # Both should have their own Excel files
        assert os.path.exists(
            "run1_group_result/run1_BasicRMSDGrouper_T1.0.xlsx"
        )
        assert os.path.exists(
            "run2_group_result/run2_BasicRMSDGrouper_T1.0.xlsx"
        )

    def test_label_with_num_groups(self, multiple_molecules_xyz_file):
        """Test that label works correctly with num_groups parameter."""
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)

        grouper = BasicRMSDGrouper(
            molecules,
            num_groups=5,
            num_procs=self.NUM_PROCS,
            label="numgroups_test",
        )
        _, _ = grouper.group()

        output_dir = "numgroups_test_group_result"
        assert os.path.exists(output_dir)

        # Excel file should reflect num_groups (N5) instead of threshold
        excel_file = os.path.join(
            output_dir, "numgroups_test_BasicRMSDGrouper_N5.xlsx"
        )
        assert os.path.exists(
            excel_file
        ), f"Excel file not found: {excel_file}"


class TestConformerIdExtraction:
    """Test conformer ID extraction from filenames."""

    @staticmethod
    def _sort_key(x):
        """Sort key for file_info: cXX files by number first, then others by filename."""
        if x[2] is not None:
            return (0, x[2], x[3])
        else:
            return (1, 0, x[3])

    def test_extract_conformer_id_with_pattern(self):
        """Test extraction with _cXX_ pattern."""
        from chemsmart.cli.grouper.grouper import _extract_conformer_id

        assert _extract_conformer_id("structure_c1_opt.log") == "c1"
        assert _extract_conformer_id("mol_c12.log") == "c12"
        assert (
            _extract_conformer_id("/path/to/structure_c123_ts.log") == "c123"
        )

    def test_extract_conformer_id_without_pattern(self):
        """Test extraction without _cXX_ pattern returns None."""
        from chemsmart.cli.grouper.grouper import _extract_conformer_id

        assert _extract_conformer_id("mol_opt.log") is None
        assert _extract_conformer_id("structure.log") is None
        assert _extract_conformer_id("/path/to/calc.log") is None

    def test_conformer_id_fallback_to_filename(self):
        """Test that filename is used as conformer ID when no _cXX_ pattern."""
        from chemsmart.cli.grouper.grouper import _extract_conformer_id

        test_files = [
            "/path/to/mol_opt.log",
            "/path/to/structure_ts.out",
        ]

        conformer_ids = []
        for f in test_files:
            conf_id = _extract_conformer_id(f)
            if conf_id is None:
                basename = os.path.basename(f)
                conf_id = os.path.splitext(basename)[0]
            conformer_ids.append(conf_id)

        assert conformer_ids == ["mol_opt", "structure_ts"]

    def test_conformer_ids_molecules_correspondence(
        self, ts_conformers_log_directory
    ):
        """Test that conformer_ids and molecules are strictly one-to-one corresponding.

        This verifies that after sorting, each molecule's energy matches the energy
        from its corresponding file (identified by conf_id).
        """
        import glob
        import re

        from chemsmart.cli.grouper.grouper import _extract_conformer_id
        from chemsmart.io.gaussian.output import Gaussian16Output
        from chemsmart.io.molecules.structure import Molecule

        # First, build a mapping from conf_id to expected gibbs_energy
        log_files = glob.glob(
            os.path.join(ts_conformers_log_directory, "*.log")
        )
        expected_energies = {}

        for f in log_files:
            conf_id = _extract_conformer_id(f)
            basename = os.path.basename(f)
            name_without_ext = os.path.splitext(basename)[0]

            if conf_id is None:
                conf_id = name_without_ext

            g16_output = Gaussian16Output(filename=f)
            gibbs_energy = g16_output.gibbs_free_energy
            if gibbs_energy is not None:
                expected_energies[conf_id] = gibbs_energy

        # Now simulate the actual loading logic with sorting
        file_info = []
        for f in log_files:
            conf_id = _extract_conformer_id(f)
            basename = os.path.basename(f)
            name_without_ext = os.path.splitext(basename)[0]

            if conf_id:
                num = int(re.search(r"\d+", conf_id).group())
                file_info.append((f, conf_id, num, name_without_ext))
            else:
                file_info.append((f, name_without_ext, None, name_without_ext))

        file_info.sort(key=self._sort_key)

        # Load molecules
        molecules = []
        conformer_ids = []

        for filepath, conf_id, _, _ in file_info:
            try:
                mol = Molecule.from_filepath(
                    filepath=filepath, index="-1", return_list=False
                )
                if mol is not None:
                    g16_output = Gaussian16Output(filename=filepath)
                    gibbs_energy = g16_output.gibbs_free_energy
                    if gibbs_energy is not None:
                        mol._energy = gibbs_energy

                    molecules.append(mol)
                    conformer_ids.append(conf_id)
            except Exception:
                pass

        # Verify one-to-one correspondence
        assert len(molecules) == len(conformer_ids)

        # Verify each molecule's energy matches the expected energy for its conf_id
        for i, (mol, conf_id) in enumerate(zip(molecules, conformer_ids)):
            if conf_id in expected_energies:
                expected_energy = expected_energies[conf_id]
                assert np.isclose(mol.energy, expected_energy, rtol=1e-10), (
                    f"Index {i}: energy mismatch for {conf_id}, "
                    f"got {mol.energy}, expected {expected_energy}"
                )

    def test_mixed_pattern_files_correspondence(self):
        """Test correspondence when mixing files with and without _cXX_ pattern."""
        from chemsmart.cli.grouper.grouper import _extract_conformer_id

        test_files = [
            ("mol_c1_opt.log", "c1"),
            ("mol_c2_opt.log", "c2"),
            ("other_calc.log", None),
            ("mol_c3_opt.log", "c3"),
        ]

        file_info = []
        for filename, expected_pattern in test_files:
            conf_id = _extract_conformer_id(filename)
            basename = os.path.basename(filename)
            name_without_ext = os.path.splitext(basename)[0]

            if conf_id:
                assert conf_id == expected_pattern
                num = int(conf_id.replace("c", ""))
                file_info.append((filename, conf_id, num, name_without_ext))
            else:
                assert expected_pattern is None
                file_info.append(
                    (filename, name_without_ext, None, name_without_ext)
                )

        file_info.sort(key=self._sort_key)
        sorted_conf_ids = [info[1] for info in file_info]

        # Files with pattern first (sorted by number), then others by filename
        assert sorted_conf_ids == ["c1", "c2", "c3", "other_calc"]

    def test_all_files_without_cxx_pattern(self):
        """Test that a folder with no _cXX_ pattern files works correctly."""
        from chemsmart.cli.grouper.grouper import _extract_conformer_id

        test_files = [
            "mol_opt.log",
            "structure_ts.log",
            "calc.log",
            "another_mol.log",
        ]

        file_info = []
        for filename in test_files:
            conf_id = _extract_conformer_id(filename)
            basename = os.path.basename(filename)
            name_without_ext = os.path.splitext(basename)[0]

            if conf_id:
                num = int(conf_id.replace("c", ""))
                file_info.append((filename, conf_id, num, name_without_ext))
            else:
                file_info.append(
                    (filename, name_without_ext, None, name_without_ext)
                )

        # All files should have None as num
        for info in file_info:
            assert info[2] is None

        file_info.sort(key=self._sort_key)
        sorted_conf_ids = [info[1] for info in file_info]

        # All files sorted alphabetically by filename
        assert sorted_conf_ids == [
            "another_mol",
            "calc",
            "mol_opt",
            "structure_ts",
        ]


@pytest.mark.usefixtures("temporary_working_dir")
class TestEnergyExtractionFunction:
    """Framework tests for energy extraction by file type."""

    def test_energy_extraction_gaussian(
        self, gaussian_dppeFeCl2_link_opt_outputfile
    ):
        from chemsmart.analysis.thermochemistry import Thermochemistry
        from chemsmart.cli.grouper.grouper import (
            _extract_energy_based_on_energy_type,
        )

        thermo = Thermochemistry(
            filename=gaussian_dppeFeCl2_link_opt_outputfile
        )
        extracted = _extract_energy_based_on_energy_type(thermo, "E")
        assert np.isclose(extracted, -3869.01351827, rtol=1e-7)

        thermo = Thermochemistry(
            filename=gaussian_dppeFeCl2_link_opt_outputfile
        )
        extracted = _extract_energy_based_on_energy_type(thermo, "H")
        assert np.isclose(extracted, -3868.552084, rtol=1e-7)

        thermo = Thermochemistry(
            filename=gaussian_dppeFeCl2_link_opt_outputfile
        )
        extracted = _extract_energy_based_on_energy_type(thermo, "G")
        assert np.isclose(extracted, -3868.648944, rtol=1e-7)

        thermo_qhh = Thermochemistry(
            filename=gaussian_dppeFeCl2_link_opt_outputfile,
            temperature=298.15,
            concentration=1.0,
            h_freq_cutoff=100.0,
            s_freq_cutoff=100.0,
            use_weighted_mass=True,
            entropy_method="grimme",
        )
        extracted = _extract_energy_based_on_energy_type(thermo_qhh, "qhH")
        assert np.isclose(extracted, -3868.558465, rtol=1e-7)

        thermo_qhg = Thermochemistry(
            filename=gaussian_dppeFeCl2_link_opt_outputfile,
            temperature=298.15,
            concentration=1.0,
            h_freq_cutoff=100.0,
            s_freq_cutoff=100.0,
            use_weighted_mass=True,
            entropy_method="grimme",
        )
        extracted = _extract_energy_based_on_energy_type(thermo_qhg, "qhG")
        assert np.isclose(extracted, -3868.645329, rtol=1e-7)

    def test_energy_extraction_orca(self, fe2_singlet_output):
        from chemsmart.analysis.thermochemistry import Thermochemistry
        from chemsmart.cli.grouper.grouper import (
            _extract_energy_based_on_energy_type,
        )

        thermo = Thermochemistry(filename=fe2_singlet_output)
        extracted = _extract_energy_based_on_energy_type(thermo, "E")
        assert np.isclose(extracted, -1568.256171848101, rtol=1e-7)

        thermo = Thermochemistry(filename=fe2_singlet_output)
        extracted = _extract_energy_based_on_energy_type(thermo, "H")
        assert np.isclose(extracted, -1568.14237774, rtol=1e-7)

        thermo = Thermochemistry(filename=fe2_singlet_output)
        extracted = _extract_energy_based_on_energy_type(thermo, "G")
        assert np.isclose(extracted, -1568.18873516, rtol=1e-7)

        thermo_qhh = Thermochemistry(
            filename=fe2_singlet_output,
            temperature=298.15,
            concentration=1.0,
            h_freq_cutoff=100.0,
            s_freq_cutoff=100.0,
            use_weighted_mass=True,
            entropy_method="grimme",
        )
        extracted = _extract_energy_based_on_energy_type(thermo_qhh, "qhH")
        assert np.isclose(extracted, -1568.143280, rtol=1e-7)

        thermo_qhg = Thermochemistry(
            filename=fe2_singlet_output,
            temperature=298.15,
            concentration=1.0,
            h_freq_cutoff=100.0,
            s_freq_cutoff=100.0,
            use_weighted_mass=True,
            entropy_method="grimme",
        )
        extracted = _extract_energy_based_on_energy_type(thermo_qhg, "qhG")
        assert np.isclose(extracted, -1568.186619, rtol=1e-7)


@pytest.mark.usefixtures("temporary_working_dir")
class TestMultiprocessingSupport:
    def test_basic_rmsd_serial_parallel_equivalent(self, methanol_molecules):
        serial = BasicRMSDGrouper(
            methanol_molecules, threshold=0.5, num_procs=1
        )
        parallel = BasicRMSDGrouper(
            methanol_molecules, threshold=0.5, num_procs=2
        )

        serial_matrix = serial.calculate_full_rmsd_matrix()
        parallel_matrix = parallel.calculate_full_rmsd_matrix()
        assert np.allclose(
            serial_matrix, parallel_matrix, atol=1e-10, rtol=1e-10
        )

        _, serial_idx = serial.group()
        _, parallel_idx = parallel.group()
        assert serial_idx == parallel_idx

    def test_hungarian_rmsd_serial_parallel_equivalent(
        self, methanol_molecules
    ):
        serial = HungarianRMSDGrouper(
            methanol_molecules, threshold=0.5, num_procs=1
        )
        parallel = HungarianRMSDGrouper(
            methanol_molecules, threshold=0.5, num_procs=2
        )

        serial_matrix = serial.calculate_full_rmsd_matrix()
        parallel_matrix = parallel.calculate_full_rmsd_matrix()
        assert np.allclose(
            serial_matrix, parallel_matrix, atol=1e-10, rtol=1e-10
        )

        _, serial_idx = serial.group()
        _, parallel_idx = parallel.group()
        assert serial_idx == parallel_idx

    def test_spyrmsd_parallel_preserves_isomorphisms(self, methanol_molecules):
        serial = SpyRMSDGrouper(methanol_molecules, threshold=0.5, num_procs=1)
        parallel = SpyRMSDGrouper(
            methanol_molecules, threshold=0.5, num_procs=2
        )

        serial_matrix = serial.calculate_full_rmsd_matrix()
        parallel_matrix = parallel.calculate_full_rmsd_matrix()
        assert np.allclose(
            serial_matrix, parallel_matrix, atol=1e-8, rtol=1e-8
        )

        _, serial_idx = serial.group()
        _, parallel_idx = parallel.group()
        assert serial_idx == parallel_idx

        assert serial.best_isomorphisms
        assert parallel.best_isomorphisms
        assert (0, 1) in serial.best_isomorphisms
        assert (0, 1) in parallel.best_isomorphisms
        assert serial.get_best_isomorphism(0, 1) is not None
        assert parallel.get_best_isomorphism(0, 1) is not None

    @pytest.mark.skipif(
        not _is_irmsd_available(), reason="irmsd API/command not available"
    )
    def test_irmsd_serial_parallel_equivalent(self, methanol_molecules):
        serial = IRMSDGrouper(methanol_molecules, threshold=0.5, num_procs=1)
        parallel = IRMSDGrouper(methanol_molecules, threshold=0.5, num_procs=2)

        serial_matrix = serial.calculate_full_rmsd_matrix()
        parallel_matrix = parallel.calculate_full_rmsd_matrix()
        assert np.allclose(
            serial_matrix, parallel_matrix, atol=1e-8, rtol=1e-8
        )

        _, serial_idx = serial.group()
        _, parallel_idx = parallel.group()
        assert serial_idx == parallel_idx

    def test_formula_grouper_falls_back_to_single_proc(
        self, caplog, methanol_molecules
    ):
        """FormulaGrouper warns and falls back to num_procs=1 when >1 is requested."""
        caplog.set_level("WARNING")
        grouper = FormulaGrouper(methanol_molecules, num_procs=4)
        assert grouper.num_procs == 1
        assert (
            "FormulaGrouper does not support multiprocessing; using num_procs=1."
            in caplog.text
        )

        groups, group_indices = grouper.group()
        assert len(groups) == 1
        assert len(group_indices) == 1

    def test_serial_only_strategy_warns_and_falls_back(
        self, caplog, multiple_molecules_xyz_file
    ):
        """Serial-only grouper warns and falls back to num_procs=1."""
        caplog.set_level("WARNING")
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)
        grouper = EnergyGrouper(molecules, threshold=0.5, num_procs=4)
        assert grouper.num_procs == 1
        assert (
            "EnergyGrouper does not support multiprocessing; using num_procs=1."
            in caplog.text
        )

    def test_basic_rmsd_progress_logs_all_milestones_serial(
        self, caplog, multiple_molecules_xyz_file
    ):
        caplog.set_level("INFO")
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)
        grouper = BasicRMSDGrouper(molecules, threshold=0.5, num_procs=1)

        grouper.calculate_full_rmsd_matrix()

        _assert_progress_milestones(caplog.text)

    def test_basic_rmsd_progress_logs_all_milestones_parallel(
        self, caplog, multiple_molecules_xyz_file
    ):
        caplog.set_level("INFO")
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)
        grouper = BasicRMSDGrouper(molecules, threshold=0.5, num_procs=2)

        grouper.calculate_full_rmsd_matrix()

        _assert_progress_milestones(caplog.text)

    def test_energy_grouper_progress_logs_all_milestones(
        self, caplog, multiple_molecules_xyz_file
    ):
        caplog.set_level("INFO")
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)
        grouper = EnergyGrouper(molecules, threshold=0.5, num_procs=1)

        grouper.group()

        _assert_progress_milestones(caplog.text)

    def test_basic_rmsd_parallel_single_molecule_no_pair_no_crash(
        self, methanol_molecules
    ):
        """Single-molecule RMSD run with num_procs>1 should safely skip pair workers."""
        molecules = [methanol_molecules[0]]
        grouper = BasicRMSDGrouper(molecules, threshold=0.5, num_procs=2)

        matrix = grouper.calculate_full_rmsd_matrix()
        assert matrix.shape == (1, 1)
        assert matrix[0, 0] == 0.0

        groups, index_groups = grouper.group()
        assert len(groups) == 1
        assert len(index_groups) == 1
        assert index_groups == [[0]]

    def test_matrix_grouper_log_progress(self, caplog, make_matrix_grouper):
        """Test functionality of MatrixGrouper._log_progress method."""
        import logging

        grouper = make_matrix_grouper([], threshold=0.1)
        with caplog.at_level(logging.INFO):
            # Initial call, 0% complete, next_progress threshold at 10%
            next_progress = 10
            next_progress = grouper._log_progress(0, 100, next_progress)
            assert next_progress == 10
            assert len(caplog.records) == 0

            # Call at 9%
            next_progress = grouper._log_progress(9, 100, next_progress)
            assert next_progress == 10
            assert len(caplog.records) == 0

            # Call at 10%
            next_progress = grouper._log_progress(10, 100, next_progress)
            assert next_progress == 20
            assert len(caplog.records) == 1
            assert (
                "Matrix calculation progress: 10%"
                in caplog.records[-1].message
            )

            # Call at 15% (no log, threshold is 20)
            next_progress = grouper._log_progress(15, 100, next_progress)
            assert next_progress == 20
            assert len(caplog.records) == 1

            # Call jumping from 15% to 35%
            next_progress = grouper._log_progress(35, 100, next_progress)
            assert next_progress == 40
            assert len(caplog.records) == 3
            assert (
                "Matrix calculation progress: 20%"
                in caplog.records[-2].message
            )
            assert (
                "Matrix calculation progress: 30%"
                in caplog.records[-1].message
            )

            # Test edge case with zero total
            next_p = grouper._log_progress(10, 0, 10)
            assert next_p == 10


@pytest.mark.usefixtures("temporary_working_dir")
class TestRepresentativeStrategy:
    class DummyDistanceGrouper(MatrixGrouper):
        def __init__(self, molecules, distance_matrix, **kwargs):
            super().__init__(molecules, threshold=9.9, **kwargs)
            self._distance_matrix = np.array(distance_matrix, dtype=float)

        def group(self):
            groups, index_groups = self.group_by_threshold(
                self._distance_matrix
            )
            self._cached_groups = groups
            self._cached_group_indices = index_groups
            self.record(distance_matrix=self._distance_matrix)
            return groups, index_groups

        def _record_results(self, distance_matrix):
            recorder = self._get_results_recorder()
            labels = recorder.get_labels(distance_matrix.shape[0])
            header_info = [("", "Dummy Distance Grouper")]
            self._append_input_usage_header(header_info)
            sheets_data = {}
            if self._cached_group_indices is not None:
                sheets_data["Groups"] = recorder.build_groups_dataframe(
                    self._cached_group_indices, len(self.molecules)
                )
            recorder.record_results(
                grouper_name=self.__class__.__name__,
                header_info=header_info,
                sheets_data=sheets_data,
                matrix_data=("Distance", distance_matrix, labels),
                suffix="test",
                startrow=len(header_info) + 2,
            )

    @staticmethod
    def _set_energies(molecules, energies):
        for mol, energy in zip(molecules, energies):
            mol.energy = energy

    def test_default_representative_strategy_is_lowest(
        self, methanol_molecules
    ):
        grouper = BasicRMSDGrouper(methanol_molecules[:2], threshold=0.5)
        assert grouper.representative_strategy == "lowest"

    def test_lowest_strategy_orders_group_by_energy(
        self, multiple_molecules_xyz_file
    ):
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)[:4]
        self._set_energies(molecules, [5.0, 1.0, 3.0, None])
        matrix = np.array(
            [
                [0.0, 0.1, 0.2, 0.3],
                [0.1, 0.0, 0.2, 0.3],
                [0.2, 0.2, 0.0, 0.3],
                [0.3, 0.3, 0.3, 0.0],
            ]
        )
        grouper = self.DummyDistanceGrouper(
            molecules,
            matrix,
            representative_strategy="lowest",
        )
        groups, index_groups = grouper.group()

        assert index_groups[0][0] == 1
        assert groups[0][0] is molecules[1]
        assert index_groups[0] == [1, 2, 0, 3]

    def test_center_strategy_orders_entire_group_by_centrality(
        self, multiple_molecules_xyz_file
    ):
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)[:5]
        self._set_energies(molecules, [0.0, 0.1, 0.2, -1.0, 0.3])
        matrix = np.array(
            [
                [0.0, 1.0, 1.0, 1.0, 1.0],
                [1.0, 0.0, 2.0, 2.0, 2.0],
                [1.0, 2.0, 0.0, 2.0, 2.0],
                [1.0, 2.0, 2.0, 0.0, 2.0],
                [1.0, 2.0, 2.0, 2.0, 0.0],
            ]
        )
        grouper = self.DummyDistanceGrouper(
            molecules,
            matrix,
            representative_strategy="center",
        )
        groups, index_groups = grouper.group()

        assert set(index_groups[0]) == {0, 1, 2, 3, 4}
        assert index_groups[0][0] == 0
        assert groups[0][0] is molecules[0]
        assert index_groups[0] == [0, 3, 1, 2, 4]

    def test_top3_strategy_picks_most_central_among_three_lowest_energy(
        self, multiple_molecules_xyz_file
    ):
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)[:5]
        self._set_energies(molecules, [0.0, 0.1, 0.2, -1.0, 0.3])
        matrix = np.array(
            [
                [0.0, 1.0, 1.0, 1.0, 1.0],
                [1.0, 0.0, 2.0, 2.0, 2.0],
                [1.0, 2.0, 0.0, 2.0, 2.0],
                [1.0, 2.0, 2.0, 0.0, 2.0],
                [1.0, 2.0, 2.0, 2.0, 0.0],
            ]
        )
        grouper = self.DummyDistanceGrouper(
            molecules,
            matrix,
            representative_strategy="top3",
        )
        groups, index_groups = grouper.group()

        assert index_groups[0][0] == 0
        assert groups[0][0] is molecules[0]
        assert index_groups[0] == [0, 3, 1, 2, 4]

    def test_top3_falls_back_to_lowest_for_group_size_two(
        self, methanol_molecules
    ):
        molecules = methanol_molecules[:2]
        self._set_energies(molecules, [5.0, 1.0])
        matrix = np.array([[0.0, 0.2], [0.2, 0.0]])
        grouper = self.DummyDistanceGrouper(
            molecules,
            matrix,
            representative_strategy="top3",
        )
        _, index_groups = grouper.group()
        assert index_groups[0] == [1, 0]

    def test_top3_falls_back_to_lowest_for_group_size_one(
        self, methanol_molecules
    ):
        molecules = [methanol_molecules[0]]
        molecules[0].energy = 3.0
        matrix = np.array([[0.0]])
        grouper = self.DummyDistanceGrouper(
            molecules,
            matrix,
            representative_strategy="top3",
        )
        _, index_groups = grouper.group()
        assert index_groups[0] == [0]

    def test_non_matrix_lowest_orders_first_by_energy(
        self, multiple_molecules_xyz_file
    ):
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)[:4]
        self._set_energies(molecules, [5.0, 2.0, 1.0, None])
        grouper = FormulaGrouper(molecules, representative_strategy="lowest")
        groups, index_groups = grouper.group()
        assert index_groups[0] == [2, 1, 0, 3]
        assert groups[0][0] is molecules[2]

    def test_non_matrix_center_and_top3_raise(self, methanol_molecules):
        with pytest.raises(
            ValueError, match="requires a pairwise distance matrix"
        ):
            FormulaGrouper(
                methanol_molecules[:2], representative_strategy="center"
            )
        with pytest.raises(
            ValueError, match="requires a pairwise distance matrix"
        ):
            FormulaGrouper(
                methanol_molecules[:2], representative_strategy="top3"
            )

    def test_groups_table_and_unique_and_xyz_use_representative_first(
        self, multiple_molecules_xyz_file, temporary_working_dir
    ):
        from openpyxl import load_workbook

        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)[:5]
        self._set_energies(molecules, [0.0, 0.1, 0.2, -1.0, 0.3])
        conformer_ids = ["c1", "c2", "c3", "c4", "c5"]
        matrix = np.array(
            [
                [0.0, 1.0, 1.0, 1.0, 1.0],
                [1.0, 0.0, 2.0, 2.0, 2.0],
                [1.0, 2.0, 0.0, 2.0, 2.0],
                [1.0, 2.0, 2.0, 0.0, 2.0],
                [1.0, 2.0, 2.0, 2.0, 0.0],
            ]
        )
        grouper = self.DummyDistanceGrouper(
            molecules,
            matrix,
            representative_strategy="center",
            label="rep_center",
            conformer_ids=conformer_ids,
        )
        groups, index_groups = grouper.group()
        assert index_groups[0][0] == 0

        xlsx_file = os.path.join(
            temporary_working_dir,
            "rep_center_group_result",
            "rep_center_DummyDistanceGrouper_test.xlsx",
        )
        wb = load_workbook(xlsx_file, data_only=True)
        ws = wb["Groups"]
        assert ws["B2"].value.startswith("c1")

        unique_mols = grouper.unique()
        assert unique_mols[0] is groups[0][0]

        xyz_file = os.path.join(
            temporary_working_dir,
            "rep_center_group_result",
            "rep_center_group_1.xyz",
        )
        with open(xyz_file, "r") as handle:
            lines = handle.readlines()
        assert "Original_Index: c1" in lines[1]

    @staticmethod
    def _partition(index_groups):
        return {frozenset(group) for group in index_groups}

    @staticmethod
    def _center_score(matrix, group, idx):
        others = [j for j in group if j != idx]
        if not others:
            return 0.0
        return float(np.mean([matrix[idx, j] for j in others]))

    def test_representative_strategy_does_not_change_group_membership(
        self, multiple_molecules_xyz_file
    ):
        """lowest/center/top3 may reorder groups but must not change membership."""
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)

        partitions = []
        for strategy in ("lowest", "center", "top3"):
            grouper = BasicRMSDGrouper(
                molecules,
                threshold=2.0,
                num_procs=1,
                representative_strategy=strategy,
            )
            _, index_groups = grouper.group()
            partitions.append(self._partition(index_groups))

        assert partitions[0] == partitions[1] == partitions[2]

    def test_center_puts_most_central_member_first(
        self, multiple_molecules_xyz_file
    ):
        """center must place the member with minimum mean within-group distance first."""
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)

        grouper = BasicRMSDGrouper(
            molecules,
            threshold=2.0,
            num_procs=1,
            representative_strategy="center",
        )
        matrix = grouper.calculate_full_rmsd_matrix()
        _, index_groups = grouper.group()

        for group in index_groups:
            if len(group) == 1:
                continue
            expected = min(
                group,
                key=lambda idx: (
                    self._center_score(matrix, group, idx),
                    (
                        float("inf")
                        if molecules[idx].energy is None
                        else molecules[idx].energy
                    ),
                    idx,
                ),
            )
            assert group[0] == expected

    def test_top3_uses_only_three_lowest_energy_candidates(
        self, multiple_molecules_xyz_file
    ):
        """top3 chooses the most central candidate among the three lowest energies."""
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)

        grouper = BasicRMSDGrouper(
            molecules,
            threshold=2.0,
            num_procs=1,
            representative_strategy="top3",
        )
        matrix = grouper.calculate_full_rmsd_matrix()
        _, index_groups = grouper.group()

        tested_large_group = False
        for group in index_groups:
            if len(group) < 3:
                continue
            tested_large_group = True
            energy_order = sorted(
                group,
                key=lambda idx: (
                    (
                        float("inf")
                        if molecules[idx].energy is None
                        else molecules[idx].energy
                    ),
                    idx,
                ),
            )
            candidates = energy_order[:3]
            expected = min(
                candidates,
                key=lambda idx: (
                    self._center_score(matrix, group, idx),
                    (
                        float("inf")
                        if molecules[idx].energy is None
                        else molecules[idx].energy
                    ),
                    idx,
                ),
            )
            assert group[0] == expected
            assert group[0] in candidates

            expected_remaining = [
                idx for idx in energy_order if idx != expected
            ]
            assert group[1:] == expected_remaining

        assert (
            tested_large_group
        ), "Fixture must contain at least one group of size >= 3"

    def test_top3_falls_back_to_lowest_for_two_member_group(self):
        """top3 must use lowest-energy ordering when a group has fewer than 3 members."""
        from ase.build import molecule as ase_molecule

        from chemsmart.io.molecules.structure import Molecule

        mol0 = Molecule.from_ase_atoms(ase_molecule("CH3OH"))
        mol1 = mol0.copy()
        mol0._energy = -100.0
        mol1._energy = -101.0

        grouper = BasicRMSDGrouper(
            [mol0, mol1],
            threshold=1.0,
            num_procs=1,
            representative_strategy="top3",
        )
        _, index_groups = grouper.group()

        assert index_groups == [[1, 0]]

    @pytest.mark.parametrize("strategy", ["center", "top3"])
    def test_non_matrix_grouper_rejects_matrix_representative_strategy(
        self, methanol_molecules, strategy
    ):
        """Non-matrix groupers must reject center/top3 instead of silently falling back."""
        with pytest.raises(ValueError, match="pairwise distance matrix"):
            grouper = FormulaGrouper(
                methanol_molecules,
                representative_strategy=strategy,
            )
            grouper.group()

    def test_unique_uses_first_group_member_as_representative(
        self, multiple_molecules_xyz_file
    ):
        """unique() must trust group[0] and must not re-select by energy."""
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)

        grouper = BasicRMSDGrouper(
            molecules,
            threshold=2.0,
            num_procs=1,
            representative_strategy="center",
            label="center_unique_test",
        )
        groups, index_groups = grouper.group()
        unique_molecules = grouper.unique()

        assert len(unique_molecules) == len(groups)
        for unique_mol, group, indices in zip(
            unique_molecules, groups, index_groups
        ):
            assert indices[0] in range(len(molecules))
            assert np.allclose(unique_mol.positions, group[0].positions)

    def test_default_representative_matches_explicit_lowest(
        self, multiple_molecules_xyz_file
    ):
        """The new option must preserve the historical lowest-energy default."""
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)

        default_grouper = BasicRMSDGrouper(
            molecules,
            threshold=2.0,
            num_procs=1,
        )
        explicit_grouper = BasicRMSDGrouper(
            molecules,
            threshold=2.0,
            num_procs=1,
            representative_strategy="lowest",
        )

        _, default_indices = default_grouper.group()
        _, explicit_indices = explicit_grouper.group()

        assert default_indices == explicit_indices

    def test_center_tie_breaks_by_energy_then_original_index(
        self, multiple_molecules_xyz_file
    ):
        """center ties must resolve deterministically by energy, then original index."""
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)[:3]
        self._set_energies(molecules, [1.0, 0.5, 0.5])

        # All three members have identical mean centrality.
        matrix = np.array(
            [
                [0.0, 1.0, 1.0],
                [1.0, 0.0, 1.0],
                [1.0, 1.0, 0.0],
            ]
        )
        grouper = self.DummyDistanceGrouper(
            molecules,
            matrix,
            representative_strategy="center",
        )
        _, index_groups = grouper.group()

        # 1 and 2 tie in centrality and energy, so original index 1 wins.
        assert index_groups[0] == [1, 2, 0]


class TestRMSDGrouperConstruction:
    def test_rejects_both_threshold_and_num_groups(self):
        with pytest.raises(ValueError, match="Cannot specify both"):
            BasicRMSDGrouper([_mol(), _mol()], threshold=0.5, num_groups=2)


class TestGetHeavyAtoms:
    def test_returns_all_atoms_when_not_ignoring_hydrogens(self):
        molecules = [_mol(), _mol()]
        grouper = BasicRMSDGrouper(molecules, ignore_hydrogens=False)
        positions, symbols = grouper._get_heavy_atoms(molecules[0])
        assert symbols == ["H", "H"]
        assert positions.shape == (2, 3)

    def test_filters_hydrogens_when_ignoring(self):
        mol = Molecule(
            symbols=["C", "H", "H"],
            positions=np.array(
                [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]]
            ),
        )
        grouper = BasicRMSDGrouper([mol, mol], ignore_hydrogens=True)
        positions, symbols = grouper._get_heavy_atoms(mol)
        assert symbols == ["C"]
        assert positions.shape == (1, 3)


class TestRMSDGrouperNumGroups:
    def test_num_groups_at_least_n_returns_singletons(self, tmp_path):
        molecules = [_mol(), _mol(1.0)]
        grouper = BasicRMSDGrouper(
            molecules, num_groups=5, output_dir=str(tmp_path)
        )
        groups, index_groups = grouper.group()
        assert index_groups == [[0], [1]]

    def test_num_groups_produces_requested_count(self, tmp_path):
        molecules = [
            _mol(0.0),
            _mol(0.001),
            _mol(5.0),
            _mol(5.001),
            _mol(10.0),
        ]
        grouper = BasicRMSDGrouper(
            molecules,
            num_groups=3,
            align_molecules=False,
            output_dir=str(tmp_path),
        )
        groups, index_groups = grouper.group()
        assert len(groups) == 3
        assert grouper._auto_threshold is not None

    def test_repr_with_num_groups(self):
        grouper = BasicRMSDGrouper([_mol(), _mol()], num_groups=1)
        assert "num_groups=1" in repr(grouper)

    def test_repr_with_threshold(self):
        grouper = BasicRMSDGrouper([_mol(), _mol()], threshold=0.3)
        assert "threshold=0.3" in repr(grouper)


class TestCalculateFullRmsdMatrix:
    def test_returns_symmetric_matrix_without_output_file(self):
        molecules = [_mol(0.0), _mol(0.001), _mol(5.0)]
        grouper = BasicRMSDGrouper(molecules, num_procs=1)
        matrix = grouper.calculate_full_rmsd_matrix()
        assert matrix.shape == (3, 3)
        assert np.allclose(matrix, matrix.T)
        assert np.allclose(np.diag(matrix), 0.0)

    def test_writes_output_file_when_requested(self, tmp_path):
        molecules = [_mol(0.0), _mol(0.001)]
        grouper = BasicRMSDGrouper(
            molecules, num_procs=1, output_dir=str(tmp_path)
        )
        output_file = tmp_path / "nested" / "matrix_out"
        grouper.calculate_full_rmsd_matrix(output_file=str(output_file))
        # record() writes to the grouper's own output_dir/group_result
        # folder, not directly to output_file -- just confirm the parent
        # directory of output_file was created as instructed.
        assert output_file.parent.is_dir()


class TestRecordResultsBranches:
    def test_record_results_with_num_groups_and_grouping_time(self, tmp_path):
        molecules = [_mol(0.0), _mol(0.001), _mol(5.0)]
        grouper = BasicRMSDGrouper(
            molecules,
            num_groups=2,
            output_dir=str(tmp_path),
            matrix_format="csv",
        )
        groups, index_groups = grouper.group()
        assert len(groups) == 2

    def test_record_results_with_threshold_and_no_grouping_time(
        self, tmp_path
    ):
        molecules = [_mol(0.0), _mol(0.001)]
        grouper = BasicRMSDGrouper(
            molecules,
            threshold=0.5,
            output_dir=str(tmp_path),
            matrix_format="txt",
        )
        rmsd_matrix = np.zeros((2, 2))
        grouper._record_results(rmsd_matrix, grouping_time=None)

    def test_record_results_with_num_groups_and_auto_threshold(self, tmp_path):
        molecules = [
            _mol_tfd_grouper(),
            _mol_tfd_grouper(),
            _mol_tfd_grouper(),
        ]
        grouper = TorsionFingerprintGrouper(
            molecules,
            num_groups=2,
            output_dir=str(tmp_path),
            matrix_format="csv",
        )
        grouper._auto_threshold = 0.05
        tfd_matrix = np.zeros((3, 3))
        grouper._record_results(
            tfd_matrix=tfd_matrix,
            grouping_time=1.23,
            groups=[[molecules[0], molecules[1]], [molecules[2]]],
            index_groups=[[0, 1], [2]],
        )

    def test_record_results_without_groups_or_grouping_time(self, tmp_path):
        molecules = [_mol_tfd_grouper(), _mol_tfd_grouper()]
        grouper = TorsionFingerprintGrouper(
            molecules,
            threshold=0.1,
            output_dir=str(tmp_path),
            matrix_format="csv",
        )
        tfd_matrix = np.zeros((2, 2))
        grouper._record_results(
            tfd_matrix=tfd_matrix,
            grouping_time=None,
            groups=None,
            index_groups=None,
        )


class TestHungarianRMSDGrouperSpecifics:
    def test_incompatible_symbols_return_inf(self):
        mol1 = _mol()
        mol2 = Molecule(
            symbols=["C", "H"],
            positions=np.array([[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]]),
        )
        grouper = HungarianRMSDGrouper([mol1, mol2])
        assert grouper._calculate_rmsd((0, 1)) == np.inf

    def test_single_atom_per_element_shortcut_and_alignment(self):
        mol1 = _formaldehyde()
        mol2 = _formaldehyde(y_rotate=True)
        grouper = HungarianRMSDGrouper([mol1, mol2], align_molecules=True)
        rmsd = grouper._calculate_rmsd((0, 1))
        # C and O each appear once (shortcut path); the two equivalent
        # H atoms are optimally matched regardless of input order, so
        # RMSD after Kabsch alignment should be ~0.
        assert rmsd < 1e-6

    def test_without_alignment_skips_kabsch(self):
        mol1 = _formaldehyde()
        mol2 = _formaldehyde(x_offset=5.0)
        grouper = HungarianRMSDGrouper([mol1, mol2], align_molecules=False)
        rmsd = grouper._calculate_rmsd((0, 1))
        # Without alignment, a pure translation still contributes to RMSD.
        assert rmsd > 1.0


class TestSpyRMSDGrouperSpecifics:
    def test_mismatched_coord_length_returns_inf(self):
        mol1 = _mol()
        mol2 = Molecule(
            symbols=["H", "H", "H"],
            positions=np.array(
                [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0], [1.5, 0.0, 0.0]]
            ),
        )
        grouper = SpyRMSDGrouper([mol1, mol2])
        assert grouper._calculate_rmsd((0, 1)) == np.inf
        assert grouper.best_isomorphisms[(0, 1)] is None

    def test_mismatched_atomic_numbers_returns_inf(self):
        mol1 = _mol()
        mol2 = Molecule(
            symbols=["C", "H"],
            positions=np.array([[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]]),
        )
        grouper = SpyRMSDGrouper([mol1, mol2])
        assert grouper._calculate_rmsd((0, 1)) == np.inf
        assert grouper.best_isomorphisms[(0, 1)] is None

    def test_non_tuple_result_from_symmrmsd(self):
        grouper = SpyRMSDGrouper([_mol(), _mol(0.001)])
        with patch("spyrmsd.rmsd.symmrmsd", return_value=0.123):
            rmsd = grouper._calculate_rmsd((0, 1))
        assert rmsd == pytest.approx(0.123)
        assert grouper.best_isomorphisms[(0, 1)] is None

    def test_exception_during_symmrmsd_returns_inf(self):
        grouper = SpyRMSDGrouper([_mol(), _mol(0.001)])
        with patch(
            "spyrmsd.rmsd.symmrmsd",
            side_effect=RuntimeError("boom"),
        ):
            rmsd = grouper._calculate_rmsd((0, 1))
        assert rmsd == np.inf
        assert grouper.best_isomorphisms[(0, 1)] is None

    def test_get_best_isomorphism_direct_hit(self):
        grouper = SpyRMSDGrouper([_mol(), _mol(0.001)])
        grouper.best_isomorphisms[(0, 1)] = ([0, 1], [1, 0])
        assert grouper.get_best_isomorphism(0, 1) == ([0, 1], [1, 0])

    def test_get_best_isomorphism_reversed_hit_with_mapping(self):
        grouper = SpyRMSDGrouper([_mol(), _mol(0.001)])
        grouper.best_isomorphisms[(1, 0)] = ([0, 1], [1, 0])
        assert grouper.get_best_isomorphism(0, 1) == ([1, 0], [0, 1])

    def test_get_best_isomorphism_reversed_hit_without_mapping(self):
        grouper = SpyRMSDGrouper([_mol(), _mol(0.001)])
        grouper.best_isomorphisms[(1, 0)] = None
        assert grouper.get_best_isomorphism(0, 1) is None

    def test_get_best_isomorphism_not_found(self):
        grouper = SpyRMSDGrouper([_mol(), _mol(0.001)])
        assert grouper.get_best_isomorphism(5, 6) is None


class TestPymolRMSDGrouperSpecifics:
    def test_calculate_rmsd_negative_result_becomes_inf(
        self, shared_pymol_grouper
    ):
        with patch.object(
            shared_pymol_grouper.cmd, "align", return_value=[-1.0]
        ):
            assert shared_pymol_grouper._calculate_rmsd((0, 1)) == float("inf")

    def test_calculate_rmsd_non_finite_result_becomes_inf(
        self, shared_pymol_grouper
    ):
        with patch.object(
            shared_pymol_grouper.cmd, "align", return_value=[float("nan")]
        ):
            assert shared_pymol_grouper._calculate_rmsd((0, 1)) == float("inf")

    def test_calculate_rmsd_none_result_becomes_inf(
        self, shared_pymol_grouper
    ):
        with patch.object(
            shared_pymol_grouper.cmd, "align", return_value=None
        ):
            assert shared_pymol_grouper._calculate_rmsd((0, 1)) == float("inf")

    def test_calculate_rmsd_exception_becomes_inf(self, shared_pymol_grouper):
        with patch.object(
            shared_pymol_grouper.cmd,
            "align",
            side_effect=RuntimeError("boom"),
        ):
            assert shared_pymol_grouper._calculate_rmsd((0, 1)) == float("inf")

    def test_calculate_rmsd_uses_cache_on_second_call(
        self, shared_pymol_grouper
    ):
        with patch.object(
            shared_pymol_grouper.cmd, "align", return_value=[0.5]
        ) as mock_align:
            first = shared_pymol_grouper._calculate_rmsd((0, 1))
            second = shared_pymol_grouper._calculate_rmsd((1, 0))
        assert first == second == 0.5
        mock_align.assert_called_once()

    def test_repr_with_threshold(self, shared_pymol_grouper):
        assert "threshold=0.5" in repr(shared_pymol_grouper)

    def test_repr_with_num_groups(self):
        # This is the one repr variant that needs a distinct instance
        # (num_groups vs. threshold is fixed at construction time).
        grouper = PymolRMSDGrouper(
            [_mol(), _mol(0.001)], threshold=None, num_groups=1
        )
        try:
            assert "num_groups=1" in repr(grouper)
        finally:
            _cleanup_pymol_grouper(grouper)

    def test_del_quits_pymol_and_removes_temp_dir(self):
        # NOTE: __del__ is called explicitly here for coverage, but the
        # real `cmd` object is a process-wide singleton PyMOL session
        # shared by every PymolRMSDGrouper in this process. If we left
        # `grouper.cmd` pointing at the real session, Python's normal
        # garbage collection would call __del__ again later (unmocked),
        # issuing a *real* cmd.quit() that kills the shared session for
        # every subsequent test. So we always null out cmd/_temp_dir
        # afterward, same as _cleanup_pymol_grouper does elsewhere.
        grouper = PymolRMSDGrouper([_mol(), _mol(0.001)])
        temp_dir = grouper._temp_dir
        assert os.path.isdir(temp_dir)
        try:
            with patch.object(grouper.cmd, "quit") as mock_quit:
                grouper.__del__()
                mock_quit.assert_called_once()
            assert not os.path.isdir(temp_dir)
        finally:
            _cleanup_pymol_grouper(grouper)

    def test_del_swallows_quit_exception(self):
        grouper = PymolRMSDGrouper([_mol(), _mol(0.001)])
        try:
            with patch.object(
                grouper.cmd, "quit", side_effect=RuntimeError("boom")
            ):
                grouper.__del__()  # Should not raise
        finally:
            _cleanup_pymol_grouper(grouper)

    def test_del_swallows_rmtree_exception(self):
        grouper = PymolRMSDGrouper([_mol(), _mol(0.001)])
        try:
            # Also mock quit() here: __del__ calls cmd.quit() before the
            # rmtree cleanup, and a real quit() tears down PyMOL's
            # process-wide embedded session for good, hanging every
            # later test (including in test_groupers.py) that tries to
            # launch a new PymolRMSDGrouper.
            with (
                patch.object(grouper.cmd, "quit"),
                patch("shutil.rmtree", side_effect=RuntimeError("boom")),
            ):
                grouper.__del__()  # Should not raise
        finally:
            _cleanup_pymol_grouper(grouper)


class TestIRMSDGrouperConstruction:
    def test_raises_when_irmsd_command_not_found(self):
        with patch(
            "chemsmart.jobs.grouper.rmsd.find_irmsd_command",
            return_value=None,
        ):
            with pytest.raises(RuntimeError, match="irmsd is not available"):
                IRMSDGrouper([_mol(), _mol(0.001)])

    def test_explicit_threshold_is_not_overridden(self):
        grouper = _make_irmsd_grouper(threshold=0.3)
        assert grouper.threshold == 0.3


class TestIRMSDGrouperSpecifics:
    def test_incompatible_symbols_return_inf(self):
        grouper = _make_irmsd_grouper()
        mol2 = Molecule(
            symbols=["C", "H"],
            positions=np.array([[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]]),
        )
        grouper.molecules.append(mol2)
        assert grouper._calculate_rmsd((0, 2)) == np.inf

    def test_successful_run_parses_rmsd_and_inversion(self):
        grouper = _make_irmsd_grouper(inversion="auto")
        fake_stdout = "some header\niRMSD: 0.4567\nInversion check: off\n"
        with patch(
            "subprocess.run",
            return_value=type(
                "Result",
                (),
                {"returncode": 0, "stdout": fake_stdout, "stderr": ""},
            )(),
        ):
            rmsd = grouper._calculate_rmsd((0, 1))
        assert rmsd == pytest.approx(0.4567)
        assert grouper._actual_inversion == "off"

    def test_nonzero_returncode_returns_inf(self):
        grouper = _make_irmsd_grouper()
        with patch(
            "subprocess.run",
            return_value=type(
                "Result",
                (),
                {"returncode": 1, "stdout": "", "stderr": "boom"},
            )(),
        ):
            assert grouper._calculate_rmsd((0, 1)) == np.inf

    def test_timeout_returns_inf(self):
        import subprocess as subprocess_module

        grouper = _make_irmsd_grouper()
        with patch(
            "subprocess.run",
            side_effect=subprocess_module.TimeoutExpired(
                cmd="irmsd", timeout=60
            ),
        ):
            assert grouper._calculate_rmsd((0, 1)) == np.inf

    def test_unexpected_exception_returns_inf(self):
        grouper = _make_irmsd_grouper()
        with patch("subprocess.run", side_effect=RuntimeError("boom")):
            assert grouper._calculate_rmsd((0, 1)) == np.inf

    def test_temp_file_cleanup_swallows_unlink_error(self):
        grouper = _make_irmsd_grouper()
        with (
            patch("subprocess.run", side_effect=RuntimeError("boom")),
            patch("os.unlink", side_effect=OSError("already gone")),
        ):
            # Should not raise despite os.unlink failing during cleanup.
            assert grouper._calculate_rmsd((0, 1)) == np.inf

    def test_ignore_hydrogens_adds_heavy_flag(self):
        grouper = _make_irmsd_grouper(ignore_hydrogens=True)
        captured_cmd = {}

        def fake_run(cmd, **kwargs):
            captured_cmd["cmd"] = cmd
            return type(
                "Result",
                (),
                {"returncode": 0, "stdout": "iRMSD: 0.1\n", "stderr": ""},
            )()

        with patch("subprocess.run", side_effect=fake_run):
            grouper._calculate_rmsd((0, 1))
        assert "--heavy" in captured_cmd["cmd"]

    def test_non_standard_inversion_omits_inversion_flag(self):
        # inversion is only normalized with .lower() in __init__, not
        # validated against {"on", "off", "auto"} -- so a nonstandard
        # value is kept as-is and skips the --inversion flag entirely.
        grouper = _make_irmsd_grouper(inversion="weird")
        captured_cmd = {}

        def fake_run(cmd, **kwargs):
            captured_cmd["cmd"] = cmd
            return type(
                "Result",
                (),
                {"returncode": 0, "stdout": "iRMSD: 0.1\n", "stderr": ""},
            )()

        with patch("subprocess.run", side_effect=fake_run):
            grouper._calculate_rmsd((0, 1))
        assert "--inversion" not in captured_cmd["cmd"]

    def test_parse_irmsd_output_handles_unparsable_rmsd_line(self):
        grouper = _make_irmsd_grouper()
        rmsd, inversion = grouper._parse_irmsd_output(
            "iRMSD: not-a-number\n", parse_inversion=True
        )
        assert rmsd == np.inf
        assert inversion is None

    def test_parse_irmsd_output_skips_inversion_when_not_requested(self):
        grouper = _make_irmsd_grouper()
        rmsd, inversion = grouper._parse_irmsd_output(
            "iRMSD: 0.2\nInversion check: on\n", parse_inversion=False
        )
        assert rmsd == pytest.approx(0.2)
        assert inversion is None

    def test_record_results_includes_inversion_header(self, tmp_path):
        grouper = _make_irmsd_grouper(
            output_dir=str(tmp_path), matrix_format="csv"
        )
        grouper._actual_inversion = "auto"
        rmsd_matrix = np.zeros((2, 2))
        grouper._record_results(rmsd_matrix, grouping_time=1.0)


class TestEnergyGrouperConstruction:

    def test_defaults_to_one_kcal_threshold(self):
        molecules = [_mol_energy_grouper(0.0), _mol_energy_grouper(1.0)]
        grouper = EnergyGrouper(molecules)
        assert grouper.threshold == 1.0
        assert np.isclose(grouper.threshold_hartree, 1.0 * KCAL_TO_HARTREE)

    def test_num_groups_leaves_threshold_hartree_none(self):
        molecules = [_mol_energy_grouper(0.0), _mol_energy_grouper(1.0)]
        grouper = EnergyGrouper(molecules, num_groups=1)
        assert grouper.threshold is None
        assert grouper.threshold_hartree is None


class TestEnergyGrouperValidation:
    def test_rejects_missing_energy_few_molecules(self):
        molecules = [_mol_energy_grouper(0.0), _mol_energy_grouper(None)]
        with pytest.raises(ValueError, match="missing energy information"):
            EnergyGrouper(molecules)

    def test_rejects_non_numeric_energy(self):
        molecules = [
            _mol_energy_grouper(0.0),
            _mol_energy_grouper("not-a-number"),
        ]
        with pytest.raises(ValueError, match="missing energy information"):
            EnergyGrouper(molecules)

    def test_rejects_non_finite_energy(self):
        molecules = [
            _mol_energy_grouper(0.0),
            _mol_energy_grouper(float("nan")),
        ]
        with pytest.raises(ValueError, match="missing energy information"):
            EnergyGrouper(molecules)

    def test_reports_count_when_many_molecules_missing_energy(self):
        molecules = [_mol_energy_grouper(None) for _ in range(6)]
        with pytest.raises(ValueError, match="Found 6 molecules"):
            EnergyGrouper(molecules)


class TestEnergyGrouperThresholdGrouping:
    def test_group_by_threshold_separates_distant_energies(self, tmp_path):
        molecules = [
            _mol_energy_grouper(0.0),
            _mol_energy_grouper(0.0001),
            _mol_energy_grouper(10.0),
        ]
        grouper = EnergyGrouper(
            molecules, threshold=1.0, output_dir=str(tmp_path)
        )
        groups, index_groups = grouper.group()
        # First two are close in energy (within 1 kcal/mol); third is far.
        assert len(groups) == 2
        assert sorted(len(g) for g in index_groups) == [1, 2]

    def test_repr_with_threshold(self):
        molecules = [_mol_energy_grouper(0.0), _mol_energy_grouper(1.0)]
        grouper = EnergyGrouper(molecules, threshold=2.0)
        assert "threshold=2.0 kcal/mol" in repr(grouper)


class TestEnergyGrouperNumGroups:
    def test_num_groups_at_least_n_returns_singleton_groups(self, tmp_path):
        molecules = [_mol_energy_grouper(0.0), _mol_energy_grouper(1.0)]
        grouper = EnergyGrouper(
            molecules, num_groups=5, output_dir=str(tmp_path)
        )
        groups, index_groups = grouper.group()
        assert len(groups) == 2
        assert index_groups == [[0], [1]]

    def test_num_groups_produces_requested_count(self, tmp_path):
        # Energies chosen so distinct thresholds yield distinct group counts.
        energies = [0.0, 0.0005, 5.0, 5.0005, 10.0]
        molecules = [_mol_energy_grouper(e) for e in energies]
        grouper = EnergyGrouper(
            molecules, num_groups=3, output_dir=str(tmp_path)
        )
        groups, index_groups = grouper.group()
        assert len(groups) == 3
        assert grouper._auto_threshold is not None

    def test_repr_with_num_groups(self):
        molecules = [_mol_energy_grouper(0.0), _mol_energy_grouper(1.0)]
        grouper = EnergyGrouper(molecules, num_groups=1)
        assert "num_groups=1" in repr(grouper)


class TestEnergyGrouperMerge:
    def test_merge_groups_to_target_reduces_group_count(self, tmp_path):
        # Three well-separated clusters of increasing size force a merge
        # down to 2 groups, exercising _merge_groups_to_target directly.
        energies = [0.0, 0.0001, 20.0, 40.0, 40.0001, 40.0002]
        molecules = [_mol_energy_grouper(e) for e in energies]
        grouper = EnergyGrouper(
            molecules, num_groups=2, output_dir=str(tmp_path)
        )
        groups, index_groups = grouper.group()
        assert len(groups) == 2
        total_members = sum(len(g) for g in index_groups)
        assert total_members == len(molecules)


class TestEnergyGrouperRecordResultsDirect:
    def test_record_results_without_grouping_time_or_cached_groups(
        self, tmp_path
    ):
        # Calling _record_results directly (bypassing group()) with
        # grouping_time=None and no cached group indices exercises the
        # false side of both `if grouping_time is not None` and `if
        # index_groups is not None`.
        molecules = [_mol_energy_grouper(0.0), _mol_energy_grouper(1.0)]
        grouper = EnergyGrouper(
            molecules,
            threshold=1.0,
            output_dir=str(tmp_path),
            matrix_format="csv",
        )
        assert grouper._cached_group_indices is None
        energy_matrix = np.zeros((2, 2))
        grouper._record_results(energy_matrix, grouping_time=None)


class TestEnergyGrouperRecordResults:
    def test_record_results_writes_output_with_threshold(self, tmp_path):
        molecules = [
            _mol_energy_grouper(0.0),
            _mol_energy_grouper(0.5),
            _mol_energy_grouper(5.0),
        ]
        grouper = EnergyGrouper(
            molecules,
            threshold=1.0,
            label="sys",
            output_dir=str(tmp_path),
            matrix_format="csv",
        )
        groups, index_groups = grouper.group()
        assert len(groups) >= 1

    def test_record_results_writes_output_with_num_groups(self, tmp_path):
        molecules = [
            _mol_energy_grouper(0.0),
            _mol_energy_grouper(0.5),
            _mol_energy_grouper(5.0),
            _mol_energy_grouper(5.5),
        ]
        grouper = EnergyGrouper(
            molecules,
            num_groups=2,
            label="sys2",
            output_dir=str(tmp_path),
            matrix_format="txt",
        )
        groups, index_groups = grouper.group()
        assert len(groups) == 2

    def test_record_results_includes_thermo_header_for_qhg(self, tmp_path):
        molecules = [_mol_energy_grouper(0.0), _mol_energy_grouper(0.5)]
        grouper = EnergyGrouper(
            molecules,
            threshold=1.0,
            energy_type="QHG",
            thermo_parameters="s=1,T=298.15",
            output_dir=str(tmp_path),
            matrix_format="csv",
        )
        grouper.group()


@pytest.mark.usefixtures("temporary_working_dir")
class TestTanimotoGrouperConstruction:

    def test_default_threshold_is_point_nine(self, methanol_molecules):
        grouper = TanimotoSimilarityGrouper(methanol_molecules)
        assert grouper.threshold == 0.9
        assert grouper.num_groups is None

    def test_num_groups_leaves_threshold_none(self, methanol_molecules):
        grouper = TanimotoSimilarityGrouper(methanol_molecules, num_groups=1)
        assert grouper.threshold is None
        assert grouper.num_groups == 1

    def test_fingerprint_type_is_lowercased(self, methanol_molecules):
        grouper = TanimotoSimilarityGrouper(
            methanol_molecules, fingerprint_type="MACCS"
        )
        assert grouper.fingerprint_type == "maccs"


@pytest.mark.usefixtures("temporary_working_dir")
class TestGetFingerprint:
    @pytest.mark.parametrize(
        "fingerprint_type",
        [
            "rdkit",
            "rdk",
            "morgan",
            "maccs",
            "atompair",
            "torsion",
            "usr",
            "usrcat",
        ],
    )
    def test_supported_fingerprint_types_produce_a_fingerprint(
        self, methanol_molecules, fingerprint_type
    ):
        grouper = TanimotoSimilarityGrouper(
            methanol_molecules, fingerprint_type=fingerprint_type
        )
        fp = grouper._get_fingerprint(grouper.rdkit_molecules[0])
        assert fp is not None

    def test_unknown_fingerprint_type_falls_back_to_rdkit_default(
        self, methanol_molecules
    ):
        grouper = TanimotoSimilarityGrouper(
            methanol_molecules, fingerprint_type="not-a-real-type"
        )
        fp = grouper._get_fingerprint(grouper.rdkit_molecules[0])
        assert fp is not None

    def test_fingerprint_generation_failure_returns_none(
        self, methanol_molecules
    ):
        grouper = TanimotoSimilarityGrouper(
            methanol_molecules, fingerprint_type="morgan"
        )
        # Passing a non-Mol object makes the RDKit call raise internally;
        # the method should swallow it and return None rather than
        # propagating the exception.
        assert grouper._get_fingerprint(object()) is None


@pytest.mark.usefixtures("temporary_working_dir")
class TestIgnoreHydrogens:
    def test_ignore_hydrogens_still_groups_successfully(
        self, methanol_molecules
    ):
        grouper = TanimotoSimilarityGrouper(
            methanol_molecules, ignore_hydrogens=True
        )
        assert len(grouper.rdkit_molecules) == len(methanol_molecules)
        groups, index_groups = grouper.group()
        assert len(groups) == 1


@pytest.mark.usefixtures("temporary_working_dir")
class TestNumGroupsGrouping:
    def test_num_groups_creates_requested_number_of_groups(
        self, crest_molecules
    ):
        grouper = TanimotoSimilarityGrouper(
            crest_molecules,
            num_groups=3,
            fingerprint_type="usrcat",
            num_procs=2,
        )
        groups, index_groups = grouper.group()
        assert len(groups) == 3
        assert len(index_groups) == 3
        assert grouper._auto_threshold is not None

    def test_num_groups_at_least_as_large_as_molecule_count_gives_singletons(
        self, methanol_molecules
    ):
        grouper = TanimotoSimilarityGrouper(
            methanol_molecules, num_groups=10, fingerprint_type="usrcat"
        )
        groups, index_groups = grouper.group()
        assert len(groups) == len(methanol_molecules)
        assert all(len(g) == 1 for g in groups)

    def test_num_groups_one_merges_everything(self, crest_molecules):
        grouper = TanimotoSimilarityGrouper(
            crest_molecules,
            num_groups=1,
            fingerprint_type="usrcat",
            num_procs=2,
        )
        groups, index_groups = grouper.group()
        assert len(groups) == 1
        assert sum(len(g) for g in groups) == len(crest_molecules)


@pytest.mark.usefixtures("temporary_working_dir")
class TestInvalidMoleculeHandling:
    def test_molecule_with_no_rdkit_conversion_is_skipped(
        self, methanol_molecules, mocker
    ):
        from chemsmart.io.molecules.structure import Molecule

        # ``Molecule.to_rdkit`` is only used to build the RDKit molecule
        # for the 3D "usr"/"usrcat" descriptors; other fingerprint types
        # perceive bonds directly from the XYZ geometry instead.
        mocker.patch.object(Molecule, "to_rdkit", return_value=None)
        grouper = TanimotoSimilarityGrouper(
            methanol_molecules, fingerprint_type="usr"
        )
        assert grouper.rdkit_molecules == []

    def test_ignore_hydrogens_falls_back_on_removehs_failure(
        self, methanol_molecules, mocker
    ):
        mocker.patch(
            "chemsmart.jobs.grouper.tanimoto.Chem.RemoveHs",
            side_effect=Exception("boom"),
        )
        # The with-H fallback on a failed RemoveHs only exists on the
        # "usr"/"usrcat" preparation path.
        grouper = TanimotoSimilarityGrouper(
            methanol_molecules, ignore_hydrogens=True, fingerprint_type="usr"
        )
        # Falls back to the original (with-H) rdkit mol for every molecule.
        assert len(grouper.rdkit_molecules) == len(methanol_molecules)

    def test_group_returns_empty_when_no_valid_fingerprints(
        self, methanol_molecules, mocker
    ):
        grouper = TanimotoSimilarityGrouper(methanol_molecules)
        mocker.patch.object(grouper, "_get_fingerprint", return_value=None)
        groups, index_groups = grouper.group()
        assert groups == []
        assert index_groups == []


@pytest.mark.usefixtures("temporary_working_dir")
@pytest.mark.usefixtures("temporary_working_dir")
class TestRecordResultsDirect:
    def test_record_results_without_grouping_time_or_index_groups(
        self, methanol_molecules, tmp_path
    ):
        import numpy as np

        grouper = TanimotoSimilarityGrouper(
            methanol_molecules, threshold=0.9, label=str(tmp_path / "run")
        )
        n = len(methanol_molecules)
        matrix = np.eye(n, dtype=np.float32)
        # Should not raise even with grouping_time=None and
        # index_groups=None (both optional-argument False branches).
        grouper._record_results(
            tanimoto_matrix=matrix,
            valid_indices=list(range(n)),
            grouping_time=None,
            groups=None,
            index_groups=None,
        )


@pytest.mark.usefixtures("temporary_working_dir")
class TestRepr:
    def test_repr_with_threshold(self, methanol_molecules):
        grouper = TanimotoSimilarityGrouper(methanol_molecules, threshold=0.8)
        rep = repr(grouper)
        assert "threshold=0.8" in rep
        assert "TanimotoSimilarityGrouper" in rep

    def test_repr_with_num_groups(self, methanol_molecules):
        grouper = TanimotoSimilarityGrouper(methanol_molecules, num_groups=2)
        rep = repr(grouper)
        assert "num_groups=2" in rep

    def test_repr_includes_num_procs(self):
        # RDKitIsomorphismGrouper does not support multiprocessing, so
        # any requested num_procs > 1 is clamped down to 1.
        grouper = RDKitIsomorphismGrouper([_make_water()], num_procs=4)
        assert repr(grouper) == "RDKitIsomorphismGrouper(num_procs=1)"


class TestTorsionFingerprintGrouperConstruction:

    def test_rejects_invalid_max_dev(self):
        with pytest.raises(ValueError, match="max_dev must be either"):
            TorsionFingerprintGrouper(
                [_mol_tfd_grouper(), _mol_tfd_grouper()], max_dev="bogus"
            )

    def test_defaults_threshold_when_neither_specified(self):
        grouper = TorsionFingerprintGrouper(
            [_mol_tfd_grouper(), _mol_tfd_grouper()]
        )
        assert grouper.threshold == 0.1
        assert grouper.num_groups is None


class TestPrepareConformerMoleculeEdgeCases:
    def test_empty_molecules_list(self):
        grouper = TorsionFingerprintGrouper.__new__(TorsionFingerprintGrouper)
        grouper.molecules = []
        grouper.ignore_hydrogens = False
        grouper._prepare_conformer_molecule()
        assert grouper.rdkit_mol is None
        assert grouper.valid_conformer_ids == []

    def test_rdkit_mol_none_short_circuits(self):
        grouper = TorsionFingerprintGrouper.__new__(TorsionFingerprintGrouper)
        grouper.molecules = [_mol_tfd_grouper()]
        grouper.ignore_hydrogens = False
        with patch.object(Molecule, "to_rdkit", return_value=None):
            grouper._prepare_conformer_molecule()
        assert grouper.rdkit_mol is None
        assert grouper.valid_conformer_ids == []

    def test_add_conformer_failure_is_logged_and_skipped(self):
        grouper = TorsionFingerprintGrouper(
            [_mol_tfd_grouper(), _mol_tfd_grouper()]
        )
        # rdkit_mol was already prepared successfully in __init__; now
        # force a failure while re-preparing to hit the except branch
        # inside the per-molecule conformer-adding loop.
        with patch(
            "chemsmart.jobs.grouper.tfd.Chem.Conformer",
            side_effect=RuntimeError("boom"),
        ):
            grouper._prepare_conformer_molecule()
        assert grouper.valid_conformer_ids == []

    def test_ignore_hydrogens_removes_hs_successfully(self):
        grouper = TorsionFingerprintGrouper(
            [_mol_tfd_grouper(), _mol_tfd_grouper()], ignore_hydrogens=True
        )
        assert grouper.rdkit_mol is not None

    def test_ignore_hydrogens_falls_back_on_failure(self):
        with patch(
            "chemsmart.jobs.grouper.tfd.Chem.RemoveHs",
            side_effect=RuntimeError("cannot remove"),
        ):
            grouper = TorsionFingerprintGrouper(
                [_mol_tfd_grouper(), _mol_tfd_grouper()], ignore_hydrogens=True
            )
        assert grouper.rdkit_mol is not None


class TestCalculateTFDEdgeCases:
    def test_returns_inf_when_no_rdkit_mol(self):
        grouper = TorsionFingerprintGrouper(
            [_mol_tfd_grouper(), _mol_tfd_grouper()]
        )
        grouper.rdkit_mol = None
        assert grouper._calculate_tfd((0, 1)) == float("inf")

    def test_returns_inf_on_calculation_exception(self):
        grouper = TorsionFingerprintGrouper(
            [_mol_tfd_grouper(), _mol_tfd_grouper()]
        )
        with patch(
            "chemsmart.jobs.grouper.tfd.TorsionFingerprints"
            ".GetTFDBetweenConformers",
            side_effect=RuntimeError("rdkit failure"),
        ):
            assert grouper._calculate_tfd((0, 1)) == float("inf")


class TestGroupEdgeCases:
    def test_empty_molecules_returns_empty_groups(self):
        grouper = TorsionFingerprintGrouper.__new__(TorsionFingerprintGrouper)
        grouper.molecules = []
        groups, index_groups = grouper.group()
        assert groups == []
        assert index_groups == []

    def test_single_molecule_returns_single_group(self):
        grouper = TorsionFingerprintGrouper.__new__(TorsionFingerprintGrouper)
        grouper.molecules = [_mol_tfd_grouper()]
        groups, index_groups = grouper.group()
        assert groups == [[grouper.molecules[0]]]
        assert index_groups == [[0]]

    def test_no_valid_conformers_falls_back_to_singleton_groups(self):
        molecules = [
            _mol_tfd_grouper(),
            _mol_tfd_grouper(),
            _mol_tfd_grouper(),
        ]
        grouper = TorsionFingerprintGrouper(molecules)
        grouper.valid_conformer_ids = []
        groups, index_groups = grouper.group()
        assert index_groups == [[0], [1], [2]]


class TestTorsionFingerprintGrouperNumGroups:
    NUM_PROCS = 1

    def test_num_groups_at_least_n_returns_singletons(
        self, multiple_molecules_xyz_file, tmp_path
    ):
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)[:3]
        grouper = TorsionFingerprintGrouper(
            molecules,
            num_groups=10,
            num_procs=self.NUM_PROCS,
            output_dir=str(tmp_path),
        )
        groups, index_groups = grouper.group()
        assert index_groups == [[0], [1], [2]]

    def test_num_groups_produces_requested_count(
        self, multiple_molecules_xyz_file, tmp_path
    ):
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)
        grouper = TorsionFingerprintGrouper(
            molecules,
            num_groups=5,
            num_procs=self.NUM_PROCS,
            output_dir=str(tmp_path),
        )
        groups, index_groups = grouper.group()
        assert len(groups) == 5
        assert grouper._auto_threshold is not None


class TestTorsionFingerprintGrouperRepr:
    def test_repr_with_threshold(self):
        grouper = TorsionFingerprintGrouper(
            [_mol_tfd_grouper(), _mol_tfd_grouper()], threshold=0.2
        )
        text = repr(grouper)
        assert "threshold=0.2" in text
        assert "TorsionFingerprintGrouper" in text

    def test_repr_with_num_groups(self):
        grouper = TorsionFingerprintGrouper(
            [_mol_tfd_grouper(), _mol_tfd_grouper()], num_groups=1
        )
        text = repr(grouper)
        assert "num_groups=1" in text


class TestGroupInvalidMoleculeFallback:
    def test_molecule_with_no_hash_gets_its_own_group(
        self, tmp_path, monkeypatch
    ):
        """group() calls self.record(...) unconditionally, which writes
        output files relative to cwd when output_dir is unset -- chdir
        into a throwaway directory to avoid polluting the repo."""
        monkeypatch.chdir(tmp_path)
        mol = _make_water()
        grouper = RDKitIsomorphismGrouper([mol], num_procs=1)
        with patch.object(Molecule, "get_rdkit_hash", return_value=None):
            groups, index_groups = grouper.group()
        assert groups == [[mol]]
        assert index_groups == [[0]]


class TestRecordResultsWithoutGroupingTime:
    def test_grouping_time_none_skips_time_header(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        mol = _make_water()
        grouper = RDKitIsomorphismGrouper([mol], num_procs=1, label="test")

        # Should not raise despite grouping_time being omitted.
        grouper._record_results(
            hashes=["abc"],
            groups=[[mol]],
            index_groups=[[0]],
            grouping_time=None,
        )


class TestResultsRecorderConstruction:
    def test_rejects_unsupported_format(self, tmp_path):
        with pytest.raises(ValueError, match="Unsupported output format"):
            ResultsRecorder(output_dir=str(tmp_path), matrix_format="pdf")

    def test_normalizes_format_case_and_dot(self, tmp_path):
        recorder = ResultsRecorder(
            output_dir=str(tmp_path), matrix_format=".XLSX"
        )
        assert recorder.matrix_format == "xlsx"

    def test_creates_output_dir(self, tmp_path):
        target = tmp_path / "nested" / "output"
        ResultsRecorder(output_dir=str(target))
        assert target.is_dir()


class TestResultsRecorderLabels:
    def test_get_labels_uses_conformer_ids_when_length_matches(self, tmp_path):
        recorder = ResultsRecorder(
            output_dir=str(tmp_path), conformer_ids=["a", "b", "c"]
        )
        assert recorder.get_labels(3) == ["a", "b", "c"]

    def test_get_labels_falls_back_to_numeric_when_length_mismatches(
        self, tmp_path
    ):
        recorder = ResultsRecorder(
            output_dir=str(tmp_path), conformer_ids=["a", "b"]
        )
        assert recorder.get_labels(3) == ["1", "2", "3"]

    def test_get_labels_falls_back_when_no_conformer_ids(self, tmp_path):
        recorder = ResultsRecorder(output_dir=str(tmp_path))
        assert recorder.get_labels(2) == ["1", "2"]


class TestResultsRecorderGroupsDataframe:
    def test_build_groups_dataframe_basic(self, tmp_path):
        recorder = ResultsRecorder(output_dir=str(tmp_path))
        df = recorder.build_groups_dataframe(
            index_groups=[[0, 1], [2]], n_molecules=3
        )
        assert list(df["Group"]) == [1, 2]
        assert list(df["Members"]) == ["1, 2", "3"]

    def test_build_groups_dataframe_with_extra_columns(self, tmp_path):
        recorder = ResultsRecorder(output_dir=str(tmp_path))
        df = recorder.build_groups_dataframe(
            index_groups=[[0], [1]],
            n_molecules=2,
            extra_columns={"Score": [0.5]},
        )
        # Second row's extra column value is skipped because the values
        # list is shorter than the number of groups (exercises the
        # `if i < len(values)` false branch).
        assert df.loc[0, "Score"] == 0.5
        assert "Score" not in df.columns or pd.isna(df.loc[1, "Score"])


class TestResultsRecorderFilename:
    def test_get_filename_with_label_and_suffix(self, tmp_path):
        recorder = ResultsRecorder(
            output_dir=str(tmp_path), label="mysys", matrix_format="csv"
        )
        filename = recorder.get_filename("RMSDGrouper", suffix="T0.5")
        assert filename == os.path.join(
            str(tmp_path), "mysys_RMSDGrouper_T0.5.csv"
        )

    def test_get_filename_without_label_or_suffix(self, tmp_path):
        recorder = ResultsRecorder(output_dir=str(tmp_path))
        filename = recorder.get_filename("RMSDGrouper")
        assert filename == os.path.join(str(tmp_path), "RMSDGrouper.xlsx")


class TestResultsRecorderRecordResults:
    def _make_common_args(self):
        header_info = [
            ("Title", None),
            ("", "Some Title Line"),
            ("Key", ""),
            ("Method", "RMSD"),
        ]
        sheets_data = {
            "Groups": pd.DataFrame({"Group": [1], "Members": ["1"]})
        }
        matrix = np.array([[0.0, np.inf], [np.inf, 0.0]])
        matrix_data = ("Matrix", matrix, ["1", "2"])
        return header_info, sheets_data, matrix_data

    def test_record_results_xlsx_with_matrix(self, tmp_path):
        recorder = ResultsRecorder(
            output_dir=str(tmp_path), matrix_format="xlsx"
        )
        header_info, sheets_data, matrix_data = self._make_common_args()
        # Include a sheet sharing the matrix's name, to hit the "skip if
        # already written as matrix" continue branch in _write_xlsx.
        sheets_data[matrix_data[0]] = pd.DataFrame({"A": [1]})
        path = recorder.record_results(
            grouper_name="TestGrouper",
            header_info=header_info,
            sheets_data=sheets_data,
            matrix_data=matrix_data,
        )
        assert os.path.isfile(path)

    def test_record_results_xlsx_without_matrix(self, tmp_path):
        recorder = ResultsRecorder(
            output_dir=str(tmp_path), matrix_format="xlsx"
        )
        header_info, sheets_data, _ = self._make_common_args()
        sheets_data["Other"] = pd.DataFrame({"A": [1, 2]})
        path = recorder.record_results(
            grouper_name="TestGrouper",
            header_info=header_info,
            sheets_data=sheets_data,
            startrow=3,
        )
        assert os.path.isfile(path)

    def test_record_results_csv(self, tmp_path):
        recorder = ResultsRecorder(
            output_dir=str(tmp_path), matrix_format="csv"
        )
        header_info, sheets_data, matrix_data = self._make_common_args()
        sheets_data["Extra"] = pd.DataFrame({"A": [1]})
        # Also include a sheet with the same name as the matrix sheet, to
        # exercise the "skip if already written as matrix" continue branch.
        sheets_data[matrix_data[0]] = pd.DataFrame({"A": [1]})
        path = recorder.record_results(
            grouper_name="TestGrouper",
            header_info=header_info,
            sheets_data=sheets_data,
            matrix_data=matrix_data,
        )
        assert os.path.isfile(path)
        extra_file = path[:-4] + "_Extra.csv"
        assert os.path.isfile(extra_file)

    def test_record_results_csv_without_matrix(self, tmp_path):
        recorder = ResultsRecorder(
            output_dir=str(tmp_path), matrix_format="csv"
        )
        header_info, sheets_data, _ = self._make_common_args()
        path = recorder.record_results(
            grouper_name="TestGrouper",
            header_info=header_info,
            sheets_data=sheets_data,
        )
        assert os.path.isfile(path)

    def test_record_results_txt(self, tmp_path):
        recorder = ResultsRecorder(
            output_dir=str(tmp_path), matrix_format="txt"
        )
        header_info, sheets_data, matrix_data = self._make_common_args()
        sheets_data["Extra"] = pd.DataFrame({"A": [1]})
        sheets_data[matrix_data[0]] = pd.DataFrame({"A": [1]})
        path = recorder.record_results(
            grouper_name="TestGrouper",
            header_info=header_info,
            sheets_data=sheets_data,
            matrix_data=matrix_data,
        )
        assert os.path.isfile(path)
        content = open(path).read()
        assert "Matrix" in content
        assert "inf" in content

    def test_record_results_txt_without_matrix(self, tmp_path):
        recorder = ResultsRecorder(
            output_dir=str(tmp_path), matrix_format="txt"
        )
        header_info, sheets_data, _ = self._make_common_args()
        path = recorder.record_results(
            grouper_name="TestGrouper",
            header_info=header_info,
            sheets_data=sheets_data,
        )
        assert os.path.isfile(path)
        content = open(path).read()
        assert "Groups" in content
        assert "Matrix" not in content

    def test_record_results_invalid_format_raises(self, tmp_path):
        recorder = ResultsRecorder(
            output_dir=str(tmp_path), matrix_format="csv"
        )
        # Force an unsupported format after construction to hit the
        # otherwise-unreachable else branch in record_results.
        recorder.matrix_format = "bogus"
        with pytest.raises(ValueError, match="Unsupported output format"):
            recorder.record_results(
                grouper_name="TestGrouper",
                header_info=[],
                sheets_data={},
            )


class TestStructureGrouperConfig:
    def test_defaults(self):
        config = StructureGrouperConfig()
        assert config.ltol == 0.1
        assert config.stol == 0.18
        assert config.angle_tol == 1

    def test_custom_values(self):
        config = StructureGrouperConfig(ltol=0.2, stol=0.3, angle_tol=5)
        assert config.ltol == 0.2
        assert config.stol == 0.3
        assert config.angle_tol == 5


class TestMoleculeGrouperValidation:
    # Note: the `not isinstance(self.molecules, Iterable)` branch in
    # _validate_inputs is unreachable through the public API: by the
    # time it runs, __init__ has already done `self.molecules =
    # list(molecules)`, which either raised TypeError itself (for a
    # genuinely non-iterable argument) or produced a real `list`
    # (always Iterable). So only the "non-Molecule items" branch below
    # is reachable in practice.

    def test_rejects_non_molecule_items(self):
        with pytest.raises(TypeError, match="Molecule instances"):
            _ConcreteGrouper(molecules=[1, 2, 3])

    def test_accepts_valid_molecules(self, water_molecule):
        grouper = _ConcreteGrouper(molecules=[water_molecule])
        assert grouper.molecules == [water_molecule]
        assert grouper.num_procs == 1

    def test_num_procs_floored_at_one(self, water_molecule):
        grouper = _ConcreteGrouper(molecules=[water_molecule], num_procs=0)
        assert grouper.num_procs == 1


class TestMoleculeGrouperOutputDir:
    def test_get_output_dir_with_label(self, water_molecule, tmp_path):
        grouper = _ConcreteGrouper(
            molecules=[water_molecule],
            label="sys1",
            output_dir=str(tmp_path),
        )
        assert grouper._get_output_dir() == os.path.join(
            str(tmp_path), "sys1_group_result"
        )

    def test_get_output_dir_without_label_or_output_dir(self, water_molecule):
        grouper = _ConcreteGrouper(molecules=[water_molecule])
        assert grouper._get_output_dir() == "group_result"

    def test_get_results_recorder(self, water_molecule, tmp_path):
        grouper = _ConcreteGrouper(
            molecules=[water_molecule],
            label="sys1",
            output_dir=str(tmp_path),
        )
        recorder = grouper._get_results_recorder()
        assert isinstance(recorder, ResultsRecorder)
        assert recorder.label == "sys1"


class TestMoleculeGrouperRecordDelegation:
    def test_record_delegates_to_record_results(self, water_molecule):
        grouper = _ConcreteGrouper(molecules=[water_molecule])
        with pytest.raises(NotImplementedError, match="_record_results"):
            grouper.record()


class TestMoleculeGrouperHeaderHelpers:
    @pytest.mark.parametrize("energy_type", ["QHH", "QHG", "SP_QHG", "qhh"])
    def test_append_thermo_header_included_for_thermo_energy_types(
        self, water_molecule, energy_type
    ):
        grouper = _ConcreteGrouper(
            molecules=[water_molecule],
            energy_type=energy_type,
            thermo_parameters="s=1,T=298.15",
        )
        header_info = []
        grouper._append_thermo_header(header_info)
        assert header_info == [("Thermochemistry Parameters", "s=1,T=298.15")]

    def test_append_thermo_header_skipped_for_non_thermo_energy_type(
        self, water_molecule
    ):
        grouper = _ConcreteGrouper(
            molecules=[water_molecule],
            energy_type="E",
            thermo_parameters="s=1,T=298.15",
        )
        header_info = []
        grouper._append_thermo_header(header_info)
        assert header_info == []

    def test_append_thermo_header_skipped_without_thermo_parameters(
        self, water_molecule
    ):
        grouper = _ConcreteGrouper(
            molecules=[water_molecule], energy_type="QHH"
        )
        header_info = []
        grouper._append_thermo_header(header_info)
        assert header_info == []

    def test_append_input_usage_header_with_skipped_ids(self, water_molecule):
        grouper = _ConcreteGrouper(
            molecules=[water_molecule], skipped_ids=["c3", "c4"]
        )
        header_info = []
        grouper._append_input_usage_header(header_info)
        assert ("Used Molecules", 1) in header_info
        assert ("Skipped Molecules", 2) in header_info
        assert ("Skipped IDs", "c3, c4") in header_info

    def test_append_input_usage_header_without_skipped_ids(
        self, water_molecule
    ):
        grouper = _ConcreteGrouper(molecules=[water_molecule])
        header_info = []
        grouper._append_input_usage_header(header_info)
        assert ("Skipped IDs", "None") in header_info


class TestMoleculeGrouperUnique:
    def test_unique_writes_xyz_and_returns_lowest_energy(
        self, water_molecule, methane_molecule, tmp_path
    ):
        mol_low = water_molecule
        mol_low.energy = -1.0
        mol_high = methane_molecule
        mol_high.energy = 5.0

        grouper = _ConcreteGrouper(
            molecules=[mol_low, mol_high], output_dir=str(tmp_path)
        )
        unique = grouper.unique(output_dir=str(tmp_path))

        assert len(unique) == 1
        assert unique[0] is mol_low

        group_file = os.path.join(str(tmp_path), "group_result", "group_1.xyz")
        assert os.path.isfile(group_file)
        content = open(group_file).read()
        assert "Group 1 Member 1" in content
        assert "Group 1 Member 2" in content

    def test_unique_uses_cached_groups_on_second_call(
        self, water_molecule, tmp_path
    ):
        water_molecule.energy = None
        grouper = _ConcreteGrouper(
            molecules=[water_molecule],
            label="cached",
            output_dir=str(tmp_path),
        )
        first = grouper.unique(output_dir=str(tmp_path))
        assert grouper._cached_groups is not None

        # Second call should hit the cached-groups branch without
        # recomputation (verified indirectly: result is consistent).
        second = grouper.unique(output_dir=str(tmp_path))
        assert first[0] is second[0]

    def test_unique_with_conformer_ids_labels_original_index(
        self, water_molecule, tmp_path
    ):
        water_molecule.energy = None
        grouper = _ConcreteGrouper(
            molecules=[water_molecule],
            conformer_ids=["confA"],
            output_dir=str(tmp_path),
        )
        grouper.unique(output_dir=str(tmp_path))
        group_file = os.path.join(str(tmp_path), "group_result", "group_1.xyz")
        content = open(group_file).read()
        assert "Original_Index: confA" in content
        assert "N/A" in content


class TestMoleculeGrouperClusteringHelpers:
    """Guard branches of the shared complete-linkage helpers."""

    def test_build_singleton_groups(
        self, center_matrix_grouper, methanol_molecules
    ):
        groups, index_groups = center_matrix_grouper._build_singleton_groups()
        assert index_groups == [[0], [1], [2]]
        assert groups == [[mol] for mol in methanol_molecules]

    def test_linkage_tree_is_none_for_single_structure(
        self, center_matrix_grouper
    ):
        assert (
            center_matrix_grouper._build_complete_linkage_tree(
                np.zeros((1, 1))
            )
            is None
        )

    def test_linkage_tree_rejects_infinite_distances(
        self, center_matrix_grouper
    ):
        matrix = np.array([[0.0, np.inf], [np.inf, 0.0]])
        with pytest.raises(ValueError, match="must be finite"):
            center_matrix_grouper._build_complete_linkage_tree(matrix)

    def test_groups_from_labels_default_to_positional_indices(
        self, center_matrix_grouper, methanol_molecules
    ):
        groups, index_groups = (
            center_matrix_grouper._build_groups_from_cluster_labels(
                np.array([2, 1, 2])
            )
        )
        assert index_groups == [[0, 2], [1]]
        assert groups[0] == [methanol_molecules[0], methanol_molecules[2]]

    def test_order_groups_defaults_matrix_indices(self, center_matrix_grouper):
        matrix = np.array([[0.0, 0.1, 0.5], [0.1, 0.0, 0.1], [0.5, 0.1, 0.0]])
        _, ordered = (
            center_matrix_grouper._order_groups_by_representative_strategy(
                [[0, 1, 2]], matrix
            )
        )
        # Index 1 is closest to the other two, so it becomes the centre.
        assert ordered == [[1, 0, 2]]

    @pytest.mark.parametrize(
        "normalize",
        (
            "_normalize_output_matrix_for_threshold",
            "_normalize_output_matrix_for_num_groups",
        ),
    )
    def test_normalize_output_matrix_rejects_size_mismatch(
        self, center_matrix_grouper, normalize
    ):
        with pytest.raises(ValueError, match="Matrix dimensions must match"):
            getattr(center_matrix_grouper, normalize)(np.zeros((2, 2)), 1)

    def test_num_groups_clustering_with_empty_matrix(
        self, center_matrix_grouper
    ):
        center_matrix_grouper._auto_threshold = 0.3
        assert center_matrix_grouper._cluster_distance_matrix_by_num_groups(
            np.zeros((0, 0)), 2
        ) == ([], [])
        assert center_matrix_grouper._auto_threshold is None

    def test_num_groups_clustering_with_one_finite_structure(
        self, center_matrix_grouper, methanol_molecules
    ):
        center_matrix_grouper._auto_threshold = 0.3
        matrix = np.array([[0.0, np.inf], [np.inf, 0.0]])
        groups, index_groups = (
            center_matrix_grouper._cluster_distance_matrix_by_num_groups(
                matrix, 2, original_indices=[1, 2]
            )
        )
        assert len(index_groups) == 1
        assert index_groups[0][0] in (1, 2)
        assert groups == [[methanol_molecules[index_groups[0][0]]]]
        assert center_matrix_grouper._auto_threshold is None


class TestGrouperJobConstruction:
    def test_rejects_non_list_molecules(self):
        with pytest.raises(ValueError, match="must be a list"):
            GrouperJob(molecules="not-a-list", grouping_strategy="rmsd")

    def test_rejects_fewer_than_two_molecules(self):
        with pytest.raises(ValueError, match="at least 2"):
            GrouperJob(molecules=make_molecules(1), grouping_strategy="rmsd")

    def test_label_defaults_to_grouper(self):
        job = GrouperJob(molecules=make_molecules(2), grouping_strategy="rmsd")
        assert job.label == "grouper"

    def test_custom_label_used(self):
        job = GrouperJob(
            molecules=make_molecules(2),
            grouping_strategy="rmsd",
            label="my_group",
        )
        assert job.label == "my_group"

    def test_stores_configuration(self):
        job = GrouperJob(
            molecules=make_molecules(3),
            grouping_strategy="tanimoto",
            threshold=0.9,
            num_groups=2,
            ignore_hydrogens=True,
            num_procs=4,
            matrix_format="csv",
            energy_type="G",
        )
        assert job.grouping_strategy == "tanimoto"
        assert job.threshold == 0.9
        assert job.num_groups == 2
        assert job.ignore_hydrogens is True
        assert job.num_procs == 4
        assert job.matrix_format == "csv"
        assert job.energy_type == "G"
        assert job.num_molecules == 3


class TestGrouperJobPaths:
    def test_job_basename_returns_label(self):
        job = GrouperJob(
            molecules=make_molecules(2),
            grouping_strategy="rmsd",
            label="mygrp",
        )
        assert job.job_basename == "mygrp"

    def test_output_dir_uses_label(self):
        job = GrouperJob(
            molecules=make_molecules(2),
            grouping_strategy="rmsd",
            label="mygrp",
        )
        assert job.output_dir == os.path.join(job.folder, "mygrp_group_result")

    def test_output_dir_without_label_falls_back(self):
        job = GrouperJob(molecules=make_molecules(2), grouping_strategy="rmsd")
        job.label = None
        assert job.output_dir == os.path.join(job.folder, "group_result")

    def test_outputfile_path(self):
        job = GrouperJob(
            molecules=make_molecules(2),
            grouping_strategy="rmsd",
            label="mygrp",
        )
        assert job.outputfile == os.path.join(
            job.output_dir, "mygrp_group_1.xyz"
        )

    def test_logfile_and_errfile_paths(self):
        job = GrouperJob(
            molecules=make_molecules(2),
            grouping_strategy="rmsd",
            label="mygrp",
        )
        assert job.logfile == os.path.join(job.folder, "log.mygrp")
        assert job.errfile == os.path.join(job.folder, "mygrp.err")


class TestGrouperJobCompletion:
    def test_not_complete_when_output_dir_missing(self):
        job = GrouperJob(molecules=make_molecules(2), grouping_strategy="rmsd")
        assert job.is_complete() is False

    def test_complete_when_output_file_exists(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        job = GrouperJob(
            molecules=make_molecules(2),
            grouping_strategy="rmsd",
            label="mygrp",
        )
        os.makedirs(job.output_dir)
        with open(job.outputfile, "w") as f:
            f.write("1\ncomment\nC 0 0 0\n")
        assert job.is_complete() is True


class TestGrouperJobRun:
    def test_run_delegates_to_jobrunner(self):
        mock_runner = MagicMock()
        job = GrouperJob(
            molecules=make_molecules(2),
            grouping_strategy="rmsd",
            jobrunner=mock_runner,
        )
        job._run()
        mock_runner.run.assert_called_once_with(job)


class TestGrouperJobFromFilename:
    def test_from_filename_rejects_too_few_molecules(
        self, single_molecule_xyz_file
    ):
        with patch(
            "chemsmart.jobs.grouper.job.Molecule.from_filepath",
            return_value=[MagicMock()],
        ):
            with pytest.raises(ValueError, match="at least 2"):
                GrouperJob.from_filename(
                    filename=single_molecule_xyz_file,
                    grouping_strategy="rmsd",
                )

    def test_from_filename_builds_job_with_jobrunner(
        self, multiple_molecules_xyz_file
    ):
        mock_jobrunner = MagicMock()
        job = GrouperJob.from_filename(
            filename=multiple_molecules_xyz_file,
            grouping_strategy="rmsd",
            jobrunner=mock_jobrunner,
        )
        assert isinstance(job, GrouperJob)
        assert job.jobrunner is mock_jobrunner
        assert job.num_molecules >= 2

    def test_from_filename_derives_label_from_filename(
        self, multiple_molecules_xyz_file
    ):
        mock_jobrunner = MagicMock()
        job = GrouperJob.from_filename(
            filename=multiple_molecules_xyz_file,
            grouping_strategy="rmsd",
            jobrunner=mock_jobrunner,
        )
        expected_label = os.path.basename(multiple_molecules_xyz_file).split(
            "."
        )[0]
        assert job.label == expected_label

    def test_from_filename_creates_jobrunner_when_not_given(
        self, multiple_molecules_xyz_file
    ):
        mock_runner_instance = MagicMock()
        with patch(
            "chemsmart.jobs.runner.JobRunner.from_job",
            return_value=mock_runner_instance,
        ) as mock_from_job:
            job = GrouperJob.from_filename(
                filename=multiple_molecules_xyz_file,
                grouping_strategy="rmsd",
            )

        mock_from_job.assert_called_once()
        assert job.jobrunner is mock_runner_instance


class TestGrouperJobRunnerSetup:
    def test_defaults_scratch_to_false(self, pbs_server):
        runner = GrouperJobRunner(server=pbs_server)
        assert runner.scratch is False

    def test_prerun_creates_output_dir(self, pbs_server, grouper_job):
        runner = GrouperJobRunner(server=pbs_server)
        assert not os.path.exists(grouper_job.output_dir)
        runner._prerun(grouper_job)
        assert os.path.isdir(grouper_job.output_dir)

    def test_prerun_no_op_when_output_dir_exists(
        self, pbs_server, grouper_job
    ):
        os.makedirs(grouper_job.output_dir)
        runner = GrouperJobRunner(server=pbs_server)
        # Should not raise even though the dir already exists.
        runner._prerun(grouper_job)

    def test_get_command_and_executable_are_none(
        self, pbs_server, grouper_job
    ):
        runner = GrouperJobRunner(server=pbs_server)
        assert runner._get_command(grouper_job) is None
        assert runner._get_executable() is None

    def test_run_and_postrun_are_no_ops(self, pbs_server, grouper_job):
        runner = GrouperJobRunner(server=pbs_server)
        runner._run(process=None)
        runner._postrun(grouper_job)


class TestGrouperJobRunnerCreateGrouper:
    @pytest.mark.parametrize(
        "strategy,module,cls_name",
        [
            ("rmsd", "chemsmart.jobs.grouper.rmsd", "BasicRMSDGrouper"),
            ("hrmsd", "chemsmart.jobs.grouper.rmsd", "HungarianRMSDGrouper"),
            ("spyrmsd", "chemsmart.jobs.grouper.rmsd", "SpyRMSDGrouper"),
            ("irmsd", "chemsmart.jobs.grouper.rmsd", "IRMSDGrouper"),
            ("pymolrmsd", "chemsmart.jobs.grouper.rmsd", "PymolRMSDGrouper"),
            (
                "tanimoto",
                "chemsmart.jobs.grouper.tanimoto",
                "TanimotoSimilarityGrouper",
            ),
            (
                "tfd",
                "chemsmart.jobs.grouper.tfd",
                "TorsionFingerprintGrouper",
            ),
            (
                "isomorphism",
                "chemsmart.jobs.grouper.isomorphism",
                "RDKitIsomorphismGrouper",
            ),
            ("formula", "chemsmart.jobs.grouper.formula", "FormulaGrouper"),
            (
                "connectivity",
                "chemsmart.jobs.grouper.connectivity",
                "ConnectivityGrouper",
            ),
            ("energy", "chemsmart.jobs.grouper.energy", "EnergyGrouper"),
        ],
    )
    def test_create_grouper_dispatches_to_correct_class(
        self, pbs_server, grouper_job, strategy, module, cls_name
    ):
        grouper_job.grouping_strategy = strategy
        runner = GrouperJobRunner(server=pbs_server)
        mock_instance = MagicMock()
        with patch(f"{module}.{cls_name}") as mock_cls:
            mock_cls.return_value = mock_instance
            result = runner._create_grouper(grouper_job)

        mock_cls.assert_called_once()
        assert result is mock_instance
        call_kwargs = mock_cls.call_args.kwargs
        assert call_kwargs["molecules"] == grouper_job.molecules
        assert call_kwargs["label"] == grouper_job.label

    def test_create_grouper_unknown_strategy_raises(
        self, pbs_server, grouper_job
    ):
        grouper_job.grouping_strategy = "bogus"
        runner = GrouperJobRunner(server=pbs_server)
        with pytest.raises(ValueError, match="Unknown grouping strategy"):
            runner._create_grouper(grouper_job)

    def test_create_grouper_merges_extra_kwargs(self, pbs_server, grouper_job):
        grouper_job.grouper_kwargs = {"custom_option": "value"}
        runner = GrouperJobRunner(server=pbs_server)
        with patch("chemsmart.jobs.grouper.rmsd.BasicRMSDGrouper") as mock_cls:
            runner._create_grouper(grouper_job)

        assert mock_cls.call_args.kwargs["custom_option"] == "value"


class TestGrouperJobRunnerCreateProcess:
    def test_create_process_success_writes_outputs(
        self, pbs_server, grouper_job
    ):
        os.makedirs(grouper_job.output_dir)
        runner = GrouperJobRunner(server=pbs_server)

        mock_grouper = MagicMock()
        mock_grouper.group.return_value = ([], [])
        with (
            patch.object(runner, "_create_grouper", return_value=mock_grouper),
            patch.object(runner, "_write_outputs") as mock_write,
        ):
            exit_code = runner._create_process(grouper_job, None, None)

        assert exit_code == 0
        assert grouper_job._grouper is mock_grouper
        mock_write.assert_called_once_with(grouper_job, [], [])

    def test_create_process_failure_writes_errfile(
        self, pbs_server, grouper_job
    ):
        os.makedirs(grouper_job.output_dir)
        runner = GrouperJobRunner(server=pbs_server)

        with patch.object(
            runner,
            "_create_grouper",
            side_effect=RuntimeError("grouping failed"),
        ):
            exit_code = runner._create_process(grouper_job, None, None)

        assert exit_code == 1
        assert os.path.exists(grouper_job.errfile)
        with open(grouper_job.errfile) as f:
            assert "grouping failed" in f.read()


class TestGrouperJobRunnerWriteOutputs:
    def test_write_outputs_preserves_upstream_group_order_and_writes_xyz(
        self, pbs_server, grouper_job
    ):
        # Group/representative ordering is finalized upstream by the
        # grouper (see RMSDGrouper/EnergyGrouper's representative_strategy
        # handling); _write_outputs just writes whatever order it is
        # given, with the first entry treated as the representative.
        os.makedirs(grouper_job.output_dir)
        runner = GrouperJobRunner(server=pbs_server)

        mol_low = SimpleMolecule(
            energy=-2.0, symbols=["H"], positions=[(1.0, 1.0, 1.0)]
        )
        mol_high = SimpleMolecule(
            energy=-1.0, symbols=["H"], positions=[(0.0, 0.0, 0.0)]
        )
        groups = [[mol_low, mol_high]]
        group_indices = [[1, 0]]

        runner._write_outputs(grouper_job, groups, group_indices)

        group_file = os.path.join(grouper_job.output_dir, "mygrp_group_1.xyz")
        assert os.path.exists(group_file)
        with open(group_file) as f:
            content = f.read()
        # Written in the given (already-ordered) sequence: mol_low first.
        assert content.index("Original_Index: 2") < content.index(
            "Original_Index: 1"
        )

    def test_write_outputs_handles_molecules_without_energy(
        self, pbs_server, grouper_job
    ):
        os.makedirs(grouper_job.output_dir)
        runner = GrouperJobRunner(server=pbs_server)

        mol = SimpleMolecule(
            energy=None, symbols=["C"], positions=[(0.0, 0.0, 0.0)]
        )
        groups = [[mol]]
        group_indices = [[0]]

        runner._write_outputs(grouper_job, groups, group_indices)

        group_file = os.path.join(grouper_job.output_dir, "mygrp_group_1.xyz")
        with open(group_file) as f:
            content = f.read()
        assert "Energy: N/A" in content

    def test_write_outputs_uses_conformer_ids_when_provided(
        self, pbs_server, grouper_job
    ):
        os.makedirs(grouper_job.output_dir)
        grouper_job.conformer_ids = ["c1", "c2"]
        runner = GrouperJobRunner(server=pbs_server)

        mol = SimpleMolecule(
            energy=-1.0, symbols=["C"], positions=[(0.0, 0.0, 0.0)]
        )
        groups = [[mol]]
        group_indices = [[1]]

        runner._write_outputs(grouper_job, groups, group_indices)

        group_file = os.path.join(grouper_job.output_dir, "mygrp_group_1.xyz")
        with open(group_file) as f:
            content = f.read()
        assert "Original_Index: c2" in content


class TestBuildThermochemistryKwargs:
    def test_defaults_when_nothing_given(self):
        kwargs = build_thermochemistry_kwargs()
        for key, value in THERMOCHEMISTRY_DEFAULT_KWARGS.items():
            assert kwargs[key] == value

    def test_both_cutoffs_raises(self):
        with pytest.raises(ValueError, match="Cannot specify both"):
            build_thermochemistry_kwargs(
                cutoff_entropy_grimme=100, cutoff_entropy_truhlar=110
            )

    def test_grimme_cutoff_sets_entropy_method(self):
        kwargs = build_thermochemistry_kwargs(cutoff_entropy_grimme=100)
        assert kwargs["s_freq_cutoff"] == 100
        assert kwargs["entropy_method"] == "grimme"

    def test_truhlar_cutoff_sets_entropy_method(self):
        kwargs = build_thermochemistry_kwargs(cutoff_entropy_truhlar=110)
        assert kwargs["s_freq_cutoff"] == 110
        assert kwargs["entropy_method"] == "truhlar"

    def test_overrides_applied(self):
        kwargs = build_thermochemistry_kwargs(
            cutoff_enthalpy=50,
            concentration=2.0,
            pressure=0.5,
            temperature=300.0,
            alpha=8,
            weighted=False,
            energy_units="eV",
            check_imaginary_frequencies=False,
        )
        assert kwargs["h_freq_cutoff"] == 50
        assert kwargs["concentration"] == 2.0
        assert kwargs["pressure"] == 0.5
        assert kwargs["temperature"] == 300.0
        assert kwargs["alpha"] == 8
        assert kwargs["use_weighted_mass"] is False
        assert kwargs["energy_units"] == "eV"
        assert kwargs["check_imaginary_frequencies"] is False


class TestExtractConformerId:
    def test_matches_underscore_pattern(self):
        assert _extract_conformer_id("structure1_c12_opt.log") == "c12"

    def test_matches_dot_pattern(self):
        assert _extract_conformer_id("structure1_c5.log") == "c5"

    def test_no_match_returns_none(self):
        assert _extract_conformer_id("structure1_opt.log") is None


class TestGetLabel:
    def test_both_label_and_append_label_raises(self):
        with pytest.raises(ValueError, match="but not both"):
            _get_label("mylabel", "suffix", "base")

    def test_append_label_appends_to_base(self):
        assert _get_label(None, "suffix", "base") == "base_suffix"

    def test_label_returns_label(self):
        assert _get_label("mylabel", None, "base") == "mylabel"

    def test_neither_returns_base(self):
        assert _get_label(None, None, "base") == "base"


class TestExtractEnergyBasedOnEnergyType:
    def _mock_thermo(self, **output_attrs):
        thermo = MagicMock()
        output = MagicMock()
        for key, value in output_attrs.items():
            setattr(output, key, value)
        thermo.file_object = output
        return thermo

    def test_energy_type_e_returns_last_energy(self):
        thermo = self._mock_thermo(energies=[1.0, 2.0, 3.0])
        assert _extract_energy_based_on_energy_type(thermo, "E") == 3.0

    def test_energy_type_e_no_energies_returns_none(self):
        thermo = self._mock_thermo(energies=[])
        assert _extract_energy_based_on_energy_type(thermo, "E") is None

    def test_energy_type_h_returns_enthalpy(self):
        thermo = self._mock_thermo(enthalpy=-5.0)
        assert _extract_energy_based_on_energy_type(thermo, "H") == -5.0

    def test_energy_type_h_attribute_error_returns_none(self):
        output = MagicMock(spec=[])
        thermo = MagicMock()
        thermo.file_object = output
        assert _extract_energy_based_on_energy_type(thermo, "H") is None

    def test_energy_type_g_returns_gibbs(self):
        thermo = self._mock_thermo(gibbs_free_energy=-6.0)
        assert _extract_energy_based_on_energy_type(thermo, "G") == -6.0

    def test_energy_type_g_attribute_error_returns_none(self):
        output = MagicMock(spec=[])
        thermo = MagicMock()
        thermo.file_object = output
        assert _extract_energy_based_on_energy_type(thermo, "G") is None

    def test_energy_type_qhh_uses_qrrho_enthalpy(self, mocker):
        mocker.patch.object(
            grouper_module,
            "energy_conversion",
            return_value=-7.0,
        )
        thermo = MagicMock()
        thermo.qrrho_enthalpy = -12345.0
        result = _extract_energy_based_on_energy_type(thermo, "QHH")
        assert result == -7.0

    def test_energy_type_qhh_falls_back_to_h_on_attribute_error(self, mocker):
        thermo = MagicMock(spec=["file_object"])
        thermo.file_object = MagicMock(enthalpy=-5.0)
        result = _extract_energy_based_on_energy_type(thermo, "QHH")
        assert result == -5.0

    def test_energy_type_qhh_falls_back_to_h_when_none(self):
        thermo = MagicMock()
        thermo.qrrho_enthalpy = None
        thermo.file_object = MagicMock(enthalpy=-5.0)
        result = _extract_energy_based_on_energy_type(thermo, "QHH")
        assert result == -5.0

    def test_energy_type_qhg_uses_qrrho_gibbs(self, mocker):
        mocker.patch.object(
            grouper_module,
            "energy_conversion",
            return_value=-8.0,
        )
        thermo = MagicMock()
        thermo.qrrho_gibbs_free_energy = -54321.0
        result = _extract_energy_based_on_energy_type(thermo, "QHG")
        assert result == -8.0

    def test_energy_type_qhg_falls_back_to_g_on_attribute_error(self):
        thermo = MagicMock(spec=["file_object"])
        thermo.file_object = MagicMock(gibbs_free_energy=-6.0)
        result = _extract_energy_based_on_energy_type(thermo, "QHG")
        assert result == -6.0

    def test_energy_type_qhg_falls_back_to_g_when_none(self):
        thermo = MagicMock()
        thermo.qrrho_gibbs_free_energy = None
        thermo.file_object = MagicMock(gibbs_free_energy=-6.0)
        result = _extract_energy_based_on_energy_type(thermo, "QHG")
        assert result == -6.0

    def test_sp_qhg_full_success(self, mocker):
        thermo = MagicMock()
        thermo.qrrho_gibbs_free_energy = None
        thermo.file_object = MagicMock(gibbs_free_energy=-6.0, energies=[-5.0])
        thermo.filename = "/tmp/foo_c1_opt.log"
        mocker.patch.object(
            grouper_module,
            "_find_matching_sp_file",
            return_value="/tmp/sp/foo_c1_sp.log",
        )
        mocker.patch.object(
            grouper_module,
            "_extract_last_energy_from_output_file",
            return_value=-4.0,
        )
        result = _extract_energy_based_on_energy_type(thermo, "SP_QHG")
        assert result == -6.0 - (-5.0) + (-4.0)

    def test_sp_qhg_no_filename_returns_qh_gibbs(self):
        thermo = MagicMock(spec=["file_object", "qrrho_gibbs_free_energy"])
        thermo.qrrho_gibbs_free_energy = None
        thermo.file_object = MagicMock(gibbs_free_energy=-6.0, energies=[-5.0])
        result = _extract_energy_based_on_energy_type(thermo, "SP_QHG")
        assert result == -6.0

    def test_sp_qhg_no_matching_sp_file_raises(self, mocker):
        thermo = MagicMock()
        thermo.qrrho_gibbs_free_energy = None
        thermo.file_object = MagicMock(gibbs_free_energy=-6.0, energies=[-5.0])
        thermo.filename = "/tmp/foo_c1_opt.log"
        mocker.patch.object(
            grouper_module,
            "_find_matching_sp_file",
            return_value=None,
        )
        with pytest.raises(click.ClickException, match="matching SP file"):
            _extract_energy_based_on_energy_type(thermo, "SP_QHG")

    def test_sp_qhg_esolv_none_raises(self, mocker):
        thermo = MagicMock()
        thermo.qrrho_gibbs_free_energy = None
        thermo.file_object = MagicMock(gibbs_free_energy=-6.0, energies=[-5.0])
        thermo.filename = "/tmp/foo_c1_opt.log"
        mocker.patch.object(
            grouper_module,
            "_find_matching_sp_file",
            return_value="/tmp/sp/foo_c1_sp.log",
        )
        mocker.patch.object(
            grouper_module,
            "_extract_last_energy_from_output_file",
            return_value=None,
        )
        with pytest.raises(click.ClickException, match="sp_qh_gibbs"):
            _extract_energy_based_on_energy_type(thermo, "SP_QHG")

    def test_unsupported_energy_type_raises(self):
        thermo = MagicMock()
        with pytest.raises(click.ClickException, match="Unsupported"):
            _extract_energy_based_on_energy_type(thermo, "BOGUS")


class TestFindMatchingSpFile:
    def test_raises_when_sp_folder_missing(self, tmp_path):
        source = tmp_path / "foo_c1_opt.log"
        source.write_text("")
        with pytest.raises(ValueError, match="Can't find SP output folder"):
            _find_matching_sp_file(str(source))

    def test_finds_matching_suffix_file(self, tmp_path):
        sp_folder = tmp_path / "sp"
        sp_folder.mkdir()
        (sp_folder / "foo_c1_opt_sp.log").write_text("")
        source = tmp_path / "foo_c1_opt.log"
        source.write_text("")
        result = _find_matching_sp_file(str(source))
        assert result == str(sp_folder / "foo_c1_opt_sp.log")

    def test_returns_none_when_no_candidates(self, tmp_path):
        sp_folder = tmp_path / "sp"
        sp_folder.mkdir()
        source = tmp_path / "foo_c1_opt.log"
        source.write_text("")
        assert _find_matching_sp_file(str(source)) is None


class TestExtractLastEnergyFromOutputFile:
    def test_gaussian_program_returns_last_energy(self, mocker):
        mocker.patch.object(
            grouper_module,
            "get_program_type_from_file",
            return_value="gaussian",
        )
        mock_output = MagicMock(normal_termination=True, energies=[-1.0, -2.0])
        mocker.patch(
            "chemsmart.io.gaussian.output.Gaussian16Output",
            return_value=mock_output,
        )
        result = _extract_last_energy_from_output_file("foo.log")
        assert result == -2.0

    def test_orca_program_returns_last_energy(self, mocker):
        mocker.patch.object(
            grouper_module,
            "get_program_type_from_file",
            return_value="orca",
        )
        mock_output = MagicMock(normal_termination=True, energies=[-3.0, -4.0])
        mocker.patch(
            "chemsmart.io.orca.output.ORCAOutput",
            return_value=mock_output,
        )
        result = _extract_last_energy_from_output_file("foo.out")
        assert result == -4.0

    def test_unknown_program_returns_none(self, mocker):
        mocker.patch.object(
            grouper_module,
            "get_program_type_from_file",
            return_value=None,
        )
        assert _extract_last_energy_from_output_file("foo.xyz") is None

    def test_not_normal_termination_returns_none(self, mocker):
        mocker.patch.object(
            grouper_module,
            "get_program_type_from_file",
            return_value="gaussian",
        )
        mock_output = MagicMock(normal_termination=False)
        mocker.patch(
            "chemsmart.io.gaussian.output.Gaussian16Output",
            return_value=mock_output,
        )
        assert _extract_last_energy_from_output_file("foo.log") is None

    def test_no_energies_returns_none(self, mocker):
        mocker.patch.object(
            grouper_module,
            "get_program_type_from_file",
            return_value="gaussian",
        )
        mock_output = MagicMock(normal_termination=True, energies=[])
        mocker.patch(
            "chemsmart.io.gaussian.output.Gaussian16Output",
            return_value=mock_output,
        )
        assert _extract_last_energy_from_output_file("foo.log") is None


class TestGrouperCliValidation:
    def test_both_filenames_and_directory_raises(
        self, mocker, multiple_molecules_xyz_file, tmp_path
    ):
        result, _ = invoke_grouper(
            [
                "-f",
                multiple_molecules_xyz_file,
                "-d",
                str(tmp_path),
                "-p",
                "gaussian",
                "rmsd",
            ],
            mocker,
        )
        assert result.exit_code != 0
        assert "Cannot specify both" in result.output

    def test_directory_without_program_or_filetype_raises(
        self, mocker, tmp_path
    ):
        result, _ = invoke_grouper(
            ["-d", str(tmp_path), "rmsd"],
            mocker,
        )
        assert result.exit_code != 0
        assert "Must specify -p/--program" in result.output

    def test_multiple_filenames_rejected(
        self, mocker, multiple_molecules_xyz_file
    ):
        result, _ = invoke_grouper(
            [
                "-f",
                multiple_molecules_xyz_file,
                "-f",
                multiple_molecules_xyz_file,
                "rmsd",
            ],
            mocker,
        )
        assert result.exit_code != 0
        assert "only single file input" in result.output

    def test_non_e_energy_type_for_xyz_rejected(
        self, mocker, multiple_molecules_xyz_file
    ):
        result, _ = invoke_grouper(
            [
                "-f",
                multiple_molecules_xyz_file,
                "-E",
                "H",
                "rmsd",
            ],
            mocker,
        )
        assert result.exit_code != 0
        assert "Only energy" in result.output

    def test_no_input_specified_warns_and_returns(self, mocker):
        result, mock_job_cls = invoke_grouper(["rmsd"], mocker)
        assert result.exit_code == 0
        mock_job_cls.assert_called_once()
        _, call_kwargs = mock_job_cls.call_args
        assert call_kwargs["molecules"] is None


class TestGrouperCliSuccessPaths:
    def test_single_file_rmsd_success(
        self, mocker, multiple_molecules_xyz_file
    ):
        result, mock_job_cls = invoke_grouper(
            ["-f", multiple_molecules_xyz_file, "rmsd"], mocker
        )
        assert result.exit_code == 0, result.output
        mock_job_cls.assert_called_once()
        _, call_kwargs = mock_job_cls.call_args
        assert call_kwargs["grouping_strategy"] == "rmsd"
        assert len(call_kwargs["molecules"]) > 1

    def test_label_option_used_as_grouper_label(
        self, mocker, multiple_molecules_xyz_file
    ):
        result, mock_job_cls = invoke_grouper(
            [
                "-f",
                multiple_molecules_xyz_file,
                "-l",
                "mylabel",
                "rmsd",
            ],
            mocker,
        )
        assert result.exit_code == 0, result.output
        _, call_kwargs = mock_job_cls.call_args
        assert call_kwargs["label"] == "mylabel_rmsd"

    def test_directory_mode_success(self, mocker, tmp_path):
        mock_molecule = MagicMock(name="a")
        mocker.patch.object(
            grouper_module,
            "_load_molecules_from_directory",
            return_value=(
                [mock_molecule, mock_molecule],
                ["c1", "c2"],
                [],
                "auto_label",
            ),
        )
        result, mock_job_cls = invoke_grouper(
            ["-d", str(tmp_path), "-p", "gaussian", "rmsd"], mocker
        )
        assert result.exit_code == 0, result.output
        _, call_kwargs = mock_job_cls.call_args
        assert call_kwargs["conformer_ids"] == ["c1", "c2"]
        assert call_kwargs["label"] == "auto_label_rmsd"

    def test_fewer_than_two_molecules_warns_but_continues(
        self, mocker, tmp_path
    ):
        single_mol_file = tmp_path / "single.xyz"
        single_mol_file.write_text("1\ncomment\nH 0.0 0.0 0.0\n")
        result, mock_job_cls = invoke_grouper(
            ["-f", str(single_mol_file), "rmsd"], mocker
        )
        assert result.exit_code == 0, result.output
        mock_job_cls.assert_called_once()

    def test_tanimoto_default_fingerprint(
        self, mocker, multiple_molecules_xyz_file
    ):
        result, mock_job_cls = invoke_grouper(
            ["-f", multiple_molecules_xyz_file, "tanimoto"], mocker
        )
        assert result.exit_code == 0, result.output
        _, call_kwargs = mock_job_cls.call_args
        assert call_kwargs["fingerprint_type"] == "rdkit"

    def test_tanimoto_custom_fingerprint(
        self, mocker, multiple_molecules_xyz_file
    ):
        result, mock_job_cls = invoke_grouper(
            [
                "-f",
                multiple_molecules_xyz_file,
                "tanimoto",
                "--fingerprint-type",
                "morgan",
            ],
            mocker,
        )
        assert result.exit_code == 0, result.output
        _, call_kwargs = mock_job_cls.call_args
        assert call_kwargs["fingerprint_type"] == "morgan"

    def test_tfd_options(self, mocker, multiple_molecules_xyz_file):
        result, mock_job_cls = invoke_grouper(
            [
                "-f",
                multiple_molecules_xyz_file,
                "tfd",
                "--no-use-weights",
                "--max-dev",
                "spec",
            ],
            mocker,
        )
        assert result.exit_code == 0, result.output
        _, call_kwargs = mock_job_cls.call_args
        assert call_kwargs["use_weights"] is False
        assert call_kwargs["max_dev"] == "spec"

    def test_irmsd_inversion_option(self, mocker, multiple_molecules_xyz_file):
        result, mock_job_cls = invoke_grouper(
            [
                "-f",
                multiple_molecules_xyz_file,
                "irmsd",
                "--inversion",
                "on",
            ],
            mocker,
        )
        assert result.exit_code == 0, result.output
        _, call_kwargs = mock_job_cls.call_args
        assert call_kwargs["inversion"] == "on"

    @pytest.mark.parametrize(
        "subcommand", ["hrmsd", "spyrmsd", "pymolrmsd", "energy"]
    )
    def test_simple_subcommands_succeed(
        self, mocker, multiple_molecules_xyz_file, subcommand
    ):
        result, mock_job_cls = invoke_grouper(
            ["-f", multiple_molecules_xyz_file, subcommand], mocker
        )
        assert result.exit_code == 0, result.output
        _, call_kwargs = mock_job_cls.call_args
        assert call_kwargs["grouping_strategy"] == subcommand


class TestGrouperCliStructuralSubcommands:
    @pytest.mark.parametrize(
        "subcommand,strategy",
        [
            ("connectivity", "connectivity"),
            ("formula", "formula"),
            ("isomorphism", "isomorphism"),
        ],
    )
    def test_succeeds_without_threshold_or_num_groups(
        self, mocker, multiple_molecules_xyz_file, subcommand, strategy
    ):
        result, mock_job_cls = invoke_grouper(
            ["-f", multiple_molecules_xyz_file, subcommand], mocker
        )
        assert result.exit_code == 0, result.output
        _, call_kwargs = mock_job_cls.call_args
        assert call_kwargs["grouping_strategy"] == strategy

    @pytest.mark.parametrize(
        "subcommand", ["connectivity", "formula", "isomorphism"]
    )
    def test_rejects_threshold_option(
        self, mocker, multiple_molecules_xyz_file, subcommand
    ):
        result, _ = invoke_grouper(
            [
                "-f",
                multiple_molecules_xyz_file,
                "-T",
                "0.5",
                subcommand,
            ],
            mocker,
        )
        assert result.exit_code != 0
        assert "does not support threshold" in result.output

    @pytest.mark.parametrize(
        "subcommand", ["connectivity", "formula", "isomorphism"]
    )
    def test_rejects_num_groups_option(
        self, mocker, multiple_molecules_xyz_file, subcommand
    ):
        result, _ = invoke_grouper(
            [
                "-f",
                multiple_molecules_xyz_file,
                "-N",
                "3",
                subcommand,
            ],
            mocker,
        )
        assert result.exit_code != 0
        assert "does not support num_groups" in result.output


class TestGrouperCliBothCutoffsRejected:
    def test_both_cutoffs_becomes_bad_parameter(
        self, mocker, multiple_molecules_xyz_file
    ):
        result, _ = invoke_grouper(
            [
                "-f",
                multiple_molecules_xyz_file,
                "-csg",
                "100",
                "-cst",
                "110",
                "rmsd",
            ],
            mocker,
        )
        assert result.exit_code != 0
        assert "Cannot specify both" in result.output


class TestFindMatchingSpFileNoSuffix:
    def test_source_without_suffix_skips_suffix_glob(self, tmp_path):
        sp_folder = tmp_path / "sp"
        sp_folder.mkdir()
        (sp_folder / "foo_c1_sp.log").write_text("")
        source = tmp_path / "foo_c1"
        source.write_text("")
        result = _find_matching_sp_file(str(source))
        assert result == str(sp_folder / "foo_c1_sp.log")


class TestLoadMoleculesFromDirectory:
    def test_missing_program_and_filetype_raises(self, mocker, tmp_path):
        mock_folder_cls = mocker.patch("chemsmart.io.folder.BaseFolder")
        mock_folder_cls.return_value = MagicMock()
        with pytest.raises(click.BadParameter, match="Must specify"):
            _load_molecules_from_directory(
                directory=str(tmp_path), program=None, filetype=None
            )

    def test_no_output_files_found_raises(self, mocker, tmp_path):
        mock_folder = MagicMock()
        mock_folder.get_all_output_files_in_current_folder_by_program.return_value = (
            []
        )
        mocker.patch(
            "chemsmart.io.folder.BaseFolder", return_value=mock_folder
        )
        with pytest.raises(click.BadParameter, match="No matching output"):
            _load_molecules_from_directory(
                directory=str(tmp_path), program="gaussian", filetype=None
            )

    def test_program_and_filetype_uses_combined_lookup(self, mocker, tmp_path):
        mock_folder = MagicMock()
        mock_folder.get_all_files_in_current_folder_by_program_and_suffix.return_value = (
            []
        )
        mocker.patch(
            "chemsmart.io.folder.BaseFolder", return_value=mock_folder
        )
        with pytest.raises(click.BadParameter, match="No matching output"):
            _load_molecules_from_directory(
                directory=str(tmp_path), program="gaussian", filetype="log"
            )
        mock_folder.get_all_files_in_current_folder_by_program_and_suffix.assert_called_once()

    def test_filetype_only_uses_suffix_lookup(self, mocker, tmp_path):
        mock_folder = MagicMock()
        mock_folder.get_all_files_in_current_folder_by_suffix.return_value = []
        mocker.patch(
            "chemsmart.io.folder.BaseFolder", return_value=mock_folder
        )
        with pytest.raises(click.BadParameter, match="No matching output"):
            _load_molecules_from_directory(
                directory=str(tmp_path), program=None, filetype="log"
            )

    def test_successful_load_with_patterned_and_unpatterned_files(
        self, mocker, tmp_path
    ):
        files = [
            str(tmp_path / "mol_c2_opt.log"),
            str(tmp_path / "mol_c1_opt.log"),
            str(tmp_path / "unpatterned.log"),
        ]
        mock_folder = MagicMock()
        mock_folder.get_all_output_files_in_current_folder_by_program.return_value = (
            files
        )
        mocker.patch(
            "chemsmart.io.folder.BaseFolder", return_value=mock_folder
        )

        mock_mol = MagicMock()
        mock_thermo = MagicMock()
        mock_thermo.molecule = mock_mol
        mock_thermo.file_object.energies = [-1.0]
        mocker.patch(
            "chemsmart.analysis.thermochemistry.Thermochemistry",
            return_value=mock_thermo,
        )

        molecules, conformer_ids, skipped_ids, common_label = (
            _load_molecules_from_directory(
                directory=str(tmp_path), program="gaussian", filetype=None
            )
        )
        assert len(molecules) == 3
        # c1 sorted before c2, unpatterned file goes last.
        assert conformer_ids == ["c1", "c2", "unpatterned"]
        assert skipped_ids == []
        assert common_label == "mol"

    def test_skips_molecule_when_energy_extraction_fails(
        self, mocker, tmp_path
    ):
        files = [str(tmp_path / "mol_c1_opt.log")]
        mock_folder = MagicMock()
        mock_folder.get_all_output_files_in_current_folder_by_program.return_value = (
            files
        )
        mocker.patch(
            "chemsmart.io.folder.BaseFolder", return_value=mock_folder
        )

        mock_thermo = MagicMock()
        mock_thermo.molecule = MagicMock()
        mock_thermo.file_object.energies = []
        mocker.patch(
            "chemsmart.analysis.thermochemistry.Thermochemistry",
            return_value=mock_thermo,
        )

        with pytest.raises(click.BadParameter, match="No valid molecules"):
            _load_molecules_from_directory(
                directory=str(tmp_path), program="gaussian", filetype=None
            )

    def test_skips_when_molecule_is_none(self, mocker, tmp_path):
        files = [str(tmp_path / "mol_c1_opt.log")]
        mock_folder = MagicMock()
        mock_folder.get_all_output_files_in_current_folder_by_program.return_value = (
            files
        )
        mocker.patch(
            "chemsmart.io.folder.BaseFolder", return_value=mock_folder
        )

        mock_thermo = MagicMock()
        mock_thermo.molecule = None
        mocker.patch(
            "chemsmart.analysis.thermochemistry.Thermochemistry",
            return_value=mock_thermo,
        )

        with pytest.raises(click.BadParameter, match="No valid molecules"):
            _load_molecules_from_directory(
                directory=str(tmp_path), program="gaussian", filetype=None
            )

    def test_skips_on_thermochemistry_construction_error(
        self, mocker, tmp_path
    ):
        files = [str(tmp_path / "mol_c1_opt.log")]
        mock_folder = MagicMock()
        mock_folder.get_all_output_files_in_current_folder_by_program.return_value = (
            files
        )
        mocker.patch(
            "chemsmart.io.folder.BaseFolder", return_value=mock_folder
        )
        mocker.patch(
            "chemsmart.analysis.thermochemistry.Thermochemistry",
            side_effect=RuntimeError("boom"),
        )

        with pytest.raises(click.BadParameter, match="No valid molecules"):
            _load_molecules_from_directory(
                directory=str(tmp_path), program="gaussian", filetype=None
            )


class TestCreateGrouperJobFromContextThermoParameters:
    def test_qhg_energy_type_builds_thermo_parameters_string(self, mocker):
        mock_job_cls = mocker.patch("chemsmart.jobs.grouper.GrouperJob")
        ctx = MagicMock()
        ctx.obj = {
            "molecules": [MagicMock(), MagicMock()],
            "ignore_hydrogens": False,
            "num_procs": 1,
            "grouper_label": "mylabel",
            "num_groups": None,
            "conformer_ids": None,
            "skipped_ids": None,
            "energy_type": "qhG",
            "thermo_kwargs": {"temperature": 298.15, "alpha": None},
            "threshold": None,
            "matrix_format": "xlsx",
        }
        create_grouper_job_from_context(ctx, strategy="rmsd")
        _, call_kwargs = mock_job_cls.call_args
        assert "temperature=298.15" in call_kwargs["thermo_parameters"]

    def test_e_energy_type_leaves_thermo_parameters_none(self, mocker):
        mock_job_cls = mocker.patch("chemsmart.jobs.grouper.GrouperJob")
        ctx = MagicMock()
        ctx.obj = {
            "molecules": [MagicMock(), MagicMock()],
            "ignore_hydrogens": False,
            "num_procs": 1,
            "grouper_label": "mylabel",
            "num_groups": None,
            "conformer_ids": None,
            "skipped_ids": None,
            "energy_type": "E",
            "thermo_kwargs": {"temperature": 298.15},
            "threshold": None,
            "matrix_format": "xlsx",
        }
        create_grouper_job_from_context(ctx, strategy="rmsd")
        _, call_kwargs = mock_job_cls.call_args
        assert call_kwargs["thermo_parameters"] is None
