"""Every geometry origin required the workspace to already have it.

A supplied file, a database record, a previous result, or a derivation,
composition or edit of one of those: that was the whole list. So a
session needing a molecule nobody handed it -- a reference couple
computed at its own level so the systematic cancels, a calibration
standard, a literature comparison -- had no route at all, while the
human CLI has carried ``-p/--pubchem`` and ``Molecule.from_pubchem`` for
the same programs all along.

OPEN-2's ino3-qwen (2026-09-07) named a same-level ferrocene/ferrocenium
pair as the route that would cancel the systematic it was reporting, and
declined it on cost. Its stream cannot say whether it knew it was
blocked; the route was blocked either way.

The hub invariant is untouched, and that is what these tests pin: the
model names an identifier and nothing else, the host owns the bytes, no
electronic state is inferred, and a lookup that fails is a typed refusal
rather than an empty molecule.
"""

from __future__ import annotations

import json

import pytest

from chemsmart.agent._contracts import ContractError, RoutedContractError
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.runtime.events import EventKind
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1
from chemsmart.io.molecules.structure import Molecule

pytestmark = pytest.mark.capability("tool:fetch_pubchem_geometry")

_TASK = "a" * 64


def _host(tmp_path):
    workspace = tmp_path / "workspace"
    workspace.mkdir(exist_ok=True)
    return CommandCompiledToolHostV1(
        event_store=RuntimeEventStore(
            tmp_path / "events.jsonl", session_id="s1"
        ),
        task_spec_sha256s=(_TASK,),
        approved_workspace=workspace,
    )


def _water():
    return Molecule(
        symbols=["O", "H", "H"],
        positions=[
            [0.0000, 0.0000, 0.1173],
            [0.0000, 0.7572, -0.4692],
            [0.0000, -0.7572, -0.4692],
        ],
    )


def _served(monkeypatch, molecule):
    monkeypatch.setattr(
        Molecule,
        "from_pubchem",
        classmethod(lambda cls, identifier, return_list=False: molecule),
    )


def test_the_model_names_an_identifier_and_the_host_owns_the_bytes(
    tmp_path, monkeypatch
):
    host = _host(tmp_path)
    _served(monkeypatch, _water())

    result = host._fetch_pubchem_geometry(
        "t1", {"artifact_id": "ref-water", "identifier": "water"}
    )

    artifact = result["artifact"]
    receipt = result["pubchem_geometry"]
    # The host wrote the file, hashed it, and registered it: the model
    # supplied no coordinates and no path.
    assert artifact.kind == "geometry_xyz"
    assert host.artifacts["ref-water"] is artifact
    written = (
        tmp_path / "workspace" / "artifacts" / "ref-water.xyz"
    ).read_text()
    assert written.splitlines()[0] == "3"
    assert "water" in written.splitlines()[1]
    assert receipt.formula == "H2O"
    assert receipt.atom_count == 3
    assert receipt.fragment_count == 1
    assert receipt.identifier_kind == "name_or_smiles"
    assert host.pubchem_geometries[artifact.sha256] is receipt

    # It binds no electronic state, and says so where the model reads it.
    assert "bind charge and multiplicity explicitly" in result["next_action"]
    assert "not a relaxed structure" in result["next_action"]

    # A numeric identifier is a CID and is recorded as one.
    other = host._fetch_pubchem_geometry(
        "t1", {"artifact_id": "ref-cid", "identifier": "962"}
    )
    assert other["pubchem_geometry"].identifier_kind == "cid"


def test_the_fetch_is_an_event_with_its_lineage(tmp_path, monkeypatch):
    host = _host(tmp_path)
    _served(monkeypatch, _water())
    host._fetch_pubchem_geometry(
        "t1", {"artifact_id": "ref-water", "identifier": "water"}
    )
    kinds = [
        json.loads(line)["kind"]
        for line in (tmp_path / "events.jsonl").read_text().splitlines()
        if line.strip()
    ]
    assert EventKind.PUBCHEM_GEOMETRY_FETCHED.value in kinds


def test_the_fetched_geometry_binds_an_identity_like_any_other(
    tmp_path, monkeypatch
):
    """Reachability: the origin reaches the rest of the plane."""

    host = _host(tmp_path)
    _served(monkeypatch, _water())
    host._fetch_pubchem_geometry(
        "t1", {"artifact_id": "ref-water", "identifier": "water"}
    )
    bound = host._bind_scientific_identity(
        "t1",
        {
            "input_artifact_id": "ref-water",
            "task_spec_sha256": _TASK,
            "charge": 0,
            "multiplicity": 1,
        },
    )
    assert bound["geometry"]["formula"] == "H2O"
    # The depositor's symmetry travels with the record, and the host
    # states it here as it does for every other origin.
    assert bound["observations"]


@pytest.mark.parametrize(
    "raised,gate",
    (
        (None, "pubchem.record_exists"),
        (ValueError("no such compound"), "pubchem.record_exists"),
        (OSError("name resolution failed"), "pubchem.record_is_reachable"),
        (
            RuntimeError(
                'Cannot convert multi-component molecule: "10219726"'
            ),
            "pubchem.record_converts_to_one_geometry",
        ),
        (ImportError("rdkit"), "pubchem.converter_is_installed"),
    ),
)
def test_each_way_a_lookup_fails_gets_its_own_route(
    tmp_path, monkeypatch, raised, gate
):
    """A cause must never be reported as a consequence.

    The live probe that earned this repair found ferrocene failing on
    the 2D-to-3D conversion, because PubChem stores it as Fe(2+) beside
    two cyclopentadienide anions. Under one flat gate that read as
    "unreachable", which would send a session to retry a lookup that can
    never succeed.
    """

    host = _host(tmp_path)

    def _serve(cls, identifier, return_list=False):
        if raised is not None:
            raise raised
        return None

    monkeypatch.setattr(Molecule, "from_pubchem", classmethod(_serve))
    with pytest.raises(RoutedContractError) as refusal:
        host._fetch_pubchem_geometry(
            "t1", {"artifact_id": "a", "identifier": "whatever"}
        )
    report = refusal.value.failure_report
    assert report["gate"] == gate
    assert report["route"]
    # No refusal leaves a half-registered artifact behind.
    assert "a" not in host.artifacts


def test_an_empty_identifier_and_a_taken_id_are_refused(tmp_path, monkeypatch):
    host = _host(tmp_path)
    _served(monkeypatch, _water())
    host._fetch_pubchem_geometry(
        "t1", {"artifact_id": "ref-water", "identifier": "water"}
    )
    with pytest.raises(ContractError):
        host._fetch_pubchem_geometry(
            "t1", {"artifact_id": "ref-water", "identifier": "water"}
        )
    with pytest.raises(ContractError):
        host._fetch_pubchem_geometry(
            "t1", {"artifact_id": "fresh", "identifier": "   "}
        )
