from __future__ import annotations

from chemsmart.agent.harness.extractors import (
    extract_cartesian_state,
    extract_gaussian_route,
    extract_orca_route,
)


def test_extract_gaussian_route_line():
    content = "%chk=ts.chk\n# b3lyp/6-31g* opt=(ts,calcfc,noeigentest)\n\n"

    assert (
        extract_gaussian_route(content)
        == "# b3lyp/6-31g* opt=(ts,calcfc,noeigentest)"
    )


def test_extract_gaussian_route_includes_continuation_lines():
    content = "%chk=ts.chk\n# b3lyp/6-31g*\nopt=(ts,calcfc,noeigentest)\n\n"

    assert (
        extract_gaussian_route(content)
        == "# b3lyp/6-31g* opt=(ts,calcfc,noeigentest)"
    )


def test_extract_gaussian_route_returns_none_for_empty_or_missing_route():
    assert extract_gaussian_route("") is None
    assert extract_gaussian_route("title\n\n0 1\n") is None
    assert extract_gaussian_route(None) is None


def test_extract_orca_route_skeleton():
    content = "! b3lyp def2-svp opt\n%pal nprocs 4 end\n! moreprint\n"

    assert extract_orca_route(content) == "! b3lyp def2-svp opt ! moreprint"


def test_extract_orca_cartesian_state():
    content = "! opt\n* xyz 0 2\nFe 0 0 0\nN 0 0 1.7\n*\n"

    assert extract_cartesian_state(content, software="orca") == {
        "charge": 0,
        "multiplicity": 2,
        "element_symbols": ["Fe", "N"],
    }


def test_extract_gaussian_cartesian_state():
    content = "# opt b3lyp/6-31g(d)\n\ntitle\n\n-1 1\nCl 0 0 0\nC 2 0 0\n\n"

    assert extract_cartesian_state(content, software="gaussian") == {
        "charge": -1,
        "multiplicity": 1,
        "element_symbols": ["Cl", "C"],
        "charge_multiplicity_pairs": [[-1, 1]],
    }


def test_extract_gaussian_oniom_state_and_atom_layers():
    content = """# opt oniom(b3lyp/6-31g(d):uff)

enzyme

-1 2 0 1 0 1
Fe 0.0 0.0 0.0 H
N  0.0 0.0 1.8 H
C  0.0 0.0 3.0 L H 2
H  0.0 1.0 3.0 L

"""

    assert extract_cartesian_state(content, software="gaussian") == {
        "charge": -1,
        "multiplicity": 2,
        "element_symbols": ["Fe", "N", "C", "H"],
        "charge_multiplicity_pairs": [[-1, 2], [0, 1], [0, 1]],
        "atom_layers": ["H", "H", "L", "L"],
        "layer_atoms": {"H": [1, 2], "M": [], "L": [3, 4]},
    }
