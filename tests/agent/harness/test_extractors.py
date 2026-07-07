from __future__ import annotations

from chemsmart.agent.harness.extractors import (
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
