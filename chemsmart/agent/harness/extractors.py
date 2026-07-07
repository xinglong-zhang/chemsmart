from __future__ import annotations


def extract_gaussian_route(content: str | None) -> str | None:
    """Return the first Gaussian route section, including continuation lines."""
    if not isinstance(content, str):
        return None
    lines = content.splitlines()
    for index, line in enumerate(lines):
        if not line.lstrip().startswith("#"):
            continue
        route_lines = [line.strip()]
        for continuation in lines[index + 1 :]:
            stripped = continuation.strip()
            if not stripped:
                break
            route_lines.append(stripped)
        return " ".join(route_lines)
    return None


def extract_orca_route(content: str | None) -> str | None:
    """Return ORCA route lines. Phase 1 does not attach hard gates to them."""
    if not isinstance(content, str):
        return None
    route_lines = [
        line.strip()
        for line in content.splitlines()
        if line.lstrip().startswith("!")
    ]
    if not route_lines:
        return None
    return " ".join(route_lines)
