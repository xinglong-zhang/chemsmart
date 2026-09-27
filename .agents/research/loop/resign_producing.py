"""Re-sign each archived goal on the code that produced it, where that code
still stands where the goal's job ran it.

    python .agents/research/loop/resign_producing.py OUT_DIR SPECS_FILE

``SPECS_FILE`` holds lines ``label|agent_dir|goal_id`` (the form
``resign.py`` reads). For each goal the job output beside it (the
``slurm-*.out`` a CUHK goal job writes) names the code it imported
("chemsmart from ...") and the tree digest it verified ("remote code tree
digest ..."). The digest of that code directory is recomputed now over the
file list the job ran (``code-files.at-run.txt`` beside the goal, else the
campaign's list), and only a directory whose digest still equals the
printed one is used: a code directory rebuilt in place after the goal ran
is not the producing code, and that goal is reported ``code_moved``.
Goals are grouped by code directory and ``resign.py`` runs once per group
with ``PYTHONPATH`` set to it; HOME is fenced by the caller.

A replay on the producing code that does not reproduce the archived word
means the harness is not faithful for that goal (R11 truth's F2): its
comparison at the pin is not evidence either way.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def goal_dir(agent: Path) -> Path:
    workspace = agent.parent
    return workspace.parent if workspace.name == "workspace" else workspace


def printed_code(folder: Path) -> tuple[str, str, str]:
    """(code dir, digest, commit) from the latest job output that says."""

    best: tuple[int, str, str, str] | None = None
    for out in folder.glob("slurm-*.out"):
        match = re.search(r"slurm-(\d+)\.out$", out.name)
        job = int(match.group(1)) if match else 0
        try:
            text = out.read_text(errors="replace")
        except OSError:
            continue
        code = re.search(r"chemsmart from (\S+)/chemsmart/__init__\.py", text)
        digest = re.search(r"code tree digest ([0-9a-f]{64})", text)
        commit = re.search(r"code commit ([0-9a-f]{40})", text)
        if code and digest and (best is None or job > best[0]):
            best = (job, code.group(1), digest.group(1), commit.group(1) if commit else "")
    if best is None:
        return "", "", ""
    return best[1], best[2], best[3]


def tree_digest(code: Path, listing: Path | None) -> str:
    if listing is not None and listing.is_file():
        names = sorted(l.strip() for l in listing.read_text().splitlines() if l.strip())
    else:
        names = []
        for folder, dirs, files in os.walk(code / "chemsmart"):
            dirs[:] = [d for d in dirs if d != "__pycache__"]
            for name in files:
                if name.endswith(".pyc") or name.startswith("._"):
                    continue
                names.append(str(Path(folder, name).relative_to(code)))
        names.append("pyproject.toml")
        names.sort()
    h = hashlib.sha256()
    for name in names:
        try:
            data = (code / name).read_bytes()
        except OSError:
            return ""
        h.update(name.encode() + b"\0" + hashlib.sha256(data).digest())
    return h.hexdigest()


def main() -> None:
    out_root = Path(sys.argv[1])
    out_root.mkdir(parents=True, exist_ok=True)
    groups: dict[str, list[str]] = {}
    status: list[dict] = []
    for line in Path(sys.argv[2]).read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        label, agent, goal_id = line.split("|")
        folder = goal_dir(Path(agent))
        code, printed, commit = printed_code(folder)
        row = {"label": label, "agent": agent, "goal_id": goal_id,
               "code": code, "printed_digest": printed, "commit": commit}
        if not code:
            row["status"] = "no_job_output"
        else:
            listing = folder / "code-files.at-run.txt"
            if not listing.is_file():
                campaign = Path(code).parent
                for candidate in (campaign / f"{Path(code).name}.files.txt",
                                  campaign / "code-files.txt"):
                    if candidate.is_file():
                        listing = candidate
                        break
                else:
                    listing = None
            now = tree_digest(Path(code), listing)
            row["digest_now"] = now
            row["status"] = "in_place" if now == printed else "code_moved"
            if now == printed:
                groups.setdefault(code, []).append(line)
        status.append(row)
    (out_root / "producing-code.json").write_text(json.dumps(status, indent=1))
    counts: dict[str, int] = {}
    for row in status:
        counts[row["status"]] = counts.get(row["status"], 0) + 1
    print("producing code:", counts, flush=True)
    for index, (code, lines) in enumerate(sorted(groups.items())):
        spec = out_root / f"specs-{index:03d}.txt"
        spec.write_text("\n".join(lines) + "\n")
        env = dict(os.environ, PYTHONPATH=code)
        target = out_root / f"group-{index:03d}"
        print(f"== group {index} {code}: {len(lines)} goal(s)", flush=True)
        subprocess.run(
            [sys.executable, str(HERE / "resign.py"), str(target), f"@{spec}"],
            env=env,
            check=False,
        )


if __name__ == "__main__":
    main()
