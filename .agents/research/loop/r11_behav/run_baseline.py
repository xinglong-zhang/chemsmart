"""R11 behav baseline runner: the pre-registered cells, interleaved, one process.

Runs under one provider lease; every sample is one invocation of the
committed instrument with the tree named by PYTHONPATH. The order is fixed
here and written to the manifest before the first sample. Stop rule: three
consecutive samples that end with no real turn after a transport failure.

    python run_baseline.py [--from INDEX]
"""

import json
import os
import subprocess
import sys
import time
from pathlib import Path

S = Path(
    "/private/tmp/claude-501/-Users-hongjiseung-developer-chemsmart/"
    "c946a1a6-a959-400c-ba03-fd5f4544cb43/scratchpad/behav"
)
W = Path(
    "/Users/hongjiseung/developer/chemsmart/.claude/worktrees/"
    "agent-a9e062387a514384b"
)
PY = "/opt/anaconda3/bin/python"
TOOL = W / ".agents/research/loop/matched_turns.py"
OUT = S / "baseline"
TREES = {
    "control": (S / "tree-9297d6ba", "tree-9297d6ba"),
    "repaired": (S / "tree-6b5adf48", "tree-6b5adf48"),
    "head": (W, "agent-a9e062387a514384b"),
}
PROFILES = {
    "deepseek": S / "profile-deepseek.yaml",
    "qwen": S / "profile-qwen-65k.yaml",
}
N = 12
CELLS = (
    ("d1", "control", "deepseek"),
    ("d1", "repaired", "deepseek"),
    ("d1", "control", "qwen"),
    ("d1", "repaired", "qwen"),
    ("d2", "head", "deepseek"),
    ("d2", "head", "qwen"),
)


def command(point, tree, model, index):
    arm = f"{point}-{tree}-{model}"
    common = [
        PY,
        "-W",
        "ignore",
        str(TOOL),
        "--home",
        str(S / f"home-{point}"),
        "--out",
        str(OUT / arm),
        "--arm",
        arm,
        "--provider-config",
        str(PROFILES[model]),
        "--profile",
        "alibaba-token-plan",
        "--first-sample",
        str(index),
        "--samples",
        "1",
        "--expect-tree",
        TREES[tree][1],
    ]
    if point == "d1":
        inputs = S / "d1-inputs"
        return arm, common + [
            "--task",
            str(inputs / "TASK.md"),
            "--inputs",
            str(inputs / "pbenzyne.xyz"),
            "--envelope",
            str(inputs / "envelope.yaml"),
            "--prefix",
            str(inputs / "prefix-transcript.json"),
            "--cut",
            "7",
            "--probe-turns",
            "3",
        ]
    inputs = S / "d2-inputs"
    return arm, common + [
        "--task",
        str(inputs / "TASK.md"),
        "--inputs",
        *(
            str(inputs / name)
            for name in (
                "benzyl-azide.xyz",
                "trifluoromethyl-ynoate.xyz",
                "triazole-ester-at-c4.xyz",
                "triazole-ester-at-c5.xyz",
            )
        ),
        "--envelope",
        str(inputs / "execution-envelope.yaml"),
        "--probe-turns",
        "4",
        "--goal",
        "goal-po3-r19",
        "--granted-by",
        "claude-researcher-behav-owner-delegated",
        "--max-revisions",
        "12",
    ]


def last_row(arm):
    path = OUT / arm / "samples.jsonl"
    lines = path.read_text().splitlines() if path.exists() else []
    return json.loads(lines[-1]) if lines else {}


def main():
    start = 0
    if "--from" in sys.argv:
        start = int(sys.argv[sys.argv.index("--from") + 1])
    # Amendment A3: a second runner may take one cell; the order of every
    # other cell is unchanged. --only ARM runs that cell's positions only;
    # --skip ARM runs every position but that cell's.
    only = sys.argv[sys.argv.index("--only") + 1] if "--only" in sys.argv else ""
    skip = sys.argv[sys.argv.index("--skip") + 1] if "--skip" in sys.argv else ""
    log_name = "runner-" + only + ".log" if only else "runner.log"
    OUT.mkdir(parents=True, exist_ok=True)
    plan = [
        command(point, tree, model, index)
        for index in range(N)
        for point, tree, model in CELLS
    ]
    manifest = OUT / "manifest.json"
    if not manifest.exists():
        manifest.write_text(
            json.dumps(
                [{"arm": arm, "argv": argv} for arm, argv in plan], indent=1
            )
        )
    log = (OUT / log_name).open("a")
    transport_streak = 0
    for position, (arm, argv) in enumerate(plan):
        if position < start:
            continue
        if (only and arm != only) or (skip and arm == skip):
            continue
        env = dict(os.environ)
        env["PYTHONPATH"] = str(TREES[arm.split("-")[1]][0])
        began = time.time()
        proc = subprocess.run(argv, env=env, capture_output=True, text=True)
        row = last_row(arm)
        line = (
            f"{time.strftime('%Y-%m-%dT%H:%M:%S')} #{position} {arm} "
            f"exit={proc.returncode} {round(time.time() - began)}s "
            f"real={row.get('real_turns')} infra={row.get('infra')} "
            f"failed_attempts={row.get('failed_attempts', 0)}"
        )
        log.write(line + "\n")
        log.flush()
        print(line, flush=True)
        if proc.returncode != 0:
            log.write(proc.stderr[-3000:] + "\n")
            log.flush()
        failed = (not row.get("real_turns")) and (
            row.get("failed_attempts") or row.get("error")
        )
        transport_streak = transport_streak + 1 if failed else 0
        if transport_streak >= 3:
            log.write("STOP: three consecutive samples lost to transport\n")
            print("STOP: three consecutive samples lost to transport")
            return
    log.write("DONE\n")
    print("DONE")


if __name__ == "__main__":
    main()
