# Agent Wizard Validation Status

Snapshot of which HPC environments the `chemsmart agent wizard` has been live-validated against, what remains untested, and how to reproduce each verdict.

Captured: 2026-05-12 against `fork/main` HEAD `8ef92261`.

## Summary

- **Live-validated schedulers**: SGE/SoGE 8.1.8, OpenPBS 23.06 (including `nodes:ppn` legacy + `select=` modern resource defaults), SLURM 22.05 (single-node + simulated multi-node).
- **Live-validated module systems**: Lmod 8.7 (tcl modulefile, lua `module show` output).
- **Audit gaps closed by live runs**: 5 of 8 P0 items.
- **Audit gaps still open**: 3 P0, 10 P1, 3 P2 — none have a live cluster yet.
- **Regression coverage**: 108 wizard tests including end-to-end PBS round-trip + 4 PBS-resource-overlay rows + 4 module-extract rows; SLURM/PBS/SGE node-overlay and Lmod setenv/path-extract paths locked in pytest.

## Validated environments

| Environment | OS / Packaging | Scheduler | Mode | Result | Evidence |
|---|---|---|---|---|---|
| Cluster #1 (external) | CentOS 7 | SoGE 8.1.8 (qconf/qstat) | A: orchestrator runs on submit host | PASS | PR #81, PR #86 |
| chemnode1 (GCP) | Rocky 8 / OpenHPC | OpenPBS 23.06.06 (queue empty defaults) | B: orchestrator → SSH login node | PASS | PR #90, locked by `tests/agent/wizard/test_e2e.py::test_mode_b_openpbs_round_trip_uses_node_overlay_and_conda_fallback` (PR #93) |
| chemnode1 (GCP) | Rocky 8 / OpenHPC | OpenPBS 23.06 with `resources_default.select = 2:ncpus=4:mem=8gb` | B | PASS — `select=` string parser ignored, node facts win | live verify 2026-05-12 |
| chemnode1 (GCP) | Rocky 8 / OpenHPC | OpenPBS 23.06 with `resources_default.nodes = 1:ppn=4` (Torque legacy) | B | PASS — auto-translated `ncpus=4` overridden by node overlay | PR #95 fixtures (`qstat_qf_workq_nodes_ppn.txt`) + live verify |
| chemnode1 (GCP) | Rocky 8 / OpenHPC | OpenPBS 23.06 with `resources_default.ncpus = 4` (direct queue default) | B | PASS — queue default overridden by node overlay | PR #95 fixtures (`qstat_qf_workq_ncpus4.txt`) + live verify |
| chemnode2 (GCP) | Rocky 9 / EPEL | SLURM 22.05.9 single-node | B | PASS | PR #91 fixtures (`scontrol_show_node_chemnode2.txt`) |
| chemnode2 (GCP) | Rocky 9 / EPEL | SLURM 22.05.9 multi-node (3 nodes via `State=FUTURE`) | B | PASS — host-node facts win over partition aggregates | cs-96 live run; partition `TotalCPUs=18 TRES=cpu=18,mem=35000M,node=3`, wizard YAML `NUM_CORES=2 MEM_GB=3` |
| chemnode1 (GCP) | Rocky 8 / OpenHPC | OpenPBS 23.06 + Lmod 8.7 with fake `gaussian/16` modulefile (`setenv g16root`, `setenv GAUSS_EXEDIR`, `prepend-path PATH`) | B | PASS — version `16` preserved in `MODULES`, `EXEFOLDER=/opt/fake/gaussian/g16` extracted from setenv, `g16root`/`GAUSS_EXEDIR` propagated to `ENVARS` | PR #97 fixtures (`lmod_module_show_gaussian_16.txt`, `lmod_module_avail_with_gaussian.txt`) + live verify |

Every "PASS" above means the full `SERVER` block validated and `CONDA_ENV` multi-line activation rendered for `GAUSSIAN`/`ORCA`/`NCIPLOT`.

## What "validated" exercised

Each live run covered:

- `SCHEDULER` family detection from probe execution (no false positives between SLURM/PBS/SGE).
- `QUEUE_NAME` extraction from native query (`scontrol show partition --oneliner`, `qstat -Qf`, `qconf -sql`).
- `NUM_CORES` sourced from **host node facts** (`scontrol show node <self>`, `pbsnodes -av` `pcpus`/`ncpus`, SGE `complex_values`), NOT from partition or queue aggregates.
- `MEM_GB` sourced from host node `RealMemory` / `resources_available.mem`, normalized to GB via `_render_server_mem_gb`.
- `NUM_HOURS` derived from queue `MaxTime` / `resources_max.walltime`.
- `CONDA_ENV` discovered through the `_WELL_KNOWN_CONDA_PATHS` fallback when `conda` is not on the SSH login `$PATH` (chemnode1, chemnode2 both meet this condition with `/opt/miniforge3`).
- `SUBMIT_COMMAND`, `HOST`, `SCRATCH_DIR`, validator clearance (no missing required keys).

## Audit gaps

The reproducibility audit from cs-95 (PR #92) enumerated wizard correctness gaps. Status as of `93f4d923`:

### Closed by live work

| ID | Gap | Closure |
|---|---|---|
| P0 #1 | SLURM partition-aggregate over-estimation on multi-node | Audited via FUTURE-node injection; node overlay always wins (cs-96) |
| P0 #2 | PBS `select=N:ncpus=M:mem=…` chunk spec on `resources_default` | Live-injected on chemnode1 2026-05-12; parser ignores `select=` string, node overlay wins. No code change required. |
| P0 #2b | PBS legacy `nodes=N:ppn=M` (Torque-style) on `resources_default` | Fixed by PR #95 (`parsers.py` + `survey.py` PBS node overlay parallel to PR #91 SLURM). PBS auto-translated `ncpus=M` from `ppn=M` no longer wins over node `pbsnodes` facts. |
| P0 #3 | SLURM `NUM_CORES` from partition `TotalCPUs` instead of host node | Fixed by PR #91 (`survey.py:112-142` node CPU override) |
| P0 #6 | `module` / `Lmod` false negative → EXEFOLDER `null` + version stripped from `module load` line | Fixed by PR #97 (`software.py` setenv/path extraction parallel to Lmod lua + tcl modulefile formats; `_normalize_module_line` rejects family-stub `gaussian/`). Closes both sub-bugs (#6a version-strip, #6b EXEFOLDER not extracted). |

### Still open (no live exposure yet)

| ID | Gap | Why it matters | Next step |
|---|---|---|---|
| P0 #4 | LSF parser bugs (no live cluster) | `bqueues -l` / `bhosts` extraction may be wrong | Defer — LSF licensing |
| P0 #5 | SGE node state `alarm` / `suspend` / `subordinate` not filtered as unusable | Wizard could pick a queue that won't dispatch | New GCE SoGE cluster, force state via `qmod -d` |
| P0 #7 | Spack-only sites — no `module`, conda discovery alone misses Gaussian | Same symptom as #6 from different cause | Future Spack chemnode |
| P1 (multiple) | Bright/xCAT vendor paths (`/cm/local/apps/*/current`), Cray PALS/ALPS, tcsh/zsh login shells, mixed-scheduler gateway nodes, remote Mode B over jumphost | Each blocks at probe execution layer | Defer until a real site appears |
| P2 | Cosmetic YAML ordering, comment style | No semantic effect | Tracked, not blocking |

## Reproduction recipes

### Mode B live wizard against chemnode1 (PBS) or chemnode2 (SLURM)

```bash
# Pre-flight: confirm editable clone matches fork/main
cd /Users/hongjiseung/developer/chemsmart
git fetch fork main
git merge-base --is-ancestor HEAD fork/main \
  || { echo "STALE — git stash -u && git checkout fork/main first"; exit 1; }

# PBS
AI_PROVIDER=openai python -m chemsmart.cli.main agent wizard chemnode1 \
  --host chemsmart@34.146.86.156

# SLURM
AI_PROVIDER=openai python -m chemsmart.cli.main agent wizard chemnode2 \
  --host chemsmart@35.189.152.165
```

Expected criteria are listed in the `tests/agent/wizard/test_e2e.py` assertions and in the [PBS](chemnode-pbs-provisioning.md) and [SLURM](chemnode-slurm-provisioning.md) runbooks.

### Multi-node SLURM regression (no extra cost)

See the `State=FUTURE` injection recipe inline below. Always back up `/etc/slurm/slurm.conf` first.

```bash
ssh chemsmart@35.189.152.165 'sudo cp /etc/slurm/slurm.conf /etc/slurm/slurm.conf.pre-multinode'
ssh chemsmart@35.189.152.165 'sudo sed -i "/^NodeName=chemnode2/a NodeName=fake01 CPUs=8 RealMemory=16000 State=FUTURE\nNodeName=fake02 CPUs=8 RealMemory=16000 State=FUTURE" /etc/slurm/slurm.conf'
ssh chemsmart@35.189.152.165 'sudo sed -i "s|^PartitionName=workq Nodes=chemnode2|PartitionName=workq Nodes=chemnode2,fake01,fake02|" /etc/slurm/slurm.conf'
ssh chemsmart@35.189.152.165 'sudo systemctl restart slurmctld'

# run wizard; expect NUM_CORES=2, MEM_GB=3

# Restore
ssh chemsmart@35.189.152.165 'sudo cp /etc/slurm/slurm.conf.pre-multinode /etc/slurm/slurm.conf && sudo systemctl restart slurmctld'
```

### Full pytest

```bash
cd /Users/hongjiseung/developer/chemsmart  # at fork/main
AI_PROVIDER=openai python -m pytest -v tests/agent/wizard/
# Expected: 108 passed
```

### PBS resource-default audit (re-runnable, 0₩)

To re-verify P0 #2 / P0 #2b closure on chemnode1:

```bash
# A1 modern select= (expect NUM_CORES=2)
ssh chemsmart@34.146.86.156 'sudo /opt/pbs/bin/qmgr -c "set queue workq resources_default.select = 2:ncpus=4:mem=8gb"'
AI_PROVIDER=openai python -m chemsmart.cli.main agent wizard scratchA1 --host chemsmart@34.146.86.156

# A2 Torque legacy nodes:ppn (expect NUM_CORES=2 after PR #95)
ssh chemsmart@34.146.86.156 'sudo /opt/pbs/bin/qmgr -c "unset queue workq resources_default.select" && sudo /opt/pbs/bin/qmgr -c "set queue workq resources_default.nodes = 1:ppn=4"'
AI_PROVIDER=openai python -m chemsmart.cli.main agent wizard scratchA2 --host chemsmart@34.146.86.156

# A3 direct ncpus (expect NUM_CORES=2 after PR #95)
ssh chemsmart@34.146.86.156 'sudo /opt/pbs/bin/qmgr -c "unset queue workq resources_default.nodes" && sudo /opt/pbs/bin/qmgr -c "unset queue workq resources_default.ncpus" && sudo /opt/pbs/bin/qmgr -c "unset queue workq resources_default.nodect" && sudo /opt/pbs/bin/qmgr -c "set queue workq resources_default.ncpus = 4"'
AI_PROVIDER=openai python -m chemsmart.cli.main agent wizard scratchA3 --host chemsmart@34.146.86.156

# Cleanup (back to empty queue defaults)
ssh chemsmart@34.146.86.156 'sudo /opt/pbs/bin/qmgr -c "unset queue workq resources_default.ncpus"'
```

## Cluster runtime references

- chemnode1: `34.146.86.156`, e2-medium, asia-northeast1-a — runbook `docs/ops/chemnode-pbs-provisioning.md`, canonical script `scripts/ops/provision-chemnode-pbs.sh`.
- chemnode2: `35.189.152.165`, e2-medium, asia-northeast1-a — runbook `docs/ops/chemnode-slurm-provisioning.md`, canonical script `scripts/ops/provision-chemnode-slurm.sh`.

Both VMs are e2-medium in `asia-northeast1-a`. Each idle hour is roughly ₩45. There is currently no GCP budget alarm on the project.

## Known caveats

- All live clusters above are **single physical node** (multi-node tested only via FUTURE simulation). Real heterogeneous compute pools remain unaudited.
- All clusters use **bash login shell**. tcsh / zsh login defaults are unexercised.
- All wizard runs are **Mode B** (orchestrator → SSH to login node). Mode A (running directly on submit host) was only exercised for the original SoGE cluster.
- The wizard's `software.py` conda discovery fallback (`_WELL_KNOWN_CONDA_PATHS`) is exercised by both GCE nodes; sites that install conda elsewhere (e.g. `~/miniconda3/`, `/usr/local/anaconda3/`) are not covered.

## Wave history

| PR | What it locked in |
|---|---|
| #81 | SoGE 8.1.8 baseline |
| #86 | Queue/host parsing + conda activation block |
| #89 | SLURM `RealMemory` fallback + `_WELL_KNOWN_CONDA_PATHS` |
| #90 | OpenPBS 23.06 live validation + integrated conda discovery |
| #91 | SLURM `NUM_CORES` from node, not partition (closes P0 #3) |
| #92 | This `docs/ops/` runbooks + `scripts/ops/` canonical provisioning |
| #93 | End-to-end OpenPBS regression test (locks PR #90 contract) |
| #94 | This document (initial snapshot) |
| #95 | PBS node overlay — `pbsnodes` per-host facts override queue `resources_default.ncpus`, closing P0 #2b |
| #96 | Doc update (this file) reflecting PRs #94/#95 |
| #97 | Module `setenv`/`prepend-path` extraction + version preservation — closes P0 #6 (Lmod lua + tcl modulefile formats) |
