#!/usr/bin/env bash
# Tier-2 wizard e2e: the real `chemsmart wizard --server --yes` binary runs
# against PATH-shim scheduler executables emitting pinned real SLURM output,
# in an isolated HOME. Proves the full onboarding flow -- detection, queue
# probe, file write, backup, EXTRA_COMMANDS carry-forward, verification
# table, and a renderable submit script -- without a scheduler installed.
set -euo pipefail

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
SHIMS="$WORK/shims"
export HOME="$WORK/home"
mkdir -p "$SHIMS" "$HOME"

cat > "$SHIMS/sinfo" <<'SHIM'
#!/usr/bin/env bash
if [ "${1:-}" = "--version" ]; then
  echo "slurm 24.05.4"
  exit 0
fi
# The wizard requests: sinfo -h -o %P|%a|%l|%c|%m|%D|%G
echo "compute*|up|2-00:00:00|64|256000|4|(null)"
echo "gpu|up|1-00:00:00|32|512000|2|gpu:a100:4"
SHIM
cat > "$SHIMS/sbatch" <<'SHIM'
#!/usr/bin/env bash
echo "Submitted batch job 1"
SHIM
chmod +x "$SHIMS"/sinfo "$SHIMS"/sbatch
export PATH="$SHIMS:$PATH"

fail() { echo "FAIL: $1" >&2; exit 1; }

echo "== first run: a fresh machine =="
chemsmart wizard --server --yes

CFG="$HOME/.chemsmart/server/SLURM.yaml"
test -f "$CFG" || fail "SLURM.yaml was not written"
grep -q "QUEUE_NAME: compute" "$CFG" || fail "default partition not chosen"
grep -q "NUM_CORES: 64" "$CFG" || fail "cores must come from the queue"
grep -q "MEM_GB: 225" "$CFG" || fail "memory must derive from the queue"
grep -q "NUM_HOURS: 24" "$CFG" || fail "hours must cap at the default"
grep -q "detected SLURM (slurm 24.05.4)" "$CFG" || fail "provenance missing"
python - "$CFG" <<'PY'
import sys, yaml
payload = yaml.safe_load(open(sys.argv[1]))
assert payload["SERVER"]["SCHEDULER"] == "SLURM"
assert isinstance(payload.get("GAUSSIAN", {}), dict), "program blocks missing"
PY

echo "== refresh run: hand-tuned EXTRA_COMMANDS must survive =="
python - "$CFG" <<'PY'
import sys
path = sys.argv[1]
text = open(path).read()
text = text.replace(
    "    EXTRA_COMMANDS: |\n        # Host commands required before execution.",
    "    EXTRA_COMMANDS: |\n        ulimit -s unlimited",
)
open(path, "w").write(text)
PY
chemsmart wizard --server --yes
ls "$HOME/.chemsmart/server/"SLURM.yaml.bak-* >/dev/null || fail "no backup"
grep -q "ulimit -s unlimited" "$CFG" || fail "EXTRA_COMMANDS was dropped"
grep -q "carried from the previous configuration" "$CFG" || fail "carry note"

echo "== the written file renders a submit script =="
RUN="$WORK/run" && mkdir -p "$RUN" && cd "$RUN"
cat > water.xyz <<'XYZ'
3
water
O        0.000000    0.000000    0.117300
H        0.000000    0.757200   -0.469200
H        0.000000   -0.757200   -0.469200
XYZ
chemsmart sub --test --fake -s SLURM xtb -p test -f water.xyz sp
SCRIPT="$RUN/water_sp/chemsmart_sub_water_sp.sh"
test -f "$SCRIPT" || fail "no submit script rendered"
grep -q -- "--ntasks-per-node=64" "$SCRIPT" || fail "queue cores not in script"
grep -q -- "--partition=compute" "$SCRIPT" || fail "partition not in script"
grep -q "ulimit -s unlimited" "$SCRIPT" || fail "extra commands not in script"

echo "== no scheduler: the wizard is honest =="
HOME2="$WORK/home2" && mkdir -p "$HOME2"
# A minimal PATH carrying only python and chemsmart -- no scheduler
# binaries -- regardless of what the invoking machine has installed.
CHEMSMART_BIN="$(dirname "$(command -v chemsmart)")"
PY_BIN="$(dirname "$(command -v python)")"
env HOME="$HOME2" PATH="$CHEMSMART_BIN:$PY_BIN:/usr/bin:/bin" \
  chemsmart wizard --server --yes
CFG2="$HOME2/.chemsmart/server/local.yaml"
test -f "$CFG2" || fail "local.yaml was not written"
grep -q "SCHEDULER: null" "$CFG2" || fail "local file must not claim a scheduler"
grep -q "no batch scheduler detected" "$CFG2" || fail "honest provenance missing"

echo "wizard shim e2e: all checks passed"
