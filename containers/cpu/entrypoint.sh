#!/usr/bin/env bash
set -euo pipefail

# xTB's documented OpenMP runtime requires a large process stack. The hard
# limit is controlled by the container host; raise the soft limit when allowed.
ulimit -s unlimited 2>/dev/null || true

export OMP_STACKSIZE="${OMP_STACKSIZE:-4G}"
export OMP_MAX_ACTIVE_LEVELS="${OMP_MAX_ACTIVE_LEVELS:-1}"

if [[ "$#" -eq 0 ]]; then
  set -- chemsmart --help
fi

exec "$@"
