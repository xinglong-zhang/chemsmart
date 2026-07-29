# ChemSmart Agent HPC Installation

This is the non-destructive installation route for the Zhang Lab cluster.
Use it only after the release notice publishes both the
`agent-hpc-rc1` tag and its full commit SHA. A tag name alone is not a
reproducibility receipt.

## Install Into A Separate Environment

Do not repair or upgrade Conda `base`, and do not modify an existing
`chemsmart` environment. Clone a separate checkout and create a dedicated
environment:

```bash
git clone https://github.com/Hongjiseung-ROK/chemsmart.git chemsmart-agent-src
cd chemsmart-agent-src
git fetch --tags --force
git checkout --detach agent-hpc-rc1
git rev-parse HEAD

conda create -n chemsmart-agent -c conda-forge \
  python=3.10 xtb=6.7.1 openbabel pip
conda activate chemsmart-agent
python -m pip install -e ".[agent,agent-tui]"
```

Compare the full `git rev-parse HEAD` value with the SHA in the release
notice before continuing. Do not infer the expected SHA from this document.

## Verify The Selected Runtime

These checks distinguish the active executable and Python package from a
stale editable installation:

```bash
command -v chemsmart
head -n 1 "$(command -v chemsmart)"
python -c 'import sys, chemsmart; print(sys.executable); print(chemsmart.__file__)'
python -m pip show chemsmart
git status --short --branch
git rev-parse HEAD
chemsmart --version
chemsmart --help
chemsmart agent --help
chemsmart agent doctor --no-ping --require-xtb
chemsmart agent doctor --tool-probe --require-xtb
```

`--no-ping` only checks local configuration. It is not proof that the
provider supports tool calls. The second doctor command forces one inert
tool call and executes no ChemSmart tool.

Provider failures are classified as follows:

- 401 or 403: account or authorization;
- 404: model configuration;
- 400 or a schema error: gateway adapter or tool-schema compatibility;
- timeout or 5xx: provider operational failure after retries.

The lab release gate requires the configured cluster
`deepseek-v4-pro` tool probe to pass.

## Diagnose An Existing Environment Without Modifying It

If an older environment lacks `chemsmart agent`, collect these read-only
facts instead of reinstalling into it:

```bash
conda info --envs
conda list -n chemsmart
conda run -n chemsmart sh -lc 'command -v chemsmart'
conda run -n chemsmart python -m pip show chemsmart
git -C "$HOME/chemsmart" status --short --branch
git -C "$HOME/chemsmart" rev-parse HEAD
conda run -n chemsmart chemsmart --version
conda run -n chemsmart chemsmart --help
```

The Agent command cannot appear when the checked-out commit predates the
Agent package, even if an editable install succeeds. Conda startup errors
mentioning `pydantic_core._pydantic_core` indicate a separate Conda
plugin/base-environment problem; this guide intentionally does not repair
Conda base.

## Capture One Failed Session

The TUI log is overwritten at every launch. Preserve it before reopening
the TUI:

```bash
cp "$HOME/.chemsmart/agent/_tui.log" \
  "$HOME/chemsmart_tui_$(date +%Y%m%d-%H%M%S).log"
chemsmart agent sessions --all
chemsmart agent debug-bundle SESSION_ID \
  --output "$HOME/chemsmart_SESSION_ID_debug.tar.gz"
```

Replace `SESSION_ID` with the explicit failed YAML session. `latest` is
rejected so a later diagnostic turn cannot silently replace the evidence.
The bundle contains only redacted decision, runtime, metadata, and
environment-manifest files. It excludes `agent.yaml`, `api.env`, and API
credentials.

Real local xTB execution remains an approval-required action. The doctor
and tool-probe commands do not start a calculation.
