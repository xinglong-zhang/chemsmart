# chemsmart agent quickstart

This page is a single-session transcript for the current preview agent layer.
All commands below were run on 2026-05-09 with `AI_PROVIDER=openai`;
session IDs and temp paths will differ on your machine.

## Setup

Install the optional TUI dependencies and export the provider before running
any agent command:

```bash
pip install -e ".[agent-tui]"
export AI_PROVIDER=openai
```

Anthropic is also supported; swap the environment variable if that is your
provider. In the current preview build used for this transcript, the API key
was loaded from `~/developer/chemsmart/api.env`, with at least:

```dotenv
ai_api_key=...
```

A healthy doctor run prints lines like:

```text
AI_PROVIDER=openai OK
api.env: OK (key length=32)
gateway: https://factchat-cloud.mindlogic.ai/v1/gateway
ping: ok (model=gpt-5.4-2026-03-05, latency=1122ms)
tools registered: 10
```

## First run — advisory question

Ask chemistry-only questions without preparing a job:

```bash
AI_PROVIDER=openai chemsmart agent ask "Recommend method/basis for a Cope rearrangement TS"
```

Sample output (abridged):

```text
╭────────────────────────────────── Request ───────────────────────────────────╮
│ Recommend method/basis for a Cope rearrangement TS                           │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─────────────────────────────────── Advice ───────────────────────────────────╮
│ For a Cope rearrangement TS, a good default is a DFT TS search with          │
│ frequency confirmation using M06-2X/6-31+G(d,p) or wB97X-D/def2-SVP, then a  │
│ single-point refinement at M06-2X/def2-TZVP or wB97X-D/def2-TZVP ...         │
╰──────────────────────────────────────────────────────────────────────────────╯
```

The important behavior is that `ask` stays advisory: it returns chemistry
reasoning directly and does not create an input file.

## First run — dry-run input generation

For the transcript below, I first created a minimal `water.xyz` in the working
directory:

```bash
cat > water.xyz <<'EOF'
3
water
O 0.000000 0.000000 0.000000
H 0.758602 0.000000 0.504284
H -0.758602 0.000000 0.504284
EOF
```

Then I ran the dry-submit path:

```bash
AI_PROVIDER=openai chemsmart agent run "single-point on water.xyz at B3LYP/6-31G(d) Gaussian" --dry-submit
```

Sample stdout (abridged):

```text
session: 20260509T064536Z-aac23479
Plan:
Rationale: This prepares a Gaussian single-point calculation on water.xyz at the exact requested B3LYP/6-31G(d) level; if you want execution, the next step would be run_local or submit_hpc after validation.
Estimated cost: Very low; seconds to minutes on a typical CPU.
1. build_molecule {"filepath": "water.xyz"}
2. build_gaussian_settings {"basis": "6-31G(d)", "functional": "B3LYP", "title": "water_sp"}
3. build_job {"kind": "gaussian.sp", "label": "water_sp", "molecule": "$step1", "settings": "$step2"}
4. dry_run_input {"job": "$step3"}
5. validate_runtime {"job": "$step3"}
critic verdict: ok
inputfile: /private/tmp/chemsmart-agent-docs.8MdMeA/water_sp.com
decision log: /Users/hongjiseung/.chemsmart/agent/sessions/20260509T064536Z-aac23479/decision_log.jsonl
```

The generated Gaussian input was:

```text
%chk=water_sp.chk
%nprocshared=12
%mem=16GB
# B3LYP 6-31G(d)
```

The command exited with code `0`.

## Resuming a session

Re-run the saved plan with the session id printed above:

```bash
AI_PROVIDER=openai chemsmart agent resume 20260509T064536Z-aac23479 --dry-submit
```

That replayed the same plan and printed the same `inputfile:` path for the
dry-run `.com` file.

## Going live (HPC submission)

`chemsmart agent run --help` documents `--dry-submit / --execute` as “Write
scripts without real remote submission, or execute submit_hpc.” In practice,
the default quickstart path is safe: it writes previews, prints the dry-run
input file, and can still surface `validate_runtime` warnings such as missing
queue/account/scratch/module fields. For a real submit, use `--execute` and
make sure your server settings under `~/.chemsmart/server/*.yaml` are ready;
under `--dry-submit`, the submit step is skipped.

## Inspecting what happened

Every run writes an audit trail under `~/.chemsmart/agent/sessions/<id>/`. For
the dry-run session above, the log lived at:

```text
~/.chemsmart/agent/sessions/20260509T064536Z-aac23479/decision_log.jsonl
```

Example event line:

```json
{"kind":"request","payload":{"request":"single-point on water.xyz at B3LYP/6-31G(d) Gaussian"},"rationale":"single-point on water.xyz at B3LYP/6-31G(d) Gaussian","ts":"2026-05-09T06:45:36.351393+00:00"}
```

Read it as: `kind=request` records the original prompt, `payload.request` is
the exact text given to the agent, and `ts` is the UTC timestamp for that
event.

## TUI mode

Running `chemsmart agent` with no subcommand opens the Textual TUI, and
`chemsmart agent --plain` keeps it inline for conservative terminals such as
HPC login nodes. In my launch test, the footer advertised `/help` and the first
`Ctrl+C` changed the footer hint to `Press Ctrl+C again within 3s to quit`;
`/jobs` is the jobs view shortcut once you are in the TUI.

## Common errors

- `AI_PROVIDER not set`: run `export AI_PROVIDER=openai` in the shell before
  any `chemsmart agent ...` command.
- Missing extras for the TUI: if `chemsmart agent` refuses to start the UI, run
  `pip install -e ".[agent-tui]"`.
- Missing input file: the actual dry-run error was
  `/private/tmp/chemsmart-agent-missing.tcJmgC/water.xyz could not be found!`;
  fix it by creating `water.xyz` in the current directory or by naming the file
  explicitly in the request.
