# `/wizard` engineering review

Audit target: `main` at `43b34012`. I read the requested wizard, agent-integration, CLI/TUI, and test files; I also cross-checked the gate expectations from `docs/research/wizard_feasibility.md` §6 via historical commit `76216b7c` because that file is not present on `main`. Validation run: `pytest -v tests/agent/wizard` → 63/63 passed in 0.23s.

## 1. Security audit

### SSH command construction

- **Good defaults exist**: remote probes use argv-based `ssh` invocation with `BatchMode=yes`, `ClearAllForwardings=yes`, `ConnectTimeout=10`, and `StrictHostKeyChecking=yes` (`chemsmart/agent/wizard/probe.py:81-97`). That correctly rejects MFA prompts and unknown host keys instead of silently accepting them.
- **The allowlist is not a real sandbox**: `_validate_ssh_command()` checks only the **first token** (`probe.py:114-123`), then the whole string is executed inside `bash -lc` (`probe.py:96`). That means shell composition is still possible whenever a caller can influence the command string.
- I verified the validator currently allows all of the following:
  - `printf %s\n $(whoami)`
  - `printf %s\n "$HOME" && rm -rf /`
  - `env python -c "print(1)"`
  - `command bash -lc "echo pwn"`
- The same design flaw exists on the local path: `_validate_local_command()` also only checks the first token (`probe.py:107-123`), so launcher commands like `env` and `command` bypass the intended read-only model.
- Bottom line: **this is a P0 remote/local command-execution risk**, not just a quoting nit.

### Allowlist coverage gaps

- `ALLOWED_COMMANDS` includes **launcher-style** or **too-powerful** entries such as `env`, `command`, `find`, `echo`, `printf`, and `module` (`probe.py:13-44`).
- `env` can execute arbitrary subprocesses; `command` can invoke arbitrary programs; `find` can be abused with `-exec`; `printf`/`echo` gain command substitution on the SSH path because the string is passed through `bash -lc`.
- The codebase duplicates free-form probe wrappers in `software.py`, `scratch.py`, `project.py`, and `survey.py`, which makes it easier for one caller to accidentally widen the surface later.
- This is an **allowlist-by-prefix**, not an allowlist-by-capability.

### Prompt injection through probe stdout

- The design has one positive property: `wizard_probe` does **not** return raw `ProbeResult.stdout` blobs to the model. It mostly returns normalized facts (`chemsmart/agent/wizard/tools.py:16-38`). That meaningfully reduces prompt-injection blast radius.
- The remaining risk is in software discovery:
  - `_find_module_candidates()` takes matching `module -t avail` lines and keeps the first whitespace-delimited token (`software.py:198-230`).
  - `render_server_yaml()` then auto-renders `module load <candidate>` from the first candidate (`render.py:205-209`).
- So a malicious or malformed module line can still flow into:
  - `module_candidates` in the tool result,
  - rendered YAML in `MODULES`,
  - any future LLM context that sees the tool output.
- I would rate this **P1, partially mitigated but not closed**.

### `write_server_yaml` path traversal

- `write_server_yaml()` directly interpolates `name` into `~/.chemsmart/server/<name>.yaml` (`write.py:8-20`).
- There is **no basename validation** in CLI, TUI, tool wrapper, or writer (`agent/cli.py:411-483`, `agent/tui/screens/chat.py:1229-1375`, `wizard/tools.py:41-58`).
- I verified that `server_name='../../../tmp/wizard_escape'` writes outside the server directory. On this host it resolved to `/private/.../tmp/wizard_escape.yaml`, not under `~/.chemsmart/server/`.
- This is a **P0 arbitrary-path write** bug.

### File permissions on written YAML

- `target.write_text(...)` relies entirely on process umask (`write.py:15-19`).
- On this machine the resulting mode was **`0644`**.
- That is too permissive for per-user cluster config by default. Even if the file usually holds no secrets, it can still expose hostnames, project/account names, scratch paths, module stacks, and future environment details.
- Production default should be **`0600`** (and ideally the parent directory should already be user-private).

## 2. Latency profile

### Sequential probe budget per turn

All probes default to a **15 s timeout** (`probe.py:67-83`). There is no overall survey deadline, no connection reuse, and almost everything is serial.

Worst-case timeout budget for a full `/wizard` turn:

- **Topology preflight** (`topology.py:30-55`): `env` + 4 local scheduler commands = **5 probes = 75 s**.
- **Scheduler survey** (`survey.py:35-172`):
  - SLURM: 2 probes
  - PBS: 2 probes
  - LSF: 2 probes
  - SGE: 1 `qconf -sql` probe **plus one `qconf -sq` per queue name**
  - Non-SGE floor: **7 probes = 105 s**.
- **Software survey** (`software.py:37-273`):
  - module detection: up to 3 probes
  - conda discovery: 2 probes
  - per program (`gaussian`, `orca`, `nciplot`): up to 3 probes each in Mode A (`command -v`, `which`, `module -t avail`) or 2 each in Mode B
  - Worst-case Mode A: **14 probes = 210 s**
  - Worst-case Mode B: **11 probes = 165 s**
- **Scratch discovery** (`scratch.py:28-125`): env probe + home scratch test + up to 4 writability tests = **up to 6 probes = 90 s**.
- **Project discovery** (`project.py:19-205`): SLURM path can do env probe + 2 `sacctmgr` variants + `groups` = **4 probes = 60 s**.
- **Extra commands probe** (`render.py:228-256`): one extra local probe in Mode A = **15 s**.

That yields:

- **Mode A worst-case (non-SGE)**: 75 + 105 + 210 + 90 + 60 + 15 = **555 s (~9.25 min)**
- **Mode B worst-case (non-SGE)**: 75 + 105 + 165 + 90 + 60 = **495 s (~8.25 min)**
- **SGE is effectively unbounded by queue count**, because `qconf -sq` is looped per queue with no cap.

This is too slow for an interactive setup command.

### SSH connection reuse opportunity

- Every remote probe is a **fresh SSH process** (`probe.py:81-104`).
- A non-SGE Mode B survey can easily hit **28 separate SSH sessions** in one turn:
  - scheduler 7
  - software 11
  - scratch 6
  - project 4
- That is the single biggest latency design flaw. On real clusters with ProxyJump, Kerberos, or high RTT, connection setup will dominate the survey.
- `ControlMaster`/`ControlPersist`/`ControlPath` or another persistent transport would likely cut end-to-end wall time by a large factor; **10x is realistic** when handshake cost dwarfs the actual read-only command.
- There is also easy intra-turn deduplication available:
  - `module -t avail` is re-run separately for each program,
  - local/remote probe helpers do not share results.

### Cache hit-rate expected behavior

- The sidecar cache is only consulted by refresh flows (`refresh.py:31-36`, `71-90`).
- `run_wizard()` itself does **not** read cache before reprobe (`orchestrator.py:32-87`).
- Therefore expected behavior today is:
  - first `/wizard`: **cache miss**
  - second `/wizard` five seconds later: **still a cache miss**
  - `/wizard-refresh` within 24 h and `--force` absent: **cache hit**
- So the cache is operationally useful for `wizard-refresh` and downstream consumers, but it provides **0% hit rate for the primary `/wizard` latency path**.

## 3. Failure-mode coverage gaps (what `test_e2e.py` doesn't cover — MFA SSH, MTU/packet loss, transient cluster outage, hostkey rotation)

`tests/agent/wizard/test_e2e.py` is a stubbed logical round-trip, not a transport/system integration test. It validates orchestration shape, but it does not exercise the real failure modes that matter in production.

Missing coverage includes:

- **MFA / keyboard-interactive SSH**
  - No test proves the user experience when `BatchMode=yes` hits MFA and the probe aborts.
  - No test covers the recommended workaround path (existing ControlMaster session / agent-backed auth / rerun on login node).
- **MTU issues, packet loss, partial stdout, slow links**
  - No test injects truncated JSON, partial `module avail` output, or mid-stream scheduler output corruption.
  - `ProbeRunner` timeout handling is unit-tested, but higher-level parser recovery and user-facing messaging are not.
- **Transient cluster outage / scheduler brownout**
  - No end-to-end test covers one scheduler command succeeding and the next failing, or JSON/text format drift during outage windows.
- **Hostkey rotation / `known_hosts` mismatch / first-seen host**
  - There is no integration test for `StrictHostKeyChecking=yes` failure modes, so the gate exists in code but not in regression coverage.
- **Real local login-node shell behavior**
  - This is the biggest blind spot: Mode A local probes are stubbed as if `module`, `type module`, and similar shell-context-dependent commands behave like normal executables. In reality, local Mode A and remote Mode B use different execution semantics.

Additional untested adversarial cases that should exist before calling this production-ready:

- malicious `server_name` path traversal,
- file mode assertion (`0600` vs ambient umask),
- malicious `module -t avail` candidate text,
- command-prefix bypass attempts (`env`, `command`, shell metacharacters),
- large SGE queue counts and bounded survey behavior.

## 4. Deployment recommendation

### pip extras

- `agent-tui` already isolates Textual-only dependencies (`pyproject.toml:54-84`).
- For production cluster setup, **that is not enough**. `/wizard` is deterministic and non-TUI, but it is still shipped inside a very broad runtime that includes AI SDKs and even test-oriented packages.
- My recommendation:
  - **Yes, introduce a minimal wizard-capable install target** — either a dedicated `[wizard]` extra or, better, a new `[agent-core]` extra that includes wizard/verify/refresh without TUI and without LLM-provider baggage.
  - Do **not** stop at a cosmetic new extra while leaving the core dependency layering untouched; the real problem is that probe-only deploys currently inherit too much unrelated runtime surface.

### conda-forge readiness

As packaged today, I would call this **not ready yet** for a clean conda-forge story:

- runtime dependencies include **`pytest`** and `pytest-mock`-adjacent expectations (`pyproject.toml:18-52`),
- dependency policy is inconsistent (many exact pins, many open-ended lower bounds, no disciplined upper bounds for fast-moving SDKs),
- wizard is coupled to a large monolithic package surface instead of a small probe-oriented one.

None of these are unsalvageable, but they are release-engineering debt.

### container approach

Recommended production approach:

- publish **two images**:
  - a minimal non-root CLI/wizard image,
  - an optional TUI image layered on top.
- mount `~/.chemsmart` at runtime; do not bake mutable per-user config into the image,
- never bake SSH private keys, agent sockets, or `known_hosts` into image layers,
- if SSH is needed from inside the container, prefer host agent forwarding or read-only bind mounts plus the existing `StrictHostKeyChecking=yes` posture,
- for clusters where login-node containers are discouraged, a **locked micromamba/venv** install is a better operational default than a heavyweight image.

### version pinning audit

Current pinning is neither solver-friendly nor cleanly reproducible:

- many scientific dependencies are **exact-pinned**,
- many fast-moving packages are **lower-bound-only** (`requests`, `rich`, `tenacity`, `anthropic`, `openai`, `pydantic`, etc.),
- test packages leak into runtime requirements.

Recommendation:

- for the application: keep a **tested constraints/lock** artifact for supported deployment environments,
- for published package metadata: relax exact pins except where ABI compatibility truly requires them,
- add explicit upper bounds for volatile SDKs,
- remove test-only packages from runtime metadata.

## 5. Code quality

### Test coverage gaps

- The wizard test suite is broad numerically (63 passing tests), but **too stub-heavy** to prove production behavior.
- CLI/TUI routing tests exist and are useful (`tests/agent/test_cli_wizard*.py`, `tests/agent/test_slash_commands_permission.py`), but they also stop at mocked boundaries.
- The suite is strongest on normalization logic and weakest on transport realism, shell behavior, and adversarial inputs.

### Duplication

There is repeated probe/plumbing logic in several places:

- `_run_probe()` duplicated in `software.py`, `scratch.py`, `project.py`, and `survey.py`,
- `_probe_env_values()` duplicated in `software.py` and `project.py`,
- server-path resolution duplicated in `verify.py` and `refresh.py`,
- `module -t avail` is re-probed per program instead of once per survey.

This is not just style debt; it increases security drift risk because quoting and transport assumptions are re-implemented in multiple files.

### Complexity hotspots

By inspection, the highest-risk hotspots are:

- `render_server_yaml()` (`render.py:53-145`) — many policy decisions in one function,
- `ProbeRunner._run()` (`probe.py:126-177`) — timeout/error/truncation semantics concentrate here,
- `refresh_cache()` + `_build_node_summary()` (`refresh.py:31-68`, `191-242`) — cache semantics and summary shaping are hand-rolled,
- `parsers.py` (427 lines) — four scheduler families in one file,
- `project._discover_sacctmgr_project()` (`project.py:79-133`) — multiple fallbacks and parsing branches.

The most important qualitative hotspot is **Mode A vs Mode B execution asymmetry**:

- remote probes run in `bash -lc`,
- local probes run as bare subprocess argv,
- so module functions, shell initialization, and environment expansion behave differently.

That asymmetry makes the code look deterministic in tests while being much less deterministic on real login nodes.

### Missing typed exceptions

Some domain-specific exceptions exist (`ProbeError`, `NoTargetError`, `AmbiguousSchedulerError`), but the error model is still too generic overall:

- `run_schedule_survey()` raises a generic `RuntimeError` for total failure (`survey.py:61-63`),
- `orchestrator.py`, `refresh.py`, `verify.py`, and `registry.py` catch broad `Exception` in multiple places,
- `_safe_parse()` swallows parse failures without surfacing which parser/payload failed (`survey.py:175-179`).

That makes telemetry, retry policy, and user messaging less precise than they should be for production operations.

## 6. P0/P1/P2 issue list with concrete mitigation per item

| Priority | Issue | Evidence | Concrete mitigation |
|---|---|---|---|
| **P0** | Probe allowlist can be bypassed into arbitrary command execution | `probe.py:114-123`, `probe.py:96`; validator accepted `env python -c ...`, `command bash -lc ...`, and `printf ... && rm ...` | Replace free-form command strings with a **closed probe catalog** (enum + fixed argv/script template). Remove launcher commands from `ALLOWED_COMMANDS`. Reject shell metacharacters entirely. |
| **P0** | `write_server_yaml()` allows path traversal and writes outside `~/.chemsmart/server` | `write.py:15-19`; verified with `../../../tmp/wizard_escape` | Validate `server_name` against a strict basename regex, resolve the target, and enforce parent containment before writing. |
| **P0** | YAML file permissions are ambient-umask driven (`0644` here) | `write.py:19`; observed mode `0o644` | Create files with explicit `0600` (e.g. `os.open(..., 0o600)` or post-write `chmod`) and keep the parent dir user-private. |
| **P1** | Full survey latency is interactive-hostile (8–9+ minute worst-case) | `probe.py:67-83`, `topology.py`, `survey.py`, `software.py`, `scratch.py`, `project.py` | Add a global survey budget, dedupe probes, parallelize independent discovery steps after scheduler detection, and cap SGE queue fan-out. |
| **P1** | Remote probing opens a new SSH session for every command | `probe.py:81-104`; ~28 SSH sessions for a typical Mode B survey | Add `ControlMaster`/`ControlPersist`/`ControlPath` or another persistent transport layer; batch multiple read-only probes per session. |
| **P1** | Local Mode A probe semantics differ from remote and will miss module-shell state | local direct subprocess vs remote `bash -lc`; `software.py:37-64`, `203-208` | Make local probing use the same shell semantics as remote (e.g. `bash -lc` locally too), or explicitly document and code a shell-wrapper abstraction. |
| **P1** | Multiple module candidates are silently auto-selected by shortest name, contrary to the feasibility gate | `software.py:212-219`, `render.py:197-209`; feasibility §6 called for warn-pause | Treat multi-candidate module discovery as an ambiguity object. Require user confirmation or a deterministic ranking rule richer than “shortest string wins.” |
| **P1** | Non-writable scratch still renders program `SCRATCH: True` | `render.py:113-117`, `179-193`; feasibility §6 recommended forcing false | If scratch writability is not confirmed, render program `SCRATCH: False` and annotate the YAML/comment path explicitly. |
| **P1** | Cache does not accelerate `/wizard` itself; hit rate is effectively 0% on repeated probes | `orchestrator.py:32-87`, `refresh.py:31-36` | Read fresh cache before reprobe, surface staleness to the user, and only force reprobe when inputs changed or TTL expired. |
| **P2** | Stub-heavy tests miss SSH/MFA/hostkey/network/security regressions | `tests/agent/wizard/test_e2e.py`, `test_probe.py`, `test_survey.py` | Add hermetic integration tests with a fake `ssh` binary, fixture `known_hosts`, malformed/truncated outputs, and malicious inputs. |
| **P2** | Generic exception handling obscures failure class and recoverability | `survey.py:61-63`, `survey.py:175-179`, broad catches in `orchestrator.py`, `refresh.py`, `verify.py`, `registry.py` | Introduce typed exceptions for transport failure, scheduler parse failure, module ambiguity, invalid server name, and cache refresh failure; preserve error codes in tool results. |
