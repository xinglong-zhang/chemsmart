# ChemSmart Agent — terminal design guideline

This file is the design authority for `chemsmart/agent/tui/`. Any change to
what a human sees or reads in the terminal lands in the same commit as its
entry here. The audience is a computational chemist at a terminal — often
over SSH, often inside tmux, sometimes color-blind, always busy.

## Identity (from the CHEMSMART logo)

The logo: a half-atom/half-gear mark (navy strokes `#003078`, cyan nucleus
`#48f0f0`, steel gear `#787890`), the wordmark **CHEM** in purple `#901890`
+ **SMART** in green `#186018`, underlined by a vibrational-spectrum
waveform running green → cyan → magenta → pink.

Terminal translation:

- **Wordmark** — `CHEM` bold magenta + `SMART` bold green + ` Agent` plain
  + dim version (`panels.wordmark()`). Magenta/green are the wordmark's
  native ANSI pair; it renders everywhere.
- **The spectrum rule** — the waveform as a one-line divider of block
  glyphs `▂▄▆█▆▄▂▁` in four segments `#51cf66 → #4dd3e8 → #c2255c →
  #ff8fab` (`panels.spectrum_rule`). Exactly two sanctioned placements:
  under the wordmark at startup, and before an "Analysis results" panel.
  Under `--plain` it degrades to a dim `─` rule. Do not add placements.
- **Busy pulse** — while the host works, the phase dot breathes through the
  spectrum colors once per second; static in `--plain`.
- **Glyph language** — `⚛` molecule/identity/program, `⚙` machinery, `✱`
  inspection/preview, `→`/`←` read/write, `Σ` analysis, `$` engine, `▣`
  model turn, `△` blocked, states `[ ] [•] [»] [✓] [✗] ⊘` (panels/humanize
  registries). A glyph always accompanies color; color is never the only
  signal.

## Palette

`theme.py` is the single semantic source. Interaction accent: cyan
(`#22b8cf`, the nucleus). Semantics: success green, warning yellow, error
red — never repurposed. The wordmark magenta/green and the spectrum colors
are brand marks, not semantic states. Navy and steel stay in the logo; the
terminal gets depth from background steps, not extra hues.

## The one state rule

Finished goes muted · active stays text · pending human decision is
warning · failure is error · validated success is green. (`_ROW_STYLE`,
`panels.state_glyph`.)

## Voice

`voice.py` is the vocabulary authority. Machine tokens are translated
exactly once (`engine_complete` → "validating result",
`waiting_for_approval` → "ready for your review", `blocked_unsupported` →
"not executable in this release"); unknown values pass through because
host-authored sentences are already written for humans. "provider-free" is
charter vocabulary and stays; a standalone "provider" meaning the LLM is
"model". Node and artifact ids (`xtb-sp-gfn2-calc`, `water-monomer`) are
scientific names — never rewritten.

## The no-hash rule

No 64-hex digest is ever rendered to a human — not in the TUI, not in the
default CLI output. Digests are provenance and live on, byte-untouched, in
the durable records (events, bundles, review files, the report `.md`) and
behind the `--json` flags. The regex `[0-9a-f]{64}` over every rendered
panel is the test guard (`assert_no_bare_hex`). Units, molecular identity
(formula · charge/multiplicity), and physical conditions are always shown;
that is the science, and science is never elided.

## Accessibility

`NO_COLOR` is honored (Textual driver); `--plain` renders inline, linear,
mouse-free, animation-free, with brand flourishes degraded to plain rules.
Wide tables scroll in place; the transcript is keyboard-scrollable
(pageup/pagedown, ctrl+home/end); every command is reachable by typing and
via ctrl+p. Red and green are never the only distinction between two
states.

## Grounding

clig.dev (human-first defaults, announce discovered configuration, never
silent); scientific-software usability (design for the named user base —
here, a computational chemist); accessible-CLI practice (linear output,
NO_COLOR). Insight taken from Claude Code, Codex CLI, and opencode; copied
from none.
