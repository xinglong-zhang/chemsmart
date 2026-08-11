# Contributing to CHEMSMART

CHEMSMART welcomes changes that improve reusable computational-chemistry
workflows while preserving one CLI and project-YAML authority across programs.

## Development environment

Python 3.10 is required. On Linux or macOS:

```bash
git clone https://github.com/Hongjiseung-ROK/chemsmart.git
cd chemsmart
make env
conda activate chemsmart
make install-dev
make pre-commit
```

For a headless x86 Ubuntu benchmark host, follow
`docs/source/installation-ubuntu-cpu-server.rst` before configuring external
programs.

## Before editing

1. Read `AGENTS.md` and the smallest relevant `.agents/skills` entry.
2. Inspect the active branch and dirty worktree.
3. Reproduce the observed behavior with the live CLI or a focused source
   probe.
4. Decide which existing layer owns the behavior: CLI, project settings,
   job/writer, runner, parser, analysis, or agent.

Do not introduce a second execution path, hand-authored native input, or a
benchmark-specific special case.

## Development workflow

Create a focused branch and keep unrelated changes out of the commit:

```bash
git switch -c feature/short-description
```

Implement the smallest general change. Add or update user documentation when
the public behavior changes. One focused test is usually sufficient while
iterating; run the broader checks before publication.

```bash
pytest tests/path/to/focused_test.py
make docs
make test
```

Formatting uses Black and isort with a 79-character Python line length. Ruff
configuration lives in `pyproject.toml`.

## Computational-chemistry changes

A program integration is complete only when the relevant path is coherent:

1. project YAML and CLI option;
2. run/sub command registration;
3. program setting and native materialisation;
4. runner and environment selection;
5. parser observation and normal termination;
6. structured quantity and units where applicable; and
7. user documentation.

Preserve geometry frame, atom order, coordinate units, charge, multiplicity,
state, constraints, method and physical conditions. Report program or method
limitations rather than silently substituting a different calculation.

## Agent and skill changes

Agent changes should be motivated by a visible behavioral failure. Run a real
model trial when provider access is authorized, read the answer as a
computational scientist, and accept creative valid routes. Tests can confirm
mechanics but do not grade scientific intelligence.

Update a project skill when a repeated general lesson would help future
research. Keep skills concise, self-improving, and free of molecule, paper,
DOI, expected-value, or private-DAG answers. Validate edited skill folders with
the skill creator's `quick_validate.py`.

## Documentation and generated artifacts

Maintained documentation belongs in `README.md` and `docs/source`. Do not
commit provider credentials, `~/.chemsmart`, engine binaries, scratch data,
generated inputs/outputs, benchmark transcripts, one-off reviews, diagnostic
reports, or agent reasoning.

## Commits and publication

Use focused conventional commits, for example:

```bash
git commit -m "fix(orca): validate semantic normal termination"
git commit -m "docs: add Ubuntu CPU benchmark setup"
```

Before pushing, inspect the complete diff and verify that relevant checks ran.
A release tag or package upload is a separate maintainer decision; `make
release` must not be used as part of an ordinary contribution.
