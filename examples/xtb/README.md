# xTB examples

These small, fixed XYZ files exercise ChemSmart's command path without asking
a model to author native xTB input. The `test` project uses GFN2-xTB and the
project loader rejects unknown sections or settings.

Safe fake local previews:

```bash
chemsmart run --fake --no-scratch xtb -p test -f examples/xtb/water.xyz sp
chemsmart run --fake --no-scratch xtb -p test -f examples/xtb/water.xyz opt
chemsmart run --fake --no-scratch xtb -p test -f examples/xtb/methane.xyz hess
```

Submission-script previews (no scheduler submission):

```bash
chemsmart sub -s configured_server --test --fake xtb -p test -f examples/xtb/water.xyz sp
chemsmart sub -s configured_server --test --fake xtb -p test -f examples/xtb/water.xyz opt
chemsmart sub -s configured_server --test --fake xtb -p test -f examples/xtb/methane.xyz hess
```

Real execution is a separate approval boundary. The runner resolves `xtb` from
the configured server block or `PATH`, uses argv rather than shell syntax, and
refuses to overwrite stale output artifacts.
