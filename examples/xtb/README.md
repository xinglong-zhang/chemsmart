# xTB Examples

Small XYZ inputs for xTB job submission smoke tests. These examples exercise
the same xTB job creation path used by `run` and `sub`; fake runs do not need
the external `xtb` executable.

Fake local generation checks:

```bash
chemsmart run --fake --no-scratch xtb -p test -f examples/xtb/water.xyz sp
chemsmart run --fake --no-scratch xtb -p test -f examples/xtb/water.xyz opt
chemsmart run --fake --no-scratch xtb -p test -f examples/xtb/methane.xyz hess
```

Submission-script generation without queue submission:

```bash
chemsmart sub -s path/to/server.yaml --test xtb -p test -f examples/xtb/water.xyz sp
chemsmart sub -s path/to/server.yaml --test xtb -p test -f examples/xtb/water.xyz opt
chemsmart sub -s path/to/server.yaml --test xtb -p test -f examples/xtb/methane.xyz hess
```

If `xtb` is available on `PATH` or configured under the server `XTB` section,
one small real smoke test can be run with:

```bash
chemsmart run --no-scratch xtb -p test -f examples/xtb/water.xyz -l water_real_sp sp
chemsmart run --no-scratch xtb -p test -f examples/xtb/water.xyz -l water_real_hess hess
```

Expected generated artifacts are per-label folders containing the copied
`<label>.xyz`, `<label>.out`, and transient `<label>.err` files for local runs;
`sub --test` also writes scheduler/run scripts without submitting to the queue.
