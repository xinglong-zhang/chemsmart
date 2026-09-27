# The broken-symmetry request, as the hub writes it

ORCA 6.1.1 on CUHK Charles, R10 Q18 oracle O1 (CUHK Slurm 2153479), written
and run through the ordinary CLI on code `5da66f9c`, project
`gas: {functional: b3lyp, basis: def2-svp, broken_symmetry: true,
ri_approximation: none, scf_convergence: tight, defgrid: defgrid3}`, the fixed
geometries of oracle O0. The hub wrote `%scf HFTyp UHF / GuessMix 45 end`
(ORCA's B3LYP/G, the Gaussian VWN3 form).

| file | geometry | printed | reading |
|---|---|---|---|
| `o_h2_074_bs_gas_phase.out` | H2 at 0.74 A | `FINAL SINGLE POINT ENERGY -1.173496798071`, `<S**2>` 0.000000 | stayed spin-symmetric: the RKS energy |
| `o_pbenzyne_bs_gas_phase.out` | p-benzyne, regular hexagon | `-230.704724039355`, `<S**2>` 0.970279 | the broken-symmetry singlet; identical to the native GuessMix input of O0 (CUHK 2153330) to the printed digit |

## Two S = 1 centres flipped from the high-spin state

ORCA 6.1.1 on CUHK Charles, R11 truth-2's Ni(II)2 oracle (CUHK Slurm
2157086): ax41 ino2's hydroxo- and chloro-bridged dinickel core with three
ammines per nickel (`dinickel-oh-cl.xyz`, sha256 `be1a5c68...`, built from
reported metrics, not optimised), charge 2, multiplicity 5 on the
coordinate line, B3LYP/G def2-SVP, SlowConv, TightSCF, maxiter 500. Run
through the ordinary CLI from a native `input_string` (an oracle outside
the typed layer); ORCA switched to HFTyp UHF itself.

| file | request | printed | reading |
|---|---|---|---|
| `ni2_fs1_gas_phase.out` | `FlipSpin 1`, `FinalMs 0.0` (ORCA counts from 0: the second Ni) | `-3890.968544659585`, `<S**2>` 1.998044 | the antiferromagnetic state, Ni spin +1.726 / -1.723 (Mulliken) |
| `ni2_bs22_gas_phase.out` | `BrokenSym 2,2` | `-3890.968534751536`, `<S**2>` 1.998045 | the same state |
