# experiments-public — two CHEMSMART Agent cases, as they actually ran

These are not demonstrations assembled after the fact. Each directory holds
one autonomous run exactly as the Agent produced it, plus a reconstruction of
what the Agent decided and what evidence it had when it decided.

| case | question | settlement |
|---|---|---|
| `ino3-r12` | Where does the +/0 couple of a square-planar Ni(II) bis(thiolate) fall versus ferrocene in acetonitrile, and does the hole sit on nickel or on sulfur? | `achieved_with_observations`, 7 engine calls |
| `po3-r19` | Which regioisomer does a thermal azide-alkyne cycloaddition favour, and by how much at 353 K? | `unreachable_from_evidence`, 11 engine calls |

`unreachable_from_evidence` is a delivered result, not a failure. The value was
computed and delivered; the host verified that the **precision the question
demanded** -- 0.5 kcal/mol, fixed by the chemist -- is not reachable in that
envelope. Three levels of theory gave two different answers about which
regioisomer wins, and saying so is the useful thing to tell someone about to
commit 300 mg of a four-step alkyne.

## What is in each directory

- `REVIEW.md` -- the case reconstructed from its own evidence: the problem as
  the chemist posed it, what the Agent decided at each stage, what evidence
  existed *before* each claim, which claims were supported, revised, rejected
  or left uncertain, and why each calculation followed the one before it. The
  Agent's own reasoning is quoted, including the approaches it **rejected**,
  which is often where the science is.
- `figures/` and `make_figures.py` -- every figure, and the script that
  regenerates each number it plots from the workspace beside it.
- `workspace/` -- the run itself: the goal ledger, the run event streams, the
  public transcripts, the project YAML the Agent wrote, the materialised
  program inputs, and the engine outputs. Scratch directories and lock files
  are removed, as are ORCA's binary restart artefacts -- 89 files and 148 MB
  of `.gbw`, `.densities`, `.opt` and `.cpcm` state, each listed with its
  SHA-256 in `EXCLUDED-ARTEFACTS.md`. Nothing in either report, and nothing in
  the typed analysis plane, reads them; a converged wavefunction scales as the
  square of the basis size, so the eight def2-TZVP single points alone wrote
  50 MB that no reader can open. What is kept is what a claim can stand on:
  every engine log, every submitted input, every Hessian, every geometry, and
  the Agent's own receipts and event streams.

`PEER-REVIEW.md` is an independent referee's audit of both reports, with its
findings and the corrections it requested. The reports were revised against
it; the referee's chemistry concerns about the *underlying work* are recorded
in each report's final section rather than corrected, because the runs are
historical.

## Reading order, if you want the short path

1. `po3-r19/REVIEW.md`, first section -- the Agent measured the ring
   connectivity of the chemist's two product files and found the file names
   transposed relative to IUPAC numbering. That cost no engine call and is the
   most immediately useful thing either case produced.
2. `po3-r19/figures/figB1-scans.png` -- two relaxed scans, and the evidence
   that their maxima were walls rather than barriers.
3. `ino3-r12/figures/figA2-spin-per-functional.png` -- where the hole sits,
   and how much that answer depends on the functional.
