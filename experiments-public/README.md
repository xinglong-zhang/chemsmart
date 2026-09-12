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
- `figures/` -- the structures, the scan profiles, the saddles and the
  level-to-level comparison, each caption naming the node its numbers came
  from.
- `workspace/` -- the run itself: the goal ledger, the run event streams, the
  public transcripts, the project YAML the Agent wrote, the materialised
  program inputs, and the engine outputs. Every claim in `REVIEW.md` is
  traceable to a file in here, and the provenance table at the end of each
  report says which. Scratch directories, lock files and ORCA's binary
  restart state (`.gbw`, `.densities`, `.opt`, `.cpcm`) are removed -- a
  converged wavefunction scales as the square of the basis size, and nothing
  in the typed analysis plane reads one. What is kept is what a claim can
  stand on: every engine log, every submitted input, every Hessian, every
  geometry, and the Agent's own receipts.

## Reading order, if you want the short path

1. `po3-r19/REVIEW.md`, first section -- the Agent measured the ring
   connectivity of the chemist's two product files and found the file names
   transposed relative to IUPAC numbering. That cost no engine call and is the
   most immediately useful thing either case produced.
2. `po3-r19/figures/figB1-scans.png` -- two relaxed scans, and the evidence
   that their maxima were walls rather than barriers.
3. `ino3-r12/figures/figA2-spin-per-functional.png` -- where the hole sits,
   and how much that answer depends on the functional.
