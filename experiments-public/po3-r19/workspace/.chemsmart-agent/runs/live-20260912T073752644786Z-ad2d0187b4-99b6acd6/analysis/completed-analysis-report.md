# Host-validated structured analysis

Completion receipt: `bcb9555f70411a0e26e3854bb8be16cdab7bc91c61065a16d48c7153e99ca961`
Toolchain plan: `ff8fb42ff010148e8c0cf95d3d8f32214c48f8497846feb4238e61e8f995cf8c`

## Validation verdicts

| Node | Rule | Predicate | Observed | Bound | Unit | Verdict |
|---|---|---|---:|---:|---|---|
| saddle-order-verdict | c4chan-exactly-one-imaginary-mode | integer_equals | `1` | `1` |  | passed |
| saddle-order-verdict | c5chan-exactly-one-imaginary-mode | integer_equals | `1` | `1` |  | passed |

## Expected versus delivered

| Observable | Expected | Delivered | Unit | Agreement | Basis |
|---|---|---:|---|---|---|
| c4-restart-imaginary-mode-count | 1.0..1.0 | `` |  | not_comparable | The cycle-2 search was cut off by the iteration limit rather than by leaving the ridge: its reached structure still shows both forming N...C contacts at 2.07 and 2.27 A and its parsed mode list already carries an imaginary mode, so continuing the same walk from that structure should terminate on the same first-order saddle of the ester-at-C4 channel. |
| ddg-activation-regio-353k | positive 0.0..3.0 | `` |  | not_comparable | Thermal azide-alkyne cycloadditions are documented to be only weakly regioselective, and here both alkyne termini carry acceptors that direct the terminal azide nitrogen in opposite senses, so the two asynchronous transition states should differ by tenths to a few kcal/mol rather than by many. |
| ddg-basis-shift-tzvp-kcal | positive 0.0..0.5 | `` |  | not_comparable | The two saddles are the same 30-atom neutral formula differing only in which alkyne terminus each azide nitrogen bonds; def2-SVP already carries polarization on every heavy atom, so the one-particle basis error is largely systematic and cancels in the difference, leaving a few tenths at most. |
| ddg-functional-shift-tzvp-kcal | positive 0.0..1.0 | `` |  | not_comparable | Global hybrids and range-separated or double-hybrid functionals are documented to differ by up to about 1 kcal/mol on differences between competing asynchronous cycloaddition barriers, where the error is set by how each treats the delocalised, charge-transfer-flavoured forming bonds. |
| ddg-regio-signed-spread-across-levels-353k | positive 0.5..1.5 | `` |  | not_comparable | The host-written workspace record already carries the three levels' signed differences (-0.456, -0.221 and +0.675 kcal/mol for the channel the workflow names esterc4 minus the one it names esterc5), so their range cannot be smaller than about 1.13 kcal/mol; this diagnostic is declared to make that level dependence a typed, scored quantity rather than prose. |
| dg-ts-esterc4-minus-esterc5-353k | positive 0.0..2.5 | `` |  | not_comparable | Frontier-orbital and electrostatic reading of the alkyne: the ester is a pi-acceptor and pushes LUMO density onto the distal CF3-bearing alkyne carbon, while CF3 acts mainly through sigma-induction and should not outweigh that resonance effect; the CF3-bearing carbon should therefore be both the larger-LUMO-coefficient and the more delta-positive terminus, so the nucleophilic terminal azide nitrogen (which becomes N3 and bonds C4) attacks it, placing CF3 at C4 and the ester at C5. This opposes the task's conventional 'ester to C4' expectation, so it is declared as a falsifiable diagnostic rather than assumed. |

An expectation is displayed, never scored: a diverging row settles nothing and means the chemistry disagreed with the reasoning, which is a result the reader owns.

Scientific decision: not recorded -- interpretation is a session act; this run executed extraction, thermochemistry, expressions, validation verdicts, and claim rendering only.
