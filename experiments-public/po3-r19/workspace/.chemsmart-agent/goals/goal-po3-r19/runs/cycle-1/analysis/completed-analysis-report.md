# Host-validated structured analysis

Completion receipt: `33a98e74270c60292725a0844f1c9104111762e668bd2c88000f7cb3340783cd`
Toolchain plan: `94b21486a0fc4d08485939d7cba17cfe9b603bec3c8a2c497e2a7208d8646683`

## Expected versus delivered

| Observable | Expected | Delivered | Unit | Agreement | Basis |
|---|---|---:|---|---|---|
| ddg-activation-regio-353k | positive 0.0..3.0 | `` |  | not_comparable | Thermal azide-alkyne cycloadditions are documented to be only weakly regioselective, and here both alkyne termini carry acceptors that direct the terminal azide nitrogen in opposite senses, so the two asynchronous transition states should differ by tenths to a few kcal/mol rather than by many. |
| dg-ts-esterc4-minus-esterc5-353k | positive 0.0..2.5 | `` |  | not_comparable | Frontier-orbital and electrostatic reading of the alkyne: the ester is a pi-acceptor and pushes LUMO density onto the distal CF3-bearing alkyne carbon, while CF3 acts mainly through sigma-induction and should not outweigh that resonance effect; the CF3-bearing carbon should therefore be both the larger-LUMO-coefficient and the more delta-positive terminus, so the nucleophilic terminal azide nitrogen (which becomes N3 and bonds C4) attacks it, placing CF3 at C4 and the ester at C5. This opposes the task's conventional 'ester to C4' expectation, so it is declared as a falsifiable diagnostic rather than assumed. |

An expectation is displayed, never scored: a diverging row settles nothing and means the chemistry disagreed with the reasoning, which is a result the reader owns.

Scientific decision: not recorded -- interpretation is a session act; this run executed extraction, thermochemistry, expressions, validation verdicts, and claim rendering only.

## Record delivery

| Record | Nodes (reached state) | Calculation | Analysis |
| --- | --- | --- | --- |
| 1 | opt-benzyl-azide (validated) | validated | none |
| 2 | opt-tfm-ynoate (validated) | validated | none |
| 3 | scan-esterc4-path (validated) | validated | executed |
| 4 | scan-esterc5-path (validated) | validated | executed |

A batch of N is N observations; each verdict above is one record's, and no aggregate quantity is rendered.
