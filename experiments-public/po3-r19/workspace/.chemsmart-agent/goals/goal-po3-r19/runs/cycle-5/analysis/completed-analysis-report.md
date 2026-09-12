# Host-validated structured analysis

Completion receipt: `9437307b1a5b5573c7c0a4b875e469029673ab17589e6be21a8aebeeb43300d9`
Toolchain plan: `cdb6bcd982aaa783954484fb81e459c1f4bbd9eaf81e972092f6ea2c306afa97`

## Thermochemical conditions

| Stage | Temperature (K) | Standard state | Entropy model | Frequency scale |
|---|---:|---|---|---:|
| thermo-c4-qrrho | `353` | 1 atm | grimme (cutoff 100 cm^-1) | `1` |
| thermo-c4-rrho | `353` | 1 atm | rrho | `1` |
| thermo-c4-unconv-rrho | `353` | 1 atm | rrho | `1` |
| thermo-c5-qrrho | `353` | 1 atm | grimme (cutoff 100 cm^-1) | `1` |
| thermo-c5-rrho | `353` | 1 atm | rrho | `1` |

## Precision requirements

| Observable | Required | Stated | Basis | Whose | Combined by | State | The host observed |
|---|---:|---:|---|---|---|---|---|
| ddg-activation-regio-353k | 0.5 kcal/mol | 1.1796354033136858 kcal/mol | measured | task | not stated | short | uncertainty.model_authored_constants[literal_value=0.75@lit-model-allowance; power_exponent=2.0@sq-conv; power_exponent=2.0@sq-ent; power_exponent=2.0@sq-method; power_exponent=2.0@sq-model], uncertainty.source_receipts=12, uncertainty.arithmetic_nodes=6 conventions=none, uncertainty.upstream_scope=receipt |

Only `met` discharges a requirement, and it states that the host resolved a stated uncertainty inside a declared tolerance -- never that the estimate is adequate, which is the session's claim to defend.

Claim record: `200fbae3a0747966e7adc39417383e571dd3a5267e36b1d9c02185670000ff38`

## Host-rendered numerical claims

| Claim | Value | Unit | Source receipt |
|---|---:|---|---|
| ddg-activation-regio-353k | `0.6752142797220688` | kcal/mol | `f25315d1f2c91abaa56ae578bdc86a5b29714c14efa68f3dbf5c1678981d3a4b` |

Claim record: `dd40b003296a6bb05034280a843cef24c2cba69ba51bc3b09245fe8dd6c14ba8`

## Host-rendered numerical claims

| Claim | Value | Unit | Source receipt |
|---|---:|---|---|
| ddg-activation-regio-b3lyp-tzvp-353k | `0.22133609714186211` | kcal/mol | `47909c18c633487f70559e36737dfc5577a5b633ee82da0bd5edb4e053c416bc` |
| ddg-activation-regio-dsd-tzvp-353k | `0.6752142797220688` | kcal/mol | `47909c18c633487f70559e36737dfc5577a5b633ee82da0bd5edb4e053c416bc` |
| ddg-activation-regio-svp-protocol-353k | `0.45559011342754085` | kcal/mol | `f19a9cd19adedb481237f7d60f6f59486dfa205c33154e5c5c855b01a64085dd` |
| ddg-basis-shift-tzvp-kcal | `0.23425401628567874` | kcal/mol | `47909c18c633487f70559e36737dfc5577a5b633ee82da0bd5edb4e053c416bc` |
| ddg-electronic-regio-dsd-tzvp-353k | `1.456543588317599` | kcal/mol | `47909c18c633487f70559e36737dfc5577a5b633ee82da0bd5edb4e053c416bc` |
| ddg-functional-shift-tzvp-kcal | `0.8965503768639309` | kcal/mol | `47909c18c633487f70559e36737dfc5577a5b633ee82da0bd5edb4e053c416bc` |
| ddg-minus-uncertainty | `-0.5044211235916171` | kcal/mol | `f25315d1f2c91abaa56ae578bdc86a5b29714c14efa68f3dbf5c1678981d3a4b` |
| ddg-thermal-contribution-regio-svp-353k | `-0.7813293085955302` | kcal/mol | `f19a9cd19adedb481237f7d60f6f59486dfa205c33154e5c5c855b01a64085dd` |
| ddg-uncertainty-estimator | `1.1796354033136858` | kcal/mol | `f25315d1f2c91abaa56ae578bdc86a5b29714c14efa68f3dbf5c1678981d3a4b` |
| dg-ts-esterc4-minus-esterc5-353k | `0.6752142797220688` | kcal/mol | `47909c18c633487f70559e36737dfc5577a5b633ee82da0bd5edb4e053c416bc` |
| dg-ts-esterc4-minus-esterc5-qrrho-353k | `-0.45559011342754085` | kcal/mol | `f19a9cd19adedb481237f7d60f6f59486dfa205c33154e5c5c855b01a64085dd` |
| regio-ratio-major-to-minor-353k | `2.6183712129579` | 1 | `f25315d1f2c91abaa56ae578bdc86a5b29714c14efa68f3dbf5c1678981d3a4b` |
| thermal-ddg-margin-over-threshold | `0.17521427972206877` | kcal/mol | `f25315d1f2c91abaa56ae578bdc86a5b29714c14efa68f3dbf5c1678981d3a4b` |
| u-convergence-kcal | `0.15886191014899` | kcal/mol | `f19a9cd19adedb481237f7d60f6f59486dfa205c33154e5c5c855b01a64085dd` |
| u-entropy-model-kcal | `0.0` | kcal/mol | `f19a9cd19adedb481237f7d60f6f59486dfa205c33154e5c5c855b01a64085dd` |
| u-method-kcal | `0.8965503768639309` | kcal/mol | `47909c18c633487f70559e36737dfc5577a5b633ee82da0bd5edb4e053c416bc` |

## Expected versus delivered

| Observable | Expected | Delivered | Unit | Agreement | Basis |
|---|---|---:|---|---|---|
| c4-restart-imaginary-mode-count | 1.0..1.0 | `` |  | not_comparable | The cycle-2 search was cut off by the iteration limit rather than by leaving the ridge: its reached structure still shows both forming N...C contacts at 2.07 and 2.27 A and its parsed mode list already carries an imaginary mode, so continuing the same walk from that structure should terminate on the same first-order saddle of the ester-at-C4 channel. |
| ddg-activation-regio-353k | positive 0.0..3.0 | `0.6752142797220688` | kcal/mol | agreed | Thermal azide-alkyne cycloadditions are documented to be only weakly regioselective, and here both alkyne termini carry acceptors that direct the terminal azide nitrogen in opposite senses, so the two asynchronous transition states should differ by tenths to a few kcal/mol rather than by many. |
| ddg-basis-shift-tzvp-kcal | positive 0.0..0.5 | `0.23425401628567874` | kcal/mol | agreed | The two saddles are the same 30-atom neutral formula differing only in which alkyne terminus each azide nitrogen bonds; def2-SVP already carries polarization on every heavy atom, so the one-particle basis error is largely systematic and cancels in the difference, leaving a few tenths at most. |
| ddg-functional-shift-tzvp-kcal | positive 0.0..1.0 | `0.8965503768639309` | kcal/mol | agreed | Global hybrids and range-separated or double-hybrid functionals are documented to differ by up to about 1 kcal/mol on differences between competing asynchronous cycloaddition barriers, where the error is set by how each treats the delocalised, charge-transfer-flavoured forming bonds. |
| dg-ts-esterc4-minus-esterc5-353k | positive 0.0..2.5 | `0.6752142797220688` | kcal/mol | indeterminate | Frontier-orbital and electrostatic reading of the alkyne: the ester is a pi-acceptor and pushes LUMO density onto the distal CF3-bearing alkyne carbon, while CF3 acts mainly through sigma-induction and should not outweigh that resonance effect; the CF3-bearing carbon should therefore be both the larger-LUMO-coefficient and the more delta-positive terminus, so the nucleophilic terminal azide nitrogen (which becomes N3 and bonds C4) attacks it, placing CF3 at C4 and the ester at C5. This opposes the task's conventional 'ester to C4' expectation, so it is declared as a falsifiable diagnostic rather than assumed. |

An expectation is displayed, never scored: a diverging row settles nothing and means the chemistry disagreed with the reasoning, which is a result the reader owns.

Scientific decision: not recorded -- interpretation is a session act; this run executed extraction, thermochemistry, expressions, validation verdicts, and claim rendering only.

## Record delivery

| Record | Nodes (reached state) | Calculation | Analysis |
| --- | --- | --- | --- |
| 1 | sp-c4-b3lyp-tzvp (validated) | validated | executed |
| 2 | sp-c5-b3lyp-tzvp (validated) | validated | executed |
| 3 | sp-c4-dsd-tzvp (validated) | validated | executed |
| 4 | sp-c5-dsd-tzvp (validated) | validated | executed |
| shared | claim-diagnostics-and-ladder (executed), claim-requested-observable (executed), ddg-at-task-tolerance (blocked_unsupported), expr-delivered-budget (executed), expr-level-sensitivity (executed), expr-svp-baseline (executed), extr-c4-svp-energy (executed), extr-c5-svp-energy (executed), thermo-c4-qrrho (executed), thermo-c4-rrho (executed), thermo-c4-unconv-rrho (executed), thermo-c5-qrrho (executed), thermo-c5-rrho (executed) | - | partial |

A batch of N is N observations; each verdict above is one record's, and no aggregate quantity is rendered.
