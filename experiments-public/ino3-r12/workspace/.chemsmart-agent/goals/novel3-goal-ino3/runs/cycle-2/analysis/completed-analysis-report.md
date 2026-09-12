# Host-validated structured analysis

Completion receipt: `52045219035629621ba9b86f1adbcdc4319dae3cb1a956c9964866e6ddcc1edc`
Toolchain plan: `dabbdf5b2852c05cfa10d4169ca849fc543e95c107ac37915a009779376e2fe5`

## Thermochemical conditions

| Stage | Temperature (K) | Standard state | Entropy model | Frequency scale |
|---|---:|---|---|---:|
| thermo-ph3-cat | `298.15` | 1 mol/L | rrho | `1` |
| thermo-ph3-n0 | `298.15` | 1 mol/L | rrho | `1` |
| thermo-pme3-cat | `298.15` | 1 mol/L | rrho | `1` |
| thermo-pme3-n0 | `298.15` | 1 mol/L | rrho | `1` |

Claim record: `3efd3377f93799671325d8eb00ad5a9f0e8c6c1615482ca5af27b8a52aa47303`

## Host-rendered numerical claims

| Claim | Value | Unit | Source receipt |
|---|---:|---|---|
| d-ni-s1-neut-pme3 | `2.239839845617762` | angstrom | `60b181f3167c37a76170fbce335434731da22d50cd76809ded6b1cf7e8937f7c` |
| d-ni-s2-neut-pme3 | `2.2407249755840186` | angstrom | `60b181f3167c37a76170fbce335434731da22d50cd76809ded6b1cf7e8937f7c` |

Claim record: `52918d77e011a12a5c7a8013426ffc368ca51a278a3c2fd3c0108096180e26d4`

## Host-rendered numerical claims

| Claim | Value | Unit | Source receipt |
|---|---:|---|---|
| spin-square-ph3 | `0.76777` | 1 | `4e9f4118555ae0b8bc7d23c2e37fd5f159cf17c273a8605cd4dfee07f844fd56` |
| spin-square-qrt-ph3 | `3.762335` | 1 | `4e386a2c42bae746b323679f7ecbc67ff2e28bcdf94264b5cbfdbeefa72c6bfc` |
| spin-square-qrt-pme3 | `3.764735` | 1 | `780cc81462b5ea488cd7eaccfd18e10605a7252cb2b08128d7147409df236ac2` |
| spin_square | `0.76777` | 1 | `4e9f4118555ae0b8bc7d23c2e37fd5f159cf17c273a8605cd4dfee07f844fd56` |

Claim record: `9701e152b3ea4badc55f61cacdc73fefb0aed5534e434609357c5ad770eb4448`

## Host-rendered numerical claims

| Claim | Value | Unit | Source receipt |
|---|---:|---|---|
| quartet-doublet-gap | `0.3501978769796836` | eV | `44684eb91d9dae2659cf35c7aceadab437597f0155efea34b3cea9a73fabaf35` |
| quartet-doublet-gap-ph3 | `0.43655518626541767` | eV | `8bd11de0d105507b9101d1dd776e5cf7d19c9474e79300abb75189d0cbf1ef8a` |

Claim record: `db219011db7dc4c3055ceed35f9e97f1799436b4b56b41a2192c71d621312cac`

## Host-rendered numerical claims

| Claim | Value | Unit | Source receipt |
|---|---:|---|---|
| dq-ni-ph3 | `0.042186` | e | `59b55d92e40c26448c3b3ff7180eff046f6e823d6b2fa79b587c13c197176c0d` |
| dq-ni-pme3 | `0.076837` | e | `235b419975c62fdd9f3cc8cd66160cbb07388798fdb981b4df2a19b8dd8ea954` |
| dq-s1-ph3 | `0.197205` | e | `59b55d92e40c26448c3b3ff7180eff046f6e823d6b2fa79b587c13c197176c0d` |
| dq-s1-pme3 | `0.149781` | e | `235b419975c62fdd9f3cc8cd66160cbb07388798fdb981b4df2a19b8dd8ea954` |
| dq-s2-ph3 | `0.19717900000000002` | e | `59b55d92e40c26448c3b3ff7180eff046f6e823d6b2fa79b587c13c197176c0d` |
| dq-s2-pme3 | `0.149739` | e | `235b419975c62fdd9f3cc8cd66160cbb07388798fdb981b4df2a19b8dd8ea954` |

Claim record: `f88549a17beef2803b74233e208619539c860b88221cf2b86709790dcc4e70b1`

## Host-rendered numerical claims

| Claim | Value | Unit | Source receipt |
|---|---:|---|---|
| couple-shift-ph3-minus-pme3 | `0.19630264228458832` | V | `6f4afbd489d878667135a79ed038e07dd96efec6d624119e9ab8709724e7b19a` |

## Validation verdicts

| Node | Rule | Predicate | Observed | Bound | Unit | Verdict |
|---|---|---|---:|---:|---|---|
| val-ph3-n0 | min-order | minimum_greater_equal | `25.59` | `-20` | cm^-1 | passed |
| val-ph3-cat | min-order | minimum_greater_equal | `30.2` | `-20` | cm^-1 | passed |
| val-ph3-qrt | min-order | minimum_greater_equal | `50.7` | `-20` | cm^-1 | passed |
| val-pme3-n0 | min-order | minimum_greater_equal | `20.57` | `-20` | cm^-1 | passed |
| val-pme3-qrt | min-order | minimum_greater_equal | `29.35` | `-20` | cm^-1 | passed |

## Expected versus delivered

| Observable | Expected | Delivered | Unit | Agreement | Basis |
|---|---|---:|---|---|---|
| quartet-doublet-gap | positive 0.05..1.5 | `0.3501978769796836` | eV | agreed | Square-planar-derived Ni thiolate cations are normally low-spin doublets; a competitive quartet would require large exchange stabilization from metal-localized d electrons. |
| quartet-doublet-gap-ph3 | positive 0.05..1.5 | `0.43655518626541767` | eV | agreed | Same low-spin expectation as for the PMe3 model; PH3 is a worse donor, which should if anything stabilize the low-spin form further. |
| redox-potential-vs-fc | positive 0.0..0.7 | `` |  | not_comparable | The user's shelf oxidants (Fc+ 0.00 V, acetylferrocenium +0.27 V, tris(4-bromophenyl)aminium +0.70 V) were chosen to bracket the likely couple; Ni(II) bis(thiolate) bis(phosphine) complexes typically oxidize within this window in acetonitrile. |
| redox-potential-vs-fc-ph3 | positive 0.0..1.0 | `` |  | not_comparable | Weaker phosphine donors destabilize the oxidized (cationic) form, so the PH3 model should oxidize at higher (more positive) potential than the PMe3 model; the magnitude of the shift is being measured. |

An expectation is displayed, never scored: a diverging row settles nothing and means the chemistry disagreed with the reasoning, which is a result the reader owns.

Scientific decision: not recorded -- interpretation is a session act; this run executed extraction, thermochemistry, expressions, validation verdicts, and claim rendering only.

## Record delivery

| Record | Nodes (reached state) | Calculation | Analysis |
| --- | --- | --- | --- |
| 1 | pme3-n0-opt2 (validated) | validated | executed |
| 2 | ph3-n0-opt2 (validated) | validated | executed |
| 3 | ph3-cat-opt2 (validated) | validated | executed |
| 4 | pme3-qrt-opt2 (validated) | validated | executed |
| 5 | ph3-qrt-opt2 (validated) | validated | executed |
| shared | claims-dq (executed), claims-gaps (executed), claims-redox-block (blocked_unsupported), claims-shift (executed), claims-spin-block (blocked_unsupported), claims-spinsq (executed), expr-dq-ph3 (executed), expr-gap-ph3 (executed), expr-shift (executed), ext-pme3-cat-existing (executed), thermo-pme3-cat (executed) | - | partial |

A batch of N is N observations; each verdict above is one record's, and no aggregate quantity is rendered.
