# Host-validated structured analysis

Completion receipt: `f05ff9a828d347c57b8b5e646bc10927c801b2b182bb1cfe8432f590508b1c0d`
Toolchain plan: `3c52e344f836244cb437d6a94af89df5deeca8f07c30a179dfb597154db66da7`

Partial analysis: 6 finding(s) over 27 analysis nodes; 1 claim record(s) rendered.

## Analysis nodes that did not execute

- expr-modes-spinsq (quantity_expression): skipped -- the result carries no value for: ex-n0-tzvp.spinsq-n0-tzvp
- claims-redox (claim_rendering): failed -- uncertainty_reference '6f96fd10e19fd5b30f349c5999d9b6b70febdd5581c2badcb6a7394e002a9b05:uncertainty-e-couple-pme3' does not resolve: the expression's own receipt names a number you supplied in it -- a literal, a scale factor or an exponent -- so its value is yours rather than the host's. Derive that number instead of typing it where the host already holds it: an electron count is the difference of the two states' own charge selectors. Otherwise cite a receipt no number of yours entered, or state the basis as 'asserted'. If a number of yours entered the chain, derive it instead of typing it -- an electron count is the difference of the two states' own charge selectors, a threshold the host owns is a convention rather than your value -- and the whole chain stays the host's. Otherwise name a receipt this host minted or a registered constant, or state the basis as 'asserted'
- val-minima (scientific_validation): skipped -- upstream analysis did not execute: expr-modes-spinsq
- val-spin-state (scientific_validation): skipped -- upstream analysis did not execute: expr-modes-spinsq
- val-svp-replication (scientific_validation): failed -- scientific validation input is not typed evidence from its planned producer
- claims-verdicts (claim_rendering): skipped -- upstream analysis did not execute: expr-modes-spinsq, val-minima, val-spin-state, val-svp-replication

## Quantities the results do not carry

| Quantity | Selector | Why the result does not carry it |
|---|---|---|
| spinsq-n0-tzvp | spin_square | this result prints no <S^2> expectation value; a spin-restricted closed-shell calculation is an eigenfunction of S^2 by construction, so read 'multiplicity' for its electronic state. This result resolves: alpha_homo, alpha_lumo, auxiliary_basis, auxiliary_basis_role, basis, beta_homo, beta_lumo, charge, connectivity, dipole_moment, dipole_moment_magnitude, dispersion_energy, energies, energy, functional, gap, homo, loewdin_atomic_charges, lumo, mulliken_atomic_charges, multiplicity, positions, reference_energy, scf_energy, solvation_cavity_surface_area, solvation_electrostatic_energy, solvation_model, solvent, spin_square_target, symbols |

## Thermochemical conditions

| Stage | Temperature (K) | Standard state | Entropy model | Frequency scale |
|---|---:|---|---|---:|
| th-cat-dbl-rrho | `298.15` | 1 atm | rrho | `1` |
| th-cat-qrt-rrho | `298.15` | 1 atm | rrho | `1` |
| th-n0-rrho | `298.15` | 1 atm | rrho | `1` |

## Literature constants

| Constant | Value | Unit | Family | Convention |
|---|---:|---|---|---|
| ferrocene_absolute_reduction_potential_acetonitrile_namazian2010 | `4.988` | V | namazian2010_g3mp2rad_cosmors_acetonitrile | absolute reduction potential of Fc+/Fc in acetonitrile at 298.15 K, G3(MP2)-RAD-Full-TZ gas-phase energies with COSMO-RS solvation, electron at rest in the gas phase; the source states a likely accuracy of 0.05-0.1 V |

Claim record: `708bb5f45f72eb1f12d50cf748cbd9e37ca618f53e8f0c5db02110bc57dd8951`

## Host-rendered numerical claims

| Claim | Value | Unit | Source receipt |
|---|---:|---|---|
| hole-on-ni-margin-pme3 | `0.573835` | 1 | `b193f6e8d4fab4372b09f398c3950ac43baee4d4000c390cdd40b048684c28e0` |
| lowd-hole-on-ni-margin-pme3 | `0.547119` | 1 | `b193f6e8d4fab4372b09f398c3950ac43baee4d4000c390cdd40b048684c28e0` |
| lowd-spin-ni-pme3 | `0.787047` | 1 | `b193f6e8d4fab4372b09f398c3950ac43baee4d4000c390cdd40b048684c28e0` |
| lowd-spin-s1-pme3 | `0.12003` | 1 | `b193f6e8d4fab4372b09f398c3950ac43baee4d4000c390cdd40b048684c28e0` |
| lowd-spin-s2-pme3 | `0.119898` | 1 | `b193f6e8d4fab4372b09f398c3950ac43baee4d4000c390cdd40b048684c28e0` |
| spin-ni-pme3 | `0.80934` | 1 | `b193f6e8d4fab4372b09f398c3950ac43baee4d4000c390cdd40b048684c28e0` |
| spin-s-sum-pme3 | `0.235505` | 1 | `b193f6e8d4fab4372b09f398c3950ac43baee4d4000c390cdd40b048684c28e0` |
| spin-s1-pme3 | `0.117821` | 1 | `b193f6e8d4fab4372b09f398c3950ac43baee4d4000c390cdd40b048684c28e0` |
| spin-s2-pme3 | `0.117684` | 1 | `b193f6e8d4fab4372b09f398c3950ac43baee4d4000c390cdd40b048684c28e0` |
| spin-total-pme3 | `0.9999989999999997` | 1 | `b193f6e8d4fab4372b09f398c3950ac43baee4d4000c390cdd40b048684c28e0` |

## Surviving receipts (evidence, not claims)

| Quantity | Value | Unit | Evidence rung | Source receipt |
|---|---:|---|---|---|
| oxidant-margin-fcplus-pme3 | `0.014520857309921031` | hartree e^-1 | derived | `1493aa78454bdb785be0a64e25f0bc9b010a8fd3e1f04fe6b053882ef9a6d15f` |
| oxidant-margin-acfcplus-pme3 | `0.02444317437815932` | hartree e^-1 | derived | `1493aa78454bdb785be0a64e25f0bc9b010a8fd3e1f04fe6b053882ef9a6d15f` |
| oxidant-margin-tbpa-pme3 | `0.04024538304239067` | hartree e^-1 | derived | `1493aa78454bdb785be0a64e25f0bc9b010a8fd3e1f04fe6b053882ef9a6d15f` |
| electronic_energy | `-3304.8362513197017` | hartree | derived | `1f0ccfc8310bfe032eaa523c3b24ebe29b584c11b3d0acfd0201fbb5cc6d233d` |
| zero_point_energy | `0.3074070545932503` | hartree | derived | `1f0ccfc8310bfe032eaa523c3b24ebe29b584c11b3d0acfd0201fbb5cc6d233d` |
| internal_energy | `-3304.5049815862358` | hartree | derived | `1f0ccfc8310bfe032eaa523c3b24ebe29b584c11b3d0acfd0201fbb5cc6d233d` |
| enthalpy | `-3304.5040374016808` | hartree | derived | `1f0ccfc8310bfe032eaa523c3b24ebe29b584c11b3d0acfd0201fbb5cc6d233d` |
| entropy_times_temperature | `0.07721135175544325` | hartree | derived | `1f0ccfc8310bfe032eaa523c3b24ebe29b584c11b3d0acfd0201fbb5cc6d233d` |
| gibbs_free_energy | `-3304.5812487534367` | hartree | derived | `1f0ccfc8310bfe032eaa523c3b24ebe29b584c11b3d0acfd0201fbb5cc6d233d` |
| thermal_internal_energy_correction | `0.33126973346651556` | hartree | derived | `1f0ccfc8310bfe032eaa523c3b24ebe29b584c11b3d0acfd0201fbb5cc6d233d` |
| thermal_enthalpy_correction | `0.3322139180208895` | hartree | derived | `1f0ccfc8310bfe032eaa523c3b24ebe29b584c11b3d0acfd0201fbb5cc6d233d` |
| enthalpy_increment_above_zero_point | `0.024806863427639256` | hartree | derived | `1f0ccfc8310bfe032eaa523c3b24ebe29b584c11b3d0acfd0201fbb5cc6d233d` |
| thermal_gibbs_correction | `0.2550025662651487` | hartree | derived | `1f0ccfc8310bfe032eaa523c3b24ebe29b584c11b3d0acfd0201fbb5cc6d233d` |
| entropy | `0.0002589681427316561` | hartree K^-1 | derived | `1f0ccfc8310bfe032eaa523c3b24ebe29b584c11b3d0acfd0201fbb5cc6d233d` |
| temperature | `298.15` | K | derived | `1f0ccfc8310bfe032eaa523c3b24ebe29b584c11b3d0acfd0201fbb5cc6d233d` |
| pressure | `1.0` | atm | derived | `1f0ccfc8310bfe032eaa523c3b24ebe29b584c11b3d0acfd0201fbb5cc6d233d` |
| near_zero_mode_count | `0` | 1 | derived | `1f0ccfc8310bfe032eaa523c3b24ebe29b584c11b3d0acfd0201fbb5cc6d233d` |
| heat_capacity_cv | `0.00013275930777489225` | hartree K^-1 | derived | `1f0ccfc8310bfe032eaa523c3b24ebe29b584c11b3d0acfd0201fbb5cc6d233d` |
| e-cat-tpssh | `-3306.005113987887` (= `-3306.005113987887` Eh as printed) | hartree | parsed | `23c889a0f74a55c1ef759d183588110141339aa31dad0263d934fe32c5056415` |
| e-n0-tzvp | `-3305.909969834469` (= `-3305.909969834469` Eh as printed) | hartree | parsed | `2ec5600588aba4189c0a32c97ed54466b0fa16b272df7fd7bc52c744283a6f04` |
| charge-n0 | `0` | 1 | parsed | `312ff2eb184c0b60e520bd68eb398f1782843a0a9b5291af39cefd0378f9d365` |
| e-n0-svp-reg | `-3304.836251319702` (= `-3304.836251319702` Eh as printed) | hartree | parsed | `312ff2eb184c0b60e520bd68eb398f1782843a0a9b5291af39cefd0378f9d365` |
| freqs-n0 | `[20.57,39.73,44.69,91.49,113.12,113.25,119.88,128.4,147.98,163.05,164.27,171.65,174.76,174.98,176.0,183.48,185.44,197.89,210.41,213.69,218.16,226.37,230.63,241.25,258.45,259.39,267.56,270.89,279.76,293.87,345.1,361.79,381.74,682.73,684.52,731.18,731.68,741.69,742.69,747.58,748.53,777.95,778.05,837.57,839.93,844.63,845.21,938.15,940.02,941.43,941.8,946.24,949.76,955.41,956.0,956.45,964.41,1265.43,1265.59,1274.35,1275.01,1288.44,1290.12,1324.92,1326.31,1401.28,1401.61,1403.11,1403.24,1408.77,1408.87,1413.72,1413.75,1421.17,1421.48,1425.67,1428.03,1433.97,1435.09,1446.42,1446.8,3040.35,3040.77,3057.75,3057.84,3060.05,3060.09,3064.81,3064.87,3138.65,3138.73,3149.47,3149.56,3172.45,3172.51,3173.34,3173.39,3178.49,3178.53,3185.72,3185.8,3192.82,3192.85,3195.38,3195.43]` | cm^-1 | parsed | `312ff2eb184c0b60e520bd68eb398f1782843a0a9b5291af39cefd0378f9d365` |
| dion-pbe0 | `0.16571781044785894` | hartree | derived | `40a2f58bf0a0803380533be29b1d305a0097f478739a3c09e6826503c4972c94` |
| dion-b3lyp | `0.16097547749086516` | hartree | derived | `40a2f58bf0a0803380533be29b1d305a0097f478739a3c09e6826503c4972c94` |
| dion-tpssh | `0.16031715323288154` | hartree | derived | `40a2f58bf0a0803380533be29b1d305a0097f478739a3c09e6826503c4972c94` |
| functional-spread-v | `0.0054006572149774` | hartree e^-1 | derived | `40a2f58bf0a0803380533be29b1d305a0097f478739a3c09e6826503c4972c94` |
| shift-b3lyp-v | `-0.004742332956993778` | hartree e^-1 | derived | `40a2f58bf0a0803380533be29b1d305a0097f478739a3c09e6826503c4972c94` |
| shift-tpssh-v | `-0.0054006572149774` | hartree e^-1 | derived | `40a2f58bf0a0803380533be29b1d305a0097f478739a3c09e6826503c4972c94` |
| repro-diff-ev | `3.6137901133770356e-05` (= `0.000983362377861444` eV as requested) | hartree | derived | `40a2f58bf0a0803380533be29b1d305a0097f478739a3c09e6826503c4972c94` |
| repro-diff-v | `3.6137901133770356e-05` | hartree e^-1 | derived | `40a2f58bf0a0803380533be29b1d305a0097f478739a3c09e6826503c4972c94` |
| e-cat-b3lyp | `-3305.383092593323` (= `-3305.383092593323` Eh as printed) | hartree | parsed | `45d976da7e14b0ba487aef89ef6eef7fee0cf354e3557d86f6443e50fe1939dc` |
| e-cat-tzvp | `-3305.741363794941` (= `-3305.741363794941` Eh as printed) | hartree | parsed | `4b1f913760eb9a527203720b02f5aa91bcff9cd421b18eecc8f4429d2a9de7a4` |
| n-elec | `1.0` | 1 | derived | `5a01f1888a20d3156b3a34460612827c27840e06be1b507e0dd87715569d1234` |
| e-abs-svp | `0.16593267201642448` | hartree e^-1 | derived | `5a01f1888a20d3156b3a34460612827c27840e06be1b507e0dd87715569d1234` |
| e-couple-svp | `-0.017372948488659223` | hartree e^-1 | derived | `5a01f1888a20d3156b3a34460612827c27840e06be1b507e0dd87715569d1234` |
| e-couple-elec-svp | `-0.01755167215609099` | hartree e^-1 | derived | `5a01f1888a20d3156b3a34460612827c27840e06be1b507e0dd87715569d1234` |
| spread-treatment | `0.00017872366743176826` | hartree e^-1 | derived | `5a01f1888a20d3156b3a34460612827c27840e06be1b507e0dd87715569d1234` |
| basis-shift-v | `0.002852091178738192` | hartree e^-1 | derived | `5a01f1888a20d3156b3a34460612827c27840e06be1b507e0dd87715569d1234` |
| e-couple-tzvp | `-0.014520857309921031` | hartree e^-1 | derived | `5a01f1888a20d3156b3a34460612827c27840e06be1b507e0dd87715569d1234` |
| de-ion-svp | `0.1657539483489927` | hartree | derived | `5a01f1888a20d3156b3a34460612827c27840e06be1b507e0dd87715569d1234` |
| de-ion-svp-ev | `0.1657539483489927` (= `4.510394673587471` eV as requested) | hartree | derived | `5a01f1888a20d3156b3a34460612827c27840e06be1b507e0dd87715569d1234` |
| de-ion-tzvp-ev | `0.1686060395277309` (= `4.588004027628808` eV as requested) | hartree | derived | `5a01f1888a20d3156b3a34460612827c27840e06be1b507e0dd87715569d1234` |
| charge-cat | `1` | 1 | parsed | `6051f6d94baa8339364100e80f92da687090db36f6a7ba632a2744b162c4df99` |
| e-cat-svp-reg | `-3304.670497371353` (= `-3304.670497371353` Eh as printed) | hartree | parsed | `6051f6d94baa8339364100e80f92da687090db36f6a7ba632a2744b162c4df99` |
| freqs-cat-dbl | `[28.3,47.79,52.71,82.29,93.84,111.57,114.71,136.76,137.16,143.58,149.46,158.14,165.2,169.83,171.5,180.31,186.87,196.19,198.15,199.5,205.47,208.19,217.41,218.93,251.65,268.24,272.72,281.08,291.71,308.02,336.35,375.01,406.7,685.29,685.33,734.63,735.77,754.43,754.91,769.92,770.39,793.11,793.21,843.67,847.79,847.97,850.22,934.17,940.44,942.08,946.21,947.78,954.66,954.72,955.92,956.29,962.02,1282.97,1283.0,1285.09,1285.8,1295.98,1299.24,1318.35,1321.37,1398.82,1398.89,1400.47,1400.48,1405.81,1405.89,1410.28,1410.34,1411.17,1411.37,1417.31,1417.71,1431.05,1431.19,1439.23,1439.42,3064.62,3064.64,3066.71,3066.86,3070.36,3070.41,3072.26,3072.31,3176.23,3176.47,3181.37,3181.4,3184.45,3184.46,3185.62,3185.68,3191.86,3192.18,3194.97,3195.01,3207.36,3207.67,3211.91,3211.94]` | cm^-1 | parsed | `6051f6d94baa8339364100e80f92da687090db36f6a7ba632a2744b162c4df99` |
| lowd-spin-cat | `[0.787047,0.12003,-0.003408,0.002565,-0.001372,0.006565,0.119898,-0.003403,0.002556,-0.00137,0.006563,-0.013671,-0.013598,-0.001967,0.001611,-0.001403,-0.001967,0.00163,-0.001409,-0.001174,-9.5e-05,-0.000158,-4.9e-05,-1.8e-05,-0.000112,-0.00099,0.000148,-3e-06,-0.001172,-9.3e-05,-0.000155,-4.6e-05,-1.9e-05,-0.000113,-0.000991,0.000147,-3e-06]` | 1 | parsed | `6051f6d94baa8339364100e80f92da687090db36f6a7ba632a2744b162c4df99` |
| mull-spin-cat | `[0.80934,0.117821,-0.008627,0.003261,-0.001717,0.008403,0.117684,-0.008614,0.003247,-0.001714,0.008402,-0.020774,-0.020684,-0.001556,0.001971,-0.000506,-0.001562,0.001994,-0.000514,-0.001603,2.5e-05,-1.2e-05,-5.8e-05,-7e-05,-0.000154,-0.001436,0.00028,9.5e-05,-0.0016,2.7e-05,-6e-06,-5.4e-05,-7.3e-05,-0.000156,-0.001437,0.000281,9.5e-05]` | 1 | parsed | `6051f6d94baa8339364100e80f92da687090db36f6a7ba632a2744b162c4df99` |
| spinsq-cat-dbl | `0.774059` | 1 | parsed | `6051f6d94baa8339364100e80f92da687090db36f6a7ba632a2744b162c4df99` |
| e-qrt-svp | `-3304.657627836642` (= `-3304.657627836642` Eh as printed) | hartree | parsed | `66985d3ce967fd11e36b11a04f9dc196caa97431e0512091d84b10d685d5d838` |
| freqs-cat-qrt | `[29.35,41.85,53.67,62.41,68.1,78.42,88.97,93.31,104.62,108.73,126.89,133.56,137.31,143.55,146.79,173.51,177.53,183.75,191.72,196.71,202.24,208.34,218.66,226.08,227.16,252.03,254.85,265.6,272.54,317.15,345.61,353.78,382.17,685.47,686.53,721.89,734.29,762.69,763.75,765.14,766.65,789.43,797.8,839.83,841.78,844.85,847.21,911.76,930.56,943.26,948.66,952.65,953.33,954.18,956.37,958.52,959.81,1283.36,1284.55,1288.67,1291.47,1298.11,1302.24,1311.89,1316.89,1403.07,1404.14,1404.64,1406.11,1407.15,1408.15,1410.52,1413.09,1413.99,1414.6,1417.73,1420.88,1422.64,1424.4,1427.51,1428.56,3055.15,3059.8,3060.47,3061.3,3062.72,3062.91,3063.67,3068.05,3166.21,3173.23,3178.08,3178.31,3178.89,3179.39,3180.09,3181.0,3184.28,3186.24,3186.47,3187.54,3189.37,3190.59,3191.26,3197.17]` | cm^-1 | parsed | `66985d3ce967fd11e36b11a04f9dc196caa97431e0512091d84b10d685d5d838` |
| spinsq-cat-qrt | `3.764735` | 1 | parsed | `66985d3ce967fd11e36b11a04f9dc196caa97431e0512091d84b10d685d5d838` |
| abs-basis-shift | `0.002852091178738192` | hartree e^-1 | derived | `6f96fd10e19fd5b30f349c5999d9b6b70febdd5581c2badcb6a7394e002a9b05` |
| u-basis-residual | `0.001426045589369096` | hartree e^-1 | derived | `6f96fd10e19fd5b30f349c5999d9b6b70febdd5581c2badcb6a7394e002a9b05` |
| u-solvation-reference | `0.005512398371243494` (= `0.15` V as requested) | hartree e^-1 | derived | `6f96fd10e19fd5b30f349c5999d9b6b70febdd5581c2badcb6a7394e002a9b05` |
| uncertainty-e-couple-pme3 | `0.007849871863739755` | hartree e^-1 | derived | `6f96fd10e19fd5b30f349c5999d9b6b70febdd5581c2badcb6a7394e002a9b05` |
| electronic_energy | `-3304.657627836642` | hartree | derived | `875e5b25d74b496052b128d86f7e22764a117c2eda3d0255a1c864e39fb5d58a` |
| zero_point_energy | `0.3066639618770147` | hartree | derived | `875e5b25d74b496052b128d86f7e22764a117c2eda3d0255a1c864e39fb5d58a` |
| internal_energy | `-3304.3260982371717` | hartree | derived | `875e5b25d74b496052b128d86f7e22764a117c2eda3d0255a1c864e39fb5d58a` |
| enthalpy | `-3304.325154052617` | hartree | derived | `875e5b25d74b496052b128d86f7e22764a117c2eda3d0255a1c864e39fb5d58a` |
| entropy_times_temperature | `0.08211106472974189` | hartree | derived | `875e5b25d74b496052b128d86f7e22764a117c2eda3d0255a1c864e39fb5d58a` |
| gibbs_free_energy | `-3304.407265117347` | hartree | derived | `875e5b25d74b496052b128d86f7e22764a117c2eda3d0255a1c864e39fb5d58a` |
| thermal_internal_energy_correction | `0.3315295994699664` | hartree | derived | `875e5b25d74b496052b128d86f7e22764a117c2eda3d0255a1c864e39fb5d58a` |
| thermal_enthalpy_correction | `0.3324737840246737` | hartree | derived | `875e5b25d74b496052b128d86f7e22764a117c2eda3d0255a1c864e39fb5d58a` |
| enthalpy_increment_above_zero_point | `0.025809822147658926` | hartree | derived | `875e5b25d74b496052b128d86f7e22764a117c2eda3d0255a1c864e39fb5d58a` |
| thermal_gibbs_correction | `0.25036271929501697` | hartree | derived | `875e5b25d74b496052b128d86f7e22764a117c2eda3d0255a1c864e39fb5d58a` |
| entropy | `0.00027540186057267116` | hartree K^-1 | derived | `875e5b25d74b496052b128d86f7e22764a117c2eda3d0255a1c864e39fb5d58a` |
| temperature | `298.15` | K | derived | `875e5b25d74b496052b128d86f7e22764a117c2eda3d0255a1c864e39fb5d58a` |
| pressure | `1.0` | atm | derived | `875e5b25d74b496052b128d86f7e22764a117c2eda3d0255a1c864e39fb5d58a` |
| near_zero_mode_count | `0` | 1 | derived | `875e5b25d74b496052b128d86f7e22764a117c2eda3d0255a1c864e39fb5d58a` |
| heat_capacity_cv | `0.0001335521099667706` | hartree K^-1 | derived | `875e5b25d74b496052b128d86f7e22764a117c2eda3d0255a1c864e39fb5d58a` |
| spin-ni-pme3 | `0.80934` | 1 | derived | `b193f6e8d4fab4372b09f398c3950ac43baee4d4000c390cdd40b048684c28e0` |
| spin-s1-pme3 | `0.117821` | 1 | derived | `b193f6e8d4fab4372b09f398c3950ac43baee4d4000c390cdd40b048684c28e0` |
| spin-s2-pme3 | `0.117684` | 1 | derived | `b193f6e8d4fab4372b09f398c3950ac43baee4d4000c390cdd40b048684c28e0` |
| spin-s-sum-pme3 | `0.235505` | 1 | derived | `b193f6e8d4fab4372b09f398c3950ac43baee4d4000c390cdd40b048684c28e0` |
| spin-total-pme3 | `0.9999989999999997` | 1 | derived | `b193f6e8d4fab4372b09f398c3950ac43baee4d4000c390cdd40b048684c28e0` |
| hole-on-ni-margin-pme3 | `0.573835` | 1 | derived | `b193f6e8d4fab4372b09f398c3950ac43baee4d4000c390cdd40b048684c28e0` |
| lowd-spin-ni-pme3 | `0.787047` | 1 | derived | `b193f6e8d4fab4372b09f398c3950ac43baee4d4000c390cdd40b048684c28e0` |
| lowd-spin-s1-pme3 | `0.12003` | 1 | derived | `b193f6e8d4fab4372b09f398c3950ac43baee4d4000c390cdd40b048684c28e0` |
| lowd-spin-s2-pme3 | `0.119898` | 1 | derived | `b193f6e8d4fab4372b09f398c3950ac43baee4d4000c390cdd40b048684c28e0` |
| lowd-spin-s-sum-pme3 | `0.239928` | 1 | derived | `b193f6e8d4fab4372b09f398c3950ac43baee4d4000c390cdd40b048684c28e0` |
| lowd-hole-on-ni-margin-pme3 | `0.547119` | 1 | derived | `b193f6e8d4fab4372b09f398c3950ac43baee4d4000c390cdd40b048684c28e0` |
| e-n0-tpssh | `-3306.16543114112` (= `-3306.16543114112` Eh as printed) | hartree | parsed | `bfca6f36f3cf3897c34fb69bdf24474ce4761079c0b94391829394173bd1d859` |
| e-n0-svp | `-3304.836279223488` (= `-3304.836279223488` Eh as printed) | hartree | parsed | `c283df95263d888c17939acdb0ca8a4405d7092014906eaf9279c4852f7431ad` |
| e-n0-b3lyp | `-3305.544068070814` (= `-3305.544068070814` Eh as printed) | hartree | parsed | `c9282792ceea9f73b31e63d3bc705e9851b899dc54816fc173b751398179516b` |
| dg-qrt-dbl | `0.008050964073390787` | hartree | derived | `cea600268e4a490b84a1c04ed7863648c7648933d8c789ec881e947de9ca10f8` |
| gap-qrt-dbl-pme3 | `0.008050964073390787` (= `0.21907789126934887` eV as requested) | hartree | derived | `cea600268e4a490b84a1c04ed7863648c7648933d8c789ec881e947de9ca10f8` |
| electronic_energy | `-3304.6704973713527` | hartree | derived | `ee33b7c88bd74cd5b1908a9318c117ab7bce358b87537a02555e3959e183c739` |
| zero_point_energy | `0.30833688870978565` | hartree | derived | `ee33b7c88bd74cd5b1908a9318c117ab7bce358b87537a02555e3959e183c739` |
| internal_energy | `-3304.3381460534756` | hartree | derived | `ee33b7c88bd74cd5b1908a9318c117ab7bce358b87537a02555e3959e183c739` |
| enthalpy | `-3304.337201868921` | hartree | derived | `ee33b7c88bd74cd5b1908a9318c117ab7bce358b87537a02555e3959e183c739` |
| entropy_times_temperature | `0.07811421249931648` | hartree | derived | `ee33b7c88bd74cd5b1908a9318c117ab7bce358b87537a02555e3959e183c739` |
| gibbs_free_energy | `-3304.4153160814203` | hartree | derived | `ee33b7c88bd74cd5b1908a9318c117ab7bce358b87537a02555e3959e183c739` |
| thermal_internal_energy_correction | `0.33235131787696354` | hartree | derived | `ee33b7c88bd74cd5b1908a9318c117ab7bce358b87537a02555e3959e183c739` |
| thermal_enthalpy_correction | `0.3332955024317726` | hartree | derived | `ee33b7c88bd74cd5b1908a9318c117ab7bce358b87537a02555e3959e183c739` |
| enthalpy_increment_above_zero_point | `0.024958613721986912` | hartree | derived | `ee33b7c88bd74cd5b1908a9318c117ab7bce358b87537a02555e3959e183c739` |
| thermal_gibbs_correction | `0.2551812899324843` | hartree | derived | `ee33b7c88bd74cd5b1908a9318c117ab7bce358b87537a02555e3959e183c739` |
| entropy | `0.00026199635250483475` | hartree K^-1 | derived | `ee33b7c88bd74cd5b1908a9318c117ab7bce358b87537a02555e3959e183c739` |
| temperature | `298.15` | K | derived | `ee33b7c88bd74cd5b1908a9318c117ab7bce358b87537a02555e3959e183c739` |
| pressure | `1.0` | atm | derived | `ee33b7c88bd74cd5b1908a9318c117ab7bce358b87537a02555e3959e183c739` |
| near_zero_mode_count | `0` | 1 | derived | `ee33b7c88bd74cd5b1908a9318c117ab7bce358b87537a02555e3959e183c739` |
| heat_capacity_cv | `0.00013254266460511536` | hartree K^-1 | derived | `ee33b7c88bd74cd5b1908a9318c117ab7bce358b87537a02555e3959e183c739` |
| e-cat-svp | `-3304.67056141304` (= `-3304.67056141304` Eh as printed) | hartree | parsed | `f597f41e758ad8cdc7d9a5ef24daa99758f6da3160a817597b630b765eb8f338` |

## Expected versus delivered

| Observable | Expected | Delivered | Unit | Agreement | Basis |
|---|---|---:|---|---|---|
| basis-shift-pme3-v | -0.05..0.05 | `` |  | not_comparable | Both partners are the same molecule differing by one electron on a metal-dominated orbital; the def2-SVP to def2-TZVP differential on such a difference is usually a few tens of meV, and the prior goal named exactly this as the unmeasured component of its 0.25 V asserted uncertainty. |
| e-couple-vs-fc-pme3 | negative -0.75..-0.2 | `` |  | not_comparable | This workspace has already delivered -0.4727 V for the same couple at PBE0-D3BJ/def2-SVP + CPCM(MeCN) from Gibbs free energies of the validated neutral and doublet-cation minima; a thiolate-rich Ni(II)/Ni(III) pair with two alkyl phosphines sits well below Fc+/Fc. |
| gap-qrt-dbl-pme3 | positive 0.1..0.4 | `` |  | not_comparable | This workspace delivered 0.2191 eV from the two validated def2-SVP/CPCM minima; a strong-field S2P2 donor set should keep the low-spin doublet below the quartet, but the quartet is not far enough away to dismiss without computing it. |
| spin-ni-pme3 | 0.6..0.95 | `0.80934` | 1 | agreed | A d7 Ni(III) centre in a square-planar S2P2 field puts most of the single hole in a metal-centred orbital, with delocalisation onto the two thiolates; this workspace has already delivered 0.809 (Mulliken) and 0.787 (Loewdin) for the same result. |
| spin-s1-pme3 | 0.05..0.25 | `0.117821` | 1 | agreed | Thiolate S 3p orbitals mix strongly with the Ni d manifold, so each sulfur carries a fraction of the hole; the same result previously gave 0.118. |
| spin-s2-pme3 | 0.05..0.25 | `0.117684` | 1 | agreed | The two thiolates are symmetry-equivalent in the trans square-planar arrangement, so the second sulfur should match the first to within relaxation asymmetry; the same result previously gave 0.118. |

An expectation is displayed, never scored: a diverging row settles nothing and means the chemistry disagreed with the reasoning, which is a result the reader owns.

Scientific decision: not recorded -- interpretation is a session act; this run executed extraction, thermochemistry, expressions, validation verdicts, and claim rendering only.

Recovery: the engine outputs and every receipt above are unchanged and remain readable. A later explicit analysis request -- a new session over this workspace planning an analysis-only toolchain on the registered results -- may re-extract, validate, and claim them without re-running any engine.

## Record delivery

| Record | Nodes (reached state) | Calculation | Analysis |
| --- | --- | --- | --- |
| 1 | n0-sp-tzvp (validated) | validated | partial |
| 2 | n0-sp-svp-pbe0 (validated) | validated | executed |
| 3 | cat-sp-svp-pbe0 (validated) | validated | executed |
| 4 | n0-sp-svp-b3lyp (validated) | validated | executed |
| 5 | cat-sp-svp-b3lyp (validated) | validated | executed |
| 6 | n0-sp-svp-tpssh (validated) | validated | executed |
| 7 | cat-sp-svp-tpssh (validated) | validated | executed |
| shared | claims-redox (failed), claims-spin (executed), claims-verdicts (skipped), ex-cat-svp-reg (executed), ex-cat-tzvp-reg (executed), ex-n0-svp-reg (executed), ex-qrt-svp-reg (executed), expr-functional-spread (executed), expr-quartet-gap (executed), expr-spin-populations (executed), expr-uncertainty (executed), th-cat-dbl-rrho (executed), th-cat-qrt-rrho (executed), th-n0-rrho (executed), val-svp-replication (failed) | - | partial |

A batch of N is N observations; each verdict above is one record's, and no aggregate quantity is rendered.
