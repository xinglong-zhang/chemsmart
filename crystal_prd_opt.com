%chk=crystal_prd_opt.chk
%nprocshared=12
%mem=40GB
# opt freq m062x def2svp scrf=(smd,solvent=dichloroethane)

Gaussian job with default settings

0 1
C        0.0000000000    0.0000000000    0.0000000000
H        1.0000000000    0.0000000000    0.0000000000
H        0.0000000000    1.0000000000    0.0000000000
H        0.0000000000    0.0000000000    1.0000000000
H       -1.0000000000    0.0000000000    0.0000000000

