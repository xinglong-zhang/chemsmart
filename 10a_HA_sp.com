%chk=10a_HA_sp.chk
%nprocshared=12
%mem=40GB
# m062x def2svp scrf=(smd,solvent=dichloroethane)

Gaussian pKa calculation job

0 1
H        0.0000000000    0.0000000000   -1.5742100000
C        0.0000000000    0.0000000000   -0.5021190000
N        0.0000000000    0.0000000000    0.6552740000

