%chk=ammonia_fake.chk
%nprocshared=12
%mem=40GB
# opt freq m062x def2svp scrf=(smd,solvent=dichloroethane)

Gaussian job with default settings

0 1
N        0.0000000000    0.0000000000    0.0000000000
H       -0.4417000000    0.2906000000    0.8711000000
H        0.7256000000    0.6896000000   -0.1907000000
H        0.4875000000   -0.8701000000    0.2089000000

