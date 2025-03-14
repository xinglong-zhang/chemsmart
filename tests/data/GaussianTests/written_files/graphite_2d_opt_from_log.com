%chk=graphite_2d_opt_from_log.chk
%nprocshared=64
%mem=400GB
# opt freq b3lyp empiricaldispersion=gd3bj def2svp

Job prepared from Gaussian file graphite_2d_opt.log

0 1
C       -0.0017240000   -0.7146210000   -0.0000000000
C        0.0017240000    0.7146210000    0.0000000000
TV       2.4755330000   -0.0060580000   -0.0000000000
TV      -1.2326020000    2.1468790000   -0.0000000000
