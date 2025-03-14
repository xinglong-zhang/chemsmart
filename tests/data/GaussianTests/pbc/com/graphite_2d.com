%chk=graphite_2d.chk
%nprocshared=24
%mem=60GB
# PBEPBE/6-31g(d,p)/Auto SCF=Tight

graphite_2d

0 1
C                  0.000000    0.000000    0.000000
C                  0.000000    1.429118    0.000000
TV                 2.475315    0.000000    0.000000
TV                -1.219952    2.133447    0.000000
