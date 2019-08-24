In this example we search for stable crystalline configurations in the well known Kob-Anderson Binary Lennard-Jones system, at the A:B 80:20 composition.

host:1.8 cjp10$ cat AB.cell
#VARVOL=15
#SPECIES=A%NUM=4,B%NUM=1
#NFORM=2
#MINSEP=1.5

host:1.8 cjp10$ cat AB.pp
2 12 6 5
A B
# Epsilon
1.00 1.50
0.50
# Sigma
2.00 1.60
1.76
host:1.8 cjp10$ airss.pl -pp3 -max 100 -seed AB
host:1.8 cjp10$ ca -r -t 10
AB-85154-3384-55       0.00    30.784     -45.812   2 BA4          C2/m       1
AB-85154-3384-43       0.00    30.784       0.000   2 BA4          C2/m       1
AB-85154-3384-13      -0.00    30.836       0.045   2 BA4          P21/m      1
AB-85154-3384-64      -0.00    30.836       0.045   2 BA4          P21/m      1
AB-85154-3384-96       0.00    30.836       0.045   2 BA4          P21/m      1
AB-85154-3384-99       0.00    30.836       0.045   2 BA4          P21/m      1
AB-85154-3384-77      -0.00    30.748       0.361   2 BA4          I4/mmm     1
AB-85154-3384-14      -0.00    30.748       0.361   2 BA4          I4/mmm     1
AB-85154-3384-56      -0.00    30.748       0.361   2 BA4          I4/mmm     1
AB-85154-3384-17       0.00    31.273       0.716   2 BA4          Cmmm       1
host:1.8 cjp10$
