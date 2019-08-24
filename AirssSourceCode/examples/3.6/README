In this example we search for small SiO2 clusters using a Vashishta potential.

host:3.6 cjp10$ cat SiO2.cell
%BLOCK LATTICE_CART
20 0 0
0 20 0
0 0 20
#FIX
%ENDBLOCK LATTICE_CART

#CLUSTER
#ELLIPSOID=2.5 0.75

#SPECIES=Si%NUM=1,O%NUM=2
#NFORM=4

#MINSEP=1.0 Si-Si=3.00 Si-O=1.60 O-O=2.58

A random cluster of 4SiO2 is built, subject to the specified minimum distances, and confined within an ellipsoid of volume 4/3.Pi.r**3, where r=2.5. The second number, 0.75 in this case, describes how spherical the ellipsoids will be on average. A value of 0.0 enforces sphericity.

host:3.6 cjp10$ airss.pl -gulp -cluster -max 10 -seed SiO2

The -cluster setting enforces a constant volume optimisation, and the point group analysis of the resulting structures.

host:3.6 cjp10$ ca -r
SiO2-33704-5555-6        0.00  2000.000     -28.907   4 SiO2         Cs         1
SiO2-33704-5555-10       0.00  2000.000       0.002   4 SiO2         Cs         1
SiO2-33704-5555-7        0.00  2000.000       0.051   4 SiO2         C1         1
SiO2-33704-5555-5        0.00  2000.000       0.100   4 SiO2         C2v        1
SiO2-33704-5555-4        0.00  2000.000       0.100   4 SiO2         C2v        1
SiO2-33704-5555-8        0.00  2000.000       1.880   4 SiO2         C1         1
SiO2-33704-5555-3        0.00  2000.000       2.698   4 SiO2         C1         1
SiO2-33704-5555-9        0.00  2000.000       4.079   4 SiO2         C1         1
SiO2-33704-5555-2        0.00  2000.000       4.595   4 SiO2         C1         1
SiO2-33704-5555-1        0.00  2000.000      11.250   4 SiO2         C1         1
