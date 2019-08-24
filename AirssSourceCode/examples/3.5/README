In this example we will perform a search for 2D silica structures using GULP and the Vashishta potential. First take a look at the input file for buildcell:

host:3.5 cjp10$ cat SiO2.cell
%BLOCK LATTICE_CART
1 0 0
0 1 0
0 0 15
#CFIX
%ENDBLOCK LATTICE_CART

#TARGVOL=50

#SLAB
#WIDTH=4.0
#SHIFT=7.5

#SPECIES=Si%NUM=1,O%NUM=2
#NFORM=4

#MINSEP=1.0 Si-Si=3.00 Si-O=1.60 O-O=2.58

#COMPACT

The lattice vectors are specified so as to provide a fixed length for the c-axis (15 Ang). The a- and b-axis lengths (and directions) will be chosen by buildcell so that the overall unit cell has a volume of 50 Ang^3 (set by TARGVOL). A larger value for TARGVOL will lead to thinner layers (for a given number of atoms).

The SLAB setting instructs buildcell to use internal settings appropriate for this geometry. Random structures are built of 4 formula units of SiO2, so that the MINSEP distances are satisfied. They are confined by hard wall planes perpendicular to the c-axis, and separated by 4 Ang. The centre of the resulting slab is shifted by 7.5 Ang, so that it sits in the centre of the unit cell.

host:3.5 cjp10$ airss.pl -gulp -slab -max 100 -seed SiO2

The -slab flag is set to force gulp to perform a fixed volume calculation. If you were using CASTEP, the cell constraints would be set through the seed.cell file.

host:3.5 cjp10$  ca -r -t 10
SiO2-1291-404-74         0.00    75.000     -31.962   4 SiO2         P-1        1
SiO2-1291-404-49         0.00    74.996       0.637   4 SiO2         P1         1
SiO2-1291-404-48         0.00    75.001       0.742   4 SiO2         Cm         1
SiO2-1291-404-91         0.00    75.005       0.789   4 SiO2         Cm         1
SiO2-1291-404-9          0.00    74.996       0.817   4 SiO2         P1         1
SiO2-1291-404-72         0.00    75.012       0.831   4 SiO2         Pc         1
SiO2-1291-404-96         0.00    74.995       0.868   4 SiO2         Cm         1
SiO2-1291-404-59         0.00    75.003       0.956   4 SiO2         P-1        1
SiO2-1291-404-54         0.00    75.002       0.965   4 SiO2         P1         1
SiO2-1291-404-83         0.00    75.007       0.971   4 SiO2         P-1        1

Visualise your results. The lower energy structures that you find should be some variant of bilayer hexagonal silica. It is clear that 100 samples is not exhaustive, since the are no repeated low energy structures.
