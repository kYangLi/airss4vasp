In this example we search for well packed methane (CH4) molecular crystals. The Tersoff potential is used to determine the interatomic interactions, and it is expected to provide only a coarse approximation.

host:3.7 cjp10$ cat CH4.cell
#TARGVOL=7.69

%BLOCK POSITIONS_ABS
    C   0.86267   0.87843   0.69644 # 1-Td % NUM=1
    H  -0.15687   1.20067   0.63115 # 1-Td
    H   1.19106   0.93720   1.71422 # 1-Td
    H   0.93994  -0.13355   0.35439 # 1-Td
    H   1.47618   1.50890   0.08615 # 1-Td
%ENDBLOCK POSITIONS_ABS

#SYMMOPS=4
#SGRANK=20
#SLACK=0.25
#OVERLAP=1
##SYSTEM=Cubi
#MINSEP=1.0 C-C=2.10 C-H=1.45 H-H=1.00

Using this input file, buildcell will attempt to construct a random unit cell, choosing a space group randomly with four symmetry operators from the top 20 most common space groups for molecular crystals. The CH4 tetrahedral unit (or molecule) is specified in absolute coordinates. All atoms are labelled in the same way (1-Td) to tie the atoms together into a unit. The miniumum distance constraints are first applied by shifting the molecules, using a SLACK parameter of 0.25 (the MINSEP distances are multiplied by 0.75). Once these relaxed constaints are satisfied an optimisation using hard sphere potentials for the specified MINSEP distances is applied, and rotations of the unit are permitted. A final OVERLAP of 1 is tolerated. This may be reduced to 0.1 or lower, at a greater computational cost.

host:3.7 cjp10$ airss.pl -gulp -press 1 -max 10 -seed CH4

A pressure of 1 GPa is applied to favour well packed solutions.

host:3.7 cjp10$ ca -r
CH4-32238-381-10         1.00     8.111     -20.738   4 CH4          P212121    1
CH4-32238-381-8          1.00     8.390       0.002   4 CH4          P212121    1
CH4-32238-381-5          1.00     8.516       0.003   4 CH4          P21        1
CH4-32238-381-9          1.00     8.606       0.003   4 CH4          Pca21      1
CH4-32238-381-2          1.00     8.538       0.003   4 CH4          Pca21      1
CH4-32238-381-3          1.00     8.660       0.003   4 CH4          Pca21      1
CH4-32238-381-4          1.00     8.881       0.005   4 CH4          P21212     1
CH4-32238-381-1          1.00     8.902       0.005   4 CH4          P21        1
CH4-32238-381-6          1.00     9.625       0.009   4 CH4          P21212     1
CH4-32238-381-7          1.00     8.980       0.730   4 CH4          P1         1

The ten samples above are not exhuastive. A low energy solution is known to be cubic with space group P213, and cubic solutions can be encouraged by uncommenting the SYSTEM=Cubi line in the CH4.cell. 