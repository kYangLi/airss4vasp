In this example we use coordination contraints to rapidly locate the icosahedral C20 fullerene.

host:3.3 cjp10$ cat C.cell
%BLOCK LATTICE_CART
20 0 0
0 20 0
0 0 20
#FIX
%ENDBLOCK LATTICE_CART

%BLOCK POSITIONS_FRAC
C 0.0 0.0 0.0 # C1 % NUM=20
%ENDBLOCK POSITIONS_FRAC

#SYMMOPS=1
#NFORM=1

#SPHERE=2.1
#POSAMP=2.5

#CLUSTER

#MINSEP=1.44
#COORD=3
#MINBANGLE=91

FIX_ALL_CELL : true
host:3.3 cjp10$

In the above, 20 carbon atoms are initially randomly placed in a sphere of radius 2.5 Ang centered on the origin.
They are subsequently confined in a sphere of radius 2.1 Ang, such that the minimum separation is 1.44 Ang, and
all the atoms are connected to three others, with a minimum bond angle of 91 degrees (so as to exclude 4-fold rings).

host:3.3 cjp10$ airss.pl -gulp -max 10 -cluster -seed C
host:3.3 cjp10$ ca -r
C-68766-2151-4           0.00   400.000      -5.779  20 C            Ih         1
C-68766-2151-10          0.00   400.000       0.000  20 C            Ih         1
C-68766-2151-3           0.00   400.000       0.072  20 C            C2v        1
C-68766-2151-1           0.00   400.000       0.072  20 C            C2v        1
C-68766-2151-7           0.00   400.000       0.072  20 C            C2v        1
C-68766-2151-2           0.00   400.000       0.072  20 C            C2v        1
C-68766-2151-8           0.00   400.000       0.072  20 C            C2v        1
C-68766-2151-5           0.00   400.000       0.072  20 C            C2v        1
C-68766-2151-6           0.00   400.000       0.072  20 C            C2v        1
C-68766-2151-9           0.00   400.000       0.109  20 C            Cs         1
host:3.3 cjp10$
