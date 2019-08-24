In this example we will explore the energy landscape of Carbon at 100GPa (or 1MBar). For comparison, the pressure at the centre of the Earth is about 350GPa.

host:4.1 cjp10$ cat C2.cell 
%BLOCK LATTICE_CART
1.709975 0 0
0 1.709975 0
0 0 1.709975
%ENDBLOCK LATTICE_CART

#VARVOL=5

%BLOCK POSITIONS_FRAC
C 0.0 0.0 0.0 # C1 % NUM=1
C 0.0 0.0 0.0 # C2 % NUM=1
%ENDBLOCK POSITIONS_FRAC

#MINSEP=1.3

The only AIRSS related command in the cell file is "#MINSEP=1.3". This is not essential for the light elements, but for transition metals it is important to avoid core overlap to prevent poor convergence of the electronic structure. 

The search is a repeat of the performed in example 2.1, except this time the search will be performed using the VASP code. First create a C2.POTCAR files. The airss.pl script is used in the usual way, except the -vasp flag it set.

host:2.1 cjp10$ airss.pl -vasp -press 100 -max 20 -seed C2

The resuts of the search can be ranked using the "ca" wrapper to the cryan code.

host:2.1 cjp10$ ca -r
C2-13385-9319-6        100.00     4.545      -6.999   2 C            Fd-3m      1
C2-13385-9319-10        99.99     5.070       0.862   2 C            C2/m       1
C2-13385-9319-7        100.02     4.990       0.902   2 C            C2/m       1
C2-13385-9319-5        100.01     4.565       2.196   2 C            P2/m       1
C2-13385-9319-3        100.04     4.270       2.386   2 C            P2/m       1
C2-13385-9319-9         99.97     4.245       2.447   2 C            Pm-3m      1
C2-13385-9319-8        100.01     4.555       2.653   2 C            P-1        1
C2-13385-9319-1         99.92     4.500       2.789   2 C            P-1        1
C2-13385-9319-2        100.02     4.575       3.259   2 C            Cmcm       1
C2-13385-9319-4        100.01     4.500       3.653   2 C            P-1        1
