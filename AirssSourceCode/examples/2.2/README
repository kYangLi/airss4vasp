In this example we explore hydrogen at 100 GPa. First we perform a free search

host:2.2 cjp10$ cat H.cell
#VARVOL=2.5
#SPECIES=H
#NATOM=8

#SLACK=0.25
#OVERLAP=0.1
#MINSEP=1 H-H=0.7
#COMPACT

KPOINTS_MP_SPACING 0.07

SYMMETRY_GENERATE
SNAP_TO_SYMMETRY

%BLOCK SPECIES_POT
QC5
%ENDBLOCK SPECIES_POT
host:2.2 cjp10$ airss.pl -press 100 -max 10 -seed H
host:2.2 cjp10$ ca -r
H-99432-8459-3       100.01     2.301     -13.670   8 H            P21/c      1
H-99432-8459-10      100.04     2.301       0.001   8 H            Pnma       1
H-99432-8459-6       100.00     2.305       0.002   8 H            C2         1
H-99432-8459-9       100.02     2.309       0.013   8 H            P1         1
H-99432-8459-5        99.95     2.334       0.015   8 H            Cmce       1
H-99432-8459-1       100.01     2.335       0.016   8 H            I41/acd    1
H-99432-8459-8        99.99     2.336       0.016   8 H            I41/acd    1
H-99432-8459-7        99.97     2.341       0.022   8 H            C2/c       1
H-99432-8459-4        99.98     2.243       0.028   8 H            Fmmm       1
H-99432-8459-2       100.00     2.172       0.098   8 H            P-1        1
host:2.2 cjp10$

Then we extract the molecular units that emerge, using the cryan code directly to analyse a single structure:

host:2.2 cjp10$ cryan -g -bs 1 < H-99432-8459-3.res > H2.cell

Degrees of freedom:    14  0.737

Number of units:        4

host:2.2 cjp10$ cat H2.cell
%BLOCK LATTICE_CART
   1.89795   0.00000   0.00000
  -0.00132   2.97936   0.00000
   0.00189   0.00094   3.25581
%ENDBLOCK LATTICE_CART

#TARGVOL=4.60

%BLOCK POSITIONS_ABS
    H   1.43077   1.95640   3.12157 # 1-D(inf)h % NUM=1
    H   0.91703   2.16131   2.62420 # 1-D(inf)h
%ENDBLOCK POSITIONS_ABS

#SYMMOPS=1
##SGRANK=20
#NFORM=4
#SLACK=0.25
#OVERLAP=0.1
#COMPACT
#MINSEP=1.0 H-H=1.44

host:2.2 cjp10$

Add the lines "KPOINTS_MP_SPACING 0.07" and below from H.cell to your new H2.cell, and copy H.param to H2.param.

host:2.2 cjp10$ airss.pl -press 100 -max 10 -seed H2
host:2.2 cjp10$ ca -r
H-99432-8459-3       100.01     2.301     -13.670   8 H            P21/c      1
H2-3639-9598-1       100.01     2.301       0.000   8 H            P21/c      1
H2-3639-9598-9        99.99     2.305       0.000   8 H            Pca21      1
H-99432-8459-10      100.04     2.301       0.001   8 H            Pnma       1
H2-3639-9598-7       100.04     2.304       0.002   8 H            C2         1
H-99432-8459-6       100.00     2.305       0.002   8 H            C2         1
H2-3639-9598-3        99.93     2.306       0.004   8 H            P-1        1
H-99432-8459-9       100.02     2.309       0.013   8 H            P1         1
H-99432-8459-5        99.95     2.334       0.015   8 H            Cmce       1
H2-3639-9598-8        99.98     2.334       0.016   8 H            Cmce       1
H2-3639-9598-2       100.02     2.334       0.016   8 H            Cmce       1
H2-3639-9598-5       100.01     2.335       0.016   8 H            Pnnm       1
H-99432-8459-1       100.01     2.335       0.016   8 H            I41/acd    1
H-99432-8459-8        99.99     2.336       0.016   8 H            I41/acd    1
H2-3639-9598-10       99.99     2.341       0.020   8 H            Cc         1
H-99432-8459-7        99.97     2.341       0.022   8 H            C2/c       1
H2-3639-9598-4       100.03     2.263       0.025   8 H            C2/m       1
H-99432-8459-4        99.98     2.243       0.028   8 H            Fmmm       1
H2-3639-9598-6        99.96     2.193       0.097   8 H            P-1        1
H-99432-8459-2       100.00     2.172       0.098   8 H            P-1        1
host:2.2 cjp10$
