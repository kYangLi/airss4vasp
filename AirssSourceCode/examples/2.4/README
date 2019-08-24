We can analyse the Pmmn gamma-B structure identified in example 2.3 using a modularity detection algorithm as
applied to the complex network of atomic contacts.

See the following for more details:

Revealing and exploiting hierarchical material structure through complex atomic networks
SE Ahnert, WP Grant, CJ Pickard
npj Computational Materials 3 (1), 35, 2017

http://www.nature.com/articles/s41524-017-0035-x

First we perform the modular decomposition.

host:2.4 cjp20$ cryan -g -bs -1.1 < ../2.3/B-69670-3073-1.res
Iteration:   1  0
 }   0.3260869565    2
Maximum q:   2 0.32609
Iteration:   1  0
 }   0.3260869565    2
Iteration:   2  1
 }   0.3260869565    2
Iteration:   3  2
 }   0.3260869565    2

Degrees of freedom:    12  0.152

Number of units:        2

%BLOCK LATTICE_CART
   5.05440   0.00000   0.00000
   0.00000   5.61990   0.00000
   0.00000   0.00000   6.98730
%ENDBLOCK LATTICE_CART

#TARGVOL=99.24

%BLOCK POSITIONS_ABS
    B   3.91409  -0.15113   0.95482 # 1-D2h % NUM=1
    B   3.91410  -0.15061   2.71893 # 1-D2h
    B   2.28386   2.31961   0.95492 # 1-D2h
    B   1.39268   1.04217   1.83717 # 1-D2h
    B   2.31622   0.56462   3.30280 # 1-D2h
    B   4.76984  -1.62220   1.83839 # 1-D2h
    B   2.31571   0.56411   0.37147 # 1-D2h
    B   3.88218   1.60450   0.37142 # 1-D2h
    B   4.80548   1.12636   1.83718 # 1-D2h
    B   3.88189   1.60431   3.30284 # 1-D2h
    B   3.85415   2.60738   1.83745 # 1-D2h
    B   2.28415   2.31959   2.71899 # 1-D2h
    B   2.34405  -0.43858   1.83738 # 1-D2h
    B   1.42838   3.79048   1.83843 # 1-D2h
%ENDBLOCK POSITIONS_ABS

#SYMMOPS=1
##SGRANK=20
#NFORM=2
#SLACK=0.25
#OVERLAP=0.1
#COMPACT
#MINSEP=1.0 B-B=1.64

Using this decompositon of the structure, we can rapidly rediscover it by randomly placing
these D2h units into the same fixed cell.

host:2.4 cjp20$ cat B.cell
%BLOCK LATTICE_ABC
5.0544 5.6199 6.9873
90 90 90
#FIX
%ENDBLOCK LATTICE_ABC

%BLOCK POSITIONS_ABS
    B   2.58544   4.77680   4.32147 # 1-D2h % NUM=2
    B   1.69796   8.00288   2.85360 # 1-D2h
    B   2.61373   3.77359   2.85525 # 1-D2h
    B   4.18411   4.06190   3.73764 # 1-D2h
    B   5.07510   5.33975   2.85470 # 1-D2h
    B   4.18329   4.06211   1.97364 # 1-D2h
    B   2.55309   6.53193   1.97357 # 1-D2h
    B   2.58485   4.77663   1.39006 # 1-D2h
    B   2.55305   6.53194   3.73757 # 1-D2h
    B   5.03988   2.59079   2.85345 # 1-D2h
    B   4.15056   5.81724   1.39010 # 1-D2h
    B   1.66210   5.25490   2.85543 # 1-D2h
    B   4.12347   6.82062   2.85477 # 1-D2h
    B   4.15140   5.81769   4.32009 # 1-D2h
%ENDBLOCK POSITIONS_ABS

#SLACK=0.25
#OVERLAP=0.1
#COMPACT
#MINSEP=1.64

KPOINTS_MP_SPACING : 0.1

FIX_ALL_CELL : true

%BLOCK SPECIES_POT
B 2|1.6|7|7|9|20:21(qc=4)
%ENDBLOCK SPECIES_POT

host:2.4 cjp20$ spawn airss.pl -mpinp 4 -steps 0 -seed B
host:2.4 cjp20$ despawn
host:2.4 cjp20$ tidy.pl
Files will be removed - <ENTER> to continue
host:2.4 cjp20$ ca -r -t 10
B-31951-340-1          -22.43     7.088     -79.005  28 B            Pnnm       1
B-50523-9074-1         -22.43     7.088       0.000  28 B            Pnnm       1
B-55899-7485-2         -22.45     7.088       0.000  28 B            Pnnm       1
B-48183-3899-1         -19.27     7.088       0.131  28 B            P1         1
B-77413-4166-1         -19.68     7.088       0.133  28 B            P1         1
B-48179-3887-2         -21.96     7.088       0.144  28 B            P1         1
B-74323-1468-1         -18.67     7.088       0.165  28 B            P1         1
B-55889-7368-1         -19.42     7.088       0.166  28 B            P1         1
B-53810-2553-1         -19.66     7.088       0.173  28 B            P1         1
B-31957-225-1          -17.15     7.088       0.176  28 B            P1         1
host:2.4 cjp20$ ca -s
B-31951-340-1          -22.43     7.088     -79.005  28 B            Pnnm      1   187      ~
Number of structures   :    187
Number of compositions :      1

This time the gamma-B structure was located three times out of 187 attempts.
