In this example gulp and a Tersoff potential are used to relax structures generated with coordination constraints.
Random 'sensible' structures with 8 Carbon atoms are constructed, such that each atom is four-fold coordinated, and
the minimum bondangle encountered is 91 degrees. This choice excludes the generation of structure with 4-rings.

host:3.2 cjp10$ cat C.cell
#TARGVOL=5.6
#SPECIES=C
#NATOM=8

#SYMMOPS=1
#NFORM=1

#COMPACT
#MINSEP=1.5
#COORD=4
#MINBANGLE=91
#MAXBANGLE=180

Run the airss.pl, using gulp. Similar structures are unified on-the-fly.

host:3.2 cjp10$ airss.pl -gulp -max 100 -sim 0.1 -seed C
host:3.2 cjp10$ ca -r
C-17050-333-94           0.00     6.408      -7.371   8 C            C2/m       1
C-17050-333-1            0.00     5.667       0.001   8 C            P63/mmc   11
C-17050-333-2            0.00     5.667       0.001   8 C            Fd-3m     57
C-17050-333-12           0.00     5.853       0.002   8 C            C2/m       4
C-17050-333-89           0.00     6.012       0.041   8 C            Cmmm       1
C-17050-333-82           0.00     7.699       0.045   8 C            Cm         1
C-17050-333-85           0.00     6.088       0.053   8 C            C2/m       1
C-17050-333-61           0.00     6.093       0.055   8 C            P-1        1
C-17050-333-68           0.00     5.806       0.114   8 C            P2/m       1
C-17050-333-45           0.00     6.197       0.118   8 C            P-1        1
C-17050-333-96           0.00     5.822       0.124   8 C            P-1        1
C-17050-333-27           0.00     5.816       0.124   8 C            C2/m       1
C-17050-333-49           0.00     5.959       0.217   8 C            P1         1
C-17050-333-3            0.00     5.776       0.235   8 C            P-1        1
C-17050-333-77           0.00     5.641       0.235   8 C            C2/m       1
C-17050-333-25           0.00     5.889       0.238   8 C            Cmmm       1
C-17050-333-81           0.00     6.094       0.457   8 C            P-1        1
C-17050-333-14           0.00     6.020       0.457   8 C            P-1        1
C-17050-333-21           0.00     6.028       0.457   8 C            P-1        1
C-17050-333-65           0.00     6.224       0.524   8 C            P1         1
C-17050-333-80           0.00     7.437       0.546   8 C            Immm       1
C-17050-333-97           0.00     5.936       0.549   8 C            P-1        1
C-17050-333-40           0.00     5.867       0.580   8 C            C2/m       1
C-17050-333-66           0.00     5.647       0.585   8 C            P1         1
C-17050-333-16           0.00     5.840       0.627   8 C            C2         2
C-17050-333-24           0.00     6.087       0.725   8 C            P1         1
C-17050-333-78           0.00     5.997       0.860   8 C            P2         1
C-17050-333-60           0.00     5.576       0.972   8 C            R-3        1
C-17050-333-36           0.00     6.051       1.068   8 C            P1         1
C-17050-333-11           0.00     5.555       1.326   8 C            C2         1

Now we use the cryan script (through the ca wrapper) to extract all the 2-d structures
in the set. To determine the dimensionality, a contact distance of 2 Ang is chosen.
The adjancency matrix generated using this contact distance is analysed, and the layered
structures are highlighted.

host:3.2 cjp10$ ca -bl 2 -d 2 -r
C-17050-333-49           0.00     5.959      -7.154   8 C            P1         1    3 2.00
C-17050-333-77           0.00     5.641       0.235   8 C            C2/m       1    1 2.00
C-17050-333-3            0.00     5.776       0.235   8 C            P-1        1    1 2.00
C-17050-333-81           0.00     6.094       0.457   8 C            P-1        1    1 2.00
C-17050-333-14           0.00     6.020       0.457   8 C            P-1        1    1 2.00
C-17050-333-21           0.00     6.028       0.457   8 C            P-1        1    1 2.00
C-17050-333-24           0.00     6.087       0.725   8 C            P1         1    1 2.00


