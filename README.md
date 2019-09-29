# AIRSS for VASP (a4v)

## Basic Info

`Author` Liyang@cmt.Tsinghua

`Start Date` 2019.9.21

`Last Update` 2019.9.29

`Version` Alpha 0.1.0

## Description

`AIRSS` (Ab Initio Random Structure Searching) is a fantastic, efficient, and easy for parallel calculation package, but with pool supporting for `VASP`. `airss4vasp`, based on `PBS` or `NSCC(Tianhe)` job management system, is an interface program design for the better communication between AIRSS and VASP.

## Installtion

Before the `make`, you may need to check the `Makefile` to modify some setting for libs. There also are some tips inside the `Makefile` for the installtion.
```bash
vim Makefile # Makefile.GCC & Makefile.INTEL are also ready to use
``` 

To compile and install a4v, please excute the following command in current folder.
```bash
make; make install
```

After the `make install`, the last step you MUST do is adding the a4v `bin` folder to the env `PATH`. Read the output of `make install` for more details.
```bash
echo "export PATH=<a4v>/<bin>/<path>/:${PATH}" >> ~/.bashrc
```
 
## Input File
  
To enable this script, you need perpare another 3 or 4 kinds of files:

- `<seedname>.cell`
- `<seedname>.INCAR-[1-9]`
- `<seedname>.POTCAR`
- `<seedname>.KPOINTS-[1-9]` (Optional)

### `<seedname>.cell`
This is the key file for the whole structure searching task. Please first learn how to use `AIRSS` before using `a4v`.

There is one thing to be care of. 
 
During the generation of new random structures in `AIRSS`, there are two steps: random move, and push a atom. The step `push` is applied to make sure two atoms are not connected to close.

In the `<seedname>.cell`, there are two atom tags called `NOMOVE` and `FIX`. The first one designed for disable the `push` step during the generation, while, the last one designed for disbale the `push` and fix the atom during the `CASTEP` relaxzation. 

Since now we are using `VASP`, in `a4v`,  `FIX` and `NOMOVE` actually have the same effect, and if you mean to fix a atom during the relaxztion, use the tag `SD-*` (the `*` can replace with `X`, `Y`, `Z`, `XY`, `YZ`, `ZX`, `XYZ`). This tag will enable the `Selective dynamics` mode of `VASP`, and fix the atoms you mean to.

Here is an example of the `.cell` file.

```bash
%BLOCK LATTICE_CART
 0.0    2.75    2.75
 2.75   0.0     2.75
 2.75   2.75    0.0
#FIX
%ENDBLOCK LATTICE_CART
 
%BLOCK POSITIONS_FRAC
Si 0.0 0.0 0.0 # Si1 % NUM=1 NOMOVE SD-XYZ
Si 0.0 0.0 0.0 # Si2 % NUM=1
%ENDBLOCK POSITIONS_FRAC
 
#MINSEP=2.0

```

### `<seedname>.INCAR-[1-9]`

The `INCAR` is a input file for VASP relaxzation. The quantity of `INCAR` decided how many times will the structure be relaxed. 

E.g. If there is `Si.INCAR-1`, `Si.INCAR-2`, `Si.INCAR-3` in the calculation file, then the same `Si` structure will first be relaxed using `INCAR-1`, then `INCAR-2`, and at last `INCAR-3`. you can also setting KPOINTS for each INCAR with name `Si.KPOINTS-*`

### `<seedname>.POTCAR`
The order of the elements in `POTCAR` must agree with that in `<seedname>.cell`.

## Start the Search

### Submit Task

Afte the input file getting ready, input the following command to submit the job.
```bash
mass-election.pbs # Use mass-election.nscc if under NSCC system.
```

### `COACH`

There is a conception called `COACH` in the main script. Basically, the `COACH` is a parallel unit among the whole task. Each `COACH` will run independently. They will pick up the calculation task from the `POSCAR-POOL`, and push the result to the `RES-POOL`, then pick up another `POSCAR` to calculate. The POSCAR that already be calculated will be marked in the `TRAIN.record` file.

The nodes number of each `COACH` is simply calculated as *(int)(nodes_number / coach_number)*. Each `COACH` at least will use one node in current version.

Frankly speaking, the more coach there are, the faster your calculation will be.

### Build Cell Mode

There are two different mode for `build cell`. `Build All Cell First` and `Build Cell While Run`. 

The first mode will geneate all structure first, and then use `COACH` to pick up each of them to relax. You can also add or delete some coaches, nodes, and structures during searching. Under this mode, it cannot promise the `RES-POOL` has the same number of structure as `POSCAR-POOL`. And the structure generation will execute on the task submit node, some of the computational cluster may forbidden that. 

The second mode will generate structures one by one during the relaxzation in each `COACH`. This mode can promise the `POSCAR-POOL` and `RES-POOL` has the same number of structures, and structure generation will also execute on calculation nodes. But this mode are not that friendly to parallel, compare to the first mode. And you must decided how many structures, nodes, and coaches you want to use at the very begining.

### Process Check

During the calculation, after enter the main calulation folder(the folder has `PARAM.CONF`), you can use the command below to check the current processing.

```bash
prochk
```

### Result Output

After or during the searching process, you can enter the `RES-POOL` folder and use `match` to check the result. You may need to learn how to use the `cryan` in `AIRSS` first.

E.g.
```bash
match -r -u 0.01 -t 5
```

### Kill the Job

```bash
./_KILLJOB.sh
```

### Clean the foder to init.

```bash
./_CLEAN.sh
```

