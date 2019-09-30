# AIRSS for VASP (a4v)

## Basic Info

`Author` Liyang@cmt.Tsinghua

`Start Date` 2019.9.21

`Last Update` 2019.9.30

`Version` Alpha 0.1.0

## Description

`AIRSS` (Ab Initio Random Structure Searching) is a fantastic, efficient, and easy to parallel structure searching package, but with pool supporting for `VASP`. `airss4vasp`, based on `PBS` or `NSCC(Tianhe)` job management system, is an interface program design for the better communication between AIRSS and VASP.

## Installtion

Before the `make`, please check the `Makefile` to modify some setting for libs. There also are some tips for the installtion inside the `Makefile`.

```bash
vim Makefile # Makefile.GNU & Makefile.INTEL are also ready for use
```

To compile and install `a4v`, excute the following command in the main folder.

```bash
make; make install
```

The last step that MUST be done is adding the `bin` folder of `a4v` to the env `PATH`. Read the output of `make install` for more details.

```bash
echo "export PATH=<a4v>/<bin>/<path>/:${PATH}" >> ~/.bashrc
```

## Input File
  
To enable this script, you need perpare another 3 or 4 kinds of files:

- `<seedname>.cell`
- `<seedname>.INCAR-[1-9]`
- `<seedname>.KPOINTS-[1-9]` (Optional)
- `<seedname>.POTCAR`

### `<seedname>.cell`

This is the key file for the whole structure searching task. Please first learn how to use `AIRSS` before using `a4v`. 

There is one thing need to be explained more specifically.

During the generation of the new random structures in `AIRSS`, basically, there are two steps for the movement of a single atom: `random shift` and `push`. The step `push` is applied to make sure two atoms are not connected to close.

In the `<seedname>.cell`, there are two atomic tags called `NOMOVE` and `FIX`. The first one designed for disable the `push` step, while, the last one designed for disbale the `push` **and** fix the atom during the `CASTEP` relaxzation. 

Since now we are using `VASP`, in `a4v`,  `FIX` and `NOMOVE` actually have the same effect, and if you mean to fix a atom during the relaxztion, use the tag `SD-*` (where the `*` can replace with `X`, `Y`, `Z`, `XY`, `YZ`, `ZX`, `XYZ`). This tag will enable the `Selective dynamics` mode of `VASP` in `POSCAR`.

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

E.g. If there are `Si.INCAR-1`, `Si.INCAR-2`, `Si.INCAR-3` in the calculation file, then the same `Si` structure will first be relaxed using `INCAR-1`, then `INCAR-2`, and at last `INCAR-3`. You can also setting KPOINTS for each INCAR with name `Si.KPOINTS-[1-9]`

### `<seedname>.POTCAR`

The order of the elements in `POTCAR` must agree with that in `<seedname>.cell`.

## Start the Search

### Modify the Default Parameters

Before submit the job, please enter the main task submit script to setting the default value of some paramters.

```bash
vim $(which mass-election.pbs) # mass-election.nscc if under NSCC system.
```

Here is the list of the default parameters you can modify:
|Parameter Name            |Type|Descripution|
|-------------------------:|:----|:-----------|
|DEFAULT_INTEL_MODULE      |Char |Intel Module load command, if you are using Intel Lib for VASP and AIRSS.|
|DEFAULT_VASP_EXEC         |Char |Path of VASP executive program.|
|DEFAULT_NODES_NUM         |Int  |Nodes number used in task.|
|DEFAULT_CORES_NUM_PER_NODE|Int  |The number of cores of each nodes in your machine.|
|DEFAULT_STR_NUM           |Int  |Total structure number you want to search.|
|DEFAULT_SYMM_PREC         |Float|Symmetry precise used in `cellsym`.|
|DEFAULT_TIME_LIMIT        |Int  |Wall time for a single VASP relazation(one INCAR step).|
|DEFAULT_DEL_CALC_DETAILS  |Bool |Keep the calculation detail file or not.|
|DEFAULT_BUILD_CELL_BEFORE |y/n  |Select the build cell mode.|
|DEFAULT_IS_2D_MATERIAL    |y/n  |IS this a 2D material? If it is, then add '-f' flag when generate KPOINTS.|
|DEFAULT_KP_SEP            |Float|K poinnts separation in k-space.|

### Submit Task

Afte the input file getting ready, input the following command to submit the job.

```bash
mass-election.pbs # Use mass-election.nscc if under NSCC system.
```

### COACH

The `COACH` is a parallel unit among the whole task. Each `COACH` will run independently. They will pick up the structure from the `POSCAR-POOL`, push the result to the `RES-POOL`. The POSCAR that already be calculated will be marked in the `TRAIN.record` file.

The nodes number of each `COACH` is simply calculated as ***(int)(nodes_number/coach_number)*** . Each `COACH` at least will use one node in current version.

### Build Cell Mode

There are two different mode for `build cell`. `Build All Cell First` and `Build Cell While Run`.

The first mode will geneate all structure first, and then use `COACH` to pick up each of them to relax. You can add or delete some coaches, nodes, and structures during the searching. Under this mode, it cannot promise the `RES-POOL` has the same number of structure as `POSCAR-POOL`, since some of the structure may too strange to get a convergent result. And the structure generation will execute on the task submit node, some of the computational cluster may forbidden that.

The second mode will generate structures one by one during the relaxzation in each `COACH`. This mode can promise the `POSCAR-POOL` and `RES-POOL` has the same number of structures, and structure generation will also execute on calculation nodes. But, under this mode, it is not that friendly to parallel compare to the first mode. And you must decided how many structures, nodes, and coaches you want to use at the very begining.

### Process Check

During the calculation, after enter the main calulation folder(the folder has `PARAM.CONF`), you can use the command below to check the current processing.

```bash
prochk
```

### Result Output

After or during the searching process, you can enter the `RES-POOL` folder and use `match` to check the result. You may need to learn how to use the `cryan` in `AIRSS` first.

For example,

```bash
match -r -u 0.01 -t 5
```

### Kill the Job

```bash
./_KILLJOB.sh
```

### Clean the Foder to Init.

```bash
./_CLEAN.sh
```
