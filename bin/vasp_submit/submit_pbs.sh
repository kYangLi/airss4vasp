#!/bin/bash
#PBS -N __task_name__
#PBS -l nodes=__nodes_num__:ppn=__ppn__
#PBS -l walltime=__pbs_walltime__:00:00
#PBS -q __pbs_queue__
#PBS -j oe

declare -r VASP_CALC_SCRIPT=__vasp_calc_script__
# Enter the Calculate Folder
cd ${PBS_O_WORKDIR}
# Copy the Machinefile
cp ${PBS_NODEFILE} __mpi_mechinefile__

${VASP_CALC_SCRIPT}
