#!/bin/bash
#PBS -N __task_name__
#PBS -l nodes=__nodes_num__:ppn=__ppn__
#PBS -l walltime=__job_walltime__:00:00
#PBS -q __job_queue__
#PBS -j oe

declare -r  VASP_CALC_SCRIPT='__vasp_calc_script__'
declare -r  VASP_PROG='__vasp_prog__'
declare -r  MPI_MF='__mpi_machinefile__'
declare -ir CORES_PER_COACH=__cores_per_coach__
declare -ir VASP_WALLTIME_INS=__vasp_walltime__

# Enter the Calculate Folder
cd ${PBS_O_WORKDIR}
# Copy the Machinefile
uniq ${PBS_NODEFILE} > ${MPI_MF}
# Load module
__prog_module__
# MPI job timeout
export I_MPI_JOB_TIMEOUT=${VASP_WALLTIME_INS}
# MPIRUN command 
declare -r VR_MPIRUN_COMMAND="mpirun -machinefile ${MPI_MF} -np ${CORES_PER_COACH} ${VASP_PROG}"

${VASP_CALC_SCRIPT} "${VR_MPIRUN_COMMAND}"
