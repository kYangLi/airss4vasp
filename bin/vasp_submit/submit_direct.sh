#!/bin/bash
#

declare -r  TASK_NAME='__task_name__'
declare -r  VASP_CALC_SCRIPT='__vasp_calc_script__'
declare -r  VASP_PROG='__vasp_prog__'
declare -ir CORES_PER_COACH=__cores_per_coach__
declare -ir VASP_WALLTIME_INS=__vasp_walltime__

# Load module
__prog_module__
# MPI job timeout
export I_MPI_JOB_TIMEOUT=${VASP_WALLTIME_INS}
# MPIRUN command 
declare -r VR_MPIRUN_COMMAND="mpirun -np ${CORES_PER_COACH} ${VASP_PROG}"

echo "[calc] ${TASK_NAME} Start!!!"
${VASP_CALC_SCRIPT} "${VR_MPIRUN_COMMAND}" >> ${TASK_NAME}.out
