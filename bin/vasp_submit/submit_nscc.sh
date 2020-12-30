#!/bin/bash
#SBATCH --job-name=__task_name__
#SBATCH --partition=__job_queue__
#SBATCH --nodes=__nodes_num__
#SBATCH --ntasks-per-node=__cores_per_coach__
#SBATCH --exclusive

declare -r  VASP_CALC_SCRIPT='__vasp_calc_script__'
declare -r  VASP_PROG='__vasp_prog__'
declare -ir VASP_WALLTIME_INS=__vasp_walltime__
declare -ir VASP_WALLTIME_INM=$((VASP_WALLTIME_INS/60))

# Load module
__prog_module__
# Define the mpirun command
declare -r VR_MPIRUN_COMMAND="yhrun --time=${VASP_WALLTIME_INM} ${VASP_PROG}"

${VASP_CALC_SCRIPT} "${VR_MPIRUN_COMMAND}"
