#!/bin/bash
#SBATCH -J __task_name__
#SBATCH -N __nodes_num__
#SBATCH -n __cores_per_coach__
#SBATCH -p __pbs_queue__

declare -r VASP_CALC_SCRIPT=__vasp_calc_script__

${VASP_CALC_SCRIPT}
