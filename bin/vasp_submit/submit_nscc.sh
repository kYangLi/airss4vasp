#!/bin/bash
#SBATCH -J __task_name__
#SBATCH -N __nodes_num__
#SBATCH -n __total_cores__

declare -r VASP_CALC_SCRIPT=__vasp_calc_script__

${VASP_CALC_SCRIPT}
