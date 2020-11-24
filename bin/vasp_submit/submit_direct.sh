#!/bin/bash
#

declare -r VASP_CALC_SCRIPT=__vasp_calc_script__
declare -r TASK_NAME=__task_name__

echo "[calc] ${TASK_NAME} Start!!!"
${VASP_CALC_SCRIPT} >> ${TASK_NAME}.out
