#!/bin/bash
#

# +---------------+
# | Tools Defined |
# +---------------+
grep_input_file(){
  keyword="$1"
  res=$(grep "${keyword}" a4v.allparas 2>/dev/null | tail -1 | 
        awk -F '#' '{print $1}' | awk -F '=' '{print $NF}')
  echo ${res}
}

# +--------------------+
# | Parameters Setting |
# +--------------------+
declare -r  SEED_NAME=$(grep_input_file 'SEED_NAME')
declare -r  INTEL_MODULE=$(grep_input_file 'INTEL_MODULE')
declare -r  VASP_PROG=$(grep_input_file 'VASP_PROG')
declare -ri STR_NUM=$(grep_input_file 'STR_NUM')
declare -r  NODES_PER_COACH=$(grep_input_file 'NODES_PER_COACH')
declare -r  CORES_PER_COACH=$(grep_input_file 'CORES_PER_COACH')
declare -r  CORES_PER_NODE=$(grep_input_file 'CORES_PER_NODE')
declare -r  SYMM_PREC=$(grep_input_file 'SYMM_PREC')
declare -r  VASP_WALLTIME_INS=$(grep_input_file 'VASP_WALLTIME')
declare -r  VASP_WALLTIME_INM=$((VASP_WALLTIME_INS/60))
declare -r  KEEP_CALC_DETAILS=$(grep_input_file 'KEEP_CALC_DETAILS')
declare -r  MPI_MACHINEFILE=$(grep_input_file 'MPI_MACHINEFILE')
declare -r  JOB_QUEUE=$(grep_input_file 'JOB_QUEUE')
declare -r  SYS_TYPE=$(grep_input_file 'SYS_TYPE')
declare -r  A4V_PATH=$(grep_input_file 'A4V_PATH')
declare -r  RELAX4RES="${A4V_PATH}/vasp_submit/relax4res"
declare -r  CABAL="${A4V_PATH}/cabal"

# +----------------------------+
# | Prepare Before Calculation |
# +----------------------------+
# Get Current Worker Folder Name
current_coach=$(pwd | awk -F '/' '{print $NF}')
current_coach_index=${current_coach##*-}
# Parallel Calculation Parameters
case ${SYS_TYPE} in
'pbs')
  vasp_mpirun="export I_MPI_JOB_TIMEOUT=${VASP_WALLTIME_INS}; ${INTEL_MODULE}; mpirun -machinefile ${MPI_MACHINEFILE} -np ${CORES_PER_COACH} -envall ${VASP_PROG}"
  ;;
'slurm')
  vasp_mpirun="${INTEL_MODULE}; srun --time=${VASP_WALLTIME_INM} ${VASP_PROG}"
  ;;
'nscc')
  vasp_mpirun="${INTEL_MODULE}; yhrun -p ${JOB_QUEUE} -N ${NODES_PER_COACH} -n ${CORES_PER_COACH} --time=${VASP_WALLTIME_INM} ${VASP_PROG}"
  if [ "${JOB_QUEUE}" == "unset-queue" ]; then
    vasp_mpirun="${INTEL_MODULE}; yhrun -N ${NODES_PER_COACH} -n ${CORES_PER_COACH} --time=${VASP_WALLTIME_INM} ${VASP_PROG}"
  fi
  ;;
'direct')
  vasp_mpirun="export I_MPI_JOB_TIMEOUT=${VASP_WALLTIME_INS}; ${INTEL_MODULE}; mpirun -np ${CORES_PER_COACH} -envall ${VASP_PROG}"
  ;;
esac

# +----------------------+
# | Structure Relazation |
# +----------------------+
## Calculation loop
while true; do
  # Pick up POSCAR from the POSCAR POOL, if in BUILD BEFORE RUN mode.
  # Get the current POSCAR list in the POSCAR POOL
  current_poscar_list=$(find ../POSCAR-POOL/ -maxdepth 1 -name "*.vasp" | 
                        awk -F '/' '{print $NF}' | sed 's/-/ /g' | 
                        sort -k2,2 -k3,3 -k4,4 -n | sed 's/ /-/g')
  # Check if the job done, on the BUILD BEFORE RELAX mode.
  if [ -z "${current_poscar_list}" ]; then
    echo "[info] ${current_coach} JOB DONE!"
    echo '[done]'
    break
  fi
  # Pick up one POSCAR
  pick_up_poscar=$(echo ${current_poscar_list} | awk '{print $1}')
  mv ../POSCAR-POOL/${pick_up_poscar} .
  # Check if there is any I/O collision
  if [ ! -s "${pick_up_poscar}" ]; then
    echo "[fail] ${current_coach}: ${pick_up_poscar} has been picked up by others..."
    sleep 0.1
    continue
  fi
  cp ${pick_up_poscar} POSCAR
  # Record the pick up opration
  while true; do
    io_is_busy=$(echo ../IO-BUSY-*.remark | awk '{print $1}')
    if [ ! -e ${io_is_busy} ]; then
      touch ../IO-BUSY-${current_coach}.remark
      echo "[do] ${pick_up_poscar} --> ${current_coach}" \
            >> ../TRAIN.record
      column -t ../TRAIN.record > train.record.temp
      mv train.record.temp ../TRAIN.record
      sleep 1
      rm ../IO-BUSY-${current_coach}.remark
      break
    fi
    echo "[info] ${current_coach} is waitting for I/O..."
    sleep 0.2
  done
  # Get the Name Stamp for this POSCAR
  name_stamp=${pick_up_poscar%\.*}
  # Run the Relaxzation Task 
  ${RELAX4RES} "${name_stamp}" "${vasp_mpirun}" "${SYMM_PREC}"
  # Collect the result
  if [ -e *.res ]; then
    ${CABAL} res poscar < ${name_stamp}.res > ../RES-POOL/${name_stamp}.vasp
    mv ${name_stamp}.res ../RES-POOL/
    cp OUTCAR ../RES-POOL/${name_stamp}.outcar
    return_info="[done] ${pick_up_poscar} --> ${current_coach} --> ${name_stamp}.res"
  else
    return_info="[fail] ${pick_up_poscar} --> ${current_coach} xx> ${name_stamp}.res"
  fi
  # Record the collect opration
  while true; do
    io_is_busy=$(echo ../IO-BUSY-*.remark | awk '{print $1}')
    if [ ! -e ${io_is_busy} ]; then
      touch ../IO-BUSY-${current_coach}.remark
      sed -i "/${pick_up_poscar}/c${return_info}" ../TRAIN.record
      column -t ../TRAIN.record > train.record.temp
      mv train.record.temp ../TRAIN.record
      sleep 1
      rm ../IO-BUSY-${current_coach}.remark
      break
    fi
    echo "[info] ${current_coach} is waitting for I/O..."
    sleep 1
  done
  # keep calculation details if permited
  if [ "${KEEP_CALC_DETAILS}" == "T" ]; then
    mkdir save-calc-details
    mv CHG CHGCAR CONTCAR* DOSCAR EIGENVAL IBZKPT           save-calc-details/
    mv OSZICAR OUTCAR PCDAT POSCAR* INCAR KPOINTS PROCAR    save-calc-details/
    mv REPORT vasprun.xml WAVECAR XDATCAR ${pick_up_poscar} save-calc-details/
    mv save-calc-details ../RES-POOL/${name_stamp}
  else
    rm CHG CHGCAR CONTCAR* DOSCAR EIGENVAL IBZKPT OSZICAR PCDAT OUTCAR
    rm POSCAR* INCAR KPOINTS PROCAR REPORT vasprun.xml WAVECAR XDATCAR
    rm ${pick_up_poscar}
  fi
done