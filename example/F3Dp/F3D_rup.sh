#!/bin/sh
#SBATCH --exclusive
#SBATCH --constraint=my_OS
#SBATCH -N 8 
#SBATCH -n 128                  # num cores
#SBATCH -t 12:00:00             # wall time  
#SBATCH -p my_cluster           # partition name
#SBATCH -J F3D_rup              # sensible name for the job
#SBATCH --mem=64000             # this is specified in megabytes per node

# defmod root
export EXE_DIR=${HOME}/defmod-swpc/bin
export TOOL_DIR=${HOME}/defmod-swpc/src/tool
if [ -f ${PWD}/defmod ]; then
    rm defmod 
fi
ln -s ${EXE_DIR}/defmod ${PWD}/defmod

mpirun ./defmod -f F3D_tet.inp -pc_asm_overlap 2 -fd 1 
# -pc_type asm  -ksp_grmres_restart 31 -sub_ksp_type preonly -sub_pc_type lu -pc_asm_overlap 2 -ksp_monitor 
rm defmod
${TOOL_DIR}/def_sort.py F3D_tet poro slip
