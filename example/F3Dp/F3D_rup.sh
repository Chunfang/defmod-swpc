#!/bin/sh
#SBATCH --exclusive
#SBATCH --constraint=my OS
#SBATCH -N 8 
#SBATCH -n 128                  # num cores
#SBATCH -t 12:00:00             # wall time  
#SBATCH -p my part              # partition name
#SBATCH -J F3D_rup              # sensible name for the job
#SBATCH --mem=64000             # this is specified in megabytes per node

# defmod root
export EXE_DIR=${HOME}/defmod-swpc/bin
#export TOOL_DIR=${HOME}/defmod-swpc/src/tool

mpirun ${EXE_DIR}/defmod -f F3D_syn.inp -pc_asm_overlap 2 -fd 1 
#${TOOL_DIR}/def_sort.py F3D_syn poro slip
