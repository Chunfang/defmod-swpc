#!/bin/sh
#SBATCH --constraint=my os
#SBATCH --exclusive         # not sharing nodes
#SBATCH -n 192              # number of cores
#SBATCH -t 72:00:00         # wall time limit 
#SBATCH -p my cluster       # domain name
#SBATCH -J F3DX_rup         # sensible name for the job
#SBATCH --mem=64000         # megabytes memory per node

# root dir 
export EXE_DIR=${HOME}/defmod-swpc/bin
export TOOL_DIR=${HOME}/defmod-swpc/src/tool

# run 
time mpirun ${EXE_DIR}/defmod-mkl-icc -f F3DX14.inp -pc_asm_overlap 2 -fd 1
time mpirun ${EXE_DIR}/defmod-mkl-icc -f F3DX15.inp -pc_asm_overlap 2 -fd 1

# sort outcome (SCEC_sort.sh) 
${TOOL_DIR}/def_sort.py F3DX14 slip
${TOOL_DIR}/def_sort.py F3DX15 slip
