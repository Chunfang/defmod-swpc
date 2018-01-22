#!/bin/sh
#SBATCH --exclusive
#SBATCH --constraint=my OS
#SBATCH -N 8                    # num nodes
#SBATCH -n 16                   # num MPI 
#SBATCH -c 8                    # num thread per MPI  
#SBATCH -t 12:00:00             # wall time  
#SBATCH -p my part              # partition name
#SBATCH -J F3D_sw               # sensible name for the job
#SBATCH --mem=64000             # this is specified in megabytes per node
export I_MPI_PIN_DOMAIN=omp
export OMP_NUM_THREADS=8

# swpc root
export EXE_DIR=${HOME}/defmod-swpc/bin
#export TOOL_DIR=${HOME}/defmod-swpc/src/tool

# run swpc
mpirun ${EXE_DIR}/swpc_3d.x -i F3D_syn.inf -r F3D_syn -e 15 

#${TOOL_DIR}/fd_sort.py F3D_syn
#./F3D_plot.py F3D_syn
