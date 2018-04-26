#!/bin/sh
#SBATCH --constraint=my os
#SBATCH --exclusive
#SBATCH -N 12             # num nodes
#SBATCH -n 24             # num MPI 
#SBATCH -c 8              # num thread per MPI  
#SBATCH -t 72:00:00       # wall time  
#SBATCH -p my cluster     # partition name
#SBATCH -J HF3D_sw        # sensible name for the job
#SBATCH --mem=64000       # this is specified in megabytes per node
export I_MPI_PIN_DOMAIN=omp
export OMP_NUM_THREADS=8

# root dir
export EXE_DIR=${HOME}/defmod-swpc/bin
export TOOL_DIR=${HOME}/defmod-swpc/src/tool

# run spwc
mpirun ${EXE_DIR}/swpc_3d.x -i F3DX.inf -r F3DX14
mv out out14
mpirun ${EXE_DIR}/swpc_3d.x -i F3DX.inf -r F3DX15
mv out out15

# sort FD waveform
${TOOL_DIR}/fd_sort.py F3DX14 
${TOOL_DIR}/fd_sort.py F3DX15
