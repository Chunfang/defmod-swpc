#!/bin/sh
#SBATCH --constraint=my OS
#SBATCH --exclusive
#SBATCH -N 12             # num nodes
#SBATCH -n 24             # num MPI 
#SBATCH -c 8              # num thread per MPI  
#SBATCH -t 72:00:00       # wall time  
#SBATCH -p my part        # partition name
#SBATCH -J HF3D_sw        # sensible name for the job
#SBATCH --mem=64000       # this is specified in megabytes per node

export I_MPI_PIN_DOMAIN=omp
export OMP_NUM_THREADS=8

# root dir
export EXE_DIR=${HOME}/defmod-swpc/bin
export TOOL_DIR=${HOME}/defmod-swpc/src/tool

# run spwc
mpirun ${EXE_DIR}/swpc_3d.x -i input.inf -r SCEC14
mv out SCEC14_snp
mpirun ${EXE_DIR}/swpc_3d.x -i input.inf -r SCEC15 
mv out SCEC15_snp
mpirun ${EXE_DIR}/swpc_3d.x -i input.inf -r HF3D -e 5
mv out HF3D_snp

# sort FD waveform (HF3D_sort.sh)
${TOOL_DIR}/fd_sort.py SCEC14 
${TOOL_DIR}/fd_sort.py SCEC15 
${TOOL_DIR}/fd_sort.py HF3D 

# compare FE-FE against FE and others (HF3D_plot.sh)
./SCEC_spl_plot.py -m SCEC14 -p 14 -r 1 -fd 1
./SCEC_spl_plot.py -m SCEC15 -p 15 -r 1 -fd 1

# plot FD snapshots (HF3D_plot.sh)
./SCEC_fd_plot.py -p 14
./SCEC_fd_plot.py -p 15 
./HF3D_fd_plot.py
