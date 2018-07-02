#!/bin/sh
#SBATCH --exclusive
#SBATCH --constraint=centos6
#SBATCH -N 10                   # num nodes
#SBATCH -n 10                   # num MPI 
#SBATCH -c 16                   # num thread per MPI  
#SBATCH -t 12:00:00             # wall time  
#SBATCH -p sched_mit_erl        # partition name
#SBATCH -J CR3D_sw          # sensible name for the job
#SBATCH --mem=64000             # this is specified in megabytes per node
export I_MPI_PIN_DOMAIN=omp
export OMP_NUM_THREADS=16

# swpc run 
export DIR_FD=${HOME}/defmod-swpc/bin
export DIR_TOOL=${HOME}/defmod-swpc/src/tool
mpirun ${DIR_FD}/swpc_3d.x -i CR3D.inf -r CR3D -e 10

# postprocess
time ${DIR_TOOL}/def_sort.py CR3D poro #slip
time ${DIR_TOOL}/fd_sort.py CR3D
