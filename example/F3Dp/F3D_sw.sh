#!/bin/sh
#SBATCH --exclusive
#SBATCH --constraint=my_OS
#SBATCH -N 8                    # num nodes
#SBATCH -n 16                   # num MPI 
#SBATCH -c 8                    # num thread per MPI  
#SBATCH -t 12:00:00             # wall time  
#SBATCH -p my_cluster           # partition name
#SBATCH -J F3D_sw               # sensible name for the job
#SBATCH --mem=64000             # this is specified in megabytes per node
export I_MPI_PIN_DOMAIN=omp
export OMP_NUM_THREADS=8

# swpc root
export EXE_DIR=${HOME}/defmod-swpc/bin
export TOOL_DIR=${HOME}/defmod-swpc/src/tool
if [ -f ${PWD}/swpc_3d.x ]; then
    rm swpc_3d.x
fi
ln -s ${EXE_DIR}/swpc_3d.x ${PWD}/swpc_3d.x

# run swpc
mpirun ./swpc_3d.x -i input.inf -r F3D_tet -e 11 
rm swpc_3d.x

${TOOL_DIR}/fd_sort.py F3D_tet
./F3D_post.py F3D_tet
