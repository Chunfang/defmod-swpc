#!/bin/sh
#SBATCH --constraint=my_OS
#SBATCH --exclusive
#SBATCH -N 8              # num nodes
#SBATCH -n 16             # num MPI 
#SBATCH -c 8              # num thread per MPI  
#SBATCH -t 12:00:00       # wall time  
#SBATCH -p my_cluster     # partition name
#SBATCH -J SCEC_sw        # sensible name for the job
#SBATCH --mem=64000       # this is specified in megabytes per node
export I_MPI_PIN_DOMAIN=omp
export OMP_NUM_THREADS=8

# root dir
export EXE_DIR=${HOME}/defmod-swpc/bin
export TOOL_DIR=${HOME}/defmod-swpc/src/tool
if [ -f ${PWD}/swpc_3d.x ]; then
    rm swpc_3d.x
fi
ln -s ${EXE_DIR}/swpc_3d.x ${PWD}/swpc_3d.x

# run spwc
mpirun ./swpc_3d.x -i SCEC205.inf -r SCEC205
mpirun ./swpc_3d.x -i SCEC10.inf  -r SCEC10
mpirun ./swpc_3d.x -i SCEC102_dsp.inf -r SCEC102_dsp
rm swpc_3d.x

# sort FD waveform
${TOOL_DIR}/fd_sort.py SCEC205
${TOOL_DIR}/fd_sort.py SCEC10
${TOOL_DIR}/fd_sort.py SCEC102_dsp

# compare FE-FE against FE and others
./FEFD_post.py -m SCEC205 -p 205
./FEFD_post.py -m SCEC10  -p 10
./FEFD_post.py -m SCEC102_dsp -p 102 -d 1
