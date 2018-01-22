#!/bin/sh
#SBATCH --constraint=my OS # delete if uniform OS
#SBATCH --exclusive        # no sharing node 
#SBATCH -N 8               # num nodes
#SBATCH -n 16              # num MPI 
#SBATCH -c 8               # num thread per MPI  
#SBATCH -t 12:00:00        # wall time  
#SBATCH -p my part         # partition name
#SBATCH -J SCEC_sw         # job name
#SBATCH --mem=64000        # megabytes memory per node
export I_MPI_PIN_DOMAIN=omp
export OMP_NUM_THREADS=8

# root dir
export EXE_DIR=${HOME}/defmod-swpc/bin
export TOOL_DIR=${HOME}/defmod-swpc/src/tool

# run spwc
time mpirun ${EXE_DIR}/swpc_3d.x -i SCEC205.inf     -r SCEC205
time mpirun ${EXE_DIR}/swpc_3d.x -i SCEC10.inf      -r SCEC10
time mpirun ${EXE_DIR}/swpc_3d.x -i SCEC102_dsp.inf -r SCEC102_dsp

# sort in .mat (SCEC_sort.sh)
# FE waveform and slip profile
${TOOL_DIR}/def_sort.py SCEC205 slip clean
${TOOL_DIR}/def_sort.py SCEC10  slip clean
${TOOL_DIR}/def_sort.py SCEC102_dsp slip rsf clean
# FD waveform 
${TOOL_DIR}/fd_sort.py SCEC205
${TOOL_DIR}/fd_sort.py SCEC10
${TOOL_DIR}/fd_sort.py SCEC102_dsp

# compare FE, FE-FE against others (SCEC_plot.sh)
./SCEC_plot.py -m SCEC205     -p 205 -r 1 -fd 1
./SCEC_plot.py -m SCEC10      -p 10  -r 1 -fd 1
./SCEC_plot.py -m SCEC102_dsp -p 102 -r 1 -fd 1 -d 1
