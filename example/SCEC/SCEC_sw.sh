#!/bin/sh
#SBATCH -p sched_mit_erl          # partition name
#SBATCH --constraint=centos7
#SBATCH --exclusive
#SBATCH --nodelist=node[297-304]
#SBATCH -n 16                     # num MPI 
#SBATCH -c 8                      # num thread per MPI  
#SBATCH -t 12:00:00               # wall time  
#SBATCH -J SCEC_sw                # sensible name for the job
#SBATCH --mem=51200               # this is specified in megabytes per node

##SBATCH -N 8                     # num nodes
##SBATCH --ntasks-per-node 2       # number of cores per node
##SBATCH --mem=64000

export I_MPI_PIN_DOMAIN=omp
export OMP_NUM_THREADS=8

# root dir
export EXE_DIR=${HOME}/GeoMech/defmod-swpc/bin
export TOOL_DIR=${HOME}/GeoMech/defmod-swpc/src/tool

# run spwc
time mpirun ${EXE_DIR}/swpc_3d.x -i SCEC205_dsp.inf -r SCEC205_dsp
time mpirun ${EXE_DIR}/swpc_3d.x -i SCEC10_dsp.inf  -r SCEC10_dsp
time mpirun ${EXE_DIR}/swpc_3d.x -i SCEC102_dsp.inf -r SCEC102_dsp

# sort FD waveform
time ${TOOL_DIR}/fd_sort.py SCEC205_dsp
time ${TOOL_DIR}/fd_sort.py SCEC10_dsp
time ${TOOL_DIR}/fd_sort.py SCEC102_dsp

# Benchmark comparison
./SCEC_plot.py -m SCEC205_dsp -p 205 -fd 1 -r 1 -d 1
./SCEC_plot.py -m SCEC10_dsp  -p 10  -fd 1 -r 1 -d 1
./SCEC_plot.py -m SCEC102_dsp -p 102 -fd 1 -r 1 -d 1
