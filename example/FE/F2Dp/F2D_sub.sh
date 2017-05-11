#!/bin/sh
#SBATCH --constraint=my_OS
#SBATCH --exclusive
#SBATCH -n 64            # num cores
#SBATCH -t 12:00:00      # wall time  
#SBATCH -p my_cluster    # partition name
#SBATCH -J F2D_rup       # sensible name for the job
#SBATCH --mem=64000      # this is specified in megabytes per node

# defmod root
export EXE_DIR=${HOME}/defmod-swpc/bin
export TOOL_DIR=${HOME}/defmod-swpc/src/tool
if [ -f ${PWD}/defmod ]; then
    rm defmod 
fi
ln -s ${EXE_DIR}/defmod ${PWD}/defmod

mpirun ../defmod -f F2D.inp -pc_type lu -pc_factor_mat_solver_package mumps
rm defmod
${TOOL_DIR}/def_sort.py F2D rsf slip 2D poro
${TOOL_DIR}/def_plot.py F2D rsf slip
