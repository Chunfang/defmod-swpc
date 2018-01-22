#!/bin/sh
#SBATCH --constraint=my OS
#SBATCH --exclusive
#SBATCH -n 64            # num cores
#SBATCH -t 12:00:00      # wall time  
#SBATCH -p my partition  # partition name
#SBATCH -J quake         # sensible name for the job
#SBATCH --mem=64000      # this is specified in megabytes per node

# defmod root
export EXE_DIR=${HOME}/defmod-swpc/bin
export TOOL_DIR=${HOME}/defmod-swpc/src/tool
# hybrid run
time mpirun ${EXE_DIR}/defmod -f SCEC10-2D_dsp.inp        -pc_type lu -pc_factor_mat_solver_package mumps
time mpirun ${exe_dir}/defmod -f SCEC11-2D_dsp.inp        -pc_type lu -pc_factor_mat_solver_package mumps
time mpirun ${exe_dir}/defmod -f SCEC102-2D_dsp.inp       -pc_type lu -pc_factor_mat_solver_package mumps
# implicit run
time mpirun ${EXE_DIR}/defmod -f SCEC10-2D_dsp-alpha.inp  -pc_type lu -pc_factor_mat_solver_package mumps
time mpirun ${exe_dir}/defmod -f SCEC11-2D_dsp-alpha.inp  -pc_type lu -pc_factor_mat_solver_package mumps
time mpirun ${exe_dir}/defmod -f SCEC102-2D_dsp-alpha.inp -pc_type lu -pc_factor_mat_solver_package mumps
