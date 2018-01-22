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

time mpirun ${EXE_DIR}/defmod -f F2D_dt12_pp2.inp      -pc_type lu -pc_factor_mat_solver_package mumps
time mpirun ${EXE_DIR}/defmod -f F2D_dt12_pp2_dsp.inp  -pc_type lu -pc_factor_mat_solver_package mumps
time mpirun ${EXE_DIR}/defmod -f F2D_dt12_pp4.inp      -pc_type lu -pc_factor_mat_solver_package mumps
time mpirun ${EXE_DIR}/defmod -f F2D_dt12_pp4_dsp.inp  -pc_type lu -pc_factor_mat_solver_package mumps
time mpirun ${EXE_DIR}/defmod -f F2D_dt12_pp8.inp      -pc_type lu -pc_factor_mat_solver_package mumps
time mpirun ${EXE_DIR}/defmod -f F2D_dt12_pp8_dsp.inp  -pc_type lu -pc_factor_mat_solver_package mumps
time mpirun ${EXE_DIR}/defmod -f F2D_dt12_pp16.inp     -pc_type lu -pc_factor_mat_solver_package mumps
time mpirun ${EXE_DIR}/defmod -f F2D_dt12_pp16_dsp.inp -pc_type lu -pc_factor_mat_solver_package mumps

${TOOL_DIR}/def_sort.py F2D_dt12_pp2      rsf slip 2D poro
${TOOL_DIR}/def_sort.py F2D_dt12_pp4      rsf slip 2D poro
${TOOL_DIR}/def_sort.py F2D_dt12_pp8      rsf slip 2D poro
${TOOL_DIR}/def_sort.py F2D_dt12_pp16     rsf slip 2D poro

${TOOL_DIR}/def_sort.py F2D_dt12_pp2_dsp  rsf slip 2D poro
${TOOL_DIR}/def_sort.py F2D_dt12_pp4_dsp  rsf slip 2D poro
${TOOL_DIR}/def_sort.py F2D_dt12_pp8_dsp  rsf slip 2D poro
${TOOL_DIR}/def_sort.py F2D_dt12_pp16_dsp rsf slip 2D poro
