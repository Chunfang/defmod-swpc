#!/bin/sh
#SBATCH --exclusive
#SBATCH --constraint=centos6
#SBATCH -N 12 
#SBATCH -n 192                  # num cores/node
#SBATCH -t 12:00:00             # wall time  
#SBATCH -p sched_mit_erl        # partition name
#SBATCH -J CR3D_rup         # sensible name for the job
#SBATCH --mem=64000             # this is specified in megabytes per node

export DIR_FV=${HOME}/pflotran/src/pflotran
export DIR_FE=${HOME}/defmod-swpc/bin
#trelis -nojournal -nographics F3D_bd.jou
#./F3D_bd_usg.py F3D_bd.exo
time mpirun ${DIR_FV}/pflotran-mkl-icc -pflotranin CR3D.in
time mpirun ${DIR_FE}/defmod-mkl-icc -f CR3D.inp -pc_asm_overlap 2 -fv 1 -fd 1 
