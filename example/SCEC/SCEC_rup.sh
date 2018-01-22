#!/bin/sh
#SBATCH --constraint=my OS  # delete if uniform OS
#SBATCH --exclusive         # no sharing node 
#SBATCH -n 192              # number of cores
#SBATCH -t 12:00:00         # wall time limit 
#SBATCH -p my part          # partition name
#SBATCH -J SCEC_rup         # job name
#SBATCH --mem=64000         # megabytes memory per node

# root dir 
export EXE_DIR=${HOME}/defmod-swpc/bin
export TOOL_DIR=${HOME}/defmod-swpc/src/tool

# mesh (SCEC_msh.sh)
export CUBIT_DIR=my CUBIT root
if [ -f ${CUBIT_DIR}/bin/cubit ]; then
    export MSH_EXE=${CUBIT_DIR}/bin/cubit
else 
    export MSH_EXE=trelis # trelis should be in $PATH
fi
export MSH_OPT="-nojournal -nographics"
${MSH_EXE} ${MSH_OPT} SCEC205.jou
${MSH_EXE} ${MSH_OPT} SCEC10.jou 
${MSH_EXE} ${MSH_OPT} SCEC102.jou

# preprocess (SCEC_prep.sh)
./SCEC_prep.py -m SCEC205.exo -p 205
./SCEC_prep.py -m SCEC10.exo  -p 10
./SCEC_prep.py -m SCEC102.exo -p 102 -d 1

# run 
time mpirun ${EXE_DIR}/defmod -f SCEC205.inp     -pc_asm_overlap 2 -fd 1
time mpirun ${EXE_DIR}/defmod -f SCEC10.inp      -pc_asm_overlap 2 -fd 1
time mpirun ${EXE_DIR}/defmod -f SCEC102_dsp.inp -pc_asm_overlap 2 -fd 1
