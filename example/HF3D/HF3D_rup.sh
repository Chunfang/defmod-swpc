#!/bin/sh
#SBATCH --constraint=my OS
#SBATCH --exclusive         # not sharing nodes
#SBATCH -n 192              # number of cores
#SBATCH -t 72:00:00         # wall time limit 
#SBATCH -p my part          # domain name
#SBATCH -J HF3D_rup         # sensible name for the job
#SBATCH --mem=64000         # megabytes memory per node

# root dir 
export EXE_DIR=${HOME}/defmod-swpc/bin
export TOOL_DIR=${HOME}/defmod-swpc/src/tool

# preprocess (HF3D_prep.sh)
export CUBIT_DIR=my CUBIT root
if [ -f ${CUBIT_DIR}/bin/cubit ]; then
    export MSH_EXE=${CUBIT_DIR}/bin/cubit
else 
    export MSH_EXE=trelis # trelis should be in $PATH
fi
export MSH_OPT="-nojournal -nographics"
${MSH_EXE} ${MSH_OPT} SCEC_spl.jou
./SCEC_spl.py -m SCEC_spl.exo -p 14
./SCEC_spl.py -m SCEC_spl.exo -p 15
./HF3D.py -m SCEC_spl.exo

# run 
time mpirun ${EXE_DIR}/defmod -f SCEC14.inp -pc_asm_overlap 2 -fd 1
time mpirun ${EXE_DIR}/defmod -f SCEC15.inp -pc_asm_overlap 2 -fd 1
time mpirun ${EXE_DIR}/defmod -f HF3D.inp   -pc_asm_overlap 2 -fd 1

# sort FE outcome (HF3D_sort.sh) 
${TOOL_DIR}/def_sort.py SCEC14 slip
${TOOL_DIR}/def_sort.py SCEC15 slip
${TOOL_DIR}/def_sort.py HF3D   slip poro
