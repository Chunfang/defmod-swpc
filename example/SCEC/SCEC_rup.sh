#!/bin/sh
#SBATCH --constraint=my_OS
#SBATCH --exclusive         # not sharing nodes
#SBATCH -n 128              # number of cores
#SBATCH -t 12:00:00         # wall time limit 
#SBATCH -p my_cluster       # domain name
#SBATCH -J SCEC_rup         # sensible name for the job
#SBATCH --mem=64000         # megabytes memory per node

# root dir 
export EXE_DIR=${HOME}/defmod-swpc/bin
export TOOL_DIR=${HOME}/defmod-swpc/src/tool
if [ -f ${PWD}/defmod ]; then
    rm defmod
fi
ln -s ${EXE_DIR}/defmod ${PWD}/defmod

# mesh (SCEC_msh.sh)
#export CUBIT_DIR=
#if [ -f ${CUBIT_DIR}/bin/cubit ]; then
#    export MSH_EXE=${CUBIT_DIR}/bin/cubit
#else 
#    export MSH_EXE=trelis # trelis should be in $PATH
#fi
#export MSH_OPT="-nojournal -nographics"
#${MSH_EXE} ${MSH_OPT} SCEC205.jou
#${MSH_EXE} ${MSH_OPT} SCEC10.jou 
#${MSH_EXE} ${MSH_OPT} SCEC102.jou

# preprocess (SCEC_prep.sh)
#./SCEC_prep.py -m SCEC205.exo -p 205
#./SCEC_prep.py -m SCEC10.exo  -p 10
#./SCEC_prep.py -m SCEC102.exo -p 102 -d 1

# run 
mpirun ./defmod -f SCEC205.inp -pc_asm_overlap 2 -fd 1
mpirun ./defmod -f SCEC10.inp -pc_asm_overlap 2 -fd 1
mpirun ./defmod -f SCEC102_dsp.inp -pc_asm_overlap 2 -fd 1
rm defmod

# sort outcome (SCEC_sort.sh) 
${TOOL_DIR}/def_sort.py SCEC205 slip
${TOOL_DIR}/def_sort.py SCEC10  slip
${TOOL_DIR}/def_sort.py SCEC102_dsp slip rsf

# compare Defmod against others (SCEC_post.sh) 
./SCEC_post.py -m SCEC205 -p 205 -r 1
./SCEC_post.py -m SCEC10  -p 10 -r 1
./SCEC_post.py -m SCEC102_dsp -p 102 -d 1
