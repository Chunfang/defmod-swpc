#!/bin/sh
#SBATCH --constraint=centos7
#SBATCH --exclusive               # not sharing nodes
#SBATCH -p sched_mit_erl          # domain name
#SBATCH --nodelist=node[297-308]  #[170-178,180-194] # [361-378] #[360-383]
#SBATCH --ntasks-per-node 16      # number of cores per node
#SBATCH -t 12:00:00               # wall time limit 
#SBATCH -J SCEC_rup               # sensible name for the job
#SBATCH --mem=51200

##SBATCH -N 7
##SBATCH --mem=64000              # megabytes memory per node
##export MKL_NUM_THREADS=1


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
#./SCEC_prep.py -m SCEC205.exo -p 205 -d 1
#./SCEC_prep.py -m SCEC10.exo  -p 10  -d 1
#./SCEC_prep.py -m SCEC102.exo -p 102 -d 1

# root dir
export EXE_DIR=${HOME}/GeoMech/defmod-swpc/bin
export TOOL_DIR=${HOME}/GeoMech/defmod-swpc/src/tool

# run defmod
time mpirun ${EXE_DIR}/defmod-mkl-icc -f SCEC205_dsp.inp -pc_asm_overlap 2 -fd 1
time mpirun ${EXE_DIR}/defmod-mkl-icc -f SCEC10_dsp.inp  -pc_asm_overlap 2 -fd 1
time mpirun ${EXE_DIR}/defmod-mkl-icc -f SCEC102_dsp.inp -pc_asm_overlap 2 -fd 1
#time mpirun ${EXE_DIR}/defmod-mkl-icc -f SCEC205_damp_dsp.inp -pc_asm_overlap 2
#time mpirun ${EXE_DIR}/defmod-mkl-icc -f SCEC10_damp_dsp.inp  -pc_asm_overlap 2
#time mpirun ${EXE_DIR}/defmod-mkl-icc -f SCEC102_damp_dsp.inp -pc_asm_overlap 2

# sort outcomes
time ${TOOL_DIR}/fe_sort.py SCEC205_dsp dyn slip
time ${TOOL_DIR}/fe_sort.py SCEC10_dsp  dyn slip
time ${TOOL_DIR}/fe_sort.py SCEC102_dsp dyn slip rsf
#time ${TOOL_DIR}/fe_sort.py SCEC205_damp_dsp dyn slip
#time ${TOOL_DIR}/fe_sort.py SCEC10_damp_dsp  dyn slip
#time ${TOOL_DIR}/fe_sort.py SCEC102_damp_dsp dyn slip rsf
