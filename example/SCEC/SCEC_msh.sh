#!/bin/sh
export CUBIT_DIR=
if [ -f ${CUBIT_DIR}/bin/cubit ]; then
    export MSH_EXE=${CUBIT_DIR}/bin/cubit
else 
    export MSH_EXE=trelis # trelis should be in $PATH
fi
export MSH_OPT="-nojournal -nographics"
${MSH_EXE} ${MSH_OPT} SCEC205.jou
${MSH_EXE} ${MSH_OPT} SCEC10.jou 
${MSH_EXE} ${MSH_OPT} SCEC102.jou
