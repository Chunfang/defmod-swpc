#!/bin/sh
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
