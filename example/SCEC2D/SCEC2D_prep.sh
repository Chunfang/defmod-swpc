#!/bin/sh
# mesh
export CUBIT_DIR=
if [ -f ${CUBIT_DIR}/bin/cubit ]; then
    export MSH_EXE=${CUBIT_DIR}/bin/cubit
else 
    export MSH_EXE=trelis # trelis should be in $PATH
fi
export MSH_OPT="-nojournal -nographics"
${MSH_EXE} ${MSH_OPT} SCEC10-2D.jou
${MSH_EXE} ${MSH_OPT} SCEC102-2D.jou
# defmod prep
./SCEC2D_prep.py -m SCEC10-2D.exo  -p 10  -d 1
./SCEC2D_prep.py -m SCEC10-2D.exo  -p 11  -d 1
./SCEC2D_prep.py -m SCEC102-2D.exo -p 102 -d 1
# implicit dynamic model
./SCEC2D_prep.py -m SCEC10-2D.exo  -p 10  -d 1 -alpha 1 
./SCEC2D_prep.py -m SCEC10-2D.exo  -p 11  -d 1 -alpha 1
./SCEC2D_prep.py -m SCEC102-2D.exo -p 102 -d 1 -alpha 1
