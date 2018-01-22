#!/bin/sh
export TOOL_DIR=${HOME}/demod-swpc/src/tool
# sort FE outcome
${TOOL_DIR}/def_sort.py SCEC14 slip
${TOOL_DIR}/def_sort.py SCEC15 slip
${TOOL_DIR}/def_sort.py HF3D   slip poro
# sort FD waveform 
${TOOL_DIR}/fd_sort.py SCEC14 
${TOOL_DIR}/fd_sort.py SCEC15 
${TOOL_DIR}/fd_sort.py HF3D 
