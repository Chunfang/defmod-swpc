#!/bin/sh
# FE waveform and slip profile
export TOOL_DIR=${HOME}/defmod-swpc/src/tool
${TOOL_DIR}/def_sort.py SCEC205 slip clean
${TOOL_DIR}/def_sort.py SCEC10  slip clean
${TOOL_DIR}/def_sort.py SCEC102_dsp slip rsf clean
# FD waveform 
${TOOL_DIR}/fd_sort.py SCEC205
${TOOL_DIR}/fd_sort.py SCEC10
${TOOL_DIR}/fd_sort.py SCEC102_dsp
