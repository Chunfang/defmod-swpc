#!/bin/sh
export TOOL_DIR=${HOME}/defmod-swpc/src/tool
${TOOL_DIR}/def_sort.py SCEC205 slip clean
${TOOL_DIR}/def_sort.py SCEC10  slip clean
${TOOL_DIR}/def_sort.py SCEC102_dsp slip rsf clean
