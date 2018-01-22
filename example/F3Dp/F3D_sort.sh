#!/bin/sh
export TOOL_DIR=${HOME}/defmod-swpc/src/tool
${TOOL_DIR}/def_sort.py F3D_syn poro slip
${TOOL_DIR}/fd_sort.py F3D_syn
