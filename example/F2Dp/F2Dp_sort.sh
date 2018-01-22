#!/bin/sh
export TOOL_DIR=${HOME}/defmod-swpc/src/tool
${TOOL_DIR}/def_sort.py F2D_dt12_pp2      rsf slip 2D poro
${TOOL_DIR}/def_sort.py F2D_dt12_pp4      rsf slip 2D poro
${TOOL_DIR}/def_sort.py F2D_dt12_pp8      rsf slip 2D poro
${TOOL_DIR}/def_sort.py F2D_dt12_pp16     rsf slip 2D poro

${TOOL_DIR}/def_sort.py F2D_dt12_pp2_dsp  rsf slip 2D poro
${TOOL_DIR}/def_sort.py F2D_dt12_pp4_dsp  rsf slip 2D poro
${TOOL_DIR}/def_sort.py F2D_dt12_pp8_dsp  rsf slip 2D poro
${TOOL_DIR}/def_sort.py F2D_dt12_pp16_dsp rsf slip 2D poro
