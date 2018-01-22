#!/bin/sh
# compare FE-FE against FE and others
./SCEC_spl_plot.py -m SCEC14 -p 14 -r 1 -fd 1
./SCEC_spl_plot.py -m SCEC15 -p 15 -r 1 -fd 1

# plot FD snapshots
./SCEC_fd_plot.py -p 14
./SCEC_fd_plot.py -p 15 
./HF3D_fd_plot.py
