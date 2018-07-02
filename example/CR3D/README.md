This example duplicates [Cappa and Rutqvist 2011](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2011GL048487) model, demonstrating the FV(pflotran)->FE(defmod)->FD(swpc) workflow.

To generate FE mesh file and FE input file, CR3D.exo and CR3D.inp, run:

trelis -nojournal -nographics CR3D.jou\
./CR3D_tet.py CR3D.exo

To run the FV->FE-FV coupled model run:

mpirun ~/pflotran/src/pflotran/pflotran -pflotranin CR3D.in\
mpirun ~/defmod-swpc/bin/defmod -f CR3D.inp -pc_asm_overlap 2 -fv 1 fd 1\
mpirun ~/defmod-swpc/bin/swpc_3d.x -i CR3D.inf -r CR3D -e 10

Or, submit CR3D_rup.sh and CR3D_sw.sh to a SLURM system. 

To sort the FE and FD output files, CR3D_fe.h5 and CR3D_fd.h5, run:

~/defmod-swpc/src/tool/fe_sort.py CR3D poro slip\
~/defmod-swpc/src/tool/fd_sort.py CR3D

To produce fault profile plot and wave snapshots, run:

./CR3D_fe_plot.py\
./CR3D_fd_plot.py
