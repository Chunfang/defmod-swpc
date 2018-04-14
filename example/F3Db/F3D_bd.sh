#!/bin/bash
trelis -nojournal -nographics F3D_bd.jou
./F3D_bd.py F3D_bd.exo
mpirun -n 8 $HOME/pflotran/src/pflotran/pflotran -pflotranin F3D_bd.in
mpirun -n 8 $HOME/defmod-swpc/bin/defmod -f F3D_bd.inp -pc_asm_overlap 2 -fd 1 -fv 1
export I_MPI_PIN_DOMAIN=omp
export OMP_NUM_THREADS=2
mpirun -n 4 $HOME/defmod-swpc/bin/swpc_3d.x -i F3D_bd.inf -r F3D_bd -e 2 
