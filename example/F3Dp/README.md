Duplicate 3D production induced seismicity model (Meng 2018), cubit/trelis, python (numpy, scipy, netCDF4), meshlab and avconv required 
1, Run F3D_syn.py. If mesh file F3D_syn.exo already exist, but need to change parameters in def_syn.inp, modify F3D_msh2inp.py, and run ./F3D_msh2inp.py F3D_syn.exo 2E-5 
2, Upload F3D_rup.sh, F3D_syn.inp and F3D_syn_fefd.cfg to a cluster folder, and cd there.
3, Modify F3D_rup.sh, and submit. 
4, After exit successfully, upload F3D_sw.sh lhm.dat stloc.xy and F3D_syn.inf to that folder, modify F3D_sw.sh (MPI size has to be consistent with inf and cfg files, or don't change) and submit.
5, Download all the resulting txt, log files and out folder to this folder, and modify/run F3D_sort.py to produce mat files. 
6, defmod-swpc/src/tool/def_plot.py F3D_syn slip 
7, ./fd_plot.py  
