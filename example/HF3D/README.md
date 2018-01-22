Duplicate splay fault examples, SCEC14, 15 and HF3D (injection trigger), in (Meng 2018).

Lazy: Cubit/Trelis, python (netCDF4, numpy scipy) and avconv required on cluster, otherwise follow "recommended". Note, this approach will waste cluster time as cubit and python parts run in serial. 
1. Copy this folder to $SCRATCH partition on cluster.
2. Modify/submit HF3D_rup.sh. 
3. Once exit successfully, modify/submit HF3D_sw.sh (MPI size has to be consistent with inf and cfg files, or don't change).

Recommended: Cubit/Trelis, python (netCDF4, numpy scipy) and avconv required on desktop.
1. Modify and run HF3D_prep.sh 
2. Upload inp, inf, stloc.xy, cfg files to cluster.
3. Modify HF3D_rup.sh (comment out all the cubit and python parts), upload and submit.
4. Once exit successfully, modify HF3D_sw.sh (comment out all the cubit and python parts, MPI size has to be consistent with inf and cfg files, or don't change), upload and submit.
5. Download all the resulting txt and log files and snp folders to this folder.
6. Modify/run HF3D_sort.sh to sort the outcomes to mat files.
7. Modify/run HF3D_plot.sh to produce plots and movies.
