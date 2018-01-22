Batch run SCEC benchmark problems 10, 102 and 205 in FE-FD hybrid on cluster, and produce waveform and rupture contour comparisons against other participants in (Meng 2017).

Lazy:
Must have Cubit/Trelis, python (numpy, scipy, netCDF4, Anaconda 2.7 tested) on cluster, otherwise follow "Recommended". Note, this approach will waste cluster hours as the Cubit/python parts run in serial.
1, Upload SCEC folder to $SCRATCH partition on cluster, and cd there.
2, Modify "SCEC_rup.sh" in executable paths, "my OS/partition" and MPI size, and then submit.
3, After exiting successfully, modify "SCEC_sw.sh" (MPI size has to be consistent with inf and cfg files, or don't change), and submit.

Recommended:
1, On desktop, cd to SCEC folder, and modify/run SCEC_msh.sh, Cubit/Trelis required.
2, Modify/run SCEC_prep.sh, python (numpy, scipy, netCDF4) required.
3, Upload all inp, inf, stloc.xy, cfg files to a $SCRATCH folder on cluster, and cd there.
4, Modify and upload SCEC_rup.sh file (comment out Cubit/python parts), submit, and wait for successful exit.
5, Modify, upload and submit SCEC_sw.sh (comment out Cubit/python parts, MPI size has to be consistent with inf and cfg files, or don't change).
6, Download all resulting txt and log files to SCEC folder on desktop, and cd there. 
7, Modify/run SCEC_sort.sh to sort the outcomes into mat files.
8, Modify/Run SCEC_plot.sh to produce plots.
