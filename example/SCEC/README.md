Batch run SCEC benchmark problems 10, 102 and 205 in FE-FD hybrid on cluster, and produce waveform and rupture contour comparisons against other participants in (Meng 2017).

Lazy: Must have Cubit/Trelis, Python (numpy, scipy, netCDF4, Anaconda 2.7/3.6 tested) on cluster, otherwise follow "Recommended". Note, this approach will waste cluster hours as the Cubit/Python parts run in serial.
1. Upload SCEC folder to $SCRATCH partition on cluster, and cd there.
2. Modify "SCEC_rup.sh" in executable paths, and then submit.
3. After exiting successfully, modify "SCEC_sw.sh" (MPI size has to be consistent with inf and cfg files, or don't change), and submit.

Recommended: Must have Cubit/Trelis, Python (numpy, scipy, netCDF4, Anaconda 2.7/3.6 tested) on desktop. 
1. On desktop, cd to SCEC folder, and modify/run the Cubit/Trelis section in the SCEC_rup.sh file.
2. Modify/run the Python section (numpy, scipy, netCDF4 required) in the SCEC_rup.sh file.
3. Upload all inp, inf, stloc.xy, cfg files and SCECXXX folders to a $SCRATCH folder on cluster, and cd there.
4. Modify the SCEC_rup.sh file (comment out Cubit/Python sections that have been executed), upload to the same folder, submit, and wait for successful exit.
6, Note, use ssh -Y XXX to logon cluster. Modify, upload and submit SCEC_sw.sh (MPI size has to be consistent with inf and cfg files, or don't change).

After successful executions, waveform and contour figures should appear in the same folder. 
