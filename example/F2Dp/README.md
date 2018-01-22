Duplicate 2D production fluctuation result in (Meng 2018): cubit/trelis required.
1. cubit -nojournal -nographics F2D.jou 
2. ./F2Dp.sh 
3. Upload inp files to a cluster folder, and cd there (optional). 
4. Run commands in F2Dp_sub.sh, or upload it to the cluster folder and modify/submit (comment out python parts if if needed).
5. If run on cluster but not the python part, download all the resulting txt and log files, and run F2D_sort.sh. If python was run on cluster, download all the resulting mat files. 
7. ./F2Dp_plot.py F2D_dt12
