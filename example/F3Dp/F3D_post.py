#!/usr/bin/env python
import numpy as np
import sys
import scipy.io as io_mat 
import argparse
from scipy.interpolate import griddata
import matplotlib
matplotlib.use('Svg')
import matplotlib.pyplot as plt
font = {'weight' : 'normal',
        'size'   : 12}
matplotlib.rc('font', **font)
name_sol = sys.argv[1]
name_fe = name_sol+'.mat'
dat_fe=np.squeeze(io_mat.loadmat(name_fe)['dat_seis'])
crd_fe=np.squeeze(io_mat.loadmat(name_fe)['crd_obs' ])
dt_fe =np.squeeze(io_mat.loadmat(name_fe)['dt_dyn'  ])

name_fd = name_sol+'_fd.mat'
dat_fd=np.squeeze(io_mat.loadmat(name_fd)['dat_obs'])
crd_fd=np.squeeze(io_mat.loadmat(name_fd)['crd_obs'])
dt_fd =np.squeeze(io_mat.loadmat(name_fd)['dt_obs' ])

# Waveform comparisons
eid=10
dat_fe=dat_fe[eid]
plt.figure(figsize=(16, 12), dpi=80)
xplt_fe=range(dat_fe.shape[1])*dt_fe 
xplt_fd=range(dat_fd.shape[1])*dt_fd
for i in range(5): 
    for j in range(3):
        plt.subplot(5,3,i*3+j+1) 
        plt.plot(xplt_fe,dat_fe[i,:,j])
        plt.plot(xplt_fd,dat_fd[i,:,j])
        plt.xlim([0,1])
        if (i>0 or j>0): plt.gca().axes.get_xaxis().set_visible(False)
plt.savefig(name_sol+'_wf.png')
