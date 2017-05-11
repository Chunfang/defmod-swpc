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
        'size'   : 11}
matplotlib.rc('font', **font)
ap=argparse.ArgumentParser()
ap.add_argument('-m') # model file (no extension)
ap.add_argument('-p') # problem number
ap.add_argument('-d') # absolute output 1
name_sol = ap.parse_args().m
idcase=int(ap.parse_args().p)
if not idcase in [205,10,102]:
    print 'Choose a problem number in (205, 10, 102)'
    sys.exit(0)
if ap.parse_args().d==None:
    dsp=0
else:
    dsp=min(1,int(ap.parse_args().d))
# Read and sort SCEC wavefroms
matfile  = name_sol+'.mat'
dat_seis = np.squeeze(io_mat.loadmat(matfile)['dat_seis'])
ocoord   = np.squeeze(io_mat.loadmat(matfile)['crd_obs' ])
dt_dyn   = np.squeeze(io_mat.loadmat(matfile)['dt_dyn'  ])
dat_ES=[]; dat_FM=[]
if idcase==205:
    nobs=4; dt=.25; FDwSkp=19; FDwSep=1; id_obs=[0,2,1,3]
    label=['FD']
    sign_df=[-1,1,-1]; 
    tol=1E-3
    FDcSkp=12; FDcSep=1
    CropY=[.1,-.1]; CropZ=[.1,-.1]; ResY=300; ResZ=150; t_plt=16.
    ang_dip=np.pi/2.
elif idcase==10:
    nobs=8; dt=-.1; FDwSkp=30; FDwSep=2; id_obs=[6,5,4,0,1,2,3,7] #FDwSkp=17
    label=['FD']
    sign_df=[-1,1,-1];
    tol=1E-3
    FDcSkp=11; FDcSep=2
    CropY=[3.,-3.]; CropZ=[3.,0.]; ResY=300; ResZ=150; t_plt=20.
    ang_dip=np.pi/3.
elif idcase==102:
    nobs=2; dt=.75; FDwSkp=18; FDwSep=2; id_obs=[0,1]
    label=['DFM']
    sign_df=[-1,-1,1];
    tol=1E-1
    FDcSkp=10; FDcSep=2
    CropY=[2.5,-2.5]; CropZ=[3.,0.]; ResY=300; ResZ=150; t_plt=16.
    ang_dip=np.pi/2.
for i in range(nobs):
    if idcase==10:
        file_FM='SCEC'+str(idcase)+'/FD_'+str(i+1)+'.txt' #SGFD or FD
    else:
        file_FM='SCEC'+str(idcase)+'/SGFD_'+str(i+1)+'.txt' #SGFD or FD
    dat_FM.append(np.loadtxt(file_FM, skiprows=FDwSkp, delimiter=" "*FDwSep, unpack=False, dtype=np.float))
dat_FM=[dat_FM[i] for i in id_obs]

# Read defmod data
if idcase in [205,10]:
    dat_DF=[dat_seis.item()[i] for i in range(nobs)]
elif idcase==102:
    dat_DF=[dat_seis.item()[i] for i in [0,5]]

# Read FD data
matfile = name_sol+'_fd.mat'
dat_obs = np.squeeze(io_mat.loadmat(matfile)['dat_obs'])
dt_obs = np.squeeze(io_mat.loadmat(matfile)['dt_obs'])
nt_obs = np.squeeze(io_mat.loadmat(matfile)['nt_obs'])

# Read FD data
if idcase==205:
    dat_FD=dat_obs[[0,1,2,3],:,:]
    scl=.9
elif idcase==10:
    scl=1.
    dat_FD=dat_obs[[0,1,2,3,4,5,6,7],:,:]
elif idcase==102:
    dat_FD=dat_obs[[0,5],:,:]
    scl=1.

# Plot waveforms
id_col=[0,6,2,4]
for i in range(nobs):
    fig=plt.figure()
    dat = dat_DF[i]
    dat_fd= dat_FD[i]
    if dsp==1: dat=np.vstack((np.zeros((1,dat.shape[1])),np.diff(dat,axis=0)))/dt_dyn
    dat_fm = dat_FM[i][:,id_col]
    ymin = min(-abs(dat).max(),dat_fd.min(),dat_fm[:,1:].min())*1.05
    ymax = max( abs(dat).max(),dat_fd.max(),dat_fm[:,1:].max())*1.05
    for j in range(3):
        plt.subplot(300 + 10 + j+1)
        plt.plot(range(len(dat))*dt_dyn+dt,sign_df[j]*dat[:,j],label='FE')
        plt.plot(range(len(dat_fd))*dt_obs+dt,sign_df[j]*dat_fd[:,j]*scl,label='FE-FD')
        if j==0:    
            plt.plot(dat_fm[:,0],dat_fm[:,1+j],'--',label=label[0])
            plt.xlabel('t [s]')
            plt.ylabel('v [m/s]')
            plt.legend(loc=1)
        else:
            plt.plot(dat_fm[:,0],dat_fm[:,1+j],'--')
        plt.xlim(0,t_plt)
        if i in [0,1,2,3,4,5] and idcase==10 and j==1: plt.ylim(ymin,ymax)
    plt.savefig(name_sol+'_wave_%d'%(i+1)+'.png')
