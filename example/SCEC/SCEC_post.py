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
        'size'   : 10}
matplotlib.rc('font', **font)
ap=argparse.ArgumentParser()
ap.add_argument('-m') # model file (no extension)
ap.add_argument('-p') # problem number
ap.add_argument('-d') # absolute output 1
ap.add_argument('-r') # plot rupture contour 1
name_sol = ap.parse_args().m
idcase=int(ap.parse_args().p)
if not idcase in [205,10,102]:
    print 'Choose a problem number in (205, 10, 102)'
    sys.exit(0)
if ap.parse_args().d==None:
    dsp=0
else:
    dsp=min(1,int(ap.parse_args().d)) 
if ap.parse_args().r==None:
    rup=0
else:
    rup=min(1,int(ap.parse_args().r))
# Read and sort SCEC wavefroms
matfile  = name_sol+'.mat'
dat_seis = np.squeeze(io_mat.loadmat(matfile)['dat_seis'])
ocoord   = np.squeeze(io_mat.loadmat(matfile)['crd_obs' ])
dt_dyn   = np.squeeze(io_mat.loadmat(matfile)['dt_dyn'  ])
dat_ES=[]; dat_FM=[]
if idcase==205:
    nobs=4; dt=.25; FDwSkp=19; FDwSep=1; id_obs=[0,2,1,3]
    label=['EqSim','Roten']
    sign_df=[-1,1,-1]; 
    tol=1E-3
    FDcSkp=12; FDcSep=1
    CropY=[.1,-.1]; CropZ=[.1,-.1]; ResY=300; ResZ=150; t_plt=16.
    ang_dip=np.pi/2.
elif idcase==10:
    nobs=8; dt=0.;  FDwSkp=17; FDwSep=2; id_obs=[6,5,4,0,1,2,3,7]
    label=['EqSim','FD']
    sign_df=[-1,1,-1];
    tol=1E-3
    FDcSkp=11; FDcSep=2
    CropY=[3.,-3.]; CropZ=[3.,0.]; ResY=300; ResZ=150; t_plt=20.
    ang_dip=np.pi/3.
elif idcase==102:
    nobs=2; dt=.75; FDwSkp=18; FDwSep=2; id_obs=[0,1]
    label=['Pylith','DFM']
    sign_df=[-1,-1,1];
    tol=1E-1
    FDcSkp=10; FDcSep=2
    CropY=[2.5,-2.5]; CropZ=[3.,0.]; ResY=300; ResZ=150; t_plt=16.
    ang_dip=np.pi/2.
for i in range(nobs):
    file_ES='SCEC'+str(idcase)+'/EqSim_'+str(i+1)+'.txt'
    file_FM='SCEC'+str(idcase)+'/SGFD_'+str(i+1)+'.txt'
    dat_ES.append(np.loadtxt(file_ES, skiprows=20, delimiter="  ", unpack=False, dtype=np.float))
    dat_FM.append(np.loadtxt(file_FM, skiprows=FDwSkp, delimiter=" "*FDwSep, unpack=False, dtype=np.float))
dat_ES=[dat_ES[i] for i in id_obs]
dat_FM=[dat_FM[i] for i in id_obs]

# Read defmod data
if idcase in [205,10]:
    dat_DF=[dat_seis.item()[i] for i in range(nobs)]
elif idcase==102:
    dat_DF=[dat_seis.item()[i] for i in [0,5]]

# Plot waveforms
id_col=[0,6,2,4]
for i in range(nobs):
    fig=plt.figure()
    dat = dat_DF[i]
    if dsp==1: dat=np.vstack((np.zeros((1,dat.shape[1])),np.diff(dat,axis=0)))/dt_dyn
    dat_es = dat_ES[i][:,id_col]
    dat_fm = dat_FM[i][:,id_col]
    ymin = min(-abs(dat).max(),dat_es[:,1:].min(),dat_fm[:,1:].min())*1.05
    ymax = max( abs(dat).max(),dat_es[:,1:].max(),dat_fm[:,1:].max())*1.05
    for j in range(3):
        plt.subplot(300 + 10 + j+1)
        plt.plot(range(len(dat))*dt_dyn+dt,sign_df[j]*dat[:,j],label='Defmod')
        if j==0:    
            plt.plot(dat_es[:,0],dat_es[:,1+j],label=label[0])
            plt.plot(dat_fm[:,0],dat_fm[:,1+j],label=label[1])
            plt.xlabel('t [s]')
            plt.ylabel('v [m/s]')
            plt.legend(loc=1)
        else:
            plt.plot(dat_es[:,0],dat_es[:,1+j])
            plt.plot(dat_fm[:,0],dat_fm[:,1+j])
        plt.xlim(0,t_plt)
        plt.ylim(ymin,ymax)
    plt.savefig(name_sol+'_wave_%d'%(i+1)+'.png')

if rup==1 or dsp==1: # Plot rupture front
    fcoord   = np.squeeze(io_mat.loadmat(matfile)['crd_flt' ]) 
    dt_slip  = np.squeeze(io_mat.loadmat(matfile)['dt_slip' ])
    dat_slip = np.squeeze(io_mat.loadmat(matfile)['dat_slip'])
    if dsp==1: dat = np.linalg.norm(dat_slip.item()[:,:,:2],axis=2)
    cnt = 1E9*np.ones(shape=dat_slip.item()[0,:,0].shape)
    for i in range(len(dat)):
        t = i*dt_slip+dt
        if dsp==1:
            slip = dat[i,:]
        else:
            # Time integral for velocity output
            dat_int = np.sum(dat_slip.item()[:i+1,:,:2],axis=0)*dt_slip
            slip = np.linalg.norm(dat_int,axis=1)
        idrup = np.squeeze(np.where(slip>tol))
        idfrn = np.squeeze(np.where(cnt<1E9))
        cnt[np.setdiff1d(idrup,idfrn)] = t
    # Load SCEC data
    file_ES='SCEC'+str(idcase)+'/EqSim_c.txt'
    file_FM='SCEC'+str(idcase)+'/SGFD_c.txt'
    dat_ES=np.loadtxt(file_ES, skiprows=13, delimiter="  ", unpack=False, dtype=np.float)
    dat_FM=np.loadtxt(file_FM, skiprows=FDcSkp, delimiter=" "*FDcSep, unpack=False, dtype=np.float)
    # Plot contour
    yi = np.linspace(min(fcoord[:,1])+CropY[0], max(fcoord[:,1]+CropY[1]),ResY)
    zi = np.linspace(min(fcoord[:,2])+CropZ[0], max(fcoord[:,2]+CropZ[1]),ResZ)
    yi,zi = np.meshgrid(yi,zi)
    dip=zi/np.sin(ang_dip)
    dat_grid = griddata((np.squeeze(fcoord[:,1]), np.squeeze(fcoord[:,2])), cnt, (yi, zi), method='nearest')
    dat_grid_ES = griddata((dat_ES[:,0], dat_ES[:,1]), dat_ES[:,2], (yi*1E3, -dip*1E3), method='nearest')
    dat_grid_FM = griddata((dat_FM[:,0], dat_FM[:,1]), dat_FM[:,2], (yi*1E3, -dip*1E3), method='nearest')
    fig=plt.figure()
    cs1=plt.contour(yi, dip, dat_grid, 10, colors='b', linestyles='--')
    cs2=plt.contour(yi, dip, dat_grid_ES, 10, colors='g', linestyles='--')
    cs3=plt.contour(yi, dip, dat_grid_FM, 10, colors='r', linestyles='--',)
    lines = [ cs1.collections[0], cs2.collections[0], cs3.collections[0]]
    labels = ['Defmod',label[0],label[1]]
    plt.legend(lines, labels)
    plt.xlabel('y [km]')
    plt.ylabel('z [km]')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.savefig(name_sol+'_rup'+'.png')
