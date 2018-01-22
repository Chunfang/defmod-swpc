#!/usr/bin/env python
import numpy as np
import sys
import scipy.io as io_mat 
import argparse
from scipy.interpolate import griddata
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec
font = {'weight' : 'normal',
        'size'   : 12}
matplotlib.rc('font', **font)
ap=argparse.ArgumentParser()
ap.add_argument('-m')     # model file (no extension)
ap.add_argument('-p')     # problem number
ap.add_argument('-d')     # absolute output 1
ap.add_argument('-r')     # plot rupture contour 1
ap.add_argument('-alpha') # plot implicit waveform/rupture 
name_sol = ap.parse_args().m
idcase=int(ap.parse_args().p)
if not idcase in [10, 11, 102]:
    print ('Choose a problem number 10, 11 or 102')
    sys.exit(0)
if ap.parse_args().d==None:
    dsp=0
else:
    dsp=min(1,int(ap.parse_args().d)) 
if ap.parse_args().r==None:
    rup=0
else:
    rup=min(1,int(ap.parse_args().r))
if ap.parse_args().alpha==None:
    alpha=0
else:
    alpha = min(1,int(ap.parse_args().alpha))

id_obs=[9,8,6,1,3,4,7,5,0,2]
# Read and sort SCEC wavefroms
matfile  = name_sol+'.mat'
dat_seis = np.squeeze(io_mat.loadmat(matfile)['dat_seis'])
ocoord   = np.squeeze(io_mat.loadmat(matfile)['crd_obs' ])
dt_dyn   = np.squeeze(io_mat.loadmat(matfile)['dt_dyn'  ])
if alpha: 
    dt_alpha  = np.squeeze(io_mat.loadmat(matfile)['dt_alpha' ])
    dat_alpha = np.squeeze(io_mat.loadmat(matfile)['dat_alpha'])

if idcase in [10, 11]:
    dat_FM=[]
    nobs=10; dt=0.; FMskp=35 
    if idcase == 11: dt = 0.075; 
    label=['Barall']
    sign_df=[-1,-1]
    t_plt=16.
    for i in range(nobs):
        file_FM='SCEC'+str(idcase)+'-2D/FM2D_'+str(i+1)+'.txt'
        dat_FM.append(np.loadtxt(file_FM, skiprows=FMskp, unpack=False, dtype=np.float))
    dat_FM=[dat_FM[i] for i in id_obs]
else:
    nobs = dat_seis.item().shape[0]; dt=0.; sign_df=[-1, -1]; t_plt=16. 

dat_DF=[dat_seis.item()[i] for i in range(nobs)]
if alpha: dat_DF_alpha = [dat_alpha[i,:,:] for i in range(nobs)]

# Plot waveforms
id_col=[0,6,4]
for i in range(nobs):
    fig=plt.figure()
    dat = dat_DF[i]
    if alpha: dat_alpha = dat_DF_alpha[i]
    if dsp==1: 
        dat = np.vstack((np.zeros((1,dat.shape[1])), np.diff(dat,axis=0)))/dt_dyn
        if alpha: dat_alpha = np.vstack((np.zeros((1,dat_alpha.shape[1])), np.diff(dat_alpha,axis=0)))/dt_alpha
    if idcase in [10, 11]:
        dat_fm = dat_FM[i][:,id_col]
        ymin = min(-abs(dat).max(), dat_fm[:,1:].min())*1.05
        ymax = max( abs(dat).max(), dat_fm[:,1:].max())*1.05
    else:
        ymin = -abs(dat).max()*1.05
        ymax =  abs(dat).max()*1.05
    for j in range(2):
        plt.subplot(2,1,j+1)
        plt.plot(range(len(dat))*dt_dyn+dt,sign_df[j]*dat[:,j],label='Defmod')
        if alpha: plt.plot(range(len(dat_alpha))*dt_alpha+dt,sign_df[j]*dat_alpha[:,j],label='alpha')
        if j==0:    
            if idcase in [10, 11]: plt.plot(dat_fm[:,0],dat_fm[:,1+j],label=label[0])
            plt.xlabel('t [s]')
            plt.ylabel('v [m/s]')
            plt.legend(loc=1)
        else:
            if idcase in [10, 11]: plt.plot(dat_fm[:,0],dat_fm[:,1+j])
        plt.xlim(0,t_plt)
        plt.ylim(ymin,ymax)
    plt.tight_layout()
    plt.savefig(name_sol+'_wave_%d'%(i+1)+'.png')

if rup:
    tmax = 5. #.8 
    dt =  .5   #0.08 
    fcoord   = np.squeeze(io_mat.loadmat(matfile)['crd_flt' ])
    dt_slip  = np.squeeze(io_mat.loadmat(matfile)['dt_slip' ])
    dat_slip = np.squeeze(io_mat.loadmat(matfile)['dat_slip']).item()
    fsort = np.argsort(fcoord[:,0]) 
    fcoord = fcoord[fsort,:]
    nmax = min(dat_slip.shape[0], int(tmax/dt_slip))
    dat_sliprate = dat_slip[:nmax+1,:,:]
    # Differentiate 
    if dsp: dat_sliprate = np.vstack((np.zeros((1,dat_sliprate.shape[1],dat_sliprate.shape[2])), np.diff(dat_sliprate,axis=0)))/dt_dyn
    nviz = int(dt/dt_slip)
    n_plt = dat_sliprate.shape[0]/nviz 
    fig = plt.figure()
    # dip and strike
    if idcase in [10 ,11]:
        x_plt = np.linalg.norm(fcoord,axis=1)
    elif idcase==102:
        x_plt = fcoord[:,0]
    for i in range(n_plt):
        dat_plt = dat_sliprate[i*nviz,fsort,0]
        if i==0:
            plt.plot(x_plt, dat_plt, label='Defmod', linestyle='-', color='royalblue')
            if idcase in [10, 11]:
                plt.xlabel('dip [km]')
            elif idcase==102:
                plt.xlabel('strike [km]')
            plt.ylabel('slip rate [m/s]')
        else:
            plt.plot(x_plt, dat_plt, linestyle='-', color='royalblue')
    if alpha:
        dt_slip_alpha  = np.squeeze(io_mat.loadmat(matfile)['dt_slip_alpha' ])
        dat_slip_alpha = np.squeeze(io_mat.loadmat(matfile)['dat_slip_alpha'])
        nmax_alpha = min(dat_slip_alpha.shape[0], int(tmax/dt_slip_alpha))
        dat_sliprate_alpha = dat_slip_alpha[:nmax_alpha+1,:,:]
        # Differentiate
        if dsp: dat_sliprate_alpha = np.vstack((np.zeros((1,dat_sliprate_alpha.shape[1],dat_sliprate_alpha.shape[2])), np.diff(dat_sliprate_alpha,axis=0)))/dt_slip_alpha
        nviz_alpha = int(dt/dt_slip_alpha)
        n_plt_alpha = dat_sliprate_alpha.shape[0]/nviz_alpha 
        for i in range(n_plt_alpha):
            dat_plt = dat_sliprate_alpha[i*nviz_alpha,fsort,0]
            if i==0:
                plt.plot(x_plt,dat_plt, label='alpha', color='coral')
            else:
                plt.plot(x_plt,dat_plt, color='coral')
    plt.legend(loc=9)
    plt.tight_layout()
    plt.savefig(name_sol+'_sliprate.png')

    # Plot displacement only when dsp==1
    if dsp:
        fig = plt.figure()
        dat_slip = dat_slip[:nmax+1,:,:]
        for i in range(n_plt):
            dat_plt = dat_slip[i*nviz,fsort,0]
            if i==0:
                plt.plot(x_plt,dat_plt, label='Defmod', linestyle='-', color='royalblue')
                if idcase in [10, 11]:
                    plt.xlabel('dip [km]')
                elif idcase==102:
                    plt.xlabel('strike [km]')
                plt.ylabel('slip [m]')
            else:
                plt.plot(x_plt,dat_plt, linestyle='-', color='royalblue')   
        if alpha:
            dat_slip_alpha = dat_slip_alpha[:nmax_alpha+1,:,:]
            for i in range(n_plt_alpha):
                dat_plt = dat_slip_alpha[i*nviz_alpha,fsort,0]
                if i==0:
                    plt.plot(x_plt,dat_plt, label='alpha', color='coral')
                else:
                    plt.plot(x_plt,dat_plt, color='coral')
        plt.legend(loc=1)
        plt.tight_layout()
        plt.savefig(name_sol+'_slip.png')
