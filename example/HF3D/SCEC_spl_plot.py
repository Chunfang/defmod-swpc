#!/usr/bin/env python
import numpy as np
import sys
import scipy.io as io_mat 
import argparse
from scipy.interpolate import griddata
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import FormatStrFormatter
font = {'weight' : 'normal',
        'size'   : 12}
matplotlib.rc('font', **font)
ap=argparse.ArgumentParser()
ap.add_argument('-m') # model file (no extension)
ap.add_argument('-p') # problem number
ap.add_argument('-d') # absolute output 1
ap.add_argument('-r') # plot rupture contour 1
ap.add_argument('-fd')# Include FD result

name_sol = ap.parse_args().m
idcase=int(ap.parse_args().p)
if not idcase in [14,15]:
    print ('Choose a problem number 14 or 15')
    sys.exit(0)
if ap.parse_args().d==None:
    dsp=0
else:
    dsp=min(1,int(ap.parse_args().d)) 
if ap.parse_args().r==None:
    rup=0
else:
    rup=min(1,int(ap.parse_args().r))
if ap.parse_args().fd==None:
    fd=0
else:
    fd=min(1,int(ap.parse_args().fd))
id_obs=[7,8,9,10,0,1,2,3,4,5,6]
# Read and sort SCEC wavefroms
matfile  = name_sol+'.mat'
dat_seis = np.squeeze(io_mat.loadmat(matfile)['dat_seis'])
ocoord   = np.squeeze(io_mat.loadmat(matfile)['crd_obs' ])
dt_dyn   = np.squeeze(io_mat.loadmat(matfile)['dt_dyn'  ])
if fd:
    name_fd = name_sol+'_fd.mat'
    dat_fd=np.squeeze(io_mat.loadmat(name_fd)['dat_obs'])
    crd_fd=np.squeeze(io_mat.loadmat(name_fd)['crd_obs'])
    dt_fd =np.squeeze(io_mat.loadmat(name_fd)['dt_obs' ])
dat_DG=[]; dat_FM=[]
nobs=11; dt=-.2; FMskp=35; DGskp=22
label=['FM','DG']
sign_df=[-1,1,-1]
if idcase==14:
    CropY1=[0.,4.5];  CropY2=[-15.8,11.8]
    CropZ1=[-14.5,-.2]; CropZ2=[-14.8,-.2]
    ResY1=60; ResY2=140
    ResZ1=100;ResZ2=100
elif idcase==15:
    CropY1=[0.,10.4];  CropY2=[-15.8,2.5]
    CropZ1=[-14.5,-.2]; CropZ2=[-14.8,-.2]
    ResY1=60; ResY2=140
    ResZ1=100;ResZ2=100

t_plt=16.
for i in range(nobs):
    file_FM='SCEC%d/FM_' %(idcase) + str(i+1)+'.txt'
    file_DG='SCEC%d/DG_' %(idcase) + str(i+1)+'.txt'
    dat_FM.append(np.loadtxt(file_FM, skiprows=FMskp, unpack=False, dtype=np.float))
    dat_DG.append(np.loadtxt(file_DG, skiprows=DGskp, unpack=False, dtype=np.float))
dat_FM=[dat_FM[i] for i in id_obs]
dat_DG=[dat_DG[i] for i in id_obs]
dat_DF=[dat_seis.item()[i] for i in range(nobs)]
ang_dip=.5*np.pi
tol=1E-3

# Plot waveforms
id_col=[0,6,2,4]
for i in range(nobs):
    fig=plt.figure()
    dat = dat_DF[i]
    if dsp==1: dat=np.vstack((np.zeros((1,dat.shape[1])),np.diff(dat,axis=0)))/dt_dyn
    dat_fm = dat_FM[i][:,id_col]
    dat_dg = dat_DG[i][:,id_col]
    if fd:
        ymin = min(-abs(dat).max(),-abs(dat_fd[i,:,:]).max(),dat_fm[:,1:].min(),dat_dg[:,1:].min())*1.05
        ymax = max( abs(dat).max(), abs(dat_fd[i,:,:]).max(),dat_fm[:,1:].max(),dat_dg[:,1:].max())*1.05
    else:
        ymin = min(-abs(dat).max(),dat_fm[:,1:].min(),dat_dg[:,1:].min())*1.05
        ymax = max( abs(dat).max(),dat_fm[:,1:].max(),dat_dg[:,1:].max())*1.05
    for j in range(3):
        ax=plt.subplot(3,1,j+1)
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        plt.plot(range(len(dat))*dt_dyn+dt,sign_df[j]*dat[:,j],label='defmod',color='b',linewidth=1.)
        if fd: plt.plot(range(dat_fd.shape[1])*dt_fd+dt,sign_df[j]*dat_fd[i,:,j],label='swpc',color='m',linewidth=1.)
        if j==0:    
            plt.plot(dat_fm[:,0],dat_fm[:,1+j],label=label[0],color='g',linewidth=1.)
            plt.plot(dat_dg[:,0],dat_dg[:,1+j],label=label[1],color='r',linewidth=1.)
            plt.xlabel('t [s]')
            plt.ylabel('v [m/s]')
            plt.legend(loc=1,fontsize=10)
        else:
            plt.plot(dat_fm[:,0],dat_fm[:,1+j],color='g',linewidth=1.)
            plt.plot(dat_dg[:,0],dat_dg[:,1+j],color='r',linewidth=1.)
        if j>0: 
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
        plt.xlim(0,t_plt)
        plt.ylim(ymin,ymax)
    plt.tight_layout()
    plt.savefig(name_sol+'_wave_%d'%(i+1)+'.png')
    plt.clf()

if rup==1 or dsp==1: # Plot rupture front
    fcoord   = np.squeeze(io_mat.loadmat(matfile)['crd_flt' ]) 
    dt_slip  = np.squeeze(io_mat.loadmat(matfile)['dt_slip' ])
    dat_slip = np.squeeze(io_mat.loadmat(matfile)['dat_slip']).item()
    id2=fcoord[:,0]<tol; id1=fcoord[:,0]>=tol
    fcoord2=fcoord[id2,:]; fcoord1=fcoord[id1,:]
    dat_slip2=dat_slip[:,id2,:]; dat_slip1=dat_slip[:,id1,:]
    dat2 = np.linalg.norm(dat_slip2[:,:,:2],axis=2)
    dat1 = np.linalg.norm(dat_slip1[:,:,:2],axis=2)
    cnt1 = 1E9*np.ones(shape=dat_slip1[0,:,0].shape)
    cnt2 = 1E9*np.ones(shape=dat_slip2[0,:,0].shape)
    for i in range(len(dat1)):
        t = i*dt_slip+dt
        if dsp==1:
            slip = dat1[i,:]
        else:
            # Time integral for velocity output
            dat_int = np.sum(dat_slip1[:i+1,:,:2],axis=0)*dt_slip
            slip = np.linalg.norm(dat_int,axis=1)
        idrup = np.squeeze(np.where(slip>tol))
        idfrn = np.squeeze(np.where(cnt1<1E9))
        cnt1[np.setdiff1d(idrup,idfrn)] = t
    for i in range(len(dat2)):
        t = i*dt_slip+dt
        if dsp==1:
            slip = dat2[i,:]
        else:
            # Time integral for velocity output
            dat_int = np.sum(dat_slip2[:i+1,:,:2],axis=0)*dt_slip
            slip = np.linalg.norm(dat_int,axis=1)
        idrup = np.squeeze(np.where(slip>tol))
        idfrn = np.squeeze(np.where(cnt2<1E9))
        cnt2[np.setdiff1d(idrup,idfrn)] = t

    # Plot main fault contour
    if idcase==14:
        gs = gridspec.GridSpec(1, 2, width_ratios=[5,.965])
    elif idcase==15:
        gs = gridspec.GridSpec(1, 2, width_ratios=[1.77,1.2])
    file_FM='SCEC%d/FM_c2.txt' %(idcase)
    file_DG='SCEC%d/DG_c2.txt' %(idcase)
    dat_FM=np.loadtxt(file_FM, skiprows=FMskp, unpack=False, dtype=np.float)
    dat_DG=np.loadtxt(file_DG, skiprows=DGskp, unpack=False, dtype=np.float)
    yi = np.linspace(CropY2[0], CropY2[1],ResY2)
    zi = np.linspace(CropZ2[0], CropZ2[1],ResZ2)
    yi,zi = np.meshgrid(yi,zi)
    dip=zi/np.sin(ang_dip)
    dat_grid = griddata((np.squeeze(fcoord2[:,1]), np.squeeze(fcoord2[:,2])), cnt2, (yi, zi), method='nearest')
    dat_grid_FM = griddata((dat_FM[:,0], dat_FM[:,1]), dat_FM[:,2], (yi*1E3, -dip*1E3), method='nearest')
    dat_grid_DG = griddata((dat_DG[:,0], dat_DG[:,1]), dat_DG[:,2], (yi*1E3, -dip*1E3), method='linear')
    ax=plt.subplot(gs[0])
    if idcase==15:
        levels = np.arange(0,4.8,.5)
        cs1=plt.contour(yi, dip, dat_grid,    colors='b', linestyles='--',levels=levels,linewidths=1.)
        cs2=plt.contour(yi, dip, dat_grid_FM, colors='g', linestyles='--',levels=levels,linewidths=1.)
        cs3=plt.contour(yi, dip, dat_grid_DG, colors='r', linestyles='--',levels=levels,linewidths=1.)
    elif idcase==14:
        levels = np.arange(0,9.6,.5)
        cs1=plt.contour(yi, dip, dat_grid,    colors='b', linestyles='--',levels=levels,linewidths=1.)
        cs2=plt.contour(yi, dip, dat_grid_FM, colors='g', linestyles='--',levels=levels,linewidths=1.)
        cs3=plt.contour(yi, dip, dat_grid_DG, colors='r', linestyles='--',levels=levels,linewidths=1.)
    lines = [cs1.collections[0], cs2.collections[0], cs3.collections[0]]
    labels = ['defmod',label[0],label[1]]
    if idcase==14:
        loc=4
    elif idcase==15:
        loc=2
    plt.legend(lines, labels,loc=loc, fontsize=10)
    plt.xlabel('strik [km]')
    plt.ylabel('dip [km]')
    plt.gca().set_aspect('equal', adjustable='box')
    #plt.savefig(name_sol+'_rup2'+'.png')
    
    # Plot splay fault contour 
    file_FM='SCEC%d/FM_c1.txt' %(idcase) 
    file_DG='SCEC%d/DG_c1.txt' %(idcase)
    dat_FM=np.loadtxt(file_FM, skiprows=FMskp, unpack=False, dtype=np.float)
    dat_DG=np.loadtxt(file_DG, skiprows=DGskp, unpack=False, dtype=np.float)
    yi = np.linspace(CropY1[0], CropY1[1],ResY1)
    zi = np.linspace(CropZ1[0], CropZ1[1],ResZ1)
    yi,zi = np.meshgrid(yi,zi)
    stk=yi/np.cos(np.pi/6.)
    dip=zi/np.sin(ang_dip)
    
    # DG data completion
    i_cnt=cnt1>=1E9
    cnt_end=cnt1[i_cnt]
    dat_add=np.hstack(((fcoord1[i_cnt,1]/np.cos(np.pi/6.)*1E3).reshape(sum(i_cnt),1),
                        -fcoord1[i_cnt,2].reshape(sum(i_cnt),1)*1E3,
                        cnt1[i_cnt].reshape(sum(i_cnt),1)))
    dat_DG=np.vstack((dat_DG, dat_add))
    
    dat_grid = griddata((np.squeeze(fcoord1[:,1]), np.squeeze(fcoord1[:,2])), cnt1, (yi, zi), method='nearest')
    dat_grid_FM = griddata((dat_FM[:,0], dat_FM[:,1]), dat_FM[:,2], (stk*1E3, -dip*1E3), method='nearest')
    dat_grid_DG = griddata((dat_DG[:,0], dat_DG[:,1]), dat_DG[:,2], (stk*1E3, -dip*1E3), method='linear')
    levels = np.arange(0,6.4,.5)     
    ax=plt.subplot(gs[1])
    cs1=plt.contour(stk, dip, dat_grid,    colors='b', linestyles='--',levels=levels,linewidths=1.)
    cs2=plt.contour(stk, dip, dat_grid_FM, colors='g', linestyles='--',levels=levels,linewidths=1.)
    cs3=plt.contour(stk, dip, dat_grid_DG, colors='r', linestyles='--',levels=levels,linewidths=1.)
    lines = [cs1.collections[0], cs2.collections[0], cs3.collections[0]]
    labels = ['defmod',label[0],label[1]]
    #plt.legend(lines, labels)
    #plt.xlabel('strike [km]')
    #plt.ylabel('dip [km]')
    plt.gca().set_yticklabels(['']*10)
    plt.xticks(np.arange(stk.min(), stk.max(), 2))
    plt.gca().set_aspect('equal', adjustable='box')
    plt.tight_layout()
    plt.savefig(name_sol+'_rup'+'.png',bbox_inches='tight')
