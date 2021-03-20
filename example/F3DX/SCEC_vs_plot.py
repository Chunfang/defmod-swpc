#!/usr/bin/env python
import numpy as np
import h5py
import sys
import scipy.io as io_mat 
import argparse
from scipy.interpolate import griddata
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import FormatStrFormatter

def load_dict_from_hdf5(filename):
    with h5py.File(filename, 'r') as h5file:
        return recursively_load_dict_contents_from_group(h5file, '/')

def recursively_load_dict_contents_from_group(h5file, path):
    ans = {}
    for key, item in h5file[path].items():
        if isinstance(item, h5py._hl.dataset.Dataset):
            ans[key] = item.value
        elif isinstance(item, h5py._hl.group.Group):
            ans[key] = recursively_load_dict_contents_from_group(h5file, path + key + '/')
    return ans


font = {'weight' : 'normal',
        'size'   : 12}
matplotlib.rc('font', **font)
ap=argparse.ArgumentParser()
ap.add_argument('-d') # absolute output 1
ap.add_argument('-r') # plot rupture contour 1
ap.add_argument('-p') # problem number 14/15


if ap.parse_args().d==None:
    dsp=0
else:
    dsp=min(1,int(ap.parse_args().d)) 
if ap.parse_args().r==None:
    rup=0
else:
    rup=min(1,int(ap.parse_args().r))
if ap.parse_args().p==None:
    idcase=14
else:
    idcase=int(ap.parse_args().p)
name_sol = 'SCEC'+str(idcase)

id_obs=[7,8,9,10,0,1,2,3,4,5,6]
# Read and sort SCEC wavefroms
matfile  = 'd100/'+name_sol+'.mat'
dat_seis = np.squeeze(io_mat.loadmat(matfile)['dat_seis'])
ocoord   = np.squeeze(io_mat.loadmat(matfile)['crd_obs' ])
dt_dyn   = np.squeeze(io_mat.loadmat(matfile)['dt_dyn'  ])

# Read F3DX14 data
dat_F3DX=load_dict_from_hdf5('d200/'+name_sol+'_fe.h5')
dat_seis2=dat_F3DX['obs_dyn']['step 1']

nobs=11; dt=-.2; sign_df=[-1,1,-1]
t_plt=16.
dat_DF=[dat_seis.item()[i] for i in range(nobs)]
ang_dip=.5*np.pi
tol=1E-3

if idcase==14:
    CropY1=[0.,4.5]; CropY2=[-15.8,11.8]
    CropY1=[0.,7.5];
    CropZ1=[-14.5,-.2]; CropZ2=[-14.8,-.2]
    ResY1=60; ResY2=140
    ResZ1=100;ResZ2=100
elif idcase==15:
    CropY1=[0.,10.4]; CropY2=[-15.8,11.8] #CropY2=[-15.8,2.5]
    CropZ1=[-14.5,-.2]; CropZ2=[-14.8,-.2]
    ResY1=60; ResY2=140
    ResZ1=100;ResZ2=100

# Plot waveforms
id_col=[0,6,2,4]
for i in range(nobs):
    fig=plt.figure()
    dat = dat_DF[i]
    dat2 = dat_seis2[i,:,:]
    if dsp==1: 
        dat=np.vstack((np.zeros((1,dat.shape[1])),np.diff(dat,axis=0)))/dt_dyn
        dat2=np.hstack((np.zeros((dat2.shape[0],1)),np.diff(dat2,axis=1)))/dt_dyn
    ymin = min(-abs(dat).max(),-abs(dat2).max())*1.05
    ymax = -ymin 
    for j in range(3):
        ax=plt.subplot(3,1,j+1)
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        plt.plot(range(len(dat))*dt_dyn+dt,sign_df[j]*dat[:,j],label='d100',color='b',linewidth=1.)
        plt.plot(range(dat2.shape[1])*dt_dyn+dt,sign_df[j]*dat2[j,:],label='d200',color='g',linewidth=1.)
        if j==0:    
            plt.xlabel('t [s]')
            plt.ylabel('v [m/s]')
            plt.legend(loc=1,fontsize=10)
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

    dat_slipX = dat_F3DX['slip_dyn']['step 1']
    fcoordX = dat_F3DX['crd_flt'][:,:3]

    id2=fcoord[:,0]<tol; id1=fcoord[:,0]>=tol
    fcoord2=fcoord[id2,:]; fcoord1=fcoord[id1,:]
    dat_slip2 = dat_slip[:,id2,:]; dat_slip1=dat_slip[:,id1,:]

    id2=abs(fcoordX[:,0])<tol; id1=abs(fcoordX[:,0])>=tol
    fcoordX2=fcoordX[id2,:]; fcoordX1=fcoordX[id1,:]
    dat_slipX2 = dat_slipX[id2,:,:]; dat_slipX1=dat_slipX[id1,:,:]

    dat2 = np.linalg.norm(dat_slip2[:,:,:2],axis=2)
    dat1 = np.linalg.norm(dat_slip1[:,:,:2],axis=2)
    cnt1 = 1E9*np.ones(shape=dat_slip1[0,:,0].shape)
    cnt2 = 1E9*np.ones(shape=dat_slip2[0,:,0].shape)

    datX2 = np.linalg.norm(dat_slipX2[:,:2,:],axis=1)
    datX1 = np.linalg.norm(dat_slipX1[:,:2,:],axis=1)
    cntX1 = 1E9*np.ones(shape=dat_slipX1[:,0,0].shape)
    cntX2 = 1E9*np.ones(shape=dat_slipX2[:,0,0].shape)

    for i in range(min(len(dat1),datX1.shape[1])):
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
        
        slipX  =  datX1[:,i]
        idrupX = np.squeeze(np.where(slipX>tol))
        idfrnX = np.squeeze(np.where(cntX1<1E9))
        cntX1[np.setdiff1d(idrupX,idfrnX)] = t    

    for i in range(min(len(dat2),datX2.shape[1])):
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
        
        slipX  =  datX2[:,i]
        idrupX = np.squeeze(np.where(slipX>tol))
        idfrnX = np.squeeze(np.where(cntX2<1E9))
        cntX2[np.setdiff1d(idrupX,idfrnX)] = t

    # Plot main fault contour
    if idcase == 14:
        gs = gridspec.GridSpec(1, 2, width_ratios=[5,1.6083]) #.965
    if idcase == 15:
        #gs = gridspec.GridSpec(1, 2, width_ratios=[1.77,1.2])
        gs = gridspec.GridSpec(1, 2, width_ratios=[1.77*1.52,1.2])
    fig = plt.figure(figsize=(7., 3.5))
    yi = np.linspace(CropY2[0], CropY2[1],ResY2)
    zi = np.linspace(CropZ2[0], CropZ2[1],ResZ2)
    yi,zi = np.meshgrid(yi,zi)
    dip=zi/np.sin(ang_dip)
    dat_grid = griddata((np.squeeze(fcoord2[:,1]), np.squeeze(fcoord2[:,2])), cnt2, (yi, zi), method='nearest')
    dat_gridX= griddata((np.squeeze(fcoordX2[:,1]), np.squeeze(fcoordX2[:,2])), cntX2, (yi, zi), method='nearest')
    ax=plt.subplot(gs[0])
    if idcase==15:
        #levels = np.arange(0,4.8,.5)
        levels = np.arange(0,9.6,.5)
        cs1=plt.contour(yi, dip, dat_grid,  10, colors='b', linestyles='--',levels=levels,linewidths=1.)
        cs2=plt.contour(yi, dip, dat_gridX, 10, colors='g', linestyles='--',levels=levels)
    elif idcase==14:
        levels = np.arange(0,9.6,.5)
        cs1=plt.contour(yi, dip, dat_grid,  10, colors='b', linestyles='--',levels=levels)
        cs2=plt.contour(yi, dip, dat_gridX, 10, colors='g', linestyles='--',levels=levels)
    lines = [cs1.collections[0],cs2.collections[0]]
    labels = ['d100','d200']
    if idcase==14:
        loc=4
    elif idcase==15:
        loc=4
    plt.legend(lines, labels, loc=loc, fontsize=10)
    plt.xlabel('strik [km]')
    plt.ylabel('dip [km]')
    plt.gca().set_aspect('equal', adjustable='box')
    
    # Plot splay fault contour 
    yi = np.linspace(CropY1[0], CropY1[1],ResY1)
    zi = np.linspace(CropZ1[0], CropZ1[1],ResZ1)
    yi,zi = np.meshgrid(yi,zi)
    stk=yi/np.cos(np.pi/6.)
    dip=zi/np.sin(ang_dip)
    
    i_cnt=cnt1>=1E9
    cnt_end=cnt1[i_cnt]
    
    dat_grid = griddata((np.squeeze(fcoord1[:,1]), np.squeeze(fcoord1[:,2])), cnt1, (yi, zi), method='nearest')
    dat_gridX= griddata((np.squeeze(fcoordX1[:,1]), np.squeeze(fcoordX1[:,2])), cntX1, (yi, zi), method='nearest')
    levels = np.arange(0,6.4,.25)     
    ax=plt.subplot(gs[1])
    cs1=plt.contour(stk, dip, dat_grid,  colors='b', linestyles='--',levels=levels)
    cs2=plt.contour(stk, dip, dat_gridX, colors='g', linestyles='--',levels=levels)
    lines = [cs1.collections[0]]
    labels = ['d100', 'd200']
    plt.gca().set_yticklabels(['']*10)
    plt.xticks(np.arange(stk.min(), stk.max(), 5.))
    plt.gca().set_aspect('equal', adjustable='box')
    plt.tight_layout()
    plt.savefig(name_sol+'_rup'+'.png',bbox_inches='tight')
