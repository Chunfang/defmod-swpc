#!/usr/bin/env python
import numpy as np
import h5py
import sys
import scipy.io as io_mat
import argparse
from scipy.interpolate import griddata
import matplotlib
#matplotlib.use('Svg')
import matplotlib.pyplot as plt

def load_dict_from_hdf5(filename):
    with h5py.File(filename, 'r') as h5file:
        return recursively_load_dict_contents_from_group(h5file, '/')

def recursively_load_dict_contents_from_group(h5file, path):
    ans = {}
    for key, item in h5file[path].items():
        if isinstance(item, h5py._hl.dataset.Dataset):
            ans[key] = item[()]
        elif isinstance(item, h5py._hl.group.Group):
            ans[key] = recursively_load_dict_contents_from_group(h5file, path + key + '/')
    return ans

font = {'weight' : 'normal',
        'size'   : 12}
matplotlib.rc('font', **font)
ap=argparse.ArgumentParser()
ap.add_argument('-m') # model file (no extension)
ap.add_argument('-p') # problem number
ap.add_argument('-d') # absolute output 1
ap.add_argument('-r') # plot rupture contour 1
ap.add_argument('-h5')# choose hdf5 input
ap.add_argument('-fd')# Include FD result

name_sol = ap.parse_args().m
idcase=int(ap.parse_args().p)
if not idcase in [205,10,102]:
    print ('Choose a problem number in (205, 10, 102)')
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
if ap.parse_args().h5==None:
    h5=1
else:
    h5=min(1,int(ap.parse_args().h5))

# Read and sort SCEC wavefroms
if h5:
    dat_h5=load_dict_from_hdf5(name_sol+'_fe.h5')
    dt_dyn  = dat_h5['dt_dyn']
    ocoord  = dat_h5['crd_obs']
    dat_seis= dat_h5['obs_dyn']['step 1']
else:
    matfile  = name_sol+'.mat'
    dat_seis = np.squeeze(io_mat.loadmat(matfile)['dat_seis'])
    ocoord   = np.squeeze(io_mat.loadmat(matfile)['crd_obs' ])
    dt_dyn   = np.squeeze(io_mat.loadmat(matfile)['dt_dyn'  ])
dat_ES=[]; dat_FM=[]
if idcase==205:
    nobs=4; dt=.25; FDwSkp=19; FDwSep=1; id_obs=[0,2,1,3]
    label=['EqSim','FD']
    sign_df=[-1,1,-1];
    tol=1E-3
    FDcSkp=12; FDcSep=1
    CropY=[.1,-.1]; CropZ=[.1,-.1]; ResY=300; ResZ=150; t_plt=16.
    ang_dip=np.pi/2.
elif idcase==10:
    nobs=8; dt=0.; FDwSkp=30; FDwSep=2; id_obs=[6,5,4,0,1,2,3,7]
    label=['EqSim','FD']
    sign_df=[-1,1,-1];
    tol=1E-3
    FDcSkp=11; FDcSep=2
    CropY=[3.,-3.]; CropZ=[3.,0.]; ResY=300; ResZ=150; t_plt=20.
    ang_dip=np.pi/3.
elif idcase==102:
    nobs=2; dt=.625;
    FDwSkp=18; FDwSep=2; id_obs=[0,1]
    label=['Pylith','DFM']
    sign_df=[-1,-1,1];
    tol=1E-1
    FDcSkp=10; FDcSep=2
    CropY=[2.5,-2.5]; CropZ=[3.,0.]; ResY=300; ResZ=150; t_plt=16.
    ang_dip=np.pi/2.
for i in range(nobs):
    file_ES='SCEC'+str(idcase)+'/EqSim_'+str(i+1)+'.txt'
    if idcase==10:
        file_FM='SCEC'+str(idcase)+'/FD_'+str(i+1)+'.txt' #SGFD or FD
    else:
        file_FM='SCEC'+str(idcase)+'/SGFD_'+str(i+1)+'.txt' #SGFD or FD
    dat_ES.append(np.loadtxt(file_ES, skiprows=20, delimiter="  ", unpack=False, dtype=np.float))
    dat_FM.append(np.loadtxt(file_FM, skiprows=FDwSkp, delimiter=" "*FDwSep, unpack=False, dtype=np.float))
dat_ES=[dat_ES[i] for i in id_obs]
dat_FM=[dat_FM[i] for i in id_obs]

# Read defmod data
if idcase in [205,10]:
    if h5:
        dat_DF=dat_seis
    else:
        dat_DF=[dat_seis.item()[i] for i in range(nobs)]
elif idcase==102:
    if h5:
        dat_DF=dat_seis[[0,5],:,:]
    else:
        dat_DF=[dat_seis.item()[i] for i in [0,5]]

# Read FD data
if fd:
    if h5:
        dat_FD  = load_dict_from_hdf5(name_sol+'_fd.h5')
        dat_obs = dat_FD['dat_obs_fd']
        dt_obs  = dat_FD['dt_obs_fd']
        nt_obs  = dat_FD['nt_obs_fd']
    else:
        fileFD = name_sol+'_fd.mat'
        dat_obs = np.squeeze(io_mat.loadmat(fileFD)['dat_obs'])
        dt_obs  = np.squeeze(io_mat.loadmat(fileFD)['dt_obs'])
        nt_obs  = np.squeeze(io_mat.loadmat(fileFD)['nt_obs'])
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
    if h5:
        dat = dat_DF[i,:,:]
        if fd: dat_fd = dat_FD[i,:,:]
    else:
        dat = dat_DF[i]
        if fd: dat_fd= dat_FD[i]
    # Time derivative for displacement output
    if dsp==1:
        if h5:
            dat=np.hstack((np.zeros((dat.shape[0],1)),np.diff(dat,axis=1)))/dt_dyn
        else:
            dat=np.vstack((np.zeros((1,dat.shape[1])),np.diff(dat,axis=0)))/dt_dyn
    dat_es = dat_ES[i][:,id_col]
    dat_fm = dat_FM[i][:,id_col]
    ymin=0.; ymax=0.
    for j in range(3):
        if fd:
            ymin_i = min((sign_df[j]*dat[j,:]).min(),(sign_df[j]*dat_fd[j,:]*scl).min(),dat_es[:,1:].min(),dat_fm[:,1:].min())*1.05
            ymax_i = max((sign_df[j]*dat[j,:]).max(),(sign_df[j]*dat_fd[j,:]*scl).max(),dat_es[:,1:].max(),dat_fm[:,1:].max())*1.05
        else:
            ymin_i = min((sign_df[j]*dat[j,:]).min(),dat_fm[:,1:].min())*1.05
            ymax_i = max((sign_df[j]*dat[j,:]).min(),dat_fm[:,1:].max())*1.05
        ymin=min(ymin,ymin_i)
        ymax=max(ymax,ymax_i)
        ax = plt.subplot(3,1,j+1)
        plt.xlim(0,t_plt)
        plt.ylim(ymin,ymax)
        if h5:
            plt.plot(np.asarray(range(dat.shape[1]))*dt_dyn+dt,sign_df[j]*dat[j,:],label='defmod',linewidth=1.0,color='b')
        else:
            plt.plot(range(len(dat))*dt_dyn+dt,sign_df[j]*dat[:,j],label='defmod',linewidth=1.0,color='b')
        if fd: plt.plot(range(dat_fd.shape[1])*dt_obs+dt,sign_df[j]*dat_fd[j,:]*scl,label='swpc',linewidth=1.0,color='m')
        if j==0:
            plt.plot(dat_es[:,0],dat_es[:,1+j],label=label[0],linewidth=1.0,color='g')
            plt.plot(dat_fm[:,0],dat_fm[:,1+j],label=label[1],linewidth=1.0,color='r')
            plt.xlabel('t [s]')
            plt.ylabel('v [m/s]')
            plt.legend(loc=1,fontsize=10.)
        else:
            plt.plot(dat_es[:,0],dat_es[:,1+j],linewidth=1.0,color='g')
            plt.plot(dat_fm[:,0],dat_fm[:,1+j],linewidth=1.0,color='r')
        if j>0:
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
        plt.ylim(ymin,ymax)
    plt.tight_layout()
    plt.savefig(name_sol+'_wave_%d'%(i+1)+'.png')
    plt.clf()

if rup==1 or dsp==1: # Plot rupture front
    if h5:
        fcoord   = dat_h5['crd_flt']
        dt_slip  = dat_h5['dt_slip']
        dat_slip = dat_h5['slip_dyn']['step 1']
    else:
        fcoord   = np.squeeze(io_mat.loadmat(matfile)['crd_flt' ])
        dt_slip  = np.squeeze(io_mat.loadmat(matfile)['dt_slip' ])
        dat_slip = np.squeeze(io_mat.loadmat(matfile)['dat_slip'])
    if h5:
        if dsp==1: dat = np.linalg.norm(dat_slip[:,:2,:],axis=1)
        cnt = 1E9*np.ones(shape=dat_slip[:,0,0].shape)
        len_tot=dat.shape[1]
    else:
        if dsp==1: dat = np.linalg.norm(dat_slip.item()[:,:,:2],axis=2)
        cnt = 1E9*np.ones(shape=dat_slip.item()[0,:,0].shape)
        len_tot=len(dat)
    for i in range(len_tot):
        t = i*dt_slip+dt
        if dsp==1:
            if h5:
                slip = dat[:,i]
            else:
                slip = dat[i,:]
        else:
            # Time integral for velocity output
            if h5:
                dat_int = np.sum(dat_slip[:,:2,:i+1],axis=2)*dt_slip
            else:
                dat_int = np.sum(dat_slip.item()[:i+1,:,:2],axis=0)*dt_slip
            slip = np.linalg.norm(dat_int,axis=1)
        idrup = np.squeeze(np.where(slip>tol))
        idfrn = np.squeeze(np.where(cnt<1E9))
        set_rup = np.setdiff1d(idrup,idfrn)
        cnt[set_rup] = t
        if i>int(len_tot*0.15) and set_rup.size==0:
            t_rup=t
            break
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
    dat_grid = griddata((np.squeeze(fcoord[:,1]), np.squeeze(fcoord[:,2])), cnt, (yi, zi), method='linear')
    dat_grid_ES = griddata((dat_ES[:,0], dat_ES[:,1]), dat_ES[:,2], (yi*1E3, -dip*1E3), method='linear')
    dat_grid_FM = griddata((dat_FM[:,0], dat_FM[:,1]), dat_FM[:,2], (yi*1E3, -dip*1E3), method='linear')
    fig, ax = plt.subplots()
    #t_rup=2.
    dt_rup=1.
    levels = np.arange(0.,t_rup+dt_rup,dt_rup)
    cs1=plt.contour(yi, dip, dat_grid,    colors='b', linestyles='--',linewidths=1.,levels=levels)
    cs2=plt.contour(yi, dip, dat_grid_ES, colors='g', linestyles='--',linewidths=1.,levels=levels)
    cs3=plt.contour(yi, dip, dat_grid_FM, colors='r', linestyles='--',linewidths=1.,levels=levels)
    ax.clabel(cs1, inline=1, fmt = '%1.1f',fontsize=10)
    ax.clabel(cs2, inline=1, fmt = '%1.1f',fontsize=10)
    ax.clabel(cs3, inline=1, fmt = '%1.1f',fontsize=10)
    lines = [cs1.collections[0], cs2.collections[0], cs3.collections[0]]
    labels = ['defmod',label[0],label[1]]
    plt.legend(lines, labels, loc=4, fontsize=10.)
    plt.xlabel('strike [km]')
    plt.ylabel('dip [km]')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.title('SCEC%d, dt = %0.1f [s]'%(idcase,dt_rup))
    plt.tight_layout()
    plt.savefig(name_sol+'_rup'+'.png')
    #plt.show()
