#!/usr/bin/env python 
import numpy as np
import h5py
from scipy import interpolate
import matplotlib
import matplotlib.pyplot as plt
font = {'weight' : 'normal',
        'size'   : 12}
matplotlib.rc('font', **font)
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
nplot=250; eps=.0025
dat_fe=load_dict_from_hdf5('CR3D_fe.h5')
slip_sta=dat_fe['dat_slip_sta']
trac_sta=dat_fe['dat_fqs']
crd_flt=dat_fe['crd_flt']
dt=dat_fe['dt']; dat_log=dat_fe['dat_log']
ycnt=(crd_flt[:,1].max()+crd_flt[:,1].min())/2.
id_dip=np.where(np.logical_and(crd_flt[:,1]>ycnt-eps,crd_flt[:,1]<ycnt+eps))
x_dip,y_dip,z_dip=np.squeeze(crd_flt[id_dip,0]),np.squeeze(crd_flt[id_dip,1]),np.squeeze(crd_flt[id_dip,2])
min_trac=trac_sta[:,1,:].min(); max_trac=trac_sta[:,1,:].max()
min_sigm=trac_sta[:,2,:].min(); max_sigm=trac_sta[:,2,:].max()
min_pres=trac_sta[:,3,:].min(); max_pres=trac_sta[:,3,:].max()
range_st=max([max_trac-min_trac,max_sigm-min_sigm,max_pres-min_pres]) 

k=0
for i in [1,10,20,26,28]: #range(np.shape(slip_sta)[-1])[17:31:4]:
    dat_slip=(slip_sta[id_dip,1,i]).squeeze() #-slip_sta[id_dip,1,0]
    dat_trac=(trac_sta[id_dip,1,i]).squeeze() #-trac_sta[id_dip,1,0]
    dat_sigm=(trac_sta[id_dip,2,i]).squeeze() #-trac_sta[id_dip,2,0]
    dat_pres=(trac_sta[id_dip,3,i]).squeeze() #-trac_sta[id_dip,3,0]
    z_plot=np.linspace(crd_flt[:,2].min(),crd_flt[:,2].max(),num=nplot)
    slip_plot=interpolate.interp1d(z_dip,dat_slip,fill_value='extrapolate')(z_plot)
    trac_plot=interpolate.interp1d(z_dip,dat_trac,fill_value='extrapolate')(z_plot)
    sigm_plot=interpolate.interp1d(z_dip,dat_sigm,fill_value='extrapolate')(z_plot)
    pres_plot=interpolate.interp1d(z_dip,dat_pres,fill_value='extrapolate')(z_plot)
    #plt.plot(-z_plot,trac_plot)
    #plt.plot(-z_plot,sigm_plot)
    #plt.plot(-z_plot,trac_plot/(pres_plot-sigm_plot))
    #plt.plot(-z_dip,dat_slip,'.')
    ax=plt.subplot(1,5,1)
    if k==0:
        plt.plot(trac_plot,z_plot,label='day '+str(i))
    else:
        plt.plot(trac_plot,z_plot,label=str(i))
    plt.legend(loc='upper right')
    plt.xlim(min_trac,min_trac+range_st)
    plt.ylim(z_plot.min(),z_plot.max())
    if k==0:
        ax.set_ylabel(r'$z$ [km]')
        ax.set_xlabel(r'$\tau$ [Pa]')
        ax.set_yticks([-1.7, -1.55, -1.45, -1.3])
    ax=plt.subplot(1,5,2)
    plt.plot(-sigm_plot,z_plot)
    plt.xlim(-max_sigm,-max_sigm+range_st)
    plt.ylim(z_plot.min(),z_plot.max())
    ax.set_yticks([]) 
    if k==0: ax.set_xlabel(r'$\sigma_n$ [Pa]')
    ax=plt.subplot(1,5,3)
    plt.plot(pres_plot,z_plot)
    plt.xlim(min_pres,min_pres+range_st)
    plt.ylim(z_plot.min(),z_plot.max())
    ax.set_yticks([])
    if k==0: ax.set_xlabel(r'$p$ [Pa]')
    ax=plt.subplot(1,5,4)
    plt.plot(-trac_plot/sigm_plot,z_plot)
    plt.ylim(z_plot.min(),z_plot.max())
    ax.set_yticks([])
    if k==0:
        ax.set_xlabel(r'$\tau/\sigma_n$')
        ax.set_xticks([0.2, 0.6])
    ax=plt.subplot(1,5,5)
    plt.plot(slip_plot,z_plot)
    plt.ylim(z_plot.min(),z_plot.max())
    if k==0: 
       ax.set_xlabel('slip [m]')
       ax.set_xticks([0.01, 0.06])
    ax.set_yticks([])
    k+=1
#plt.tight_layout()
plt.show()
