#!/usr/bin/env python
import numpy as np
import sys
import scipy.io as io_mat 
from subprocess import call
import os
import matplotlib
matplotlib.use('Svg')
import matplotlib.pyplot as plt
font = {'weight' : 'normal',
        'size'   : 12}
matplotlib.rc('font', **font)

name_sol = sys.argv[1]
matfile = [name_sol+     '_2.mat',\
           name_sol+     '_4.mat',\
           name_sol+     '_8.mat',\
           name_sol+    '_16.mat',\
           name_sol+ '_2_dsp.mat',\
           name_sol+ '_4_dsp.mat',\
           name_sol+ '_8_dsp.mat',\
           name_sol+'_16_dsp.mat']

dat_slip=[]; fcoord=[]; ocoord=[]; dat_seis=[]; dt_dyn=[]; dt_slip=[]
dat_log=[]; dat_fqs=[]; dt=[]; dt_rsf=[]; dat_log_rsf=[]
for f in matfile:
    log_tmp=np.squeeze(io_mat.loadmat(f)['dat_log'])
    fcoord.append(np.squeeze(io_mat.loadmat(f)['crd_flt']))
    dat_slip.append(np.squeeze(io_mat.loadmat(f)['dat_slip']))
    dt_slip.append(np.squeeze(io_mat.loadmat(f)['dt_slip']))
    ocoord.append(np.squeeze(io_mat.loadmat(f)['crd_obs']))
    dat_seis.append(np.squeeze(io_mat.loadmat(f)['dat_seis']))
    dt_dyn.append(np.squeeze(io_mat.loadmat(f)['dt_dyn']))
    dat_fqs.append(np.squeeze(io_mat.loadmat(f)['dat_fqs']))
    dt.append(np.squeeze(io_mat.loadmat(f)['dt']))
    dt_rsf.append(np.squeeze(io_mat.loadmat(f)['dt_rsf']))
    dat_log_rsf.append(np.squeeze(io_mat.loadmat(f)['dat_log_rsf']))
    dat_log.append(log_tmp)


# Calculate event magnitude, location
mag=[]; xloc=[]; tevt=[]
vtol=1E-3
for k in [4,5,6,7]: # Use absolute output
    flt=fcoord[k][:,1]
    idx=np.argsort(flt)
    xflt=fcoord[k][idx,:]
    if len(dat_log[k].shape)==1:
        event=np.array([dat_log[k][:]]) 
    else:
        event=dat_log[k]
    tmpmag=[]; tmploc=[]; tmpt=[]
    for j,i in zip(event[:,0],range(len(event))): 
        if event[i,1]==1:
            # Magnitude  
            sdrp=dat_fqs[k][j+1,idx,0] - dat_fqs[k][j,idx,0]
            try:
                slip=dat_slip[k][i][1500,idx,0]
            except:
                slip=dat_slip[k].item()[1500,idx,0]
            idslp=np.where(sdrp<0.)
            tmpmag.append(np.log10(-sum(slip[idslp]*sdrp[idslp])))
            # Location 
            try: 
                v=dat_slip[k-4][i][5,idx,0]
            except:
                v=dat_slip[k-4].item()[5,idx,0]
            idtmp=np.where(v>=vtol)
            z=xflt[sum(idtmp[0])/len(idtmp[0]),:]
            x=xflt[0,0]+(xflt[-1,0]-xflt[0,0])*(z-xflt[0,1])/(xflt[-1,1]-xflt[0,1])
            tmploc.append(np.array([x,z]))
            # Event time
            tmpt.append((dat_log_rsf[k][event[i,0]-2]*dt_rsf[k]+dt[k])/3600)
    mag.append(tmpmag)
    xloc.append(tmploc)
    tevt.append(tmpt)

plt.figure()
color = ['b','c','m','k']
plt.plot(xflt[:,0],xflt[:,1])
for k in [0,1,2,3]:
    ax=plt.subplot(2,2, k+1)
    plt.plot(xflt[[0,-1],0],xflt[[0,-1],1])
    scatplt=[]; scatlab=[]
    for i in range(len(mag[k])):
        scatplt.append(plt.scatter(xloc[k][i][0],xloc[k][i][1],marker='o',s=mag[k][i]*8,c=color[i]))    
        if k==0 and i==0:
            lab = 't= %0.1f' %(tevt[k][i])+' hr'
            plt.ylabel('z [km]')
            tit = 'wavelength = %d' %(2**(k+1))+' days'
        else:
            if k==2: plt.xlabel('x [km]')
            lab = '%0.1f' %(tevt[k][i])
            tit = '%d' %(2**(k+1))
            if k==0: tit = 'wavelength = %d' %(2**(k+1))+' days'
        scatlab.append(lab)
        ax.set_xlim([-.05,.05])
        ax.set_ylim([-2.08,-2.0])
        plt.gca().set_aspect('equal', adjustable='box')
        plt.legend(scatplt,
               scatlab,
               scatterpoints=1,
               loc='upper right',
               ncol=1,
               fontsize=12)
        plt.title(tit)
plt.savefig(name_sol+'_evt.png')

# Stress ratio plot
plt.figure()
for k in [0,1,2,3]:
    flt=fcoord[k][:,1]
    idx=np.argsort(flt)
    dat_plt=dat_fqs[k][:,idx,:]
    n_plt=len(dat_plt[:,0,0])
    id_plt=int(len(dat_plt[0,:,0])/2)
    flt=flt[idx]
    depth=flt[id_plt]
    yplt=-dat_plt[:,id_plt,0]/dat_plt[:,id_plt,1]  
    xplt=dat_log_rsf[k]*dt_rsf[k]/3600.
    if k==0:
        lab='wavelength '+'%d' %(2**(k+1))+' days'
    else:
        lab='%d' %(2**(k+1))
    plt.plot(xplt,yplt[1:-1],label=lab)
    plt.legend(loc=3,prop={'size':11})
plt.title('depth = '+'%0.1f' %(depth) +' [km]')
ax.set_xlim([0,700])
plt.xlabel('time [hr]')
plt.ylabel(r'$\tau/\sigma_n$')
plt.savefig(name_sol+'_mu.png')
