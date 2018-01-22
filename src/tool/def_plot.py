#!/usr/bin/env python
import numpy as np
import sys
import scipy.io as io_mat 
from subprocess import call
import os
name_sol = sys.argv[1]
if 'slip' in sys.argv:
    slip=True
else:
    slip=False
if 'rsf' in sys.argv:
    rsf=True
else:
    rsf=False
matfile = name_sol+'.mat'
dat_qs      = np.squeeze(io_mat.loadmat(matfile)['dat_qs'     ])
dat_seis    = np.squeeze(io_mat.loadmat(matfile)['dat_seis'   ]) 
ocoord      = np.squeeze(io_mat.loadmat(matfile)['crd_obs'    ]) 
dt          = np.squeeze(io_mat.loadmat(matfile)['dt'         ]) 
dat_log     = np.squeeze(io_mat.loadmat(matfile)['dat_log'    ]) 
dt_dyn      = np.squeeze(io_mat.loadmat(matfile)['dt_dyn'     ]) 
dat_log_dyn = np.squeeze(io_mat.loadmat(matfile)['dat_log_dyn']) 
import matplotlib
matplotlib.use('Svg')
import matplotlib.pyplot as plt
font = {'weight' : 'normal',
        'size'   : 11}
matplotlib.rc('font', **font)
if slip: 
    try:
        fcoord       = np.squeeze(io_mat.loadmat(matfile)['crd_flt'     ]) 
        dt_slip      = np.squeeze(io_mat.loadmat(matfile)['dt_slip'     ])
        dat_log_slip = np.squeeze(io_mat.loadmat(matfile)['dat_log_slip'])
        dat_slip     = np.squeeze(io_mat.loadmat(matfile)['dat_slip'    ]) 
        dat_fqs      = np.squeeze(io_mat.loadmat(matfile)['dat_fqs'     ])
        plt.ion()
    except:
        print 'Slip data not sorted!'
        sys.exit(0)
if slip and rsf:
    dt_rsf      = np.squeeze(io_mat.loadmat(matfile)['dt_rsf'     ])
    dat_log_rsf = np.squeeze(io_mat.loadmat(matfile)['dat_log_rsf'])
    dat_rsf     = np.squeeze(io_mat.loadmat(matfile)['dat_rsf'    ])
# Plot seismic waveforms for the largest event at observations 1
nobs=len(ocoord)
try:
    n_event=len(dat_log)
    amp=[]
    for i in range(dat_seis.shape[0]):
        amp.append(max(np.linalg.norm(dat_seis[i][0],axis=1)))
    eid = np.argmax(amp)
except:
    n_event=1
    eid=0
# xyz component to plot
dmn = ocoord.shape[1]
oid=range(nobs);
if len(oid)>6: oid=oid[:6] 
ymax=np.empty((len(oid),),dtype=float)
ymin=np.empty((len(oid),),dtype=float)
for i in oid:
    if n_event>1:
        dat = dat_seis[eid][i]
    else:
        dat = dat_seis.item()[i]
    for j in range(dmn):
        plt.subplot(len(oid), dmn, j+1+i*dmn)
        plt.plot(range(len(dat))*dt_dyn,dat[:,j])
        if j>0:
            plt.gca().axes.get_xaxis().set_visible(False)
            plt.gca().axes.get_yaxis().set_visible(False)
        elif i>0:
            plt.gca().axes.get_xaxis().set_visible(False)
    ymax[i] = dat.max(); ymin[i] = dat.min()
for i in oid:
    for j in range(dmn):
        ax=plt.subplot((len(oid)), dmn, j+1+i*dmn)
        ylim=max(abs(ymin[i]),abs(ymax[i]))*1.04
        ax.set_ylim([-ylim,ylim])
        ax.set_xlim([0.,len(dat)*dt_dyn])
plt.savefig(name_sol+'_wave.png')
if slip:
    if dmn==3:
        call(["mkdir", "fig"])
        from scipy.interpolate import griddata
        yi = np.linspace(min(fcoord[:,1]), max(fcoord[:,1]),150)
        zi = np.linspace(min(fcoord[:,2]), max(fcoord[:,2]),75)
        yi,zi = np.meshgrid(yi,zi)
        theta=np.pi/6  # np.pi/6 (60 dg dip)

        # Plot 2D fault slip for given event
        if n_event>1:
            dat = dat_slip[eid] 
        else:
            dat= dat_slip.item()
        i_scale = np.shape(dat)[0]/20 
        dat_grid = griddata((np.squeeze(fcoord[:,1]), np.squeeze(fcoord[:,2])), np.linalg.norm(dat[i_scale,:,:dmn-1], axis=1), (yi, zi), method='nearest')
        vmax=dat_grid.max(); vmin=0.
        fig=plt.figure()
        plt.contourf(yi, zi/np.cos(theta), dat_grid, 20, cmap=plt.cm.rainbow, vmax=vmax, vmin=vmin)
        plt.colorbar(orientation='horizontal')
        for i in range(len(dat)):
            dat_grid = griddata((np.squeeze(fcoord[:,1]), np.squeeze(fcoord[:,2])), np.linalg.norm(dat[i,:,:dmn-1], axis=1), (yi, zi), method='linear')
            plt.contourf(yi, zi/np.cos(theta), dat_grid, 20, cmap=plt.cm.rainbow, vmax=vmax, vmin=vmin)
            plt.xlabel('length [km]')
            plt.ylabel('height [km]')
            plt.gca().set_aspect('equal', adjustable='box')
            plt.title(r't = %1.2f [s]' % (i*dt_slip))
            plt.savefig('fig/sr_%d.png' %(i))
        os.chdir('fig')
        call(['avconv', '-i', 'sr_%d.png', name_sol+'_sr.mov'])
        call(['mv', name_sol+'_sr.mov', '../'])
        call(['mv', 'sr_0.png', '../'+name_sol+'_sr.png'])
        os.chdir('../')

        # Plot QS stress
        plt.figure()
        dat_fsq_sort=[]
        vmax=0.6; vmin=0.
        for i in range(len(dat_fqs)):
            dat_grid = griddata((np.squeeze(fcoord[:,1]), np.squeeze(fcoord[:,2])), np.linalg.norm(dat_fqs[i][:,:2],axis=1)/(-dat_fqs[i][:,2]), (yi, zi), method='linear')
            plt.contourf(yi, zi/np.cos(theta), dat_grid, 20, cmap=plt.cm.rainbow, vmax=vmax, vmin=vmin)
            plt.xlabel('length [km]'); plt.ylabel('height [km]');
            plt.gca().set_aspect('equal', adjustable='box')
            plt.title(r'$\tau/\sigma_n$, day %i' % (i))
            if (i==0):
                #plt.tight_layout()
                plt.colorbar(orientation='horizontal')
            plt.savefig('fig/fqs_%d.png' %(i))
        os.chdir('fig')
        call(['avconv', '-i', 'fqs_%d.png', name_sol+'_fqs.mov'])
        call(['mv', name_sol+'_fqs.mov', '../'])
        call(['mv', 'fqs_0.png', '../'+name_sol+'_fqs.png'])
        os.chdir('../')
        #import shutil; shutil.rmtree('fig')
    
    else:
        tol_v=1E-3 # Minimum slip rate to be seismic
        flt_plot=fcoord[:,:]
        dat_flt=[]
        if n_event>1:
            for j in range(len(dat_slip)):
                dat = dat_slip[j]
                dat_flt.append(dat[:,:,:])
        else:
            dat = dat_slip.item()
            dat_flt.append(dat[:,:,:])
        if rsf: 
            dat_rsf_flt=[]
            for i in range(len(dat_log_rsf)):
                dat=dat_rsf[i]
                dat_rsf_flt.append(dat[:,:,:])
        id_seis=[]
        mslip = np.empty((0,1),dtype=np.float64)
        for j in range(len(dat_flt)):
            if dat_flt[j][:,:,0].shape[1]>0:
                idmx = np.argmax(np.max(abs(dat_flt[j][:,:,0]),axis=0))
                mslip = np.vstack((mslip,np.max(abs(dat_flt[j][:,:,0]),axis=0)[idmx]))
            else:
                idmx = 0 
                mslip = np.vstack((mslip,0.)) 
        mslip=np.squeeze(mslip)
        if len(dat_log.shape)>1:
            pick=(mslip>tol_v)*(dat_log[:,1]>0)
            mslip=mslip[pick]
            id_seis.append(pick)
        else:
            pick=(mslip>tol_v)*(dat_log[1]>0)
            id_seis.append(pick)
        
        # Plot decaying slip rate
        plt.figure()
        ax=plt.subplot()
        for j in [eid]: 
            if len(dat_log.shape)>1:
                seis=np.squeeze(id_seis)[j]
            else:
                seis=id_seis
            if seis:
                for k in range(dat_flt[j].shape[1]):
                    if max(abs(dat_flt[j][:,k,0]))>tol_v:
                        plt.semilogy(range(len(dat_flt[j][:,k,0]))*dt_slip+dt_slip,dat_flt[j][:,k,0])
                        #plt.plot(range(len(dat_flt[j][:,k,0]))*dt_slip,dat_flt[j][:,k,0])
        ax.set_ylim([1E-6,15.])
        ax.set_xlim([0.,10.])
        plt.xlabel('dynamic time [s]'); plt.ylabel('slip rate [m/s]')
        plt.savefig(name_sol+'_slip.png')

        plt.figure()
        ax=plt.subplot()
        for j in [eid]: #range(len(id_seis[i])): 
            if len(dat_log.shape)>1:
                seis=np.squeeze(id_seis)[j]
            else:
                seis=id_seis
            if seis:
                for k in range(dat_flt[j].shape[1]):
                    plt.plot(range(len(dat_flt[j][:,k,0]))*dt_slip,dat_flt[j][:,k,1])
        ax.set_xlim([0.,10.])
        plt.xlabel('dynamic time [s]'); plt.ylabel(r'$\mu$')
        plt.savefig(name_sol+'_mu.png')
        
        if rsf==1:
            plt.figure()
            ax=plt.subplot()
            for j in [eid]: #range(len(id_seis[i])): 
                if len(dat_log.shape)>1:
                    seis=np.squeeze(id_seis)[j]
                else:
                    seis=id_seis
                if seis:
                    for k in range(dat_flt[j].shape[1]):
                        plt.semilogy(range(len(dat_flt[j][:,k,0]))*dt_slip,dat_flt[j][:,k,-1])
            ax.set_xlim([0.,10.])
            plt.xlabel('dynamic time [s]'); plt.ylabel(r'$\theta$')
            plt.savefig(name_sol+'_theta.png')

        # Line plot the velocity profile along the fault
        plt.figure()
        ax=plt.subplot()
        for j in [eid]: #range(len(id_seis[i])): 
            if len(dat_log.shape)>1:
                seis=np.squeeze(id_seis)[j]
            else:
                seis=id_seis
            if seis:
                idx=np.argsort(flt_plot[:,1])
                for k in range(dat_flt[j].shape[0]): # dynamic time idx
                    #if max(dat_flt[j][k,:,0])>1E-2:
                    if k%100==0 and k<1000:
                        xplt=flt_plot[idx,1]
                        yplt=dat_flt[j][k,idx,0]
                        idtmp=np.where(yplt>=tol_v)
                        if len(idtmp[0])>0 and len(idtmp[0])<len(xplt) and k>0:
                            xplt=np.hstack((xplt[0:idtmp[0][0]],xplt[idtmp],xplt[idtmp[0][-1]+1:-1]))
                            yplt=np.hstack((yplt[0:idtmp[0][0]],yplt[idtmp],yplt[idtmp[0][-1]+1:-1]))
                        plt.plot(xplt,yplt)
                    ax.set_xlim([-2.3,-1.8])
                    #ax.set_ylim([0,.6])
        plt.xlabel('depth [km]'); plt.ylabel('slip [m]')
        plt.savefig(name_sol+'_disp.png')


        if rsf: # Plot RSF pseudo velocity
            dt_rsf=dt_rsf/3600.
            plt.figure()
            ax=plt.subplot()
            shift=dt/3600.
            for j in range(len(dat_log_rsf)): # truncated pseudo time steps
                for k in range(dat_rsf_flt[j].shape[1]): # fault node index 
                    if dat_rsf_flt[j][:,k,0].max()>=1E-12:
                        plt.semilogy(np.array(range(len(dat_rsf_flt[j][:,k,0])))*dt_rsf+shift,dat_rsf_flt[j][:,k,0])
                try:
                    shift=shift+(len(dat_rsf_flt[j][:,k,0])-1)*dt_rsf
                except: pass
            plt.xlabel('quasi-static time [hr]'); plt.ylabel('pseudo slip rate [m/s]')
            #ax.set_xlim([0,720])
            ax.set_ylim([1E-12,1.])
            plt.savefig(name_sol+'_rsf.png')

        # Plot stress stress drop 
        plt.figure()
        idx=np.argsort(flt_plot[:,1])
        xplt=flt_plot[idx,1]
        if  len(dat_log.shape)==1: dat_log=np.array([dat_log[:]])
        for j,i in zip(dat_log[:,0],range(len(dat_log))): #range(len(dat_fqs[:,0,0])): 
            if dat_log[i,1]==1:
                yplt=dat_fqs[j+1,idx,0] - dat_fqs[j,idx,0]
                plt.plot(xplt,yplt/1E6)
            plt.xlabel('depth [km]'); plt.ylabel('shear stress drop [MPa]');
            plt.savefig(name_sol+'_sd.png')
