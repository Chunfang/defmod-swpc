#!/usr/bin/env python
import numpy as np
import os, sys, netCDF4
from subprocess import call
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patches as patches 
from matplotlib import gridspec
font = {'weight' : 'normal',
        'size'   : 12}
matplotlib.rc('font', **font)
#fin='out/swpc.fs.v.nc'
#nc=netCDF4.Dataset(fin)
#x=nc.variables['x'][:]
#y=nc.variables['y'][:]
#vx=nc.variables['Vx'][:]
#vy=nc.variables['Vy'][:]
#vz=nc.variables['Vz'][:]
#fin='out/swpc.fs.u.nc'
#nc=netCDF4.Dataset(fin)
#ux=nc.variables['Ux'][:]
#uy=nc.variables['Uy'][:]
#uz=nc.variables['Uz'][:]
#ax=plt.subplot(gs[0])
#plt.contourf(x,y,-vz[it,:,:],20, cmap=plt.cm.rainbow, vmax=1E-2, vmin=-1E-2)
#ax=plt.subplot(gs[1])
#plt.contourf(x,y,-uz[it,:,:],20, cmap=plt.cm.rainbow, vmax=1E-3, vmin=-1E-3)
path='HF3D_snp'
fin=path+'/swpc.xy.v.nc'
nc=netCDF4.Dataset(fin)
x=nc.variables['x'][:]
y=nc.variables['y'][:]
t=nc.variables['t'][:]
vx_xy=nc.variables['Vx'][:]
vy_xy=nc.variables['Vy'][:]
vz_xy=nc.variables['Vz'][:]
fin=path+'/swpc.xy.u.nc'
nc=netCDF4.Dataset(fin)
ux_xy=nc.variables['Ux'][:]
uy_xy=nc.variables['Uy'][:]
uz_xy=nc.variables['Uz'][:]
fin=path+'/swpc.xz.v.nc'
nc=netCDF4.Dataset(fin)
z=nc.variables['z'][:]
vx_xz=nc.variables['Vx'][:]
vy_xz=nc.variables['Vy'][:]
vz_xz=nc.variables['Vz'][:]
fin=path+'/swpc.xz.u.nc'
nc=netCDF4.Dataset(fin)
ux_xz=nc.variables['Ux'][:]
uy_xz=nc.variables['Uy'][:]
uz_xz=nc.variables['Uz'][:]

# Aspect ratios
rangex=x.max()-x.min()
rangey=y.max()-y.min()
rangez=z.max()-z.min()
gs = gridspec.GridSpec(2, 2, width_ratios=[1,1],height_ratios=[rangey,rangez])
w=rangex+rangez
h=2*rangex
width=6.4
height=width*h/w
fig=plt.figure()
vmax=1.; vmin=1E-5
umax=2.; umin=1E-5
i_scale=5
nplot=vx_xy.shape[0];t_play=15
for it in range(nplot):
    mappablev = cm.ScalarMappable(cmap=plt.get_cmap('gray'), norm=matplotlib.colors.Normalize(vmin=0, vmax=vmax))
    mappableu = cm.ScalarMappable(cmap=plt.get_cmap('gray'), norm=matplotlib.colors.Normalize(vmin=0, vmax=umax))
    mappablev.set_array(np.sqrt(vx_xy[i_scale,:,:]**2+vy_xy[i_scale,:,:]**2+vz_xy[i_scale,:,:]**2)-vmin)
    mappableu.set_array(np.sqrt(ux_xy[i_scale,:,:]**2+uy_xy[i_scale,:,:]**2+uz_xy[i_scale,:,:]**2)-vmin)
    fig.set_size_inches(height,width, forward=True)
    ax=plt.subplot(gs[0])
    plt.xlabel('x [km]')
    plt.ylabel('y [km]')
    plt.contourf(x,y,np.sqrt(vx_xy[it,:,:]**2+vy_xy[it,:,:]**2+vz_xy[it,:,:]**2)-vmin,40, cmap=plt.get_cmap('gray'), vmax=vmax, vmin=0.)
    ax.xaxis.set_ticks_position('none')
    ax.add_patch(patches.Rectangle((-24., -32.), 48., 64., fill=False,edgecolor="white",linestyle='dotted'))
    if True:
        cbax1 = fig.add_axes([0.08, 1.0, 0.4, 0.02])
        plt.colorbar(mappablev,orientation='horizontal',cax=cbax1,format='%1.2g')
    else:
        plt.axis('off') 
    #ax.add_patch(patches.Rectangle((-1., -.6), 2.0, 1.2, fill=False,edgecolor="white",linestyle='dotted'))
    plt.title(r'$||\mathbf{v}||$ [m/s], t=%1.2f [s]'%(t[it]))
    ax=plt.subplot(gs[1])
    #plt.xlabel('x [km]')
    #plt.ylabel('y [km]')
    plt.contourf(x,y,np.sqrt(ux_xy[it,:,:]**2+uy_xy[it,:,:]**2+uz_xy[it,:,:]**2)-umin,40, cmap=plt.get_cmap('gray'), vmax=umax, vmin=0.)
    if True:
        cbax2 = fig.add_axes([.56, 1.0, 0.4, 0.02])
        plt.colorbar(mappableu,orientation='horizontal',cax=cbax2,format='%1.2g')
    else:
        plt.axis('off')
    ax.xaxis.set_ticks_position('none')
    plt.title(r'$||\mathbf{u}||$ [m], t=%1.2f [s]'%(t[it]))
    ax=plt.subplot(gs[2])
    plt.xlabel('x [km]')
    plt.ylabel('z [km]')
    if False:
        plt.axis('off')
    plt.contourf(x,-z,np.sqrt(vx_xz[it,:,:]**2+vy_xz[it,:,:]**2+vz_xz[it,:,:]**2)-vmin,40, cmap=plt.get_cmap('gray'), vmax=vmax, vmin=0.)
    ax.xaxis.set_ticks_position('none')
    ax.add_patch(patches.Rectangle((-24., -32.), 48., 33., fill=False,edgecolor="white",linestyle='dotted'))
    #plt.colorbar(mappablev,orientation='vertical')
    #ax.add_patch(patches.Rectangle((-1., -3.7), 2.0, .7, fill=False,edgecolor="white",linestyle='dotted'))
    ax=plt.subplot(gs[3])
    #plt.xlabel('x [km]')
    #plt.ylabel('z [km]')
    if False:
        plt.axis('off')
    plt.contourf(x,-z,np.sqrt(ux_xz[it,:,:]**2+uy_xz[it,:,:]**2+uz_xz[it,:,:]**2)-umin,40, cmap=plt.get_cmap('gray'), vmax=umax, vmin=0.)
    ax.xaxis.set_ticks_position('none')
    #plt.colorbar(mappableu,orientation='vertical')
    if it==0: plt.tight_layout()
    fig.savefig(path+'/wav_%d.png'%(it),bbox_inches='tight')
    plt.clf()
os.chdir(path)
call(['avconv','-r',str(nplot/t_play),'-i','wav_%d.png','-vf', 'scale=trunc(iw/2)*2:trunc(ih/2)*2','fd_wav.mov'])
call(['mv', 'fd_wav.mov', '../HF3D_fd_wav.mov'])
call(['mv', 'wav_0.png', '../'+'HF3D_fd_wav.png'])
os.chdir('../')
