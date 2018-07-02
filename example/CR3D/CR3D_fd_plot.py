#!/usr/bin/env python
import numpy as np
import h5py
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
wave = "wave" in sys.argv
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
path='snp_fd'
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

#import scipy.io as io_mat
#matfile = 'CR3D_fd.mat'
#dat_fd  = np.squeeze(io_mat.loadmat(matfile)['dat_obs'])
#crd_fd  = np.squeeze(io_mat.loadmat(matfile)['crd_obs'])
#dt_fd   = np.squeeze(io_mat.loadmat(matfile)['dt_obs' ])

dic_fd = load_dict_from_hdf5('CR3D_fd.h5')
dat_fd = dic_fd['dat_obs_fd'] 
crd_fd = dic_fd['crd_obs_fd'] 
dt_fd  = dic_fd['dt_obs_fd' ] 

stid=range(27,36,1)
ymin = dat_fd[stid,2,:].min()*1.02
ymax = dat_fd[stid,2,:].max()*1.02

# Aspect ratios
rangex=x.max()-x.min()
rangey=y.max()-y.min()
rangez=z.max()-z.min()
gs = gridspec.GridSpec(2, 2, width_ratios=[1,1],height_ratios=[rangey,rangez])
outer= gridspec.GridSpec(1, 2, width_ratios=[1,1])
inner = gridspec.GridSpecFromSubplotSpec(2, 1,subplot_spec=outer[0],height_ratios=[rangey,rangez])
w=rangex+rangez
h=2*rangex
width=6.4
height=width*h/w
fig=plt.figure()
vmax=2E-2; vmin=1E-6
umax=2E-2; umin=1E-6
i_scale=5
t_plot=6.4
nplot=min(vx_xy.shape[0],t_plot/dt_fd);t_play=15
for it in range(nplot):
    mappablev = cm.ScalarMappable(cmap=plt.get_cmap('gray'), norm=matplotlib.colors.Normalize(vmin=0, vmax=vmax))
    mappableu = cm.ScalarMappable(cmap=plt.get_cmap('gray'), norm=matplotlib.colors.Normalize(vmin=0, vmax=umax))
    mappablev.set_array(np.sqrt(vx_xy[i_scale,:,:]**2+vy_xy[i_scale,:,:]**2+vz_xy[i_scale,:,:]**2)-vmin)
    mappableu.set_array(np.sqrt(ux_xy[i_scale,:,:]**2+uy_xy[i_scale,:,:]**2+uz_xy[i_scale,:,:]**2)-vmin)
    fig.set_size_inches(height,width, forward=True)

    if wave:
        ax=plt.Subplot(fig,inner[0])
    else:
        ax=plt.subplot(gs[0])
    ax.set_xlabel('x [km]')
    ax.set_ylabel('y [km]')
    ax.contourf(x,y,np.sqrt(vx_xy[it,:,:]**2+vy_xy[it,:,:]**2+vz_xy[it,:,:]**2)-vmin,40, cmap=plt.get_cmap('gray'), vmax=vmax, vmin=0.)
    ax.xaxis.set_ticks_position('none')
    ax.add_patch(patches.Rectangle((-1., -1.), 2., 2., fill=False,edgecolor="white",linestyle='dotted'))
    j=0
    for i in stid:
        j+=1
        if wave:
            ax.plot(crd_fd[i,0],crd_fd[i,1],marker='v',color='red')
            ax.text(crd_fd[i,0]-1.2E-1,crd_fd[i,1]+1.6E-1,str(j),color='white')
    if True:
        cbax1 = fig.add_axes([0.08, 1.0, 0.4, 0.02])
        plt.colorbar(mappablev,orientation='horizontal',cax=cbax1,format='%1.2g',ticks=[0,vmax/2,vmax])
    else:
        plt.axis('off') 
    #ax.add_patch(patches.Rectangle((-1., -.6), 2.0, 1.2, fill=False,edgecolor="white",linestyle='dotted'))
    plt.title(r'$||\mathbf{v}||$ [m/s], t=%1.2f [s]'%(t[it]))
    fig.add_subplot(ax)

    if wave:
        ax=plt.Subplot(fig,outer[1])
        j=0
        for i in stid:
            if j==0:
                ax_st = fig.add_axes([0.56, .96-j*1./len(stid), 0.44, 1./len(stid)], xticklabels=[],xlim=(0,10), ylim=(ymin,ymax))
                ax_st.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
                ax_st.set_title(r'$v_z$ [m/s]')
            elif j<len(stid)-1:
                ax_st = fig.add_axes([0.56, .96-j*1./len(stid), 0.44, 1./len(stid)], xticklabels=[],yticklabels=[], xlim=(0,10), ylim=(ymin,ymax))
            else:
                ax_st = fig.add_axes([0.56, .96-j*1./len(stid), 0.44, 1./len(stid)], yticklabels=[], xlim=(0,10), ylim=(ymin,ymax))
                ax_st.set_xlabel('t [s]')
            ax_st.plot(np.array(range(dat_fd.shape[-1]))*dt_fd, dat_fd[i,2,:],linewidth=1.)
            j+=1
            ax_st.text(9,ymin+(ymax-ymin)/10.,str(j))
            ax_st.plot([t[it],t[it]],[ymin,ymax],linewidth=1.)
    else:
        ax=plt.subplot(gs[1])
        plt.xlabel('x [km]')
        plt.ylabel('y [km]')
        plt.contourf(x,y,np.sqrt(ux_xy[it,:,:]**2+uy_xy[it,:,:]**2+uz_xy[it,:,:]**2)-umin,40, cmap=plt.get_cmap('gray'), vmax=umax, vmin=0.)
        if True:
            cbax2 = fig.add_axes([.56, 1.0, 0.4, 0.02])
            plt.colorbar(mappableu,orientation='horizontal',cax=cbax2,format='%1.2g',ticks=[0,umax/2,umax])
        else:
            plt.axis('off')
        ax.xaxis.set_ticks_position('none')
        plt.title(r'$||\mathbf{u}||$ [m], t=%1.2f [s]'%(t[it]))

    if wave:
        ax=plt.Subplot(fig,inner[1])
    else:
        ax=plt.subplot(gs[2])
    ax.set_xlabel('x [km]')
    ax.set_ylabel('z [km]')
    if False:
        plt.axis('off')
    ax.contourf(x,-z,np.sqrt(vx_xz[it,:,:]**2+vy_xz[it,:,:]**2+vz_xz[it,:,:]**2)-vmin,40, cmap=plt.get_cmap('gray'), vmax=vmax, vmin=0.)
    ax.xaxis.set_ticks_position('none')
    ax.add_patch(patches.Rectangle((-1., -2.5), 2., 2., fill=False,edgecolor="white",linestyle='dotted'))
    #plt.colorbar(mappablev,orientation='vertical')
    #ax.add_patch(patches.Rectangle((-1., -3.7), 2.0, .7, fill=False,edgecolor="white",linestyle='dotted'))
    if wave:
        for i in stid:
            ax.plot(crd_fd[i,0],crd_fd[i,2]-1E-1,marker='v',color='red')
        fig.add_subplot(ax)
    else:
        ax=plt.subplot(gs[3])
        plt.xlabel('x [km]')
        plt.ylabel('z [km]')
        if False:
            plt.axis('off')
        plt.contourf(x,-z,np.sqrt(ux_xz[it,:,:]**2+uy_xz[it,:,:]**2+uz_xz[it,:,:]**2)-umin,40, cmap=plt.get_cmap('gray'), vmax=umax, vmin=0.)
        ax.xaxis.set_ticks_position('none')
        #plt.colorbar(mappableu,orientation='vertical')

    if it==0: plt.tight_layout()
    fig.savefig(path+'/wav_%d.png'%(it),bbox_inches='tight')
    plt.clf()
os.chdir(path)
call(['avconv','-r',str(nplot/t_play),'-i','wav_%d.png','-vf', 'scale=trunc(iw/2)*2:trunc(ih/2)*2', '-pix_fmt', 'yuv420p', 'fd_wav.mov'])
call(['mv', 'fd_wav.mov', '../CR3D_fd.mov'])
call(['cp', 'wav_0.png', '../'+'CR3D_fd.png'])
os.chdir('../')
