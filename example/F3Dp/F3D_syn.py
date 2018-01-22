#!/usr/bin/env python 
import numpy as np
import os,sys
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt 
import argparse

ap=argparse.ArgumentParser()
ap.add_argument('-vis')    # 1 plot cropped point cloud
ap.add_argument('-refine') # 1 refine mesh
ap.add_argument('-clean')  # 1 remove tmp files

if ap.parse_args().vis==None:
    vis=0
else:
    vis=int(ap.parse_args().vis)
if ap.parse_args().refine==None:
    refine=0
else:
    refine=int(ap.parse_args().refine)
if ap.parse_args().clean==None:
    clean=0
else:
    clean=int(ap.parse_args().clean)

# Synthetic fault pixels
z=np.linspace(.2, -.8, num=100)
y=np.linspace(-.625,.625, num=120)
grid=np.meshgrid(y,z)
x=np.zeros((len(z)*len(y),1),dtype=np.float)
dat_vert=np.hstack((x,grid[0].reshape(x.shape),grid[1].reshape(x.shape)))

# weak
wl=np.linspace(.12,.18,num=8); amp=.03125*np.sqrt(wl)
e=1.025; r=-.2 
dip=70.; zcnt=-.35
omg=[ 0.82976173,  0.89624834,  0.03829284, -0.50016345, -1.06606012,  1.40505898, -1.24256034,  1.28623393]
#omg=(np.random.rand(wl.shape[0])-.5)*np.pi
L=dat_vert[1,:].max()-dat_vert[1,:].min()
zmax=z.max(); zmin=z.min()

for i in range(len(wl)):
    phs=dat_vert[:,1]/wl[i]*np.pi+omg[i]
    dat_vert[:,0]=dat_vert[:,0]+amp[i]*np.cos(phs)*(e*zmax-dat_vert[:,2])/(e*zmax-zmin)*np.exp(r*abs(phs)/np.pi)
dat_vert[:,0]=dat_vert[:,0]+(zcnt-dat_vert[:,2])*np.tan((90.-dip)/180.*np.pi)

# ridge patch
def flt_patch(dat_vert,slope1,slope2,trunc1,trunc2,hlw,hup):
    b1=-slope1*trunc1-.7
    b2=-slope2*trunc2-.7
    in_id=np.where(np.logical_and(dat_vert[:,2]-slope1*dat_vert[:,1]<b1, dat_vert[:,2]-slope2*dat_vert[:,1]<b2))[0]
    out_id=np.setdiff1d(np.array(range(len(dat_vert)),dtype=np.int32),in_id) 
    x_shift=dat_vert[in_id,0]
    # ridge patch
    k=0
    zup=dat_vert[:,2].max()
    zlw=dat_vert[:,2].min()
    for i in in_id:
        r=abs(dat_vert[i,1]-.5*(trunc1+trunc2))
        R=.5*((dat_vert[i,2]-b2)/slope2-(dat_vert[i,2]-b1)/slope1)
        h=hlw+(dat_vert[i,2]-zlw)/(zup-zlw)*(hup-hlw)
        x_shift[k]=x_shift[k]+np.cos(r/R*np.pi/2.)*h
        k+=1
    dat_vert=np.vstack((dat_vert[out_id,:],
        np.hstack((x_shift.reshape(len(in_id),1),
        dat_vert[in_id,1].reshape(len(in_id),1),
        dat_vert[in_id,2].reshape(len(in_id),1)))))
    return dat_vert
slope1=10.;slope2=-10.
trunc1=.1;trunc2=.6
hup=0.;hlw=.08
#dat_vert=flt_patch(dat_vert,slope1,slope2,trunc1,trunc2,hlw,hup)
print omg

fout='F3D_syn.xyz'
f=open(fout,'w+')
np.savetxt(f,dat_vert,delimiter=' ', fmt='%.6f '*3)
f.close()
from subprocess import call
fin=fout
fout=fout.rsplit('.')[0]+'.stl'
mxl='xyz2stl.mlx' 
call(['meshlabserver', '-i',fin,'-o',fout,'-s',mxl])
if clean==1: os.remove(fin)

# Mesh
fin=fout
if refine==1:
    fout=fout.rsplit('.')[0]+'_dns.exo'
else:
    fout=fout.rsplit('.')[0]+'.exo'
jou='F3D_tet.jou'
txt_jou=open(jou,'r')
txt_jou_tmp=open('tmp.jou','w+')
hf=0.0025 # fault  grid length (0.0025 for ~100 m tet model, 0.003 for ~40 m)
hm=0.0075 # matrix grid length (0.0075 for ~100 m tet model, 0.010 for ~40 m)
for line in txt_jou:
    line=line.strip('\r\n')
    if 'import' in line.lower():
        line='import stl "'+fin+'"'
    if 'export' in line.lower():
        line='export mesh "'+fout+'" dimension 3 overwrite'
    if 'surface 46 94 95 97 size' in line.lower():
        line='surface 46 94 95 97 size %0.6f' %(2*hf)
    if 'volume all size' in line.lower():
        line='volume all size %0.6f' %(2*hm) 
    txt_jou_tmp.write(line+'\n')
    if 'mesh volume all' in line.lower() and refine==1:
        txt_jou_tmp.write('refine volume all\n')
txt_jou.close();txt_jou_tmp.close()
call(['trelis','-nojournal','-nographics','tmp.jou'])
if clean==1: os.remove('tmp.jou')

# Preprocessing msh=>inp
dt_dyn=2E-5 #1E-5 for dns 100 m tet model, 8E-5 for 40 m tet, 8E-4 for ~1 m tet 
import F3D_msh2inp 
_=F3D_msh2inp.msh2inp(fout,dt_dyn)

# Fault plot
if vis==1:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(dat_vert[:,0], dat_vert[:,1], dat_vert[:,2], c='b', marker='.')
    # Create cubic bounding box to simulate equal aspect ratio
    max_range = np.array([np.max(dat_vert[:,0])-np.min(dat_vert[:,0]),np.max(dat_vert[:,1])\
            -np.min(dat_vert[:,1]), np.max(dat_vert[:,2])-np.min(dat_vert[:,2])]).max()
    Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten()
    Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten()
    Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten()
    for xb, yb, zb in zip(Xb, Yb, Zb):
       ax.plot([xb], [yb], [zb], 'w',)
    plt.title('fault [km]')
    plt.grid()
    plt.show()
