#!/usr/bin/env python
import numpy as np
import os, sys, netCDF4
import argparse
ap=argparse.ArgumentParser()
ap.add_argument('-m') # mesh file
ap.add_argument('-p') # problem number
ap.add_argument('-d') # absolute output 1

# Exodus mesh file
fin= ap.parse_args().m
# SCEC problem index
idcase=int(ap.parse_args().p)
if not idcase in [14, 15]:
    print "-p arg has to be 14 or 15"
    sys.exit(0)
# output type
if ap.parse_args().d==None:
    dsp_hyb=0
else:
    dsp_hyb=min(1,int(ap.parse_args().d))

# headers
dim=3
t=.05; dt=.05; nviz=1
t_dyn=0.1; t_lim=10.; dsp=1; Xflt=2
bod_frc=0; hyb=1; nviz_dyn=100; init=0; rsf=0
alpha=0.; beta=0.0025
dt_dyn=1E-3; nviz_wave=16; nviz_slip=40
line1 = ["fault tet 29"]
line3 = np.array([t,dt,nviz,dsp]).reshape(1,4)
line4 = np.array([t_dyn,dt_dyn,nviz_dyn,t_lim,dsp_hyb,Xflt,bod_frc,hyb,rsf,init]).reshape(1,10)
line5 = np.array([nviz_wave,nviz_slip]).reshape(1,2)
line6 = np.array([alpha,beta]).reshape(1,2)

# read mesh
nc = netCDF4.Dataset(fin)
print 'Extracting mesh...'
coord = np.hstack((nc.variables['coordx'][:].\
        reshape(len(nc.variables['coordx']),1),
        nc.variables['coordy'][:].\
        reshape(len(nc.variables['coordy']),1),
        nc.variables['coordz'][:].\
        reshape(len(nc.variables['coordz']),1)))
tet_node = np.empty(shape=[0, 4], dtype=np.uint32)
nnd = len(coord)

# velocity model
vp=6000.; vs=3464.; rho=2670.
E=rho*vs**2*(3*vp**2-4*vs**2)/(vp**2-vs**2)
nu=(vp**2-2*vs**2)/2/(vp**2-vs**2)
visc=1E25 # solid viscosity
r=1.0 # power law
mat = [[E,nu,visc,r,rho,E,nu],
       [E,nu,visc,r,rho,E,nu]]
mat_typ = np.empty(shape = (0,1), dtype=np.uint32)
for i in nc.variables['eb_prop1'][:]:
    cnct = nc.variables["connect"+str(i)][:]
    n_elem = len(cnct)
    tet_node = np.vstack((tet_node, cnct))
    mat_typ = np.vstack((mat_typ, i*np.ones((len(cnct),1))))
print '%d nodes, %d elements' %(nnd, len(tet_node))

# Observation locations
ogrid=np.array([[ -3,-2, 0],
                [ -3, 2, 0],
                [ -3, 5, 0],
                [ -3, 8, 0],
                [0.6, 2, 0],
                [1.4, 5, 0],
                [2.3, 8, 0],
                [  3,-2, 0],
                [4.2, 2, 0],
                [5.9, 5, 0],
                [7.6, 8, 0]])

# fixed nodes
bcx_nodes = nc.variables['node_ns9' ][:]
bcy_nodes = nc.variables['node_ns10'][:]
bcz_nodes = nc.variables['node_ns11'][:]
bc_typ = np.ones((nnd,3), dtype=np.int8)
for node in bcx_nodes:
    bc_typ[node - 1, 0] = 0
for node in bcy_nodes:
    bc_typ[node - 1, 1] = 0
for node in bcz_nodes:
    bc_typ[node - 1, 2] = 0

# Traction and abs boundaries,id 1~8 reserved for fault faces
bnd_el = []
for i in nc.variables['ss_prop1'][8:]:
    els = nc.variables['elem_ss' + str(i)][:]
    sides = nc.variables['side_ss' + str(i)][:]
    bnd_el.append(np.hstack((els.reshape(len(els),1),sides.reshape(len(sides),1))))
trac_el1=bnd_el[3]
trac_el2=bnd_el[4]
trac_el3=bnd_el[5]
abs_bc1=bnd_el[0]
abs_bc2=bnd_el[1]
abs_bc3=bnd_el[2]

#--------------uniform traction BC-------------
trac_val=[-1E3, -1E3, -1E3]
trac_bc1=np.zeros(shape=[len(trac_el1),5])
trac_bc2=np.zeros(shape=[len(trac_el2),5])
trac_bc3=np.zeros(shape=[len(trac_el3),5])
trac_bc1[:,0]=trac_val[0]; trac_bc1[:,3]=0.; trac_bc1[:,4]=0.
trac_bc2[:,1]=trac_val[1]; trac_bc2[:,3]=0.; trac_bc2[:,4]=0.
trac_bc3[:,2]=trac_val[2]; trac_bc3[:,3]=0.; trac_bc3[:,4]=0.
trac_el = np.vstack((trac_el1, trac_el2, trac_el3))
trac_bc = np.vstack((trac_bc1, trac_bc2, trac_bc3))

# absorbing BC
abs_bc1=np.hstack((abs_bc1,   np.ones((len(abs_bc1),1))))
abs_bc2=np.hstack((abs_bc2, 2*np.ones((len(abs_bc2),1))))
abs_bc3=np.hstack((abs_bc3, 3*np.ones((len(abs_bc3),1))))
abs_bc4=np.hstack((trac_el1,  np.ones((len(trac_el1),1))))
abs_bc5=np.hstack((trac_el2,2*np.ones((len(trac_el2),1))))
abs_bc6=np.hstack((trac_el3,3*np.ones((len(trac_el3),1))))
# abs_bc6 (upper bound) is usually traction free, not absorbing.
abs_bc = np.vstack((abs_bc1, abs_bc2, abs_bc3, abs_bc4, abs_bc5))

# Fault nodes, direction vector hard coded for planar fault
print 'Forming fault constraints...'
ft_pos_nodes = []
ft_neg_nodes = []
for i in range(8):
    if (i+1)%2 == 0:
        ft_neg_nodes.append(nc.variables['node_ns%d' %(i+1)][:])
    else:
        ft_pos_nodes.append(nc.variables['node_ns%d' %(i+1)][:])

# Split nodes at intersection
ft_x= [np.intersect1d(ft_pos_nodes[0],ft_pos_nodes[3]),
       np.intersect1d(ft_neg_nodes[0],ft_pos_nodes[2]),
       np.intersect1d(ft_neg_nodes[1],ft_neg_nodes[2]),
       np.intersect1d(ft_pos_nodes[1],ft_neg_nodes[3])]
tmp=np.intersect1d(ft_x[0],ft_x[2])
ft_x[0]=np.setdiff1d(ft_x[0],tmp); ft_x[2]=np.setdiff1d(ft_x[2],tmp)
tmp=np.intersect1d(ft_x[1],ft_x[3])
ft_x[1]=np.setdiff1d(ft_x[1],tmp); ft_x[3]=np.setdiff1d(ft_x[3],tmp)
# Sort intersection nodes on different side of the fault
coordx = []
ft_xall = []
for i in range(4):
    coordx.append(coord[ft_x[0]-1,:])
ft_map = np.array(np.array(np.all((coordx[0][:,None,:]==coordx[2][None,:,:]),axis=-1).nonzero()).T.tolist())
ft_x[2]=ft_x[2][ft_map[:,1]]
ft_map = np.array(np.array(np.all((coordx[1][:,None,:]==coordx[3][None,:,:]),axis=-1).nonzero()).T.tolist())
nxfnd = len(ft_x[0])+len(ft_x[1])
for i in range(4):
    ft_xall = list(set().union(ft_xall,ft_x[i]))

# Sort non-intersection nodes on different side of the fault
nfnd = 0
for i in range(4):
    # Remove intersection nodes
    ft_pos_nodes[i]=np.setdiff1d(ft_pos_nodes[i],ft_xall)
    ft_neg_nodes[i]=np.setdiff1d(ft_neg_nodes[i],ft_xall)
    coord_pos = coord[ft_pos_nodes[i]-1,:]
    coord_neg = coord[ft_neg_nodes[i]-1,:]
    ft_map = np.array(np.array(np.all((coord_pos[:,None,:]==coord_neg[None,:,:]),axis=-1).nonzero()).T.tolist())
    ft_neg_nodes[i] = ft_neg_nodes[i][ft_map[:,1]]
    id_tmp = np.argsort(ft_map[:,0])
    ft_pos_nodes[i] = ft_pos_nodes[i][id_tmp]
    ft_neg_nodes[i] = ft_neg_nodes[i][id_tmp]
    nfnd += len(ft_pos_nodes[i])-len(np.intersect1d(ft_pos_nodes[i],ft_neg_nodes[i]))
# Add non-cross links
nfnd+=nxfnd

# fault parameters
vec_fs=np.zeros((nfnd,3),dtype=float)
vec_fd=np.zeros((nfnd,3),dtype=float)
vec_fn=np.zeros((nfnd,3),dtype=float)
fc=.677*np.ones((nfnd,1),dtype=np.float)
fcd=.525*np.ones((nfnd,1),dtype=np.float)
dc=0.4*np.ones((nfnd,1),dtype=np.float)
coh=np.zeros((nfnd,1))
dcoh=np.ones((nfnd,1))
st_init=np.zeros((nfnd,3))
frc=np.zeros((nfnd,1),dtype=np.uint32)
j=0; dx=.1; eps=dx*.1
for k in range(4):
    for node_pos,node_neg in zip(ft_pos_nodes[k],ft_neg_nodes[k]):
        if node_pos!=node_neg:
            x=coord[node_pos-1,0]
            y=coord[node_pos-1,1]
            z=coord[node_pos-1,2]
            if k < 2 and abs(y+8.)<=1.5 and abs(z+7.5)<=1.5:
                if idcase == 14:
                    st_init[j,:]=[8.16E7,0.,-12E7]
                else:
                    st_init[j,:]=[-8.16E7,0.,-12E7]
                vec_fn[j,0]=-1.; vec_fn[j,1]=0.
                vec_fs[j,0]=0.; vec_fs[j,1]=1
            elif k < 2:
                if idcase == 14:
                    st_init[j,:]=[7E7,0.,-12E7]
                else:
                    st_init[j,:]=[-7E7,0.,-12E7]
                vec_fn[j,0]=-1.; vec_fn[j,1]=0.
                vec_fs[j,0]=0.; vec_fs[j,1]=1
            elif k > 1: # Auxiliary fault
                if idcase == 14:
                    st_init[j,:]=[7E7,0.,-12E7]
                else:
                    st_init[j,:]=[-7.8E7,0.,-12E7]
                vec_fn[j,0]=-np.sqrt(3)*0.5; vec_fn[j,1]=0.5
                vec_fs[j,0]=0.5; vec_fs[j,1]=np.sqrt(3)*0.5
            vec_fd[j,2]=-1.
            if y<12: frc[j] = 1
            j += 1
# Cross-link nodes at intersection [0<->2] [3<->1]
for node_pos,node_neg in zip(np.hstack((ft_x[0],ft_x[3])),np.hstack((ft_x[2],ft_x[1]))):
    vec_fn[j,0]=-1.; vec_fs[j,1]=1; vec_fd[j,2]=-1.
    frc[j]=1
    if idcase==14:
        st_init[j,:]=[7E7,0.,-12E7]
    else:
        st_init[j,:]=[-7E7,0.,-12E7]
    j += 1
if j != nfnd: print "Warning "+ str(nfnd)+" fault pair expected, "+str(j)+" given."

# Total length of constraint function
neq = dim*nfnd; print '%d constraint equations' %(neq)

# Export to Defmod .inp file
ext = str()
if fin.find('lit') !=-1: ext = '_lit'
if dsp_hyb==1: ext += '_dsp'
fout = 'F3DX%d'%(idcase)+ext+'.inp'
print 'Writing to ' + fout + '...'
if os.path.isfile(fout): os.remove(fout)
f = open(fout, 'a')
# six parameter lines
np.savetxt(f, line1, fmt='%s')
neqNCF=0 # zero nonconformal nodes
neqPIX=0 # zero fixed pressure
nfnode=0 # zero nodal force/flux
nvsrc=0  # zero volume source
line2=np.array([len(tet_node),nnd,len(mat),neq,nfnode,len(trac_el),nvsrc,len(abs_bc),nfnd,len(ogrid),neqNCF,neqPIX]).reshape(1,12)
np.savetxt(f, line2, delimiter=' ', fmt='%d '*10)
np.savetxt(f, line3, delimiter=' ', fmt='%g %g %d %d')
np.savetxt(f, line4, delimiter=' ', fmt='%g %g %d %g '+'%d '*6)
if rsf==1:
    np.savetxt(f, line5, delimiter=' ', fmt='%d %d %g')
else:
    np.savetxt(f, line5, delimiter=' ', fmt='%d %d')
np.savetxt(f, np.hstack((line6,[[nxfnd]])), delimiter=' ', fmt='%g %g %d')
# mesh
np.savetxt(f, np.column_stack((tet_node, mat_typ)), delimiter=' ', fmt='%d '*5)
np.savetxt(f, np.column_stack((coord, bc_typ)) , delimiter = ' ', fmt='%g '*3+ '%d '*3)
np.savetxt(f, mat, delimiter=' ', fmt = '%g '*7)

# fault slip: strike, dip and open, zero for hybrid model
slip = np.array([0.0, 0.0, 0.0]).reshape(3,1)
n=[2]
j=0
vecf = np.empty(shape=(0,11))
xfnd = np.empty(shape=(0,3))
for k in range(4):
    for node_pos, node_neg in zip(ft_pos_nodes[k], ft_neg_nodes[k]):
        if node_pos != node_neg:
            vec1  = [[1, 0, 0, node_pos],
                     [-1, 0, 0, node_neg]]
            vec2  = [[0, 1, 0, node_pos],
                     [0, -1, 0, node_neg]]
            vec3  = [[0, 0, 1,  node_pos],
                     [0, 0, -1, node_neg]]
            mat_ft = np.hstack((vec_fs[j,:].reshape(3,1), vec_fd[j,:].reshape(3,1), vec_fn[j,:].reshape(3,1)))
            mat_f = np.matrix.transpose(mat_ft).reshape(1,9)
            val = np.dot(mat_ft,slip)
            cval1 = np.hstack((val[0], [0.,0.])).reshape(1,3)
            cval2 = np.hstack((val[1], [0.,0.])).reshape(1,3)
            cval3 = np.hstack((val[2], [0.,0.])).reshape(1,3)
            np.savetxt(f, n, fmt = '%d')
            np.savetxt(f, vec1, delimiter = ' ', fmt = '%g %g %g %d')
            np.savetxt(f, cval1, delimiter = ' ', fmt = "%g %g %g")
            np.savetxt(f, n, fmt = '%d')
            np.savetxt(f, vec2, delimiter = ' ', fmt = '%g %g %g %d')
            np.savetxt(f, cval2, delimiter = ' ', fmt = "%g %g %g")
            np.savetxt(f, n, fmt = '%d')
            np.savetxt(f, vec3, delimiter = ' ', fmt = '%g %g %g %d')
            np.savetxt(f, cval3, delimiter = ' ', fmt = "%g %g %g")
            vecf = np.vstack((vecf,np.hstack(([[node_pos, node_neg]], mat_f))))
            xfnd = np.vstack((xfnd, coord[node_pos-1,:]))
            j+=1
# Cross-link nodes
i=0
vecxf=np.empty(shape=(nxfnd,10))
stx_init=np.empty(shape=(nxfnd,3))
for node_pos,node_neg in zip(np.hstack((ft_x[0],ft_x[3])),np.hstack((ft_x[2],ft_x[1]))):
    vec1  = [[1, 0, 0, node_pos],
             [-1, 0, 0, node_neg]]
    vec2  = [[0, 1, 0, node_pos],
             [0, -1, 0, node_neg]]
    vec3  = [[0, 0, 1,  node_pos],
             [0, 0, -1, node_neg]]
    cval = np.array([[0.,0.,0.]])
    np.savetxt(f, n, fmt = '%d')
    np.savetxt(f, vec1, delimiter = ' ', fmt = '%g %g %g %d')
    np.savetxt(f, cval, delimiter = ' ', fmt = "%g %g %g")
    np.savetxt(f, n, fmt = '%d')
    np.savetxt(f, vec2, delimiter = ' ', fmt = '%g %g %g %d')
    np.savetxt(f, cval, delimiter = ' ', fmt = "%g %g %g")
    np.savetxt(f, n, fmt = '%d')
    np.savetxt(f, vec3, delimiter = ' ', fmt = '%g %g %g %d')
    np.savetxt(f, cval, delimiter = ' ', fmt = "%g %g %g")
    mat_f = np.hstack((vec_fs[j,:], vec_fd[j,:], vec_fn[j,:])).reshape(1,9)
    vecf = np.vstack((vecf,np.hstack(([[node_pos, node_neg]], mat_f))))
    xfnd = np.vstack((xfnd, coord[node_pos-1,:]))
    j+=1 # Fortran 1 based nfnd index
    # Auxiliary fault
    if i<len(ft_x[0]):
        vecxf[i,:] = np.hstack((j,[0.5,np.sqrt(3)*0.5,0], [0,0,-1], [-np.sqrt(3)*0.5,0.5,0]))
    else: # Reversely linked node pairs [3<->1]
        vecxf[i,:] = np.hstack((j,[-0.5,-np.sqrt(3)*0.5,0], [0,0,1], [np.sqrt(3)*0.5,-0.5,0]))
    if idcase == 14:
        stx_init[i,:]=[7E7,0.,-12E7]
    else:
        stx_init[i,:]=[-7.8E7,0.,-12E7]
    i+=1

np.savetxt(f,np.hstack((vecf,fc,fcd,dc,st_init,xfnd,frc,coh,dcoh)),delimiter=' ',\
       fmt='%d '*2+'%g '*18+'%d '+'%g '*2)
# Auxiliary fault
np.savetxt(f,np.hstack((vecxf,stx_init)),delimiter=' ', fmt='%d '+'%g '*12)
# Boundary traction
np.savetxt(f,np.column_stack((trac_el,trac_bc)),delimiter=' ',fmt='%d %d '+'%g '*5)
# Observation grid
np.savetxt(f,ogrid,delimiter=' ',fmt='%g '*3)
# Absorbing boundaries
np.savetxt(f,abs_bc,delimiter=' ',fmt='%d %d %d')
f.close();
print 'Defmod file '+fout+' created'
