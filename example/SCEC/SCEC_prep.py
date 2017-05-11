#!/usr/bin/env python
import numpy as np
import os, sys, netCDF4
import argparse 
from scipy.spatial import ConvexHull
from scipy.interpolate import griddata 

ap=argparse.ArgumentParser()
ap.add_argument('-m') # mesh file
ap.add_argument('-p') # problem number
ap.add_argument('-d') # absolute output 1 
# Exodus mesh file
fin= ap.parse_args().m
# SCEC problem index
idcase=int(ap.parse_args().p)
if not idcase in [205,10,102]:
    print 'Choose a problem number in (205, 10, 102)'
    sys.exit(0)
nc = netCDF4.Dataset(fin)
dim=3

# dsp_str=1 output fault stress; bod_frc=1 body force; hyb=1 engage hybrid; rsf=1 rate state friction
t=.05; dt=.05; nviz=1
t_dyn=0.1; t_lim=10.; dsp=1; dsp_str=1
if ap.parse_args().d==None:
    dsp_hyb=0
else:
    dsp_hyb=min(1,int(ap.parse_args().d))
if idcase in [205,10]:
    rsf=0
elif idcase==102:
    rsf=1; v_bg=1E-12
bod_frc=0; hyb=1; nviz_dyn=100; init=0 
alpha=0.; beta=0.0025; rfrac=0

# Dynamic step; output frequency
if idcase==205:
    dt_dyn=0.01; nviz_wave=2; nviz_slip=5;
elif idcase==10:
    dt_dyn=0.0025; nviz_wave=8; nviz_slip=20
elif idcase==102:
    dt_dyn=0.00125; nviz_wave=16; nviz_slip=40
line1 = ["fault hex 29"]
line3 = np.array([t,dt,nviz,dsp]).reshape(1,4)
line4 = np.array([t_dyn,dt_dyn,nviz_dyn,t_lim,dsp_hyb,dsp_str,bod_frc,hyb,rsf,init]).reshape(1,10)
if rsf==1: 
    line5 = np.array([nviz_wave,nviz_slip,v_bg]).reshape(1,3)
else:
    line5 = np.array([nviz_wave,nviz_slip]).reshape(1,2)
line6 = np.array([alpha,beta,rfrac]).reshape(1,3)

# node coordinates 
print 'Extracting mesh...'
coord = np.hstack((nc.variables['coordx'][:].\
        reshape(len(nc.variables['coordx']),1),
        nc.variables['coordy'][:].\
        reshape(len(nc.variables['coordy']),1),
        nc.variables['coordz'][:].\
        reshape(len(nc.variables['coordz']),1)))
hx_node = np.empty(shape=[0, 8], dtype=np.uint32)
nnd = len(coord)

# hex, vert and type
if idcase in [205,102]:
    vp=6000.;vs=3464.;rho=2670.
elif idcase==10:
    vp=5716.;vs=3300.;rho=2700.
E=rho*vs**2*(3*vp**2-4*vs**2)/(vp**2-vs**2)
nu=(vp**2-2*vs**2)/2/(vp**2-vs**2)
visc=1E25 # solid viscosity 
r=1.0 # power law 
if idcase in [205,102]:
    mat = [[E,nu,visc,r,rho,E,nu],
           [E,nu,visc,r,rho,E,nu]]
elif idcase==10:
    mat = [[E,nu,visc,r,rho,E,nu],
           [E,nu,visc,r,rho,E,nu],
           [E,nu,visc,r,rho,E,nu]]
mat_typ = np.empty(shape = (0,1), dtype=np.uint32)
for i in nc.variables['eb_prop1'][:]:
    cnct = nc.variables["connect"+str(i)][:]
    n_elem = len(cnct)
    hx_node = np.vstack((hx_node, cnct))
    mat_typ = np.vstack((mat_typ, i*np.ones((len(cnct),1))))
print '%d nodes, %d elements' %(nnd, len(hx_node))

# Waveform locations
if idcase==205:
    ogrid = np.array([[-3.,-12.,  0.],
                      [-3., 12.,  0.],
                      [-3.,-12.,-7.5],
                      [-3., 12.,-7.5]])
elif idcase==10:
    ogrid = np.array([[-3.,  0.,  0.],
                      [-2.,  0.,  0.],
                      [-1.,  0.,  0.],
                      [ 1.,  0.,  0.],
                      [ 2.,  0.,  0.],
                      [ 3.,  0.,  0.],
                      [ 3., 12.,  0.],
                      [-3., 12.,  0.],
                      [-1.,  0., -.3],
                      [-.5,  0., -.3],
                      [ .5,  0., -.3],
                      [ 1.,  0., -.3]])
elif idcase==102:
    ogrid = np.array([[-6.,-12., 0.],
                      [-6., 12., 0.],
                      [ 6.,-12., 0.],
                      [ 6., 12., 0.],
                      [ 9.,  0., 0.],
                      [-9.,  0., 0.]])

# fixed nodes
bcx_nodes = nc.variables['node_ns3'][:]
bcy_nodes = nc.variables['node_ns4'][:]
bcz_nodes = nc.variables['node_ns5'][:]
bc_typ = np.ones((nnd,3), dtype=np.int8)
for node in bcx_nodes:
    bc_typ[node - 1, 0] = 0
for node in bcy_nodes:
    bc_typ[node - 1, 1] = 0   
for node in bcz_nodes:
    bc_typ[node - 1, 2] = 0

# Traction and abs boundaries,id 1,2 reserved for fault faces
bnd_el = []
for i in nc.variables['ss_prop1'][2:]:
    els = nc.variables['elem_ss' + str(i)][:]
    sides = nc.variables['side_ss' + str(i)][:]
    bnd_el.append(np.hstack((els.reshape(len(els),1),sides.reshape(len(sides),1))))
trac_el1 = bnd_el[3]
trac_el2 = bnd_el[4]
trac_el3 = bnd_el[5]
abs_bc1 = bnd_el[0] 
abs_bc2 = bnd_el[1]
abs_bc3 = bnd_el[2]
#--------------uniform traction BC-------------
if idcase in [205,10]:
    trac_val = [1E5, 1E5, -1E5]
elif idcase==102:
    trac_val = [1E3, 1E3, -1E3]
#----------------------------------------------
trac_bc1 = np.zeros(shape=[len(trac_el1), 5])
trac_bc2 = np.zeros(shape=[len(trac_el2), 5])
trac_bc3 = np.zeros(shape=[len(trac_el3), 5])
trac_bc1[:,0] = trac_val[0]; trac_bc1[:,3] = 0.; trac_bc1[:,4] = 0.
trac_bc2[:,1] = trac_val[1]; trac_bc2[:,3] = 0.; trac_bc2[:,4] = 0.
trac_bc3[:,2] = trac_val[2]; trac_bc3[:,3] = 0.; trac_bc3[:,4] = 0.
trac_el = np.vstack((trac_el1, trac_el2, trac_el3))
trac_bc = np.vstack((trac_bc1, trac_bc2, trac_bc3))
# absorbing BC 
abs_bc1 = np.hstack((abs_bc1, np.ones((len(abs_bc1),1))))
abs_bc2 = np.hstack((abs_bc2, 2*np.ones((len(abs_bc2),1))))
abs_bc3 = np.hstack((abs_bc3, 3*np.ones((len(abs_bc3),1))))
abs_bc4 = np.hstack((trac_el1, np.ones((len(trac_el1),1))))
abs_bc5 = np.hstack((trac_el2, 2*np.ones((len(trac_el2),1))))
abs_bc6 = np.hstack((trac_el3, 3*np.ones((len(trac_el3),1))))
# abs_bc6 (upper bound) is usually traction free, not absorbing.
abs_bc = np.vstack((abs_bc1, abs_bc2, abs_bc3, abs_bc4, abs_bc5))

# Fault nodes, direction vector hard coded for planar fault 
ft_pos_nodes = nc.variables['node_ns1'][:]
ft_neg_nodes = nc.variables['node_ns2'][:]
if idcase in [205,102]: # vertical fault
    vec_fn = np.zeros((len(ft_pos_nodes), 3), dtype=float)
    vec_fn[:,0]=-1.
    vec_fs = np.zeros((len(ft_pos_nodes), 3), dtype=float)
    vec_fs[:,1]=1.
    vec_fd = np.zeros((len(ft_pos_nodes), 3), dtype=float)
    vec_fd[:,2]=-1.
elif idcase==10: # dipping fault
    vec_fn = np.zeros((len(ft_pos_nodes), 3), dtype=float)
    vec_fn[:,0]=-.5*np.sqrt(3); vec_fn[:,2]=.5
    vec_fs = np.zeros((len(ft_pos_nodes), 3), dtype=float)
    vec_fs[:,1]=1.
    vec_fd = np.zeros((len(ft_pos_nodes), 3), dtype=float)
    vec_fd[:,0]=-.5; vec_fd[:,2]=-.5*np.sqrt(3)
nc.close() # end of mesh reading

vecf = np.empty(shape=(0,11))
xfnd = np.empty(shape=(0,3))
# Sort the nodes on different side of the fault
print 'Forming fault constraints...'
coord_pos = coord[ft_pos_nodes - 1,:]
coord_neg = coord[ft_neg_nodes - 1,:]
ft_map = np.array(np.array(np.all((coord_pos[:,None,:]==coord_neg[None,:,:]),axis=-1).nonzero()).T.tolist())
ft_neg_nodes = ft_neg_nodes[ft_map[:,1]]
id_tmp = np.argsort(ft_map[:,0])
ft_pos_nodes = ft_pos_nodes[id_tmp]
ft_neg_nodes = ft_neg_nodes[id_tmp]
nfnd = len(ft_pos_nodes) - len(np.intersect1d(ft_pos_nodes, ft_neg_nodes))

# fault parameters 
j = 0
st_init = np.zeros((nfnd,3))
frc = np.empty((nfnd,1),dtype=np.uint32)
if idcase==205:
    fc = .667*np.ones((nfnd,1),dtype=np.float)
    fcd = .525*np.ones((nfnd,1),dtype=np.float)
    dc = .4*np.ones((nfnd,1),dtype=np.float)
    coh = np.zeros((nfnd,1))
    dcoh = np.ones((nfnd,1))
elif idcase==10: 
    fc = .76*np.ones((nfnd,1),dtype=np.float)
    fcd = .448*np.ones((nfnd,1),dtype=np.float)
    dc = .5*np.ones((nfnd,1),dtype=np.float)
    coh = 2E5*np.empty((nfnd,1),dtype=np.float)
    dcoh = 1E-6*np.empty((nfnd,1),dtype=np.float)
elif idcase==102: # Rate-state parameters
    tau0=75E6;sn0=120E6
    a = .008*np.ones((nfnd,1)); b0=.6*np.ones((nfnd,1));V0=1E-6*np.ones((nfnd,1))
    dtau0=25E6*np.ones((nfnd,1)); b=.012*np.ones((nfnd,1));L=.02*np.ones((nfnd,1))
    W=15.;w=3.;da0=.008;theta_init=np.empty((nfnd,1),dtype=np.float)
    coh = np.zeros((nfnd,1))
    dcoh = np.ones((nfnd,1))
for node_pos, node_neg, i in zip(ft_pos_nodes, ft_neg_nodes, range(len(ft_pos_nodes))):
    if node_pos != node_neg:
        y = coord[node_pos-1,1]
        z = coord[node_pos-1,2]
        if idcase==205:
            if abs(y)<=1.5 and abs(z+7.5)<=1.5:
                st_init[j,:]=[8.16E7, 0., -12E7]
            elif abs(y+7.5)<=1.5 and abs(z+7.5)<=1.5:
                st_init[j,:]=[7.8E7, 0., -12E7]
            elif abs(y-7.5)<=1.5 and abs(z+7.5)<=1.5:
                st_init[j,:]=[6.2E7, 0., -12E7]
            else:
                st_init[j,:]=[7E7, 0., -12E7]
            frc[j] = 1
        elif idcase==10: 
            dip = z*2./np.sqrt(3)
            stn = 7378*dip*1E3
            if abs(y)<=1.5 and abs(dip+12.)<=1.5:
                sts = 2E5 - (.76+.0057)*stn 
                st_init[j,:]=[0., sts, stn]
            elif abs(y)<=15. and dip>=-15.:
                sts = -.55*stn 
                st_init[j,:]=[0., sts, stn]
            else:
                sts = -.55*stn 
                st_init[j,:]=[0., sts, stn]
                fc[j]=1E4
                coh[j]=1E9
        elif idcase==102:
            st_init[j,:]=[tau0,0.,-sn0]
            if abs(y)<=W:
                By=1.
            elif abs(y)>W and abs(y)<W+w:
                By=.5*(1+np.tanh(w/(abs(y)-W-w)+w/(abs(y)-W)))
            else:
                By=0.
            if abs(z+7.5)<=.5*W:
                Bz=1.
            elif abs(z+7.5)>.5*W and abs(z+7.5)<.5*W+w: 
                Bz=.5*(1+np.tanh(w/(abs(z+7.5)-.5*W-w)+w/(abs(z+7.5)-.5*W)))
            else:
                Bz=0.
            a[j]=a[j]+da0*(1.-By*Bz)
            theta_init[j]=L[j]/V0[j]*np.exp((a[j]*np.log(2.*np.sinh(tau0/a[j]/sn0))-b0[j]-a[j]*np.log(v_bg/V0[j]))/b[j])
        frc[j]=1
        j=j+1
# Total length of constraint function
neq = dim*nfnd; print '%d constraint equations' %(neq)

# Export to Defmod .inp file
if dsp_hyb==1:
    fout = fin.rsplit('.')[0] + '_dsp.inp'
else:
    fout = fin.rsplit('.')[0] + '.inp'
print 'Writing to ' + fout + '...' 
if os.path.isfile(fout): os.remove(fout)
f = open(fout, 'a')
# six parameter lines 
np.savetxt(f, line1, fmt='%s')
neqNCF=0 # zero nonconformal nodes
nfnode=0 # zero nodal force/flux 
line2=np.array([len(hx_node),nnd,len(mat),neq,nfnode,len(trac_el),len(abs_bc),nfnd,len(ogrid),neqNCF]).reshape(1,10)
np.savetxt(f, line2, delimiter=' ', fmt='%d '*10)
np.savetxt(f, line3, delimiter=' ', fmt='%g %g %d %d')
np.savetxt(f, line4, delimiter=' ', fmt='%g %g %d %g '+'%d '*6)
if rsf==1:
    np.savetxt(f, line5, delimiter=' ', fmt='%d %d %g')
else:
    np.savetxt(f, line5, delimiter=' ', fmt='%d %d')
np.savetxt(f, line6, delimiter=' ', fmt='%g %g %g')
# mesh
np.savetxt(f, np.column_stack((hx_node, mat_typ)), delimiter=' ', fmt='%d '*9)
np.savetxt(f, np.column_stack((coord, bc_typ)) , delimiter = ' ', fmt='%g %g %g %d %d %d')
np.savetxt(f, mat, delimiter=' ', fmt = '%g '*7)

# fault slip: strike, dip and open, zero for hybrid model
slip = np.array([0.0, 0.0, 0.0]).reshape(3,1)
n=[2]
for node_pos, node_neg, i in zip(ft_pos_nodes, ft_neg_nodes, range(len(ft_pos_nodes))):
    if node_pos != node_neg:
        vec1  = [[1, 0, 0, node_pos], 
                 [-1, 0, 0, node_neg]]
        vec2  = [[0, 1, 0, node_pos], 
                 [0, -1, 0, node_neg]]
        vec3  = [[0, 0, 1,  node_pos], 
                 [0, 0, -1, node_neg]]
        mat_ft = np.hstack((vec_fs[i,:].reshape(3,1), vec_fd[i,:].reshape(3,1), vec_fn[i,:].reshape(3,1)))
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
if idcase in [205,10]:
    np.savetxt(f,np.hstack((vecf,fc,fcd,dc,st_init,xfnd,frc,coh,dcoh)),delimiter=' ',\
       fmt='%d '*2+'%g '*18+'%d '+'%g '*2)
elif idcase==102:
    np.savetxt(f,np.hstack((vecf,b0,V0,dtau0,a,b,L,theta_init,st_init,xfnd,frc,coh,dcoh)),delimiter = ' ',\
       fmt='%d '*2+'%g '*22+'%d '+'%g '*2)
# Form rotated constraint matrix for fault model
for node_pos,node_neg,i in zip(ft_pos_nodes,ft_neg_nodes,range(len(ft_pos_nodes))):
    if node_pos != node_neg:
        mat_ft = np.hstack((vec_fs[i,:].reshape(3,1),vec_fd[i,:].reshape(3,1),vec_fn[i,:].reshape(3,1)))
        vec = np.array([1., 0., 0.]).reshape(3,)
        vec = np.dot(mat_ft, vec).reshape(3,)
        vec1  = [[vec[0], vec[1], vec[2], node_pos], 
                 [-vec[0], -vec[1], -vec[2], node_neg]]
        vec = np.array([0., 1., 0.]).reshape(3,1)
        vec = np.dot(mat_ft, vec).reshape(3,)
        vec2  = [[vec[0], vec[1], vec[2], node_pos], 
                 [-vec[0], -vec[1], -vec[2], node_neg]]
        vec = np.array([0., 0., 1.]).reshape(3,1)
        vec = np.dot(mat_ft, vec).reshape(3,)
        vec3  = [[vec[0], vec[1], vec[2], node_pos], 
                 [-vec[0], -vec[1], -vec[2], node_neg]]
        np.savetxt(f, vec1, delimiter = ' ', fmt = '%g %g %g %d') 
        np.savetxt(f, vec2, delimiter = ' ', fmt = '%g %g %g %d')
        np.savetxt(f, vec3, delimiter = ' ', fmt = '%g %g %g %d')
# Boundary traction
np.savetxt(f,np.column_stack((trac_el,trac_bc)),delimiter=' ',fmt='%d %d '+'%g '*5)
# Observation grid
#np.savetxt(f,np.column_stack((ogrid,obs_nlist,obs_N)),delimiter=' ',fmt='%g '*3+'%d '*8+'%g '*8)
np.savetxt(f,ogrid,delimiter=' ',fmt='%g '*3)
# Absorbing boundaries
np.savetxt(f,abs_bc,delimiter=' ',fmt='%d %d %d')
f.close(); 
print 'Defmod file '+fout+' created'
