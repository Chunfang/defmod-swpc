#!/usr/bin/env python
import numpy as np
import os, sys, netCDF4
import argparse

ap=argparse.ArgumentParser()
ap.add_argument('-m')     # mesh file
ap.add_argument('-p')     # problem number
ap.add_argument('-alpha') # implicit solver 1 
ap.add_argument('-dsp')   # absolute output 1

# Exodus mesh file
fin= ap.parse_args().m
idcase=int(ap.parse_args().p)
if not idcase in [10,11,102]:
    print 'Choose a problem number in (10, 11, 102)'
    sys.exit(0)
nc = netCDF4.Dataset(fin)
dim = 2

# Define model type
if ap.parse_args().dsp==None:
    dsp_hyb = 0
else: 
    dsp_hyb = min(1,int(ap.parse_args().dsp))
if ap.parse_args().alpha==None:
    galpha = 0
else: 
    galpha = min(1,int(ap.parse_args().alpha)) 
if galpha:
    t = 16; dt = 0.0025; dt_dyn = dt
else:
    t = 0.07; dt = 0.07; dt_dyn = 0.0025
nviz      = 16 
nviz_dyn  = nviz 
nviz_wave = nviz 
nviz_slip = nviz 
dsp=1; dsp_str=1; rsf=0; bod_frc=0; hyb=1; t_lim=10.; init=0
alpha= 0.; beta = 0.00125; rfrac = 0
if idcase==102: rsf=1; v_bg=1E-12
if galpha:
    line1 = ["alpha qua 12"]
else:
    line1 = ["fault qua 12"]
line3 = np.array([t, dt, nviz, dsp]).reshape(1,4)
line4 = np.array([dt_dyn*1E2, dt_dyn, nviz_dyn, t_lim, dsp_hyb, dsp_str, bod_frc,hyb,rsf,init]).reshape(1,10)
if rsf==0:
    line5 = np.array([nviz_wave,nviz_slip]).reshape(1,2)
else:
    line5 = np.array([nviz_wave,nviz_slip,v_bg]).reshape(1,3)
line6 = np.array([alpha, beta, rfrac]).reshape(1,3)

# node data
print 'Extracting mesh...'
coord = np.hstack((nc.variables['coordx'][:].\
        reshape(len(nc.variables['coordx']),1),\
        nc.variables['coordy'][:].\
        reshape(len(nc.variables['coordy']),1)))
qd_node = np.empty(shape=[0, 4], dtype=np.uint32)
rho=np.array([2700.])
vp =np.array([5716.])
vs =np.array([3300.])
E_dyn=rho*vs**2*(3*vp**2-4*vs**2)/(vp**2-vs**2)
nu_dyn=(vp**2-2*vs**2)/2/(vp**2-vs**2)
E=E_dyn
nu=nu_dyn
#solid viscosity; power law
visc = 1E25; r = 1.
mat = [[E[0],nu[0],visc,r,rho[0],E_dyn[0],nu_dyn[0]]]

mat_typ = np.empty(shape = (0,1), dtype=np.uint32)
for i in nc.variables['eb_prop1'][:]:
    cnct = nc.variables["connect"+str(i)][:]
    n_elem = len(cnct)
    cnct = cnct.reshape(n_elem*4)
    cnct = cnct.reshape(n_elem,4) 
    qd_node = np.vstack((qd_node, cnct))
    mat_typ = np.vstack((mat_typ, i*np.ones((len(cnct),1))))
print '%d nodes, %d elements' %(len(coord), len(qd_node))

print 'Forming fault constraints...'
id_tmp = nc.variables['ss_prop1'][0]
el_flt = nc.variables['elem_ss' + str(id_tmp)]
sd_flt = nc.variables['side_ss' + str(id_tmp)]
nd_flt = qd_node[el_flt[:]-1,:]
nd_tap = np.empty((0,2),dtype=np.uint32)
nd_flt_p = np.empty((0),dtype=np.uint32)
crd_flt_p = np.empty((0,2),dtype=float)
nd_flt_n = np.empty((0),dtype=np.uint32)
crd_flt_n = np.empty((0,2),dtype=float)
sd_flt_p = np.empty((0),dtype=np.uint32)
spair_flt = np.empty((0,2),dtype=np.uint32)
for i in range(len(el_flt)):
    el = el_flt[i]
    sd = sd_flt[i]
    nd = nd_flt[i,:]
    if   sd ==1: nd_on = nd[[0,1]];idn=[0,1];nd_off=nd[[2,3]] 
    elif sd ==2: nd_on = nd[[1,2]];idn=[1,2];nd_off=nd[[0,3]] 
    elif sd ==3: nd_on = nd[[2,3]];idn=[2,3];nd_off=nd[[0,1]] 
    elif sd ==4: nd_on = nd[[3,0]];idn=[3,0];nd_off=nd[[1,2]] 
    # negative side has elements lower right than the fault
    if sum(coord[nd_on-1,1]-coord[nd_on-1,0]) > sum(coord[nd_off-1,1]-coord[nd_off-1,0]):
        for j in range(2):
            if nd_on[j] in nd_tap[:,0]: 
                nd_add = nd_tap[nd_tap[:,0]==nd_on[j],1]
                qd_node[el-1,idn[j]] = nd_add
            else:
                nd_add = len(coord)+len(nd_tap)+1
                qd_node[el-1,idn[j]] = nd_add
                nd_tap = np.vstack((nd_tap,[nd_on[j],nd_add]))
                nd_flt_n = np.hstack((nd_flt_n, nd_add))
                crd_flt_n = np.vstack((crd_flt_n, coord[nd_on[j]-1,:])) 
    else:
        for j in range(2):
            if not (nd_on[j] in nd_flt_p):
                nd_flt_p = np.hstack((nd_flt_p, nd_on[j]))
                crd_flt_p = np.vstack((crd_flt_p, coord[nd_on[j]-1,:]))
        sd_flt_p = np.hstack((sd_flt_p,sd))        
        spair_flt = np.vstack((spair_flt,nd_on))
# remove split node at tip, except at free boundary y==0
crd_add=np.empty((0),dtype=bool)
for i in range(len(nd_tap)):
    loc = qd_node==nd_tap[i,1]
    loc_flt = nd_flt_n==nd_tap[i,1]
    if sum(sum(loc))<2 and coord[nd_tap[i,0]-1,1] < 0.:
        qd_node[loc] = nd_tap[i,0]
        nd_flt_n[loc_flt] = nd_tap[i,0]
        crd_add=np.hstack((crd_add,False))
    else:
        qd_node[loc] = qd_node[loc]-sum(~crd_add)
        nd_flt_n[loc_flt] = nd_flt_n[loc_flt]-sum(~crd_add)
        crd_add=np.hstack((crd_add,True))
coord = np.vstack((coord,crd_flt_n[crd_add,:]))
# Pair fault nodes
ft_map = np.array(np.array(np.all((crd_flt_p[:,None,:]\
       ==crd_flt_n[None,:,:]),axis=-1).nonzero()).T.tolist())
nd_flt_n = nd_flt_n[ft_map[:,1]]
nd_flt_p = nd_flt_p[ft_map[:,0]]
crd_flt_n = crd_flt_n[ft_map[:,1],:]
crd_flt_p = crd_flt_n[ft_map[:,0],:]
crd_add = crd_add[ft_map[:,1]]
nd_flt_n = nd_flt_n[crd_add]
nd_flt_p = nd_flt_p[crd_add]
crd_flt_n = crd_flt_n[crd_add,:]
crd_flt_p = crd_flt_p[crd_add,:]

# Fault's strike and normal vectors
vec_fn = np.zeros((len(nd_flt_p), 2), dtype=float)
vec_fs = np.zeros((len(nd_flt_p), 2), dtype=float)
for i in range(len(spair_flt)):
    v1 = coord[spair_flt[i,0]-1,:]
    v2 = coord[spair_flt[i,1]-1,:]
    vec1 = v2 - v1
    vec1 = -vec1*np.sign(vec1[1])
    vec2 = np.array([vec1[1], -vec1[0]])
    vec2 = vec2*np.sign(vec2[1])
    row = np.squeeze(np.hstack((np.where(nd_flt_p==spair_flt[i,0]),\
          np.where(nd_flt_p==spair_flt[i,1]))))
    vec_fs[row,:] += vec1 
    vec_fn[row,:] += vec2
    if idcase==102:
        vec_fs[row,:] += [-1.,0.] 
        vec_fn[row,:] += [0., 1.]
vec_fs /= (np.ones((2,1))*np.linalg.norm(vec_fs, axis=1)).T
vec_fn /= (np.ones((2,1))*np.linalg.norm(vec_fn, axis=1)).T
vecf = np.empty(shape=(0,6))
xfnd = np.empty(shape=(0,2))
nfnd = len(nd_flt_p)

# Define frictional parameters (static friction), and add initial stress. 
st_init = np.zeros((nfnd,2))
frc = np.empty((nfnd,1),dtype=np.uint32)
if idcase in [10, 11]:
    fc      = np.empty((nfnd,1),dtype=float)
    fcd     = np.empty((nfnd,1),dtype=float)
    dc      = np.empty((nfnd,1),dtype=float)
    coh     = np.empty((nfnd,1),dtype=float)
    dcoh    = np.empty((nfnd,1),dtype=float)
    if idcase==10:
        fc_stat = .76
    else:
        fc_stat = 0.57 # SCEC10 0.76, SCEC11 0.57 
    for node_pos, i in zip(nd_flt_p,range(len(nd_flt_p))):
        y = coord[node_pos - 1,1]
        dip = y*2./np.sqrt(3) 
        stn = 7378*dip*1E3
        if abs(dip+12.)<=1.5:
            sts          = 2E5 - (fc_stat+.0057)*stn 
            st_init[i,:] = [sts, stn]
            fc[i]        = fc_stat 
            fcd[i]       = .448
            dc[i]        = .5
            coh[i]       = 2E5
            dcoh[i]      = 1E-6
            frc[i]       = 1
        elif dip>=-15.:
            sts          = -.55*stn 
            st_init[i,:] = [sts, stn]
            fc[i]        = fc_stat 
            fcd[i]       = .448
            dc[i]        = .5
            coh[i]       = 2E5
            dcoh[i]      = 1E-6
            frc[i]       = 1
        else:
            sts          = -.55*stn 
            st_init[i,:] = [sts, stn]
            fc[i]        = 1E4
            fcd[i]       = .448
            dc[i]        = .5
            coh[i]       = 1E9
            dcoh[i]      = 1E-6
            frc[i]       = 1
else: # Rate-state parameters
    tau0=75E6;sn0=120E6
    a = .008*np.ones((nfnd,1)); b0=0.6*np.ones((nfnd,1));V0=1E-6*np.ones((nfnd,1))
    dtau0=25E6*np.ones((nfnd,1)); b=.012*np.ones((nfnd,1));L=.02*np.ones((nfnd,1))
    W=15.;w=3.;da0=.008;theta_init=np.empty((nfnd,1),dtype=np.float)
    coh = np.empty((nfnd,1),dtype=float)
    dcoh = np.ones((nfnd,1))*1E-6
    for node_pos, i in zip(nd_flt_p,range(len(nd_flt_p))):
        x = coord[node_pos - 1,0]
        st_init[i,:]=[tau0,-sn0]
        if abs(x)<=W:
            Bx=1.
            frc[i]=1
        elif abs(x)>W and abs(x)<W+w:
            Bx=.5*(1+np.tanh(w/(abs(x)-W-w)+w/(abs(x)-W)))
            frc[i]=1
        else:
            Bx=0.
            frc[i]=0
        if abs(x)<w/2:
            coh[i]=0.
        else:
            coh[i]=1.75E6
        a[i]=a[i]+da0*(1.-Bx)
        theta_init[i]=L[i]/V0[i]*np.exp((a[i]*np.log(2.*np.sinh(tau0/a[i]/sn0))-b0[i]-a[i]*np.log(v_bg/V0[i]))/b[i])

# Observation info
if idcase==102:
    ogrid = np.array([[-12., -6.],
                      [ 12., -6.],
                      [-12.,  6.],
                      [ 12.,  6.],
                      [  0.,  9.],
                      [  0., -9.]])
elif idcase in [10, 11]:
    ogrid = np.array([[-3.,  0.],
                      [-2.,  0.],
                      [-1.,  0.],
                      [ 1.,  0.],
                      [ 2.,  0.],
                      [ 3.,  0.],
                      [-1., -.3],
                      [-.5, -.3],
                      [ .5, -.3],
                      [ 1., -.3]])

#----------------------Boundary ID-------------------------
# Side/nodesets: 0 fault, 1 upper, 2 left, 3 lower, 4 right  
#----------------------------------------------------------
bnd_el = []
# Traction and abs boundaries,id 0 preserved for fault faces
for i in nc.variables['ss_prop1'][1:]:
    els = nc.variables['elem_ss' + str(i)][:]
    sides = nc.variables['side_ss' + str(i)][:]
    bnd_el.append(np.hstack((els.reshape(len(els),1),sides.reshape(len(sides),1))))
trac_el1 = bnd_el[3]
trac_el2 = bnd_el[2]
abs_bc1 = bnd_el[0]
abs_bc2 = bnd_el[1]

# fixed nodes
bcy_nodes = nc.variables['node_ns2'][:] 
bcx_nodes = nc.variables['node_ns3'][:]
bc_typ = np.ones((len(coord),2), dtype=np.int8)
for node in bcx_nodes:
    bc_typ[node - 1, 0] = 0
for node in bcy_nodes:
    bc_typ[node - 1, 1] = 0

# traction bc 
trac_val = [5E4, -5E4]
trac_bc1 = np.zeros(shape=[len(trac_el1), 4])
trac_bc2 = np.zeros(shape=[len(trac_el2), 4])
trac_bc1[:,0] = trac_val[0]; trac_bc1[:,2] = 0.; trac_bc1[:,3] = 0.
trac_bc2[:,1] = trac_val[1]; trac_bc2[:,2] = 0.; trac_bc2[:,3] = 0.
trac_el = np.vstack((trac_el1, trac_el2))
trac_bc = np.vstack((trac_bc1, trac_bc2))

# absorbing bc 
abs_bc1 = np.hstack((abs_bc1, 2*np.ones((len(abs_bc1),1))))
abs_bc2 = np.hstack((abs_bc2,   np.ones((len(abs_bc2),1))))
abs_bc3 = np.hstack((trac_el1,  np.ones((len(trac_el1),1))))
abs_bc4 = np.hstack((trac_el2,2*np.ones((len(trac_el2),1))))
if idcase==102:
    abs_bc = np.vstack((abs_bc1, abs_bc2, abs_bc3, abs_bc4))
else:
    abs_bc = np.vstack((abs_bc1, abs_bc2, abs_bc3))
#  abs_bc4 upper bound can be free boundary for top surface 

# Total length of constraint function
neqNCF=0 # Zero non-conforming constraint equations.
nfnode=0
neqFT=dim*nfnd
neq = neqNCF+neqFT 
print '%d NCF and %d fault constraint equations.' %(neqNCF,neqFT)

# Export to Defmod .inp file
fout = 'SCEC'+str(idcase)+'-2D' 
if dsp_hyb: fout += '_dsp' 
if galpha: fout += '-alpha'
fout += '.inp' 
print 'Write to ' + fout + '...'
if os.path.isfile(fout): os.remove(fout)
f = open(fout, 'a')
line2 = np.array([len(qd_node), len(coord), len(mat), neq,\
        nfnode, len(trac_el), len(abs_bc), nfnd, len(ogrid), neqNCF]).reshape(1,10)
np.savetxt(f, line1, fmt='%s')
np.savetxt(f, line2, delimiter=' ', fmt='%d '*10)
np.savetxt(f, line3, delimiter=' ', fmt='%g %g %d %d')
np.savetxt(f, line4, delimiter=' ', fmt='%g %g %d %g %d %d %d %d %d %d')
if rsf==0:
    np.savetxt(f, line5, delimiter=' ', fmt='%d %d')
else:
    np.savetxt(f, line5, delimiter=' ', fmt='%d %d %g')
np.savetxt(f, line6, delimiter=' ', fmt='%g %g %g')
np.savetxt(f, np.column_stack((qd_node, mat_typ)), delimiter=' ', fmt='%d %d %d %d %d')
np.savetxt(f, np.column_stack((coord, bc_typ)) , delimiter = ' ', fmt='%g %g %d %d')
np.savetxt(f, mat, delimiter=' ', fmt = '%g '*7)
# fault slip: strike and open
n = [2]
ft_neg_nodes_tap = []
t_slip = [0., 0.]; val=0.
cval1 = np.hstack((val, t_slip)).reshape(1,3)
cval2 = np.hstack((val, t_slip)).reshape(1,3)
for nd_p, nd_n, i in zip(nd_flt_p, nd_flt_n, range(len(nd_flt_p))):
    vec1  = [[1, 0, nd_p], 
            [-1, 0, nd_n]]
    vec2  = [[0, 1, nd_p], 
             [0,-1, nd_n]]
    mat_ft = np.hstack((vec_fs[i,:].reshape(2,1),\
                vec_fn[i,:].reshape(2,1)))
    mat_f = np.matrix.transpose(mat_ft).reshape(1,4)
    np.savetxt(f, n, fmt = '%d') 
    np.savetxt(f, vec1, delimiter = ' ', fmt = '%g %g %d') 
    np.savetxt(f, cval1, delimiter = ' ', fmt = "%1.2E %g %g") 
    np.savetxt(f, n, fmt = '%d')
    np.savetxt(f, vec2, delimiter = ' ', fmt = '%g %g %d')
    np.savetxt(f, cval2, delimiter = ' ', fmt = "%1.2E %g %g")
    vecf = np.vstack((vecf,np.hstack(([[nd_p, nd_n]], mat_f))))
    xfnd = np.vstack((xfnd,coord[nd_p-1,:]))

# Write fault orientation tensor + frictional parameters
if rsf==1:
    np.savetxt(f, np.hstack((vecf, b0, V0, dtau0, a, b, L, theta_init, st_init, xfnd, frc,coh,dcoh)), delimiter = ' ',\
       fmt = '%d '*2 + '%g '*4 + '%g '*7 +         '%g '*2 + '%g '*2 + '%d ' + '%g '*2)
else:
    np.savetxt(f, np.hstack((vecf, fc, fcd, dc, st_init, xfnd, frc,coh,dcoh)), delimiter = ' ', \
           fmt = '%d '*2 + '%g '*4 + '%g '*3 +         '%g '*2 + '%g '*2 + '%d ' + '%g '*2)


# Point force/source and boundary traction/flux
np.savetxt(f, np.column_stack((trac_el, trac_bc)), delimiter=' ',\
        fmt ='%d %d %1.2E %1.2E %g %g') 

# Observation grid
if len(ogrid)>0:
    np.savetxt(f, ogrid , delimiter = ' ', fmt='%g '*2)
# Abs boundary
np.savetxt(f, abs_bc, delimiter=' ', fmt='%d %d %d')
f.close(); 
print 'Defmod file ' + fout + ' created'
