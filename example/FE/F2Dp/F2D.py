#!/usr/bin/env python
import numpy as np
import os, sys, netCDF4
fin = sys.argv[1]
nc = netCDF4.Dataset(fin)

# Model header
dim = 2
dt_hr=12.
ndy=72
t = 3600.*24*ndy; dt = 3600.*dt_hr; nviz = 1 
t_dyn = 0.025; dt_dyn=1.6E-5; t_lim = 5.; dsp=1; dsp_hyb=0; dsp_str=1; rsf=1; v_bg=1E-12
bod_frc=0; hyb=1; nviz_dyn=960; nviz_wave=80; nviz_slip=200; init=1
alpha= 0.; beta = 0.00125; rfrac = 0
line1 = ["fault-p qua 12"]
line3 = np.array([t, dt, nviz, dsp]).reshape(1,4)
line4 = np.array([t_dyn, dt_dyn, nviz_dyn, t_lim, dsp_hyb, dsp_str, bod_frc,hyb,rsf,init]).reshape(1,10)
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
# quad data and type
rho=np.array([2700.,2700.,2700.,2700.])
vp =np.array([4200.,3800.,5000.,3000.])
vs =np.array([2400.,2200.,3000.,1800.])
E_dyn=rho*vs**2*(3*vp**2-4*vs**2)/(vp**2-vs**2)
nu_dyn=(vp**2-2*vs**2)/2/(vp**2-vs**2)
E=E_dyn/1.6
nu=nu_dyn
K=np.array([9E-20,5E-14,9E-20,9E-20])/1.5E-4
#solid viscosity; power law
visc = 1E25; r = 1.
B=0.9 #Biot coef
phi=.15 #porosity
cf=1E9 #fluid bulk modulus 
source = -1.8E5 #source/sink -2.05E5
if init==1:
    mat = [[E[0], nu[0], visc, r, rho[0], K[0], B, phi, cf,    0.,E_dyn[0],nu_dyn[0]],
           [E[1], nu[1], visc, r, rho[1], K[1], B, phi, cf,source,E_dyn[1],nu_dyn[1]],
           [E[2], nu[2], visc, r, rho[2], K[2], B, phi, cf,    0.,E_dyn[2],nu_dyn[2]],
           [E[3], nu[3], visc, r, rho[3], K[3], B, phi, cf,    0.,E_dyn[3],nu_dyn[3]]]
else:
    mat = [[3.0E10, 0.25, 1.0E25, 1.0, 3000, 1.0E-15, 0.9, 0.15, 2.2E9],
           [3.0E10, 0.25, 1.0E25, 1.0, 3000, 1.0E-11, 0.9, 0.15, 2.2E9],
           [4.0E10, 0.25, 1.0E25, 1.0, 3000, 1.0E-15, 0.9, 0.15, 2.2E9],
           [3.0E10, 0.25, 1.0E25, 1.0, 3000, 1.0E-15, 0.9, 0.15, 2.2E9]]
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
# remove split node at tip
crd_add=np.empty((0),dtype=bool)
for i in range(len(nd_tap)):
    loc = qd_node==nd_tap[i,1]
    loc_flt = nd_flt_n==nd_tap[i,1]
    if sum(sum(loc))<2:
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
vec_fs /= (np.ones((2,1))*np.linalg.norm(vec_fs, axis=1)).T
vec_fn /= (np.ones((2,1))*np.linalg.norm(vec_fn, axis=1)).T
vecf = np.empty(shape=(0,6))
xfnd = np.empty(shape=(0,2))
nfnd = len(nd_flt_p)

print 'Locating nonconformal nodes...' 
NCF_s2m = np.zeros((0,3),np.uint32) 
NCF_wt  = np.zeros((0,2),np.float)
hit=[]
pair=np.empty((0,2,),dtype=np.uint32)
for s in [6,8]: #6,8
    mst_nodes = nc.variables['node_ns'+str(s)][:]
    slv_nodes = nc.variables['node_ns'+str(s+1)][:]
    mst_nodes=np.setdiff1d(mst_nodes,nd_flt_n)
    mst_nodes=np.setdiff1d(mst_nodes,nd_flt_p)
    slv_nodes=np.setdiff1d(slv_nodes,nd_flt_n)
    slv_nodes=np.setdiff1d(slv_nodes,nd_flt_p)
    mst_coord = coord[mst_nodes-1,:]
    slv2mst = np.zeros((0,3),np.uint32) 
    wt_mst = np.zeros((0,2),np.float)
    slv_nodes = np.setdiff1d(slv_nodes,np.array(hit))
    for slv in slv_nodes:
        dis2nd = np.linalg.norm(mst_coord - np.dot(np.ones((len(mst_nodes),1),dtype=np.float),coord[slv-1,:].reshape(1,2)),axis=1)
        mst2=np.argsort(dis2nd)[:2]
        if min(dis2nd)<1E-6: # Node on node
            ndtmp = mst_nodes[mst2[0]]
            hit.append(slv)
            pair=np.vstack((pair,[slv,ndtmp]))
        else:
            onedg=False
            for el, side in zip(nc.variables['elem_ss'+str(s)][:], nc.variables['side_ss'+str(s)][:]):
                el_node = qd_node[el-1,:]
                if   side ==1: e_node = el_node[[0,1]] 
                elif side ==2: e_node = el_node[[1,2]] 
                elif side ==3: e_node = el_node[[2,3]] 
                elif side ==4: e_node = el_node[[3,0]]
                v2j=coord[e_node[0]-1,:]-coord[slv-1,:]
                v2k=coord[e_node[1]-1,:]-coord[slv-1,:]
                onedg=np.linalg.norm(np.cross(v2j,v2k))<1E-6 and np.dot(v2j,v2k)<0 # Node on edge
                if onedg:
                    dis2j=np.linalg.norm(v2j)
                    dis2k=np.linalg.norm(v2k)
                    wj = dis2k/(dis2j+dis2k); wk = dis2j/(dis2j+dis2k)
                    slv2mst= np.vstack((slv2mst,np.hstack((slv,e_node[0],e_node[1]))))
                    wt_mst=np.vstack((wt_mst,[wj,wk]))
                    hit.append(slv)
                    break
    NCF_s2m=np.vstack((NCF_s2m,slv2mst))
    NCF_wt=np.vstack((NCF_wt,wt_mst))
print '%d nonconformal nodes located' %(len(NCF_s2m))

# Make the fault permeable at certain segment; add initial stress 
perm = np.ones((nfnd,1)); st_init = np.zeros((nfnd,2))
y_min = min(coord[:,1]); y_max = max(coord[:,1])
for node_pos, i in zip(nd_flt_p,range(len(nd_flt_p))):
    y = coord[node_pos - 1,1]
    if y > y_min + (y_max - y_min)*0.6 and y < y_min + (y_max - y_min)*0.75:
        perm[i]=1 
        #st_init[i,:]=[s*2.7E7, -5E7]
        st_init[i,:]=[0., 0.]
    else:
        perm[i]=1
        #st_init[i,:]=[s*2.7E7, -5E7]
        st_init[i,:]=[0., 0.]

if rsf==1: # Rate-state parameters
    a = .008*np.ones((nfnd,1)); b0=0.6*np.ones((nfnd,1));V0=1E-6*np.ones((nfnd,1))
    dtau0=25E6*np.ones((nfnd,1)); b=.012*np.ones((nfnd,1));L=.001*np.ones((nfnd,1))
    theta_init=L/V0*np.exp((a*np.log(2.*np.sinh(0.7/a))-b0-a*np.log(v_bg/V0))/b)

# Observation info
ogrid = np.array([[-2.0, -0.2],
                  [-1.0, -0.2],
                  [ 0.0, -0.2],
                  [ 1.0, -0.2],
                  [ 2.0, -0.2]])

#----------------------Boundary ID-------------------------
# Side/nodesets: 0 fault, 1 upper, 2 left, 3 lower, 4 right  
#----------------------------------------------------------
bnd_el = []
# Traction and abs boundaries,id 0 preserved for fault faces
for i in nc.variables['ss_prop1'][1:]:
    els = nc.variables['elem_ss' + str(i)][:]
    sides = nc.variables['side_ss' + str(i)][:]
    bnd_el.append(np.hstack((els.reshape(len(els),1),\
            sides.reshape(len(sides),1))))
trac_el1 = bnd_el[3]
trac_el2 = bnd_el[2]
abs_bc1 = bnd_el[0]
abs_bc2 = bnd_el[1]

# fixed nodes
bcy_nodes = nc.variables['node_ns2'][:] 
bcx_nodes = nc.variables['node_ns3'][:]
bc_typ = np.ones((len(coord),3), dtype=np.int8)
for node in bcx_nodes:
    bc_typ[node - 1, 0] = 0
for node in bcy_nodes:
    bc_typ[node - 1, 1] = 0

# nodal force/flux 
flux=-250.*(dt_hr/24.) #-144. -288. -112.5 250.
# Constant source 
#fnodes = [[-1, -2, 0, 0, 0,   2*24*3600, t/2.],
#         [.25, -2, 0, 0, flux,2*24*3600, t/2.]]

# Varying source 
df=16
fnodes=np.empty([0,7],dtype=np.float)
for i in range(ndy):
    if (i-2)%df==0: fnodes=np.vstack((fnodes,[.25, -2, 0, 0, flux,i*24*3600, (i+df/2-1)*24*3600]))

fnode_bc = np.empty((0,6))
for fnode in fnodes:
    dis = np.linalg.norm(coord - np.dot(np.ones((len(coord), 1)),\
            np.array(fnode[:2]).reshape(1,2)), axis=1)
    row = np.argsort(dis)[0]
    fnode_bc = np.vstack((fnode_bc, np.hstack((row + 1, fnode[2:]))))

# traction bc 
trac_val = [50E6, -65E6]
trac_bc1 = np.zeros(shape=[len(trac_el1), 5])
trac_bc2 = np.zeros(shape=[len(trac_el2), 5])
trac_bc1[:,0] = trac_val[0]; trac_bc1[:,3] = 0.; trac_bc1[:,4] = 0. 
trac_bc2[:,1] = trac_val[1]; trac_bc2[:,3] = 0.; trac_bc2[:,4] = 0.
trac_el = np.vstack((trac_el1, trac_el2))
trac_bc = np.vstack((trac_bc1, trac_bc2))

# absorbing bc (abs_bc4 upper bound can be free boundary for top surface) 
abs_bc1 = np.hstack((abs_bc1, 2*np.ones((len(abs_bc1),1))))
abs_bc2 = np.hstack((abs_bc2,   np.ones((len(abs_bc2),1))))
abs_bc3 = np.hstack((trac_el1,  np.ones((len(trac_el1),1))))
abs_bc4 = np.hstack((trac_el2,2*np.ones((len(trac_el2),1))))
abs_bc = np.vstack((abs_bc1, abs_bc2, abs_bc3,abs_bc4))

# Total length of constraint function
neqNCF=(dim+1)*len(NCF_s2m)
neqFT=dim*nfnd + sum(perm)
neq = neqNCF+neqFT 
print '%d NCF and %d fault constraint equations.' %(neqNCF,neqFT)

# Eliminate doubling nodes due to NCF interface 
for i in range(len(pair)):
    slv=pair[i,0];mst=pair[i,1]
    qd_node[qd_node==slv]=mst
    NCF_s2m[NCF_s2m==slv]=mst
    nd_flt_n[nd_flt_n==slv]=mst
    nd_flt_p[nd_flt_p==slv]=mst
    fnode_bc[fnode_bc[:,0]==slv,0]=mst
    coord=np.vstack((coord[:slv-1,:],coord[slv:,:]))
    bc_typ=np.vstack((bc_typ[:slv-1,:],bc_typ[slv:,:]))
    id_shift=qd_node>slv
    qd_node[id_shift]=qd_node[id_shift]-1
    id_shift=NCF_s2m>slv
    NCF_s2m[id_shift]=NCF_s2m[id_shift]-1
    id_shift=nd_flt_p>slv
    nd_flt_p[id_shift]=nd_flt_p[id_shift]-1
    id_shift=nd_flt_n>slv
    nd_flt_n[id_shift]=nd_flt_n[id_shift]-1
    id_shift=fnode_bc[:,0]>slv
    fnode_bc[id_shift,0] = fnode_bc[id_shift,0]-1
    id_shift=pair>slv
    pair[id_shift]=pair[id_shift]-1

# Export to Defmod .inp file
if rsf==1:
    ext='_rs.inp'
else:
    ext='_sw.inp'
ext='_%d'%(df)+'.inp'
fout = fin.rsplit('.')[0]+'_dt%d' %(dt_hr) + ext 
print 'Write to ' + fout + '...'
if os.path.isfile(fout): os.remove(fout)
f = open(fout, 'a')
line2 = np.array([len(qd_node), len(coord), len(mat), neq,\
        len(fnode_bc), len(trac_el), len(abs_bc), nfnd, len(ogrid), neqNCF]).reshape(1,10)
np.savetxt(f, line1, fmt='%s')
np.savetxt(f, line2, delimiter=' ', fmt='%d '*10)
np.savetxt(f, line3, delimiter=' ', fmt='%g %g %d %d')
np.savetxt(f, line4, delimiter=' ', fmt='%g %g %d %g %d %d %d %d %d %d')
if rsf==0:
    np.savetxt(f, line5, delimiter=' ', fmt='%d %d')
else:
    np.savetxt(f, line5, delimiter=' ', fmt='%d %d %g')
np.savetxt(f, line6, delimiter=' ', fmt='%g %g %g')
np.savetxt(f, np.column_stack((qd_node, mat_typ)), delimiter=' ',\
           fmt='%d %d %d %d %d')
np.savetxt(f, np.column_stack((coord, bc_typ)) , delimiter = ' ',\
               fmt='%g %g %d %d %d')
if init==1:
    np.savetxt(f, mat, delimiter=' ', fmt = '%g '*12)
else:
    np.savetxt(f, mat, delimiter=' ', fmt = '%g '*9)

# Write NCF constraints
for s2m, wt in zip(NCF_s2m,NCF_wt):
    npt = np.count_nonzero(s2m)
    if npt==3:
        vec1  = [[   -1, 0, 0, s2m[0]], 
                 [wt[0], 0, 0, s2m[1]],
                 [wt[1], 0, 0, s2m[2]]]
        vec2  = [[0,    -1, 0, s2m[0]], 
                 [0, wt[0], 0, s2m[1]],
                 [0, wt[1], 0, s2m[2]]]
        vec3  = [[0, 0,    -1, s2m[0]], 
                 [0, 0, wt[0], s2m[1]],
                 [0, 0, wt[1], s2m[2]]]
    np.savetxt(f, [npt], fmt = '%d') 
    np.savetxt(f, vec1, delimiter = ' ', fmt = '%g %g %g %d') 
    np.savetxt(f, [[0.,0.,0.]], delimiter = ' ', fmt = "%1.2E %g %g") 
    np.savetxt(f, [npt], fmt = '%d')
    np.savetxt(f, vec2, delimiter = ' ', fmt = '%g %g %g %d')
    np.savetxt(f, [[0.,0.,0.]], delimiter = ' ', fmt = "%1.2E %g %g")
    np.savetxt(f, [npt], fmt = '%d')
    np.savetxt(f, vec3, delimiter = ' ', fmt = '%g %g %g %d')
    np.savetxt(f, [[0.,0.,0.]], delimiter = ' ', fmt = "%1.2E %g %g")

# fault slip: strike and open
slip = np.array([0.0, 0.0]).reshape(2,1) 
c = 0; dt_slip=0; t_rup=0
n = [2]
ft_neg_nodes_tap = []
j = 0
for nd_p, nd_n, i in zip(nd_flt_p, nd_flt_n, range(len(nd_flt_p))):
    vec1  = [[1, 0, 0, nd_p], 
            [-1, 0, 0, nd_n]]
    vec2  = [[0, 1, 0, nd_p], 
             [0,-1, 0, nd_n]]
    mat_ft = np.hstack((vec_fs[i,:].reshape(2,1),\
                vec_fn[i,:].reshape(2,1)))
    mat_f = np.matrix.transpose(mat_ft).reshape(1,4)
    val = np.dot(mat_ft,slip)
    y = np.array(coord[nd_p - 1,:][1])
    t_act =  dt_slip+(1-c)*t_rup 
    t_slip = [t_act-dt_slip/2, t_act+dt_slip/2]
    cval1 = np.hstack((c*val[0], t_slip)).reshape(1,3)
    cval2 = np.hstack((c*val[1], t_slip)).reshape(1,3)
    np.savetxt(f, n, fmt = '%d') 
    np.savetxt(f, vec1, delimiter = ' ', fmt = '%g %g %g %d') 
    np.savetxt(f, cval1, delimiter = ' ', fmt = "%1.2E %g %g") 
    np.savetxt(f, n, fmt = '%d')
    np.savetxt(f, vec2, delimiter = ' ', fmt = '%g %g %g %d')
    np.savetxt(f, cval2, delimiter = ' ', fmt = "%1.2E %g %g")
    vecf = np.vstack((vecf,np.hstack(([[nd_p, nd_n]], mat_f))))
    xfnd = np.vstack((xfnd,coord[nd_p-1,:]))
    # permeable fault 
    if perm[j] > 0: 
        vec4 = [[0, 0,  1, nd_p], 
                [0, 0, -1, nd_n]]
        cval4 =[[0, 0, 0]]
        np.savetxt(f, n, fmt = '%d')
        np.savetxt(f, vec4, delimiter = ' ', fmt = '%g %g %g %d')
        np.savetxt(f, cval4, delimiter = ' ', fmt = "%1.2E %g %g")
    j = j + 1

# Define frictional parameters (static friction)
fc = np.empty((len(vecf),1),dtype=float)
fcd = np.empty((len(vecf),1),dtype=float)
dc = np.empty((len(vecf),1),dtype=float)
frc = np.empty((len(vecf),1),dtype=np.uint32)
for node_pos, i in zip(nd_flt_p,range(len(nd_flt_p))):
    x,y = coord[node_pos - 1,:]
    if y > y_min + (y_max - y_min)*0.15 and y < y_min + (y_max - y_min)*0.85:
        fc[i] = .7
        fcd[i] = .4
        dc[i] = .01
        frc[i] = 1
    else:
        fc[i] = .7
        fcd[i] = .4
        dc[i] = .01
        frc[i] = 0 
coh = np.zeros((len(vecf),1),dtype=float)
dcoh = np.ones((len(vecf),1),dtype=float)
# Write fault orientation tensor + frictional parameters
if rsf==1:
    np.savetxt(f, np.hstack((vecf, b0, V0, dtau0, a, b, L, theta_init, perm, st_init, xfnd, frc,coh,dcoh)), delimiter = ' ',\
       fmt = '%d '*2 + '%g '*4 + '%g '*7 + '%d ' + '%g '*2 + '%g '*2 + '%d ' + '%g '*2)
else:
    np.savetxt(f, np.hstack((vecf, fc, fcd, dc, perm, st_init, xfnd, frc,coh,dcoh)), delimiter = ' ',\
       fmt = '%d '*2 + '%g '*4 + '%g '*3 + '%d ' + '%g '*2 + '%g '*2 + '%d ' + '%g '*2)

# Form rotated constraint matrix for fault model
for nd_p, nd_n, i in zip(nd_flt_p, nd_flt_n, range(len(nd_flt_p))):
    mat_ft = np.hstack((vec_fs[i,:].reshape(2,1),vec_fn[i,:].reshape(2,1)))
    vec = np.array([1., 0.]).reshape(2,1)
    vec = np.dot(mat_ft, vec).reshape(2,)
    vec1  = [[ vec[0],  vec[1], nd_p], 
             [-vec[0], -vec[1], nd_n]]
    vec = np.array([0., 1.]).reshape(2,1)
    vec = np.dot(mat_ft, vec).reshape(2,)
    vec2  = [[ vec[0], vec[1],  nd_p], 
             [-vec[0], -vec[1], nd_n]]
    np.savetxt(f, vec1, delimiter = ' ', fmt = '%g %g %d') 
    np.savetxt(f, vec2, delimiter = ' ', fmt = '%g %g %d')

# Point force/source and boundary traction/flux
np.savetxt(f, fnode_bc, delimiter=' ',\
        fmt ='%d %1.2E %1.2E %1.2E %g %g')
np.savetxt(f, np.column_stack((trac_el, trac_bc)), delimiter=' ',\
        fmt ='%d %d %1.2E %1.2E %1.2E %g %g')

# Observation grid
if len(ogrid)>0:
    np.savetxt(f, ogrid , delimiter = ' ', fmt='%g '*2)
# Abs boundary
np.savetxt(f, abs_bc, delimiter=' ', fmt='%d %d %d')
f.close(); 
print 'Defmod file ' + fout + ' created'

# h* and CFL condition
lbd=rho*(vp**2-2*vs**2); mu=rho*vs**2
if rsf==1:
    L=np.mean(L); a=np.mean(a); b=np.mean(b)
    sigma_e=max(map(abs,trac_val))
    hstar=min(2*mu*L/(b-a)/np.pi/sigma_e)
    print "Critical RSF distance h*=%0.3f m" %hstar
Lc=min(dt_dyn*np.sqrt(E/rho))
print "Critical element length h=%0.3f m" %Lc
