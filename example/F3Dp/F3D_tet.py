#!/usr/bin/env python
import numpy as np
import os, sys, netCDF4
from scipy.interpolate import griddata 
from scipy.spatial import ConvexHull

fin = sys.argv[1]
nc = netCDF4.Dataset(fin)
fault = True 
poro = True 
visc = False
explicit = False
dim=3
# dsp_str=1 output fault stress; bod_frc=1 body force; hyb=1 engage hybrid
t=3600*24*30.; dt=3600*24.; nviz=1
t_dyn=0.01; dt_dyn=0.00008; t_lim=5.; dsp=1; dsp_hyb=0; dsp_str=1; rsf=0
bod_frc=0; hyb=1; nviz_dyn=120; nviz_wave=10; nviz_slip=100; init=1 
alpha=0.; beta=0.0025; rfrac=0
if  poro and visc:
    line1 = ["fault-pv tet 29"]
elif poro and not visc:
    line1 = ["fault-p tet 29"]
elif not poro and visc:
    line1 = ["fault-v tet 29"]
else:
    line1 = ["fault tet 29"]
line3 = np.array([t, dt, nviz, dsp]).reshape(1,4)
line4 = np.array([t_dyn, dt_dyn, nviz_dyn, t_lim, dsp_hyb, dsp_str, bod_frc, hyb, rsf, init]).reshape(1,10)
line5 = np.array([nviz_wave,nviz_slip]).reshape(1,2)
line6 = np.array([alpha, beta, rfrac]).reshape(1,3)

# node data
print 'Extracting mesh...'
coord = np.hstack((nc.variables['coordx'][:].\
        reshape(len(nc.variables['coordx']),1),
        nc.variables['coordy'][:].\
        reshape(len(nc.variables['coordy']),1),
        nc.variables['coordz'][:].\
        reshape(len(nc.variables['coordz']),1)))
tet_node = np.empty(shape=[0, 4], dtype=np.uint32)
nnd = len(coord)

# tet data and type
rho=np.array([2700.,2700.,2700.,2700.,2700.])
vp =np.array([4200.,3800.,4000.,5000.,3000.])
vs =np.array([2400.,2200.,2500.,3000.,1800.])
E_dyn=rho*vs**2*(3*vp**2-4*vs**2)/(vp**2-vs**2)
nu_dyn=(vp**2-2*vs**2)/2/(vp**2-vs**2)
E=E_dyn/1.6
nu=nu_dyn
K=np.array([9E-19,5E-14,9E-19,5E-14,9E-19])/1.5E-4
#solid viscosity; power law
visc = 1E25; r = 1.
B=0.9 #Biot coef
phi=.15 #porosity
cf=1E9 #fluid bulk modulus 
source = -1.75E5 #source/sink -1.75E5, -1.8E5, -1.6E5 to make ~2MPa drop

if init==1:
    mat = [[E[0], nu[0], visc, r, rho[0], K[0], B, phi, cf,    0.,E_dyn[0],nu_dyn[0]],
           [E[1], nu[1], visc, r, rho[1], K[1], B, phi, cf,source,E_dyn[1],nu_dyn[1]],
           [E[2], nu[2], visc, r, rho[2], K[2], B, phi, cf,    0.,E_dyn[2],nu_dyn[2]],
           [E[3], nu[3], visc, r, rho[3], K[3], B, phi, cf,    0.,E_dyn[3],nu_dyn[3]],
           [E[4], nu[4], visc, r, rho[4], K[4], B, phi, cf,    0.,E_dyn[4],nu_dyn[4]]]
else:
    mat = [[E[0], nu[0], visc, r, rho[0], K[0], B, phi, cf,E_dyn[0],nu_dyn[0]],
           [E[1], nu[1], visc, r, rho[1], K[1], B, phi, cf,E_dyn[1],nu_dyn[1]],
           [E[2], nu[2], visc, r, rho[2], K[2], B, phi, cf,E_dyn[2],nu_dyn[2]],
           [E[3], nu[3], visc, r, rho[3], K[3], B, phi, cf,E_dyn[3],nu_dyn[3]],
           [E[4], nu[4], visc, r, rho[4], K[4], B, phi, cf,E_dyn[4],nu_dyn[4]]]
mat_typ = np.empty(shape = (0,1), dtype=np.uint32)
for i in nc.variables['eb_prop1'][:]:
    cnct = nc.variables["connect"+str(i)][:]
    n_elem = len(cnct)
    tet_node = np.vstack((tet_node, cnct))
    mat_typ = np.vstack((mat_typ, i*np.ones((len(cnct),1))))
print '%d nodes, %d elements' %(len(coord), len(tet_node))

# If the cell hull contains a point
def pnt_in_hull(points, new_pt):
    hull = ConvexHull(points)
    new_pts = np.vstack((points,new_pt))
    new_hull = ConvexHull(new_pts)
    if list(hull.vertices) == list(new_hull.vertices): 
        return True
    else:
        return False

# Observation locations 
ogrid = np.array([[0.0, 0, -0.1],
                  [0.1, 0, -0.1],
                  [0.2, 0, -0.1],
                  [0.3, 0, -0.1],
                  [0.4, 0, -0.1]])

# boundary data
bnd_el = []
# Traction and abs boundaries,id 1,2 reserved for fault faces
for i in nc.variables['ss_prop1'][2:]:
    els = nc.variables['elem_ss' + str(i)][:]
    sides = nc.variables['side_ss' + str(i)][:]
    bnd_el.append(np.hstack((els.reshape(len(els),1),sides.\
            reshape(len(sides),1))))
trac_el1 = bnd_el[0]
trac_el2 = bnd_el[1]
trac_el3 = bnd_el[5]
abs_bc1 = bnd_el[3] 
abs_bc2 = bnd_el[4]
abs_bc3 = bnd_el[2]

# nodal force/flux 
flux = -1E5 #-7.5E4, -5E4
if poro:
    fnodes = [[-.1, 0, -.35, 0, 0, 0, flux, dt, 15*dt],
              [ .1, 0, -.35, 0, 0, 0,    0, dt, 15*dt]]
    fnode_bc = np.empty((0,7))
    for fnode in fnodes:
        dis = np.linalg.norm(coord - np.dot(np.ones((nnd, 1)),\
                np.array(fnode[:3]).reshape(1,3)), axis=1)
        row = np.argsort(dis)[0]
        fnode_bc = np.vstack((fnode_bc, np.hstack((row + 1, fnode[3:]))))
else:
    fnodes = [[-1, 0, -2, 0, 0, 0, 0, 75E4],
              [ 1, 0, -2, 0, 0, 0, 0, 75E4]]
    fnode_bc = np.empty((0,6))
    for fnode in fnodes:
        dis = np.linalg.norm(coord - np.dot(np.ones((nnd, 1)), \
                np.array(fnode[:3]).reshape(1,3)), axis=1)
        row = np.argsort(dis)[0]
        fnode_bc = np.vstack((fnode_bc, np.hstack((row + 1, fnode[3:]))))

# fixed nodes
bcx_nodes = nc.variables['node_ns6'][:]
bcy_nodes = nc.variables['node_ns7'][:]
bcz_nodes = nc.variables['node_ns5'][:]
if poro:
    bc_typ = np.ones((nnd,4), dtype=np.int8)
else:
    bc_typ = np.ones((nnd,3), dtype=np.int8)
for node in bcx_nodes:
    bc_typ[node - 1, 0] = 0
for node in bcy_nodes:
    bc_typ[node - 1, 1] = 0   
for node in bcz_nodes:
    bc_typ[node - 1, 2] = 0

#--------------uniform BC----------------------
trac_val = [50E6, 50E6, -65E6]
#----------------------------------------------
if poro:
    trac_bc1 = np.zeros(shape=[len(trac_el1), 6])
    trac_bc2 = np.zeros(shape=[len(trac_el2), 6])
    trac_bc3 = np.zeros(shape=[len(trac_el3), 6])
else:
    trac_bc1 = np.zeros(shape=[len(trac_el1), 5])
    trac_bc2 = np.zeros(shape=[len(trac_el2), 5])
    trac_bc3 = np.zeros(shape=[len(trac_el3), 5])

if poro:
    trac_bc1[:,0] = trac_val[0]; trac_bc1[:,4] = 0.; trac_bc1[:,5] =0 
    trac_bc2[:,1] = trac_val[1]; trac_bc2[:,4] = 0.; trac_bc2[:,5] =0 
    trac_bc3[:,2] = trac_val[2]; trac_bc3[:,4] = 0.; trac_bc3[:,5] =0  
else:
    trac_bc1[:,0] = trac_val[0]; trac_bc1[:,3] = 0.; trac_bc1[:,4] = 0
    trac_bc2[:,1] = trac_val[1]; trac_bc2[:,3] = 0.; trac_bc2[:,4] = 0
    trac_bc3[:,2] = trac_val[2]; trac_bc3[:,3] = 0.; trac_bc3[:,4] = 0
trac_el = np.vstack((trac_el1, trac_el2, trac_el3))
trac_bc = np.vstack((trac_bc1, trac_bc2, trac_bc3))

# absorbing bc 
abs_bc1 = np.hstack((abs_bc1, np.ones((len(abs_bc1),1))))
abs_bc2 = np.hstack((abs_bc2, 2*np.ones((len(abs_bc2),1))))
abs_bc3 = np.hstack((abs_bc3, 3*np.ones((len(abs_bc3),1))))
abs_bc4 = np.hstack((trac_el1, np.ones((len(trac_el1),1))))
abs_bc5 = np.hstack((trac_el2, 2*np.ones((len(trac_el2),1))))
abs_bc6 = np.hstack((trac_el3, 3*np.ones((len(trac_el3),1))))
abs_bc = np.vstack((abs_bc1, abs_bc2, abs_bc3, abs_bc4, abs_bc5, abs_bc6))

# Fault nodes
ft_pos_nodes = nc.variables['node_ns1'][:]
ft_neg_nodes = nc.variables['node_ns2'][:]
# fault constraint coefficient
ft_pos_nodes = nc.variables['node_ns1'][:]
ft_neg_nodes = nc.variables['node_ns2'][:]
vec_fn = np.zeros((len(ft_pos_nodes), 3), dtype=float)
vec_fs = np.zeros((len(ft_pos_nodes), 3), dtype=float)
vec_fd = np.zeros((len(ft_pos_nodes), 3), dtype=float)
for el, side in zip(nc.variables['elem_ss1'][:], nc.variables['side_ss1'][:]):
    el_node = tet_node[el-1,:]
    if   side == 1: t_node = el_node[[0,1,3,2]]; 
    elif side == 2: t_node = el_node[[1,2,3,0]]; 
    elif side == 3: t_node = el_node[[0,2,3,1]]; 
    elif side == 4: t_node = el_node[[0,1,2,3]]; 
    row = np.empty((3,), dtype=np.uint32)
    for node, i in zip(t_node[:3], [0, 1, 2]): row[i] = np.where(ft_pos_nodes == node)[0][0]
    v = coord[t_node-1,:]
    vec = np.cross(v[1] - v[0], v[2] - v[0])
    vecpos = v[3]-(v[0]+v[1]+v[2])/3.
    vec=vec*np.sign(np.dot(vec,vecpos))
    vec_fn[row[0],:] += vec
    vec_fn[row[1],:] += vec
    vec_fn[row[2],:] += vec
    vecs = np.cross(vec, [0, 0, 1])
    vec_fs[row[0],:] += vecs
    vec_fs[row[1],:] += vecs
    vec_fs[row[2],:] += vecs
    vec = np.cross(vec, vecs) 
    vec_fd[row[0],:] += vec
    vec_fd[row[1],:] += vec
    vec_fd[row[2],:] += vec
vec_fs /= (np.ones((3,1))*np.linalg.norm(vec_fs, axis=1)).T
vec_fd /= (np.ones((3,1))*np.linalg.norm(vec_fd, axis=1)).T
vec_fn /= (np.ones((3,1))*np.linalg.norm(vec_fn, axis=1)).T

# Nonconforming interface
vtri=np.array([[0,0],[1,0],[0,1]],dtype=np.float)
# Tri shape function for element coordinate
def ShapeFuncTri(spoint):
    N=np.empty(4,dtype=np.float)
    eta=spoint[0];nu=spoint[1] 
    N[0]=1.-eta-nu
    N[1]=eta
    N[2]=nu
    return N

print 'Locating nonconformal nodes...' 
NCF_s2m = np.zeros((0,4),np.uint32) 
NCF_wt  = np.zeros((0,3),np.float)
hit=[]
pair=np.empty((0,2,),dtype=np.uint32)
tol=1E-4
for s in []: #[9,11,13,15]:
    mst_nodes = nc.variables['node_ns'+str(s)][:]
    slv_nodes = nc.variables['node_ns'+str(s+1)][:]
    mst_coord = coord[mst_nodes-1,:]
    slv2mst = np.zeros((0,4),np.uint32) 
    wt_mst = np.zeros((0,3),np.float)
    slv_nodes = np.setdiff1d(slv_nodes,np.array(hit))
    for slv in slv_nodes:
        dis2nd = np.linalg.norm(mst_coord - np.dot(np.ones((len(mst_nodes),1),dtype=np.float),coord[slv-1,:].reshape(1,3)),axis=1)
        mst6=np.argsort(dis2nd)[:6]
        if min(dis2nd)<tol: # Node on node
            ndtmp = mst_nodes[mst6[0]]
            hit.append(slv)
            pair=np.vstack((pair,[slv,ndtmp]))
        else:
            onedg=False
            for el, side in zip(nc.variables['elem_ss'+str(s)][:], nc.variables['side_ss'+str(s)][:]):
                el_node = tet_node[el-1,:]
                if   side == 1: t_node = el_node[[0,1,3]]; off_node = el_node[[2]]
                elif side == 2: t_node = el_node[[1,2,3]]; off_node = el_node[[0]]
                elif side == 3: t_node = el_node[[0,2,3]]; off_node = el_node[[1]]
                elif side == 4: t_node = el_node[[0,1,2]]; off_node = el_node[[3]]
                if len(np.intersect1d(mst_nodes[mst6],t_node))==4:
                    for j,k in ([0,1],[1,2],[2,0]):
                        v2j=coord[t_node[j]-1,:]-coord[slv-1,:]
                        v2k=coord[t_node[k]-1,:]-coord[slv-1,:]
                        onedg=np.linalg.norm(np.cross(v2j,v2k))<tol and np.dot(v2j,v2k)<0 # Node on edge
                        if onedg: 
                            dis2j=np.linalg.norm(v2j)
                            dis2k=np.linalg.norm(v2k)
                            wj = dis2k/(dis2j+dis2k); wk = dis2j/(dis2j+dis2k)
                            slv2mst= np.vstack((slv2mst,np.hstack((slv,tet_node[j],tet_node[k],0,0))))
                            wt_mst=np.vstack((wt_mst,[wj,wk,0.,0.]))
                            hit.append(slv)
                            break
                    if onedg: 
                        break
                    else: # Node on face
                        tcoord = coord[tet_node-1,:]; scoord=coord[slv-1,:]; off_coord = coord[off_node-1,:]
                        # Move the slave toward master by a gap tolerance
                        vectmp = np.cross(tcoord[1,:] - tcoord[0,:], tcoord[2,:] - tcoord[0,:])
                        vectmp = np.sign(np.dot(off_coord[0,:]-tcoord[0,:],vectmp))*vectmp/np.linalg.norm(vectmp)
                        #scoord=scoord+vectmp*tol
                        if scoord[0]<=max(tcoord[:,0]) and scoord[0]>=min(tcoord[:,0]) and\
                           scoord[1]<=max(tcoord[:,1]) and scoord[1]>=min(tcoord[:,1]) and\
                           scoord[2]<=max(tcoord[:,2]) and scoord[2]>=min(tcoord[:,2]):
                            tetcoord= coord[el_node-1,:]
                            if pnt_in_hull(tetcoord,scoord): 
                                spoint=np.zeros(3,dtype=np.float)
                                for k in range(3): 
                                    spoint[k]=griddata(tetcoord,vtri[:,k],scoord.reshape(1,3),method='linear')
                                spoint=spoint[np.argsort(abs(spoint))[:2]]
                                wt = ShapeFuncTri(spoint)
                                dst2nd = np.linalg.norm(tcoord - np.dot(np.ones((3,1),dtype=np.float),scoord.reshape(1,3)),axis=1)
                                idtmp=np.argsort(dst2nd)
                                t_node=t_node[idtmp]
                                wt=sorted(wt,reverse=True) 
                                slv2mst=np.vstack((slv2mst, np.hstack((slv,t_node))))
                                wt_mst =np.vstack((wt_mst, wt))
                                hit.append(slv)
                                break 
    NCF_s2m=np.vstack((NCF_s2m,slv2mst))
    NCF_wt=np.vstack((NCF_wt,wt_mst))
print '%d nonconformal nodes located' %(len(NCF_s2m))

# All mesh data reading stay above.
nc.close()
# Eliminate doubling nodes 
for i in range(len(pair)):
    slv=pair[i,0];mst=pair[i,1]
    tet_node[tet_node==slv]=mst
    NCF_s2m[NCF_s2m==slv]=mst
    ft_neg_nodes[ft_neg_nodes==slv]=mst
    ft_pos_nodes[ft_pos_nodes==slv]=mst
    fnode_bc[fnode_bc[:,0]==slv,0]=mst
    obs_nlist[obs_nlist==slv]=mst
    coord=np.vstack((coord[:slv-1,:],coord[slv:,:]))
    bc_typ=np.vstack((bc_typ[:slv-1,:],bc_typ[slv:,:]))
    id_shift=tet_node>slv
    tet_node[id_shift]=tet_node[id_shift]-1
    id_shift=NCF_s2m>slv
    NCF_s2m[id_shift]=NCF_s2m[id_shift]-1
    id_shift=ft_pos_nodes>slv
    ft_pos_nodes[id_shift]=ft_pos_nodes[id_shift]-1
    id_shift=ft_neg_nodes>slv
    ft_neg_nodes[id_shift]=ft_neg_nodes[id_shift]-1
    id_shift=obs_nlist>slv
    obs_nlist[id_shift]=obs_nlist[id_shift]-1
    id_shift=fnode_bc[:,0]>slv
    fnode_bc[id_shift,0] = fnode_bc[id_shift,0]-1
    id_shift=pair>slv
    pair[id_shift]=pair[id_shift]-1
nnd = len(coord)

# Sort the nodes on different side of the fault
vecf = np.empty(shape=(0,11))
xfnd = np.empty(shape=(0,3))
print 'Forming fault constraints...'
coord_pos = coord[ft_pos_nodes - 1,:]
coord_neg = coord[ft_neg_nodes - 1,:]
ft_map = np.array(np.array(np.all((coord_pos[:,None,:]==coord_neg[None,:,:]),axis=-1).nonzero()).T.tolist())
ft_neg_nodes = ft_neg_nodes[ft_map[:,1]]
id_tmp = np.argsort(ft_map[:,0])
ft_pos_nodes = ft_pos_nodes[id_tmp]
ft_neg_nodes = ft_neg_nodes[id_tmp]
nfnd = len(ft_pos_nodes) - len(np.intersect1d(ft_pos_nodes, ft_neg_nodes))

if fault:
    # Heterogeneous initial state and permeability on the fault
    j = 0
    st_init = np.zeros((nfnd,3))
    frc = np.empty((nfnd,1),dtype=np.uint32)
    fc  = np.empty((nfnd,1),dtype=np.float)
    fcd = np.empty((nfnd,1),dtype=np.float)
    dc  = np.empty((nfnd,1),dtype=np.float)
    coh = np.zeros((nfnd,1),dtype=np.float)
    dcoh = np.ones((nfnd,1),dtype=np.float)
    perm = np.zeros((nfnd,1),dtype=np.uint32)
    for node_pos, node_neg, i in zip(ft_pos_nodes, ft_neg_nodes,\
        range(len(ft_pos_nodes))):
        if node_pos != node_neg:
            x = coord[node_pos-1,0]
            y = coord[node_pos-1,1]
            z = coord[node_pos-1,2]
            if y < -.5 or y > .5:
                fc[j] = .6
                fcd[j] = .3
                dc[j] = .01
                perm[j] = 1
                frc[j] = 0 
            else:
                fc[j] = .6
                fcd[j] = .3
                dc[j] = .01
                perm[j] = 1
                frc[j] = 1
            j = j+1

# Total length of constraint function
if poro:
    neqNCF=(dim+1)*len(NCF_s2m)
    neqFT=dim*nfnd + sum(perm)
else:
    neqNCF=dim*len(NCF_s2m)
    neqFT=dim*nfnd
neq = neqNCF+neqFT 
print '%d NCF and %d fault constraint equations.' %(neqNCF,neqFT)

# Export to Defmod .inp file
fout = fin.rsplit('.')[0] + '.inp'
print 'Writing to ' + fout + '...' 
if os.path.isfile(fout): os.remove(fout)
f = open(fout, 'a')
if explicit:
    line2 = np.array([len(tet_node), nnd, len(mat), neq,\
            len(fnode_bc), len(trac_el), len(abs_bc)]).reshape(1,7)
    np.savetxt(f, line1, fmt='%s')
    np.savetxt(f, line2, delimiter=' ', fmt='%d %d %d %d %d %d %d')
    np.savetxt(f, line3, delimiter=' ', fmt='%.4f %.4f %d %d')
    np.savetxt(f, line4, delimiter=' ', fmt='%g %g %g')
elif fault:
    np.savetxt(f, line1, fmt='%s')
    line2 = np.array([len(tet_node), nnd, len(mat), neq,\
           len(fnode_bc), len(trac_el), len(abs_bc), nfnd, len(ogrid), neqNCF]).reshape(1,10)
    np.savetxt(f, line2, delimiter=' ', fmt='%d '*10)
    np.savetxt(f, line3, delimiter=' ', fmt='%.4f %.4f %d %d')
    np.savetxt(f, line4, delimiter=' ', fmt='%g %g %d %g %d %d %d %d %d %d')
    np.savetxt(f, line5, delimiter=' ', fmt='%d %d')
    np.savetxt(f, line6, delimiter=' ', fmt='%g %g %g')
else:
    line2 = np.array([len(tet_node), nnd, len(mat), neq,\
            len(fnode_bc), len(trac_el), 0, nfnd]).reshape(1,8)
    np.savetxt(f, line1, fmt='%s')
    np.savetxt(f, line2, delimiter=' ', fmt='%d %d %d %d %d %d %d %d')
    np.savetxt(f, line3, delimiter=' ', fmt='%.1f %.1f %d %d')
np.savetxt(f, np.column_stack((tet_node, mat_typ)), delimiter=' ',\
           fmt='%d %d %d %d %d')
if poro:
    np.savetxt(f, np.column_stack((coord, bc_typ)) , delimiter = ' ',\
               fmt='%g %g %g %d %d %d %d')
else:
    np.savetxt(f, np.column_stack((coord, bc_typ)) , delimiter = ' ',\
               fmt='%g %g %g %d %d %d')
if init==1:
    np.savetxt(f, mat, delimiter=' ', fmt = '%g '*12)
else:
    np.savetxt(f, mat, delimiter=' ', fmt = '%g '*11) 

# Write NCF constraints
for s2m, wt in zip(NCF_s2m,NCF_wt):
    npt = np.count_nonzero(s2m)
    if npt==3:
        if explicit or not poro:
            vec1  = [[   -1, 0, 0, s2m[0]], 
                     [wt[0], 0, 0, s2m[1]],
                     [wt[1], 0, 0, s2m[2]]]
            vec2  = [[0,    -1, 0, s2m[0]], 
                     [0, wt[0], 0, s2m[1]],
                     [0, wt[1], 0, s2m[2]]]
            vec3  = [[0, 0,    -1, s2m[0]], 
                     [0, 0, wt[0], s2m[1]],
                     [0, 0, wt[1], s2m[2]]]
        elif (fault or implicit) and poro:
            vec1  = [[   -1, 0, 0, 0, s2m[0]], 
                     [wt[0], 0, 0, 0, s2m[1]],
                     [wt[1], 0, 0, 0, s2m[2]]]
            vec2  = [[0,    -1, 0, 0, s2m[0]], 
                     [0, wt[0], 0, 0, s2m[1]],
                     [0, wt[1], 0, 0, s2m[2]]]
            vec3  = [[0, 0,    -1, 0, s2m[0]], 
                     [0, 0, wt[0], 0, s2m[1]],
                     [0, 0, wt[1], 0, s2m[2]]]
            vec4  = [[0, 0, 0,    -1, s2m[0]], 
                     [0, 0, 0, wt[0], s2m[1]],
                     [0, 0, 0, wt[1], s2m[2]]]
    elif npt==4:
        if explicit or not poro:
            vec1  = [[   -1, 0, 0, s2m[0]], 
                     [wt[0], 0, 0, s2m[1]],
                     [wt[1], 0, 0, s2m[2]],
                     [wt[2], 0, 0, s2m[3]]]
            vec2  = [[0,    -1, 0, s2m[0]], 
                     [0, wt[0], 0, s2m[1]],
                     [0, wt[1], 0, s2m[2]],
                     [0, wt[2], 0, s2m[3]]]
            vec3  = [[0, 0,    -1, s2m[0]], 
                     [0, 0, wt[0], s2m[1]],
                     [0, 0, wt[1], s2m[2]],
                     [0, 0, wt[2], s2m[3]]]
        elif (fault or implicit) and poro:
            vec1  = [[   -1, 0, 0, 0, s2m[0]], 
                     [wt[0], 0, 0, 0, s2m[1]],
                     [wt[1], 0, 0, 0, s2m[2]],
                     [wt[2], 0, 0, 0, s2m[3]]]
            vec2  = [[0,    -1, 0, 0, s2m[0]], 
                     [0, wt[0], 0, 0, s2m[1]],
                     [0, wt[1], 0, 0, s2m[2]],
                     [0, wt[2], 0, 0, s2m[3]]]
            vec3  = [[0, 0,    -1, 0, s2m[0]], 
                     [0, 0, wt[0], 0, s2m[1]],
                     [0, 0, wt[1], 0, s2m[2]],
                     [0, 0, wt[2], 0, s2m[3]]]
            vec4  = [[0, 0, 0,    -1, s2m[0]], 
                     [0, 0, 0, wt[0], s2m[1]],
                     [0, 0, 0, wt[1], s2m[2]],
                     [0, 0, 0, wt[2], s2m[3]]]
    if (explicit or not poro):
        np.savetxt(f, [npt], fmt = '%d') 
        np.savetxt(f, vec1, delimiter = ' ', fmt = '%g %g %g %d') 
        np.savetxt(f, [[0.,0.,0.]], delimiter = ' ', fmt = "%1.2E %g %g") 
        np.savetxt(f, [npt], fmt = '%d')
        np.savetxt(f, vec2, delimiter = ' ', fmt = '%g %g %g %d')
        np.savetxt(f, [[0.,0.,0.]], delimiter = ' ', fmt = "%1.2E %g %g")
        np.savetxt(f, [npt], fmt = '%d')
        np.savetxt(f, vec3, delimiter = ' ', fmt = '%g %g %g %d')
        np.savetxt(f, [[0.,0.,0.]], delimiter = ' ', fmt = "%1.2E %g %g")
    elif (fault or implicit) and poro:
        np.savetxt(f, [npt], fmt = '%d') 
        np.savetxt(f, vec1, delimiter = ' ', fmt = '%g %g %g %g %d') 
        np.savetxt(f, [[0.,0.,0.]], delimiter = ' ', fmt = "%1.2E %g %g") 
        np.savetxt(f, [npt], fmt = '%d')
        np.savetxt(f, vec2, delimiter = ' ', fmt = '%g %g %g %g %d')
        np.savetxt(f, [[0.,0.,0.]], delimiter = ' ', fmt = "%1.2E %g %g")
        np.savetxt(f, [npt], fmt = '%d')
        np.savetxt(f, vec3, delimiter = ' ', fmt = '%g %g %g %g %d')
        np.savetxt(f, [[0.,0.,0.]], delimiter = ' ', fmt = "%1.2E %g %g")
        np.savetxt(f, [npt], fmt = '%d')
        np.savetxt(f, vec4, delimiter = ' ', fmt = '%g %g %g %g %d')
        np.savetxt(f, [[0.,0.,0.]], delimiter = ' ', fmt = "%1.2E %g %g")

# fault slip: strike, dip and open, zero for non explicit model
if not explicit:
    slip = np.array([0.0, 0.0, 0.0]).reshape(3,1) 
    c = 0; dt_slip=0; t_rup=0
n = [2]
ft_neg_nodes_tap = []
j = 0
# Write fault constraint
for node_pos, node_neg, i in zip(ft_pos_nodes, ft_neg_nodes,\
    range(len(ft_pos_nodes))):
    if node_pos != node_neg:
       ft_neg_nodes_tap.append(node_neg)
       if explicit or not poro:
           vec1  = [[1, 0, 0, node_pos], 
                    [-1, 0, 0, node_neg]]
           vec2  = [[0, 1, 0, node_pos], 
                    [0, -1, 0, node_neg]]
           vec3  = [[0, 0, 1,  node_pos], 
                    [0, 0, -1, node_neg]]
       elif (fault or implicit) and poro:
           vec1  = [[1, 0, 0, 0, node_pos], 
                    [-1, 0, 0, 0, node_neg]]
           vec2  = [[0, 1, 0, 0, node_pos], 
                    [0, -1, 0, 0, node_neg]]
           vec3  = [[0, 0, 1, 0, node_pos], 
                    [0, 0, -1, 0, node_neg]]
       mat_ft = np.hstack((vec_fs[i,:].reshape(3,1),\
                vec_fd[i,:].reshape(3,1), vec_fn[i,:].reshape(3,1)))
       mat_f = np.matrix.transpose(mat_ft).reshape(1,9)
       val = np.dot(mat_ft,slip)
       y, z = np.array(coord[node_pos - 1,:][[1,2]])
       if explicit: 
          c = np.sqrt(1 - 0.36*((z+2.)**2 + y**2))
       t_act =  dt_slip+(1-c)*t_rup 
       t_slip = [t_act-dt_slip/2, t_act+dt_slip/2]
       cval1 = np.hstack((c*val[0], t_slip)).reshape(1,3)
       cval2 = np.hstack((c*val[1], t_slip)).reshape(1,3)
       cval3 = np.hstack((c*val[2], t_slip)).reshape(1,3)
       if explicit or not poro:
           np.savetxt(f, n, fmt = '%d') 
           np.savetxt(f, vec1, delimiter = ' ', fmt = '%g %g %g %d') 
           np.savetxt(f, cval1, delimiter = ' ', fmt = "%1.2E %g %g") 
           np.savetxt(f, n, fmt = '%d')
           np.savetxt(f, vec2, delimiter = ' ', fmt = '%g %g %g %d')
           np.savetxt(f, cval2, delimiter = ' ', fmt = "%1.2E %g %g")
           np.savetxt(f, n, fmt = '%d')
           np.savetxt(f, vec3, delimiter = ' ', fmt = '%g %g %g %d')
           np.savetxt(f, cval3, delimiter = ' ', fmt = "%1.2E %g %g")
       elif (fault or implicit) and poro:
           np.savetxt(f, n, fmt = '%d') 
           np.savetxt(f, vec1, delimiter = ' ', fmt = '%g %g %g %g %d') 
           np.savetxt(f, cval1, delimiter = ' ', fmt = "%1.2E %g %g") 
           np.savetxt(f, n, fmt = '%d')
           np.savetxt(f, vec2, delimiter = ' ', fmt = '%g %g %g %g %d')
           np.savetxt(f, cval2, delimiter = ' ', fmt = "%1.2E %g %g")
           np.savetxt(f, n, fmt = '%d')
           np.savetxt(f, vec3, delimiter = ' ', fmt = '%g %g %g %g %d')
           np.savetxt(f, cval3, delimiter = ' ', fmt = "%1.2E %g %g")
       vecf = np.vstack((vecf,np.hstack(([[node_pos, node_neg]], mat_f))))
       xfnd = np.vstack((xfnd, coord[node_pos-1,:])) 
       # permeable fault 
       if poro and perm[j] > 0 and not explicit: 
           vec4 = [[0, 0, 0, 1, node_pos], 
                   [0, 0, 0, -1, node_neg]]
           cval4 = [[0, 0, 0]]
           np.savetxt(f, n, fmt = '%d')
           np.savetxt(f, vec4, delimiter = ' ', fmt = '%g %g %g %g %d')
           np.savetxt(f, cval4, delimiter = ' ', fmt = "%1.2E %g %g")
       j = j + 1

if fault:
    # Write fault orientation tensor + frictional parameters
    if poro:
        np.savetxt(f, np.hstack((vecf, fc, fcd, dc, perm, st_init, xfnd,frc,coh,dcoh)), delimiter = ' ',\
           fmt = '%d '*2 + '%g '*9 + '%g '*3 + '%d ' + '%g '*3 + '%g '*3 + '%d '+ '%g '*2)
    else:
        np.savetxt(f, np.hstack((vecf, fc, fcd, dc, st_init, xfnd,frc,coh,dcoh)), delimiter = ' ',\
           fmt = '%d '*2 + '%g '*9 + '%g '*3 +         '%g '*3 + '%g '*3 + '%d '+ '%g '*2)
# Form rotated constraint matrix for fault model
if fault:
    for node_pos, node_neg, i in zip(ft_pos_nodes, ft_neg_nodes,\
        range(len(ft_pos_nodes))):
        if node_pos != node_neg:
            ft_neg_nodes_tap.append(node_neg)
            mat_ft = np.hstack((vec_fs[i,:].reshape(3,1), vec_fd[i,:].\
                    reshape(3,1), vec_fn[i,:].reshape(3,1)))
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

# Point force/source and boundary traction/flux
if poro:
    np.savetxt(f, fnode_bc, delimiter=' ',\
            fmt ='%d %g %g %g %g %g %g')
    np.savetxt(f, np.column_stack((trac_el, trac_bc)), delimiter=' ',\
            fmt ='%d %d %1.2E %1.2E %1.2E %1.2E %g %g')
else:
    np.savetxt(f, fnode_bc, delimiter=' ', fmt ='%d %1.2E %1.2E %1.2E %g %g')
    np.savetxt(f, np.column_stack((trac_el, trac_bc)), delimiter=' ',\
            fmt ='%d %d %1.2E %1.2E %1.2E %g %g') 
# Observation grid
if len(ogrid)>0:
    np.savetxt(f, ogrid, delimiter = ' ', fmt='%g '*3)
# Absorbing boundaries 
if (explicit or fault): np.savetxt(f, abs_bc, delimiter=' ', fmt='%d %d %d')
f.close(); 
print 'Defmod file ' + fout + ' created'
