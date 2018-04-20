#!/usr/bin/env python
import numpy as np
import os, sys, netCDF4
fin = sys.argv[1]
nc = netCDF4.Dataset(fin)

# Define model type
explicit = False  
implicit = False
fault = True 
# One and only one of aboves must be True
poro = True 
visc = False
grad_BC = False 
dim = 3
if explicit: 
    t = 7.5; dt = 0.001; nviz = 25; dsp=0
    max_slip = np.array([0., .1, 0.]).reshape(3,1)
    dt_slip = .02
    slip = max_slip/(dt_slip/dt)
    t_rup = .6 
    alpha= 0.; beta = 0.0125; rfrac = 0
    line1 = ["explicit hex 51"]
    line3 = np.array([t, dt, nviz, dsp]).reshape(1,4)
    line4 = np.array([alpha, beta, rfrac]).reshape(1,3)
elif fault: 
    # dsp_str=1 output fault stress; bod_frc=1 body force; hyb=1 engage hybrid
    t = 30*24*3600.; dt = 24*3600.; nviz=1
    t_dyn=0.1; dt_dyn=0.002; t_lim=6.; dsp=1; dsp_hyb=0; dsp_str=1; rsf=0
    bod_frc=0; hyb=1; nviz_dyn=20; nviz_wave=1; nviz_slip=1; init=0
    alpha=0.; beta=0.0025; rfrac=0
    if  poro and visc:
        line1 = ["fault-pv hex 51"]
    elif poro and not visc:
        line1 = ["fault-p hex 51"]
    elif not poro and visc:
        line1 = ["fault-v hex 51"]
    else:
        line1 = ["fault hex 51"]
    line3 = np.array([t, dt, nviz, dsp]).reshape(1,4)
    line4 = np.array([t_dyn, dt_dyn, nviz_dyn, t_lim, dsp_hyb, dsp_str, bod_frc\
            , hyb,rsf,init]).reshape(1,10)
    line5 = np.array([nviz_wave,nviz_slip]).reshape(1,2)
    line6 = np.array([alpha, beta, rfrac]).reshape(1,3)
else:
    t = 1500000.; dt = 10000.; nviz = 5
    if poro and visc:
        line1 =["implicit-pv hex 51"]
    elif poro and not visc:
        line1 =["implicit-p hex 51"]
    elif not poro and visc:
        line1 =["implicit-v hex 51"]
    else:
        line1 =["implicit hex 51"]
    line3 = np.array([t, dt, nviz, 1]).reshape(1,4)

# node data
print 'Extracting mesh...'
coord = np.hstack((nc.variables['coordx'][:].\
        reshape(len(nc.variables['coordx']),1),
        nc.variables['coordy'][:].\
        reshape(len(nc.variables['coordy']),1),
        nc.variables['coordz'][:].\
        reshape(len(nc.variables['coordz']),1)))
nnd = len(coord)
hx_node = np.empty(shape=[0, 8], dtype=np.uint32)

# hex data and type
vp=6000.;vs=3464.;rho=2670
E=rho*vs**2*(3*vp**2-4*vs**2)/(vp**2-vs**2)
nu=(vp**2-2*vs**2)/2/(vp**2-vs**2)
visc=1E25 # solid viscosity 
K0=1E-15;K1=1E-11 # fluid mobility
r=1.0 # power law 
B=0.9 # Biot's
phi0=.15;phi1=.25# porosity
cf=2.2E9 # fluid compressibility
if poro:
    mat = [[E, nu, visc, r, rho, K1, B, phi1, cf,E,nu],
           [E, nu, visc, r, rho, K1, B, phi1, cf,E,nu],
           [E, nu, visc, r, rho, K1, B, phi1, cf,E,nu]]
else:
    mat = [[E, nu, visc, r, rho,E,nu],
           [E, nu, visc, r, rho,E,nu],
           [E, nu, visc, r, rho,E,nu]]
mat_typ = np.empty(shape = (0,1), dtype=np.uint8)

n_elem = []
for i in nc.variables['eb_prop1'][:]:
    cnct = nc.variables["connect"+str(i)][:]
    n_elem.append(len(cnct))
    hx_node = np.vstack((hx_node, cnct))
    mat_typ = np.vstack((mat_typ, i*np.ones((len(cnct),1))))
print '%d nodes, %d elements' %(nnd, len(hx_node))

# Observation locations 
ogrid = np.array([[-0.0, 0., -1.],
                  [-0.5, 0., -1.],
                  [-1.0, 0., -1.],
                  [-1.5, 0., -1.],
                  [-2.0, 0., -1.],
                  [-2.5, 0., -1.]])

# boundary data
bnd_el = []
# Traction and abs boundaries,id 1,2 reserved for fault faces
for i in nc.variables['ss_prop1'][2:]:
    els = nc.variables['elem_ss' + str(i)][:]
    sides = nc.variables['side_ss' + str(i)][:]
    bnd_el.append(np.hstack((els.reshape(len(els),1),sides.\
            reshape(len(sides),1))))
trac_el1 = bnd_el[3]
trac_el2 = bnd_el[4]
trac_el3 = bnd_el[5]
abs_bc1 = bnd_el[0] 
abs_bc2 = bnd_el[1]
abs_bc3 = bnd_el[2]

# Bc type space 1 free, 0 fix, -1 Winkler, 2 sync to FV
if poro:
    bc_typ = np.ones((nnd,4), dtype=np.int8)
else:
    bc_typ = np.ones((nnd,3), dtype=np.int8)
#bd = np.zeros((nnd,1) ,dtype=np.int8)
# 2nodal force/flux 
if poro:
    #fnodes = [[-.05, 0, -2, 0, 0, 0, 0, 2*dt, 6*dt],# 6E3 to trigger
    #          [ .5, 0, -2, 0, 0, 0, 0, 2*dt, 7*dt]] # 1E4 to trigger
    fnodes = [[-.1, 0, -3, 0, 0, 0, 0, 2*dt, 6*dt]] 
    fnode_bc = np.empty((0,7))
    for fnode in fnodes:
        dis = np.linalg.norm(coord - np.dot(np.ones((nnd, 1)),\
                np.array(fnode[:3]).reshape(1,3)), axis=1)
        row = np.argsort(dis)[0]
        bc_typ[dis<0.5,3]=2
        fnode_bc = np.vstack((fnode_bc, np.hstack((row + 1, fnode[3:]))))
else:
    fnodes = [[-1, 0, -2, 0, 0, 0, 0, 75E4],
             [1, 0, -2, 0, 0, 0, 0, 75E4]]
    fnode_bc = np.empty((0,6))
    for fnode in fnodes:
        dis = np.linalg.norm(coord - np.dot(np.ones((nnd, 1)), \
                np.array(fnode[:3]).reshape(1,3)), axis=1)
        row = np.argsort(dis)[0]
        fnode_bc = np.vstack((fnode_bc, np.hstack((row + 1, fnode[3:]))))

# fixed nodes
bcx_nodes = nc.variables['node_ns3'][:] 
bcy_nodes = nc.variables['node_ns4'][:]
bcz_nodes = nc.variables['node_ns5'][:]

for node in bcx_nodes:
    bc_typ[node - 1, 0] = 0
    bc_typ[node - 1,3] = 2 
for node in bcy_nodes:
    bc_typ[node - 1, 1] = 0   
    bc_typ[node - 1,3] = 2 
for node in bcz_nodes:
    bc_typ[node - 1, 2] = 0
    bc_typ[node - 1,3] = 2 
bcx_nodes = nc.variables['node_ns6'][:] 
bcy_nodes = nc.variables['node_ns7'][:]
bcz_nodes = nc.variables['node_ns8'][:]
for node in bcx_nodes:
    bc_typ[node - 1,3] = 2 
for node in bcy_nodes:
    bc_typ[node - 1,3] = 2
for node in bcz_nodes:
    bc_typ[node - 1,3] = 2

# fix all side walls when have gravity
if grad_BC:
    bcx_nodes = nc.variables['node_ns6'][:]
    bcy_nodes = nc.variables['node_ns7'][:]
    for node in bcx_nodes:
        bc_typ[node - 1, 0] = 0
    for node in bcy_nodes:
        bc_typ[node - 1, 1] = 0
#--------------gradient traction BC------------
dep = 0.; rho = np.array(mat)[0,4]; g = 9.80665
# Horizontal vertical stress ratio
#Kx = 0.748; Ky=0.795
Kx = 0.; Ky=0.
#--------------uniform BC----------------------
#trac_val = [20E4, 20E4, -20E4]
trac_val = [16E6, 16E6, -42E6] #-42E6 to trigger
#----------------------------------------------
if poro:
    trac_bc1 = np.zeros(shape=[len(trac_el1), 6])
    trac_bc2 = np.zeros(shape=[len(trac_el2), 6])
    trac_bc3 = np.zeros(shape=[len(trac_el3), 6])
else:
    trac_bc1 = np.zeros(shape=[len(trac_el1), 5])
    trac_bc2 = np.zeros(shape=[len(trac_el2), 5])
    trac_bc3 = np.zeros(shape=[len(trac_el3), 5])
if grad_BC:
    # uniform vertical traction on the top
    trac_bc3[:,2] = -dep*1E3*rho*g;  
    # gradient traction on the sides
    for el, side, i in zip(nc.variables['elem_ss6'][:], nc.variables['side_ss6'][:], range(len(trac_bc1))):
        el_node = hx_node[el-1,:] 
        if   side == 1: q_node = el_node[[0,1,5,4]]; 
        elif side == 2: q_node = el_node[[1,2,6,5]]; 
        elif side == 3: q_node = el_node[[2,3,7,6]]; 
        elif side == 4: q_node = el_node[[3,0,4,7]]; 
        elif side == 5: q_node = el_node[[0,1,2,3]]; 
        elif side == 6: q_node = el_node[[4,5,6,7]]; 
        qcoord = coord[q_node-1,:]
        zcnt = np.mean(qcoord,0)[2]
        trac_bc1[i,0]=Kx*g*1E3*(dep-zcnt)*rho
    for el, side, i in zip(nc.variables['elem_ss7'][:], nc.variables['side_ss7'][:], range(len(trac_bc2))):
        el_node = hx_node[el-1,:]
        if   side == 1: q_node = el_node[[0,1,5,4]]; 
        elif side == 2: q_node = el_node[[1,2,6,5]]; 
        elif side == 3: q_node = el_node[[2,3,7,6]]; 
        elif side == 4: q_node = el_node[[3,0,4,7]]; 
        elif side == 5: q_node = el_node[[0,1,2,3]]; 
        elif side == 6: q_node = el_node[[4,5,6,7]];
        qcoord = coord[q_node-1,:]
        zcnt = np.mean(qcoord,0)[2]
        trac_bc2[i,1]=Ky*g*1E3*(dep-zcnt)*rho
    if poro:
        trac_bc3[:,2] = trac_val[2]; trac_bc3[:,4] = 0.; trac_bc3[:,5] = t
    else:
        trac_bc3[:,2] = trac_val[2]; trac_bc3[:,3] = 0.; trac_bc3[:,4] = t
    trac_el = np.vstack((trac_el3))
    trac_bc = np.vstack((trac_bc3))
else:
    if explicit: trac_val = [0., 0., 0.]
    if poro:
        trac_bc1[:,0] = trac_val[0]; trac_bc1[:,4] = 0.; trac_bc1[:,5] = 0. #t/4
        trac_bc2[:,1] = trac_val[1]; trac_bc2[:,4] = 0.; trac_bc2[:,5] = 0. #t/4
        trac_bc3[:,2] = trac_val[2]; trac_bc3[:,4] = 0.; trac_bc3[:,5] = 0. #t*2/3
    else:
        trac_bc1[:,0] = trac_val[0]; trac_bc1[:,3] = 0.; trac_bc1[:,4] = t/4
        trac_bc2[:,1] = trac_val[1]; trac_bc2[:,3] = 0.; trac_bc2[:,4] = t/4
        trac_bc3[:,2] = trac_val[2]; trac_bc3[:,3] = 0.; trac_bc3[:,4] = t*9/10
    trac_el = np.vstack((trac_el1, trac_el2, trac_el3))
    trac_bc = np.vstack((trac_bc1, trac_bc2, trac_bc3))

# absorbing bc 
abs_bc1 = np.hstack((abs_bc1, np.ones((len(abs_bc1),1))))
abs_bc2 = np.hstack((abs_bc2, 2*np.ones((len(abs_bc2),1))))
abs_bc3 = np.hstack((abs_bc3, 3*np.ones((len(abs_bc3),1))))
abs_bc4 = np.hstack((trac_el1, np.ones((len(trac_el1),1))))
abs_bc5 = np.hstack((trac_el2, 2*np.ones((len(trac_el2),1))))
abs_bc6 = np.hstack((trac_el3, 3*np.ones((len(trac_el3),1))))
abs_bc = np.vstack((abs_bc1, abs_bc2, abs_bc3, abs_bc4, abs_bc5))
# abs_bc6 upper bound is traction free not absorbing

# fault constraint coefficient
print 'Forming fault constraints...'
ft_pos_nodes = nc.variables['node_ns1'][:]
ft_neg_nodes = nc.variables['node_ns2'][:]
vec_fn = np.zeros((len(ft_pos_nodes), 3), dtype=float)
vec_fs = np.zeros((len(ft_pos_nodes), 3), dtype=float)
vec_fd = np.zeros((len(ft_pos_nodes), 3), dtype=float)
for el, side in zip(nc.variables['elem_ss1'][:],\
    nc.variables['side_ss1'][:]):
    el_node = hx_node[el-1,:]
    if   side == 1: q_node = el_node[[0,1,5,4,2,3,6,7]]; 
    elif side == 2: q_node = el_node[[1,2,6,5,0,3,4,7]]; 
    elif side == 3: q_node = el_node[[2,3,7,6,0,1,4,5]]; 
    elif side == 4: q_node = el_node[[3,0,4,7,1,2,5,6]]; 
    elif side == 5: q_node = el_node[[0,1,2,3,4,5,6,7]]; 
    elif side == 6: q_node = el_node[[4,5,6,7,0,1,2,3]]; 
    row = np.empty((4,), dtype=np.uint32)
    for node, i in zip(q_node[:4], [0, 1, 2, 3]): row[i] =\
        np.where(ft_pos_nodes == node)[0][0]
    v = coord[q_node-1,:]
    vecpos = v[4]+v[4]+v[6]+v[7]-v[0]-v[1]-v[2]-v[3]
    vec1 = np.cross(v[3] - v[0], v[1] - v[0])
    vec2 = np.cross(v[0] - v[1], v[2] - v[1])
    vec3 = np.cross(v[1] - v[2], v[3] - v[2])
    vec4 = np.cross(v[2] - v[3], v[0] - v[3])
    vec1 = vec1*np.sign(np.dot(vec1,vecpos));
    vec2 = vec1*np.sign(np.dot(vec2,vecpos));
    vec3 = vec1*np.sign(np.dot(vec3,vecpos));
    vec4 = vec1*np.sign(np.dot(vec4,vecpos));
    vec_fn[row[0],:] += vec1 
    vec_fn[row[1],:] += vec2
    vec_fn[row[2],:] += vec3 
    vec_fn[row[3],:] += vec4
    vec1s = np.cross(vec1, [0, 0, 1])
    vec2s = np.cross(vec2, [0, 0, 1])
    vec3s = np.cross(vec3, [0, 0, 1])
    vec4s = np.cross(vec4, [0, 0, 1])
    vec_fs[row[0],:] += vec1s 
    vec_fs[row[1],:] += vec2s
    vec_fs[row[2],:] += vec3s 
    vec_fs[row[3],:] += vec4s
    vec1 = np.cross(vec1, vec1s) 
    vec2 = np.cross(vec2, vec2s)
    vec3 = np.cross(vec3, vec3s)
    vec4 = np.cross(vec4, vec4s)
    vec_fd[row[0],:] = vec1 
    vec_fd[row[1],:] = vec2
    vec_fd[row[2],:] = vec3 
    vec_fd[row[3],:] = vec4
vec_fs /= (np.ones((3,1))*np.linalg.norm(vec_fs, axis=1)).T
vec_fd /= (np.ones((3,1))*np.linalg.norm(vec_fd, axis=1)).T
vec_fn /= (np.ones((3,1))*np.linalg.norm(vec_fn, axis=1)).T
nc.close()
# All mesh data reading stay above.
vecf = np.empty(shape=(0,11))
xfnd = np.empty(shape=(0,3))

# Sort the nodes on different side of the fault
coord_pos = coord[ft_pos_nodes - 1,:]
coord_neg = coord[ft_neg_nodes - 1,:]
ft_map = np.array(np.array(np.all((coord_pos[:,None,:]==\
        coord_neg[None,:,:]),axis=-1).nonzero()).T.tolist())
ft_neg_nodes = ft_neg_nodes[ft_map[:,1]]
id_tmp = np.argsort(ft_map[:,0])
ft_pos_nodes = ft_pos_nodes[id_tmp]
ft_neg_nodes = ft_neg_nodes[id_tmp]
nfnd = len(ft_pos_nodes) - len(np.intersect1d(ft_pos_nodes, ft_neg_nodes))

# Make the fault permeable at certain depth; add initial stress
perm = np.ones((nfnd,1))
st_init = np.zeros((nfnd,3))
frc = np.ones((nfnd,1))

# permeable at certain depth
y_min = min(coord[:,1]); y_max = max(coord[:,1])
z_min = min(coord[:,2]); z_max = max(coord[:,2])
j = 0
for node_pos, node_neg, i in zip(ft_pos_nodes, ft_neg_nodes,\
    range(len(ft_pos_nodes))):
    if node_pos != node_neg:
        y = coord[node_pos - 1,1]
        z = coord[node_pos - 1,2]
        if z > z_min + (z_max - z_min)*0.44 and \
           z < z_min + (z_max - z_min)*0.55 and \
           y > y_min + (y_max - y_min)*0.44 and \
           y < y_min + (y_max - y_min)*0.55:
            perm[j] = 1 
            #st_init[j,:]=[0., 2.745E7, -5E7]
            frc[j] = 1 
        else:
            perm[j] = 1 
            #st_init[j,:]=[0., 2.325E7 , -5E7]
            frc[j] = 1 
        j = j+1

# Total length of constraint function
if poro:
    neq = dim*nfnd + sum(perm)
else:
    neq = dim*nfnd
print '%d constraint equations' %(neq)                

# Export to Defmod .inp file
fout = fin.rsplit('.')[0] + '.inp'
if os.path.isfile(fout): os.remove(fout)
f = open(fout, 'a')
print 'Writing defmod input to ' + fout + '...' 
if explicit:
    line2 = np.array([len(hx_node), nnd, len(mat), neq,\
            len(fnode_bc), len(trac_el), len(abs_bc)]).reshape(1,7)
    np.savetxt(f, line1, fmt='%s')
    np.savetxt(f, line2, delimiter=' ', fmt='%d '*7)
    np.savetxt(f, line3, delimiter=' ', fmt='%g %g %d %d')
    np.savetxt(f, line4, delimiter=' ', fmt='%g %g %g')
elif fault:
    np.savetxt(f, line1, fmt='%s')
    neqNCF=0 # No nonconformal nodes
    line2 = np.array([len(hx_node), nnd, len(mat), neq,\
           len(fnode_bc), len(trac_el), len(abs_bc), nfnd, len(ogrid), neqNCF]).reshape(1,10)
    np.savetxt(f, line2, delimiter=' ', fmt='%d '*10)
    np.savetxt(f, line3, delimiter=' ', fmt='%g %g %d %d')
    np.savetxt(f, line4, delimiter=' ', fmt='%g %g %d %g ' +  '%d '*6)
    np.savetxt(f, line5, delimiter=' ', fmt='%d %d')
    np.savetxt(f, line6, delimiter=' ', fmt='%g %g %g')
else:
    line2 = np.array([len(hx_node), nnd, len(mat), neq,\
            len(fnode_bc), len(trac_el), 0, nfnd]).reshape(1,8)
    np.savetxt(f, line1, fmt='%s')
    np.savetxt(f, line2, delimiter=' ', fmt='%d '*8)
    np.savetxt(f, line3, delimiter=' ', fmt='%g '*2 + '%d %d')
np.savetxt(f, np.column_stack((hx_node, mat_typ)), delimiter=' ',\
           fmt='%d '*8 + '%d')
if poro:
    np.savetxt(f, np.column_stack((coord, bc_typ)) , delimiter = ' ',\
               fmt='%g '*3 + '%d '*4)
    np.savetxt(f, mat, delimiter=' ', fmt = '%g '*11)
else:
    np.savetxt(f, np.column_stack((coord, bc_typ)) , delimiter = ' ',\
               fmt='%g '*3 + '%d '*3)
    np.savetxt(f, mat, delimiter=' ', fmt = '%g '*7) 
# fault slip: strike, dip and open, zero for non explicit model
if not explicit:
    slip = np.array([0.0, 0.0, 0.0]).reshape(3,1) 
    c = 0; dt_slip=0; t_rup=0
n = [2]
ft_neg_nodes_tap = []
j = 0
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
           np.savetxt(f, cval1, delimiter = ' ', fmt = "%g %g %g") 
           np.savetxt(f, n, fmt = '%d')
           np.savetxt(f, vec2, delimiter = ' ', fmt = '%g %g %g %d')
           np.savetxt(f, cval2, delimiter = ' ', fmt = "%g %g %g")
           np.savetxt(f, n, fmt = '%d')
           np.savetxt(f, vec3, delimiter = ' ', fmt = '%g %g %g %d')
           np.savetxt(f, cval3, delimiter = ' ', fmt = "%g %g %g")
       elif (fault or implicit) and poro:
           np.savetxt(f, n, fmt = '%d') 
           np.savetxt(f, vec1, delimiter = ' ', fmt = '%g %g %g %g %d') 
           np.savetxt(f, cval1, delimiter = ' ', fmt = "%g %g %g") 
           np.savetxt(f, n, fmt = '%d')
           np.savetxt(f, vec2, delimiter = ' ', fmt = '%g %g %g %g %d')
           np.savetxt(f, cval2, delimiter = ' ', fmt = "%g %g %g")
           np.savetxt(f, n, fmt = '%d')
           np.savetxt(f, vec3, delimiter = ' ', fmt = '%g %g %g %g %d')
           np.savetxt(f, cval3, delimiter = ' ', fmt = "%g %g %g")
       vecf = np.vstack((vecf,np.hstack(([[node_pos, node_neg]], mat_f))))
       xfnd = np.vstack((xfnd, coord[node_pos-1,:])) 
       # permeable fault 
       if poro and perm[j] > 0 and not explicit: 
           vec4 = [[0, 0, 0, 1, node_pos], 
                   [0, 0, 0, -1, node_neg]]
           cval4 = [[0, 0, 0]]
           np.savetxt(f, n, fmt = '%d')
           np.savetxt(f, vec4, delimiter = ' ', fmt = '%g %g %g %g %d')
           np.savetxt(f, cval4, delimiter = ' ', fmt = "%g %g %g")
       j = j + 1
if fault:
    # Define frictional parameters (static friction)
    fc = .55*np.ones((len(vecf),1))
    fcd = .45*np.ones((len(vecf),1))
    dc = .01*np.ones((len(vecf),1))
    # Cohesive stress and idstance 
    coh= 1E4*np.ones((len(vecf),1))
    dcoh= .01*np.ones((len(vecf),1))
    biot = 1.*np.ones((len(vecf),1))
    # Write fault orientation tensor + frictional parameters
    if poro:
        np.savetxt(f, np.hstack((vecf, fc, fcd, dc, perm, st_init, xfnd, frc,coh,dcoh, biot)), delimiter = ' ',\
           fmt = '%d '*2 + '%g '*9 + '%g '*3 + '%d ' + '%g '*3 + '%g '*3 + '%d '+ '%g '*3)
    else:
        np.savetxt(f, np.hstack((vecf, fc, fcd, dc, st_init, xfnd, frc,coh,dcoh)), delimiter = ' ',\
           fmt = '%d '*2 + '%g '*9 + '%g '*3 +         '%g '*3 + '%g '*3 + '%d '+ '%g '*2)

# Point force/source and boundary traction/flux
if poro:
    np.savetxt(f, fnode_bc, delimiter=' ',\
            fmt ='%d ' + '%g '*4 + '%g %g')
    np.savetxt(f, np.column_stack((trac_el, trac_bc)), delimiter=' ',\
            fmt ='%d %d %1.2E %1.2E %1.2E %1.2E %g %g')
else:
    np.savetxt(f, fnode_bc, delimiter=' ', fmt ='%d %1.2E %1.2E %1.2E %g %g')
    np.savetxt(f, np.column_stack((trac_el, trac_bc)), delimiter=' ',\
            fmt ='%d %d %1.2E %1.2E %1.2E %g %g') 
# Observation grid
if len(ogrid)>0:
    np.savetxt(f, ogrid, delimiter = ' ', fmt='%g '*3)
# Absorbing boundary
if (explicit or fault): np.savetxt(f, abs_bc, delimiter=' ', fmt='%d %d %d') 
f.close(); 
print 'Defmod imput ' + fout + ' created' 

# Heal fault for flow model
for i in range(len(ft_pos_nodes)):
    slv = ft_pos_nodes[i]; mst = ft_neg_nodes[i]
    if slv != mst:
        hx_node[hx_node==slv]=mst
        coord=np.vstack((coord[:slv-1,:],coord[slv:,:]))
        id_shift=hx_node>slv
        hx_node[id_shift]=hx_node[id_shift]-1
        id_shift=ft_pos_nodes>slv
        ft_pos_nodes[id_shift]=ft_pos_nodes[id_shift]-1
        id_shift=ft_neg_nodes>slv
        ft_neg_nodes[id_shift]=ft_neg_nodes[id_shift]-1
nnd = len(coord)

# h5 file for pflotran
import h5py
col1 = 8*np.ones((len(hx_node),1),dtype=np.int32)
id_cell = np.array(range(len(hx_node)))+1
fout = fin.rsplit('.')[0] + '_h5.h5'
f = h5py.File(fout, "w")
print 'Outputting pflotran file ' + fout + '...'
grp = f.create_group('Domain')
grp.create_dataset('Cells', data=np.hstack((col1,hx_node)), dtype=np.int32)
grp.create_dataset('Vertices', data=coord*1E3, dtype=np.float64)
grp_reg = f.create_group('Regions')
grp = grp_reg.create_group('All')
grp.create_dataset('Cell Ids', data=id_cell, dtype=np.int32)
grp = grp_reg.create_group('top_layer')
grp.create_dataset('Cell Ids', data=id_cell[:n_elem[0]], dtype=np.int32)
grp = grp_reg.create_group('mid_layer')
grp.create_dataset('Cell Ids', data=id_cell[n_elem[0]:sum(n_elem[:2])], dtype=np.int32)
grp = grp_reg.create_group('bottom_layer')
grp.create_dataset('Cell Ids', data=id_cell[sum(n_elem[:2]):sum(n_elem)], dtype=np.int32)

# Absorbing sides
ss = np.empty((0,5), dtype=np.int32)
for el, side in zip(abs_bc6[:,0],abs_bc1[:,1]):
    el_node = hx_node[int(el)-1,:]
    if   side == 1: q_node = el_node[[0,1,5,4]]; 
    elif side == 2: q_node = el_node[[1,2,6,5]]; 
    elif side == 3: q_node = el_node[[2,3,7,6]]; 
    elif side == 4: q_node = el_node[[3,0,4,7]]; 
    elif side == 5: q_node = el_node[[0,1,2,3]]; 
    elif side == 6: q_node = el_node[[4,5,6,7]]; 
    ss = np.vstack((ss, np.hstack((4,q_node))));
grp = grp_reg.create_group('East_sidesets')
grp.create_dataset('Vertex Ids', data=ss, dtype=np.int32)

ss = np.empty((0,5), dtype=np.int32)
for el, side in zip(abs_bc6[:,0],abs_bc2[:,1]):
    el_node = hx_node[int(el)-1,:]
    if   side == 1: q_node = el_node[[0,1,5,4]]; 
    elif side == 2: q_node = el_node[[1,2,6,5]]; 
    elif side == 3: q_node = el_node[[2,3,7,6]]; 
    elif side == 4: q_node = el_node[[3,0,4,7]]; 
    elif side == 5: q_node = el_node[[0,1,2,3]]; 
    elif side == 6: q_node = el_node[[4,5,6,7]]; 
    ss = np.vstack((ss, np.hstack((4,q_node))));
grp = grp_reg.create_group('North_sidesets')
grp.create_dataset('Vertex Ids', data=ss, dtype=np.int32)

ss = np.empty((0,5), dtype=np.int32)
for el, side in zip(abs_bc6[:,0],abs_bc3[:,1]):
    el_node = hx_node[int(el)-1,:]
    if   side == 1: q_node = el_node[[0,1,5,4]]; 
    elif side == 2: q_node = el_node[[1,2,6,5]]; 
    elif side == 3: q_node = el_node[[2,3,7,6]]; 
    elif side == 4: q_node = el_node[[3,0,4,7]]; 
    elif side == 5: q_node = el_node[[0,1,2,3]]; 
    elif side == 6: q_node = el_node[[4,5,6,7]]; 
    ss = np.vstack((ss, np.hstack((4,q_node))));
grp = grp_reg.create_group('Bottom_sidesets')
grp.create_dataset('Vertex Ids', data=ss, dtype=np.int32)

ss = np.empty((0,5), dtype=np.int32)
for el, side in zip(abs_bc6[:,0],abs_bc4[:,1]):
    el_node = hx_node[int(el)-1,:]
    if   side == 1: q_node = el_node[[0,1,5,4]]; 
    elif side == 2: q_node = el_node[[1,2,6,5]]; 
    elif side == 3: q_node = el_node[[2,3,7,6]]; 
    elif side == 4: q_node = el_node[[3,0,4,7]]; 
    elif side == 5: q_node = el_node[[0,1,2,3]]; 
    elif side == 6: q_node = el_node[[4,5,6,7]]; 
    ss = np.vstack((ss, np.hstack((4,q_node))));
grp = grp_reg.create_group('West_sidesets')
grp.create_dataset('Vertex Ids', data=ss, dtype=np.int32)

ss = np.empty((0,5), dtype=np.int32)
for el, side in zip(abs_bc6[:,0],abs_bc5[:,1]):
    el_node = hx_node[int(el)-1,:]
    if   side == 1: q_node = el_node[[0,1,5,4]]; 
    elif side == 2: q_node = el_node[[1,2,6,5]]; 
    elif side == 3: q_node = el_node[[2,3,7,6]]; 
    elif side == 4: q_node = el_node[[3,0,4,7]]; 
    elif side == 5: q_node = el_node[[0,1,2,3]]; 
    elif side == 6: q_node = el_node[[4,5,6,7]]; 
    ss = np.vstack((ss, np.hstack((4,q_node))));
grp = grp_reg.create_group('South_sidesets')
grp.create_dataset('Vertex Ids', data=ss, dtype=np.int32)

ss = np.empty((0,5), dtype=np.int32)
for el, side in zip(abs_bc6[:,0],abs_bc6[:,1]):
    el_node = hx_node[int(el)-1,:]
    if   side == 1: q_node = el_node[[0,1,5,4]]; 
    elif side == 2: q_node = el_node[[1,2,6,5]]; 
    elif side == 3: q_node = el_node[[2,3,7,6]]; 
    elif side == 4: q_node = el_node[[3,0,4,7]]; 
    elif side == 5: q_node = el_node[[0,1,2,3]]; 
    elif side == 6: q_node = el_node[[4,5,6,7]]; 
    ss = np.vstack((ss, np.hstack((4,q_node))));
grp = grp_reg.create_group('Top_sidesets')
grp.create_dataset('Vertex Ids', data=ss, dtype=np.int32)
f.close()
print fout + ' created.'
