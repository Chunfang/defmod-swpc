#!/usr/bin/env python
import numpy as np
import cubit, os, sys
cubit.init([''])
cubit.cmd("reset")
cubit.cmd("playback 'Mandel3D_msh.jou'")
cubit.cmd("reset")
cubit.cmd("import mesh geometry './Mandel3D.exo' block all use nodeset sideset")
surf_list = cubit.get_sideset_surfaces(1)
quad_list = [cubit.get_surface_quads(surf_list[i]) for i in range(len(surf_list))]
vol_list = cubit.get_block_volumes(1)
hex_list = [cubit.get_volume_hexes(vol_list[i]) for i in range(len(vol_list))] 

trac_el = []
mat_typ = np.ones(shape = (cubit.get_hex_count(),1), dtype=np.uint8)
for quadset, hexset in zip(quad_list, hex_list):
    for quad in quadset:
        nodes = cubit.get_connectivity("quad", quad)
        for hx in hexset:
            nodes_hx = cubit.get_connectivity("hex", hx) 
            if set(nodes).issubset(set(nodes_hx)):
                if set(nodes).issubset(set(np.array(nodes_hx)[[0,1,5,4]])):
                    side = 1
                elif set(nodes).issubset(set(np.array(nodes_hx)[[1,2,6,6]])):
                    side = 2
                elif set(nodes).issubset(set(np.array(nodes_hx)[[2,3,7,7]])):
                    side = 3
                elif set(nodes).issubset(set(np.array(nodes_hx)[[3,2,6,7]])):
                    side = 4
                elif set(nodes).issubset(set(np.array(nodes_hx)[[0,1,2,3]])):
                    side = 5 
                elif set(nodes).issubset(set(np.array(nodes_hx)[[4,5,6,7]])):
                    side = 6 
                trac_el.append([hx, side])
                mat_typ[hx - 1] = 2

mat = [[3.0E10, 0.25, 1.0E25, 1.0, 3000, 1.0E-12, 0.9664429530201342, 0.0, 2.2E9],
       [3.0E10, 0.25, 1.0E25, 1.0, 3000, 1.0E-12, 0.9664429530201342, 0.0, 2.2E9]]

trac_bc = np.zeros(shape=[len(trac_el),6])
trac_bc[:,2] = -1E6
coord = np.empty(shape=[0, 3], dtype=np.float16)
hx_node = np.empty(shape=[0, 8], dtype=np.uint16)
for node in range(cubit.get_node_count()):
    coord = np.vstack((coord, cubit.get_nodal_coordinates(node + 1)))

bcy_nodes = cubit.get_nodeset_nodes_inclusive(1)
bcx_nodes = cubit.get_nodeset_nodes_inclusive(2)
bcz_nodes = cubit.get_nodeset_nodes_inclusive(3)
bcp_nodes = cubit.get_nodeset_nodes_inclusive(4)

bc_typ = np.ones((cubit.get_node_count(),4), dtype=np.int8)
for node in bcx_nodes:
    bc_typ[node - 1, 0] = 0
for node in bcy_nodes:
    bc_typ[node - 1, 1] = 0   
for node in bcz_nodes:
    bc_typ[node - 1, 2] = 0
for node in bcp_nodes:
    bc_typ[node - 1, 3] = 0
for hx in range(cubit.get_hex_count()):
    hx_node = np.vstack((hx_node, cubit.get_connectivity("hex", hx + 1))) 

line1 = ["implicit-p hex 51"]
line2 = np.array([cubit.get_hex_count(), cubit.get_node_count(), len(mat), 0, 0, len(trac_el), 0, 0]).reshape(1,8) 
line3 = np.array([120000000.0, 80000.0, 10, 1]).reshape(1,4)
fname = './Mandel3D.inp'
if os.path.isfile(fname): os.remove(fname)
f = open(fname, 'a')
np.savetxt(f, line1, fmt='%s')
np.savetxt(f, line2, delimiter=' ', fmt='%d %d %d %d %d %d %d %d')
np.savetxt(f, line3, delimiter=' ', fmt='%.1f %.1f %d %d')
np.savetxt(f, np.column_stack((hx_node, mat_typ)), delimiter=' ', fmt='%d %d %d %d %d %d %d %d %d')
np.savetxt(f, np.column_stack((coord, bc_typ)) , delimiter = ' ', fmt='%g %g %g %d %d %d %d')
np.savetxt(f, mat, delimiter=' ', fmt = '%1.2E %g %1.2E %g %g %1.2E %g %g %1.2E') 
np.savetxt(f, np.column_stack((trac_el, trac_bc)), delimiter=' ', fmt ='%d %d %1.2E %1.2E %1.2E %1.2E %g %g') 
f.close()
