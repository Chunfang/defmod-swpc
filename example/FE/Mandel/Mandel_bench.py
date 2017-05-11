#!/usr/bin/env python
from cmath import *
from  numpy import ones, sinh, cosh, sqrt, pi, array, exp, vstack, hstack, shape, empty, mgrid, linalg, where, squeeze, logical_and, argsort
from scipy.special import iv
from vtk import vtkStructuredPointsReader, vtkUnstructuredGridReader 
from vtk.util import numpy_support as vtk2num 
from scipy.interpolate import griddata
import re, os, fnmatch

# load 2D defmod result
step = []; part = []
r_dfm2d = empty((0, 1), dtype=float) 
p_dfm2d = empty((0, 1), dtype=float)
for f in os.listdir('./Mandel2D'):
    if fnmatch.fnmatch(f,'*Mandel2D_*.vtk'):
        reader = vtkUnstructuredGridReader()
        reader.SetFileName('./Mandel2D/'+f)
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update()
        data = reader.GetOutput()
        disp = vtk2num.vtk_to_numpy(data.GetPointData().GetArray("displacements"))
        pres =  vtk2num.vtk_to_numpy(data.GetPointData().GetArray("pressure"))
        coord = vtk2num.vtk_to_numpy(data.GetPoints().GetData())
        crd = empty((0, 3), dtype=float) 
        prs = empty((0, 1), dtype=float)
        for i in range(len(coord)):
            if coord[i,0] > 4.9:
                crd = vstack((crd, coord[i,:]))
                prs = vstack((prs, pres[i]))
        if crd.size:
            r_dfm2d = vstack((r_dfm2d, 1. + crd[:,1].reshape(len(crd), 1))) 
            p_dfm2d = vstack((p_dfm2d, prs.reshape(len(prs), 1)))
            step.append(int(re.search("Mandel2D_(.*).vtk", f).group(1)))
r_dfm2d = r_dfm2d.reshape(len(step), len(r_dfm2d)/len(step))
p_dfm2d = p_dfm2d.reshape(len(step), len(p_dfm2d)/len(step))
idx = argsort(step)
r_dfm2d = r_dfm2d[idx,:]
p_dfm2d = p_dfm2d[idx,:]

step = []
grid_x = mgrid[0:5/sqrt(2):50j]
grid_y = grid_x
h = -9.875
grid_z = h*ones(len(grid_x))
r_dfm3d = linalg.norm(hstack((grid_x.reshape(len(grid_x),1), grid_y.reshape(len(grid_y),1))), axis=1)
p_dfm3d = [] 
for f in os.listdir('./Mandel3D'):
    if fnmatch.fnmatch(f,'*Mandel3D_*.vtk'):
        reader = vtkUnstructuredGridReader()
        reader.SetFileName('./Mandel3D/'+f)
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update()
        data = reader.GetOutput()
        coord = vtk2num.vtk_to_numpy(data.GetPoints().GetData())
        idx = where(logical_and(coord[:,2] < h + .128, coord[:,2] > h - .128))
        if len(squeeze(idx)) > 0 : 
            step.append(int(re.search("Mandel3D_(.*).vtk", f).group(1)))
            coord = squeeze(coord[idx,:])
            disp = vtk2num.vtk_to_numpy(data.GetPointData().GetArray("displacements"))
            disp = disp[idx,:]
            pres =  vtk2num.vtk_to_numpy(data.GetPointData().GetArray("pressure"))
            pres = pres[idx]
            p_dfm3d.append(griddata(coord, pres, hstack((grid_x.reshape(len(grid_x),1), grid_y.reshape(len(grid_z),1), grid_z.reshape(len(grid_z),1))), method='linear'))
p_dfm3d = vstack(p_dfm3d)
r_dfm3d = r_dfm3d/r_dfm3d[-1]
idx = argsort(step)
p_dfm3d = p_dfm3d[idx,:]
step = array(step)[idx]
def F(s):
    return 1.0/(s+1.0)
def cot(phi):
    return 1.0/tan(phi)
def csc(phi):
    return 1.0/sin(phi)
def talbot(t,N=10, F=F):
    h = 2*pi/N;
    shift = 0.0;
    ans =   0.0;
    if t == 0:
        print "ERROR:   Inverse transform can not be calculated for t=0"
        return ("Error");
    for k in range(0,N):
        theta = -pi + (k+1./2)*h;
        z = shift + N/t*(0.5017*theta*cot(0.6407*theta) - 0.6122 + 0.2645j*theta); 
        dz = N/t*(-0.5017*0.6407*theta*(csc(0.6407*theta)**2)+0.5017*cot(0.6407*theta)+0.2645j);
        ans = ans + exp(z*t)*F(z)*dz;
    return ((h/(2j*pi))*ans).real        

def p_Mandel2D(a, c, B, nu_u, nu, chi, tau):
    def D(xi):
        return (1 - nu) - (nu_u -nu)*sinh(sqrt(xi))/sqrt(xi)/cosh(sqrt(xi))
    def p_lap(xi):
        return B*2/3.*(1. + nu_u)*(1. - nu)/xi/D(xi)*(1. - cosh(sqrt(xi)*chi)/cosh(sqrt(xi)))
    return talbot(tau, F = p_lap, N = 16)

def p_Mandel3D(a, c, B, nu_u, nu, chi, tau):
    def D(xi):
        return (1. - nu)*(1. + nu_u) - 4.*(nu_u -nu)*iv(1, sqrt(xi))/sqrt(xi)/iv(0, sqrt(xi))
    def p_lap(xi):
        return B*(1. + nu_u)*(1. - nu)/xi/D(xi)*(1. - iv(0, sqrt(xi)*chi)/iv(0, sqrt(xi)))
    return talbot(tau, F = p_lap, N = 16)

# Mandel parameters
#Vs, Vp = 1600., 3000,
rhos = 3000. 
#rhof = 1E3
#G = rhos*Vs**2
#nu = (Vp**2 - 2*Vs**2)/2/(Vp**2 - Vs**2)
#E = 2*G*(1 + nu)
E = 3E10
nu = .25; nu_u = .49; B = 1. # hard code Poisson ratio
G = E/2/(1 + nu) 
mu    = E / (2.0*(1.0 + nu))
lmbda = E*nu / ((1.0 + nu)*(1.0 - 2.0*nu))
alphaB = 3*(nu_u - nu)/B/(1 + nu_u)/(1 - 2*nu) 
K = 1E-12 
#storage rate:
S = alphaB**2*(1 - 2*nu)**2*(1 - nu_u)/2/G/(nu_u - nu)/(1 - nu)
#visc = 1E-3
# diffusivity
c = K/S
# dimension
a = 1.
dt = 200000.0
t = array(step)*dt
tdl= (t*c/(a**2))[1:]*1E-6
dt = 800000.0
t = step*dt
a = 5.
tdl3D = (t*c/(a**2))[1:]*1E-6

Nx = 50
place, p2D, p3D = [], [], [] 
for tau, tau3D in zip(tdl, tdl3D):
    # tau = 1E-3*10**(5.*i/Nt)
    p2, p3, r = [], [], [] 
    for k in range(Nx+1):
        #chi = 1E-3*10**(3*k/Nx)
        chi = 1.*k/Nx
        p2.append(p_Mandel2D(a, c, B, nu_u, nu, chi, tau)) 
        p3.append(p_Mandel3D(a, c, B, nu_u, nu, chi, tau3D))
        # t.append(tau)
        r.append(chi)
    p2D.append(p2)
    p3D.append(p3)
    # time.append(t)
    place.append(r) 
# time = np.vstack(time)
place= vstack(place)
p2D =  vstack(p2D)
p3D =  vstack(p3D)
p_dfm2d = p_dfm2d/p_dfm2d[1,0]*p2D[0,0]
p_dfm3d = p_dfm3d*p3D[0,0]/p_dfm3d[1,0]

import matplotlib.pyplot as plt
plt.semilogx(tdl[:], p2D[:,0], "c-", label ="2D Mandel")
plt.semilogx(tdl[:], p_dfm2d[1:,0], "g-",label = "2D Defmod")
plt.semilogx(tdl3D[:], p3D[:,0], "b-", label ="3D Mandel")
plt.semilogx(tdl3D[:], p_dfm3d[1:,0], "m-",label = "3D Defmod")
plt.xlabel("normalized time")
plt.ylabel("normalized pressure")
plt.title(r"normalized pressure at origin $(r = 0)$")
plt.legend()
plt.show()

plt.plot(place[0,:], p2D[0,:], "c-", label ="2D Mandel")
plt.plot(r_dfm2d[0,:], p_dfm2d[1,:], "g-",label = "2D Defmod")
plt.plot(place[5,:], p2D[5,:], "c-")
plt.plot(r_dfm2d[10,:], p_dfm2d[6,:], "g-")
plt.plot(place[10,:], p2D[10,:], "c-")
plt.plot(r_dfm2d[10,:], p_dfm2d[11,:], "g-")
plt.plot(place[20,:], p2D[20,:], "c-")
plt.plot(r_dfm2d[20,:], p_dfm2d[21,:], "g-")
plt.plot(place[20,:], p2D[40,:], "c-")
plt.plot(r_dfm2d[20,:], p_dfm2d[41,:], "g-")
plt.title(r"normalized pressure at different time")
plt.xlabel(r"$r/a$")
plt.ylabel("normalized pressure")
plt.legend(loc = 3)
plt.show()

plt.plot(place[0,:], p3D[0,:], "b-", label ="3D Mandel")
plt.plot(r_dfm3d, p_dfm3d[1,:], "m-",label = "3D Defmod")
plt.plot(place[10,:], p3D[20,:], "b-")
plt.plot(r_dfm3d, p_dfm3d[21,:], "m-")
plt.plot(place[10,:], p3D[40,:], "b-")
plt.plot(r_dfm3d, p_dfm3d[41,:], "m-")
plt.plot(place[15,:], p3D[80,:], "b-")
plt.plot(r_dfm3d, p_dfm3d[81,:], "m-")
plt.plot(place[5,:], p3D[149,:], "b-")
plt.plot(r_dfm3d, p_dfm3d[150,:], "m-")
plt.xlabel(r"$r/a$")
plt.ylabel("normalized pressure")
plt.title(r"normalized pressure at different time")
plt.legend(loc = 3)
plt.show()
