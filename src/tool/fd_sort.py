#!/usr/bin/env python 
import numpy as np
import os,sys
from itertools import islice
import glob as glb
name_sol = sys.argv[1]
files_seis=glb.glob(name_sol+"_*_ofd.txt")
fopen=[]; nobs_loc=[]
for file_seis in files_seis: fopen.append(open(file_seis))
dmn=3
ocoord=np.empty(shape=[0,dmn+1],dtype=np.float)
for f in fopen:
    n,nt,dt=np.genfromtxt(islice(f,1),delimiter=" ",unpack=False,dtype=np.float)
    n=int(n); nt=int(nt)
    ocoord=np.vstack((ocoord,np.genfromtxt(islice(f,n),delimiter=" ",unpack=False,dtype=np.float)))
    nobs_loc.append(n)
oidx,opick=np.unique(ocoord[:,-1],return_index=True)
osort=np.argsort(oidx)
ocoord=ocoord[opick,:dmn]; ocoord=ocoord[osort,:dmn]
dat_seis_tmp=np.empty(shape=[0,dmn],dtype=np.float)
for i in range(nt):
    dat_tmp=np.empty(shape=[0,dmn],dtype=np.float)
    for f,j in zip(fopen,range(len(fopen))):
        dat_tmp=np.vstack((dat_tmp,np.loadtxt(islice(f,nobs_loc[j]),unpack=False,dtype=np.float)))
    dat_tmp=dat_tmp[opick,:]
    dat_seis_tmp=np.vstack((dat_seis_tmp,dat_tmp[osort,:]))
    if np.remainder(i+1,100)*np.remainder(i+1,nt)==0:print "frame " + str(i+1) +"/"+str(nt)+ " merged"
nobs=len(ocoord)
dat_seis = np.empty(shape=[nobs,nt,dmn],dtype=np.float)
for i in range(nobs): 
    dat_seis[i,:,:] = dat_seis_tmp[i::nobs,:]
import scipy.io as io_mat 
matfile = name_sol+'_fd.mat'
io_mat.savemat(matfile, mdict={'dat_obs': dat_seis,
                               'dt_obs': dt,
                               'nt_obs': nt,
                               'crd_obs':ocoord},
                               oned_as='row')
