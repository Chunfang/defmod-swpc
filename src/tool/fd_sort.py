#!/usr/bin/env python 
import numpy as np
import os,sys
from itertools import islice
import glob as glb
import h5py 
def save_dict_to_hdf5(dic, filename):
    """
    ....
    """
    with h5py.File(filename, 'w') as h5file:
        recursively_save_dict_contents_to_group(h5file, '/', dic)

def recursively_save_dict_contents_to_group(h5file, path, dic):
    """
    ....
    """
    for key, item in dic.items():
        if isinstance(item, (np.ndarray, np.int64, np.float64, str, bytes)):
            h5file[path + key] = item
        elif isinstance(item, dict):
            recursively_save_dict_contents_to_group(h5file, path + key + '/', item)
        else:
            raise ValueError('Cannot save %s type'%type(item))

def load_dict_from_hdf5(filename):
    """
    ....
    """
    with h5py.File(filename, 'r') as h5file:
        return recursively_load_dict_contents_from_group(h5file, '/')

def recursively_load_dict_contents_from_group(h5file, path):
    """
    ....
    """
    ans = {}
    for key, item in h5file[path].items():
        if isinstance(item, h5py._hl.dataset.Dataset):
            ans[key] = item.value
        elif isinstance(item, h5py._hl.group.Group):
            ans[key] = recursively_load_dict_contents_from_group(h5file, path + key + '/')
    return ans

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
nobs=len(ocoord)
dat_seis_tmp=np.empty(shape=[nt*nobs,dmn],dtype=np.float)
for i in range(nt):
    dat_tmp=np.empty(shape=[0,dmn],dtype=np.float)
    for f,j in zip(fopen,range(len(fopen))):
        dat_tmp=np.vstack((dat_tmp,np.loadtxt(islice(f,nobs_loc[j]),unpack=False,dtype=np.float)))
    dat_tmp=dat_tmp[opick,:]
    dat_seis_tmp[i*nobs:(i+1)*nobs,:]=dat_tmp[osort,:]
    if np.remainder(i+1,100)*np.remainder(i+1,nt)==0:print "frame " + str(i+1) +"/"+str(nt)+ " merged"
dat_seis = np.empty(shape=[nobs,dmn,nt],dtype=np.float)
for i in range(nobs): 
    dat_seis[i,:,:] = np.transpose(dat_seis_tmp[i::nobs,:])

mdict={}
mdict['dat_obs_fd'] = dat_seis
mdict['dt_obs_fd' ] = dt
mdict['nt_obs_fd' ] = np.array(nt, dtype=np.uint32)
mdict['crd_obs_fd'] = ocoord

h5file = name_sol+'_fd.h5'
print 'writing to '+ h5file +'...'    
save_dict_to_hdf5(mdict, h5file)
print h5file + ' created'
