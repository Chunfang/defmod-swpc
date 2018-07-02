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

# Model poroelastics, dimension and rate-state-friction 
name_sol = sys.argv[1]
if 'slip' in sys.argv: 
    slip=True 
else:
    slip=False
if 'dyn' in sys.argv:
    dyn=True
else:
    dyn=False
if 'clean' in sys.argv: 
    clean=True 
else:
    clean=False
if 'rsf' in sys.argv:
    rsf=1
else:
    rsf=0
if 'poro' in sys.argv:
    p=1
else:
    p=0
if 'alpha' in sys.argv: 
    alpha=1
else:
    alpha=0
if '2D' in sys.argv or '2d' in sys.argv:
    dmn=2
else:
    dmn=3

# File sets
log_file     = name_sol+".log"
log_dyn_file = name_sol+"_dyn.log"
files_qs = sorted(glb.glob(name_sol+"_obs_*.txt"))
files_seis=[file_qs.replace("obs","dyn_obs") for file_qs in files_qs]

print 'Merge quasi-static observation grid/data...'
fopen=[]; nobs_loc=[]
for file_qs in files_qs: fopen.append(open(file_qs))
ocoord=np.empty(shape=[0,dmn+1],dtype=np.float)
for f in fopen:
    n=int(np.genfromtxt(islice(f,1),delimiter=" ",unpack=False,dtype=np.uint32))
    ocoord=np.vstack((ocoord,np.genfromtxt(islice(f,n),delimiter=" ",unpack=False,dtype=np.float)))
    nobs_loc.append(n)
oidx,opick=np.unique(ocoord[:,-1],return_index=True)
osort=np.argsort(oidx)
ocoord=ocoord[opick,:dmn]; ocoord=ocoord[osort,:dmn]
dat_log = np.loadtxt(log_file, delimiter=" ",skiprows=1, unpack=False, dtype=np.uint32)
nframe=dat_log[-1,0]
dat_log=dat_log[:-1,:]
dat_qs_tmp=np.empty(shape=[0,dmn+p],dtype=np.float)
for i in range(nframe):
    dat_tmp=np.empty(shape=[0,dmn+p],dtype=np.float)
    for f,j in zip(fopen,range(len(fopen))):
        dat_tmp=np.vstack((dat_tmp,np.loadtxt(islice(f,nobs_loc[j]),unpack=False,dtype=np.float)))
    dat_tmp=dat_tmp[opick,:]
    dat_qs_tmp=np.vstack((dat_qs_tmp,dat_tmp[osort,:]))
    if np.remainder(i+1,10)*np.remainder(i+1,nframe)==0:print "frame " + str(i+1) +"/"+str(nframe)+ " merged"

print 'Merge seismic observation grid/data...'
fopen=[]; nobs_loc=[]
for file_seis in files_seis: fopen.append(open(file_seis))
ocoord=np.empty(shape=[0,dmn+1],dtype=np.float)
for f in fopen:
    n=int(np.genfromtxt(islice(f,1),delimiter=" ",unpack=False,dtype=np.uint32))
    ocoord=np.vstack((ocoord,np.genfromtxt(islice(f,n),delimiter=" ",unpack=False,dtype=np.float)))
    nobs_loc.append(n)
oidx,opick=np.unique(ocoord[:,-1],return_index=True)
osort=np.argsort(oidx)
ocoord=ocoord[opick,:dmn]; ocoord=ocoord[osort,:dmn]
dat_log_dyn = np.loadtxt(log_dyn_file, delimiter=" ",skiprows=1, unpack=False, dtype=np.uint32)
try:
    nframe=dat_log_dyn[-1]
except:
    nframe=dat_log_dyn.item()
dat_seis_tmp=np.empty(shape=[nframe*len(ocoord),dmn],dtype=np.float)
for i in range(nframe):
    dat_tmp=np.empty(shape=[0,dmn],dtype=np.float)
    for f,j in zip(fopen,range(len(fopen))):
        dat_tmp=np.vstack((dat_tmp,np.loadtxt(islice(f,nobs_loc[j]),unpack=False,dtype=np.float)))
    dat_tmp=dat_tmp[opick,:]
    dat_seis_tmp[i*len(osort):(i+1)*len(osort),:]=dat_tmp[osort,:]
    if np.remainder(i+1,100)*np.remainder(i+1,nframe)==0:print "frame " + str(i+1) +"/"+str(nframe)+ " merged"
dt = np.loadtxt(log_file, delimiter=" ", usecols=[0], unpack=False, dtype=np.float)[0]
dt_dyn = np.loadtxt(log_dyn_file, dtype=np.float)[0]

# For equivalent implicit dynamic solver
if alpha:
    log_dyn_file_alpha = name_sol+"-alpha_dyn.log"
    dt_alpha=np.loadtxt(log_dyn_file_alpha, dtype=np.float)[0]
    dat_log_dyn_alpha = np.loadtxt(log_dyn_file_alpha, delimiter=" ",skiprows=1, unpack=False, dtype=np.uint32)
    try: 
        nframe=dat_log_dyn_alpha[-1]
    except:
        nframe=dat_log_dyn_alpha.item()
    files_alpha=sorted(glb.glob(name_sol+"-alpha_dyn_obs_*.txt"))
    fopen=[]
    for file_alpha in files_alpha: fopen.append(open(file_alpha))
    for f,j in zip(fopen, range(len(fopen))): # Remove headers
        n = np.genfromtxt(islice(f,1),delimiter=" ",unpack=False,dtype=np.uint32)
        _ = np.genfromtxt(islice(f,n),delimiter=" ",unpack=False,dtype=np.float)
    dat_seis_alpha_tmp = np.empty(shape=[nframe*len(ocoord),dmn],dtype=np.float)
    for i in range(nframe):
        dat_tmp=np.empty(shape=[nframe*len(ocoord),dmn],dtype=np.float)
        for f,j in zip(fopen,range(len(fopen))):
            dat_tmp=np.vstack((dat_tmp,np.loadtxt(islice(f,nobs_loc[j]),unpack=False,dtype=np.float)))
        dat_tmp=dat_tmp[opick,:]
        dat_seis_alpha_tmp[i*len(osort):(i+1)*len(osort),:]=dat_tmp[osort,:]
        if np.remainder(i+1,100)*np.remainder(i+1,nframe)==0:print "alpha frame " + str(i+1) +"/"+str(nframe)+ " merged"
    
if slip:
    # File sets 
    log_slip_file = name_sol+"_slip.log"
    files_slip = sorted(glb.glob(name_sol+"_fe_*.h5")) 
    dt_slip = np.loadtxt(log_slip_file, dtype=np.float)[0]
    dat_log_slip = np.loadtxt(log_slip_file, dtype=np.uint32)[1:]
    h5 = h5py.File(files_slip[0], 'r')
    nframe_qs = np.shape(h5['slip_sta'])[-1]
    print 'Merge seismic fault slip grid/data...'
    nframe=dat_log_slip[-1]
    nfnd_all=0
    for f in files_slip:
        h5 = h5py.File(f, 'r')
        nfnd_all+=np.shape(h5['slip_sta'])[0]
    fcoord=np.empty(shape=[nfnd_all,dmn],dtype=np.float)
    if dyn: dat_slip_tmp=np.empty(shape=[nfnd_all,dmn,nframe],dtype=np.float)
    h5 = h5py.File(files_slip[0], 'r')
    dat_slip_sta=np.empty(shape=[nfnd_all,dmn,nframe_qs],dtype=np.float)
    dat_trac_sta=np.empty(shape=[nfnd_all,dmn+p,nframe_qs],dtype=np.float)
    fopen=[];n_lmnd=[]
    j0=0
    for f in files_slip:
        h5 = h5py.File(f, 'r')
        j1=j0+len(h5['fault_x'])
        fcoord[j0:j1,:]=h5['fault_x']    
        dat_slip_sta[j0:j1,:,:]=h5['slip_sta']
        dat_trac_sta[j0:j1,:,:]=h5['trac_sta']
        if (dyn): dat_slip_tmp[j0:j1,:,:]=h5['slip_dyn']
        j0=j1
        print 'fault patch ' + str(j0) +"/"+str(nfnd_all)+ " merged"
    nfnd=len(fcoord)

    # Equivalent g-alpha solver
    if alpha:
        log_slip_file_alpha = name_sol+"-alpha_slip.log"
        dat_log_slip_alpha = np.loadtxt(log_slip_file_alpha, dtype=np.uint32)[1:]
        dt_slip_alpha = np.loadtxt(log_slip_file_alpha, dtype=np.float)[0]
        nframe=dat_log_slip_alpha[-1]
        dat_slip_alpha_tmp=np.empty(shape=[nfnd_all,dmn,nframe],dtype=np.float)
        files_slip_alpha = sorted(glb.glob(name_sol+"-alpha_fe_*.h5"))
        j0=0
        for f in files_slip_alpha:
            h5 = h5py.File(f, 'r')
            j1=j0+len(h5['fault_x'])
            dat_slip_alpha_tmp[j0:j1,:,:]=h5['slip_dyn']
            j0=j1

    if rsf==1:
        log_rsf_file=name_sol+"_rsf.log"
        try:
            dt_rsf=np.loadtxt(log_rsf_file,dtype=np.float)[0]
            pseudo=True
        except:
            pseudo=False
        if pseudo: 
            print 'Merge pseudo time fault grid/data...' 
            dat_log_rsf=np.loadtxt(log_rsf_file,dtype=np.uint32)[1:]
            files_rsf=[file_slip.replace('slip','rsf') for file_slip in files_slip]
            fcoord=np.empty(shape=[0,dmn+1],dtype=np.float)
            dat_rsf_tmp=np.empty(shape=[0,2],dtype=np.float)
            fopen=[];n_lmnd=[]
            for file_rsf in files_rsf: fopen.append(open(file_rsf))
            for f in fopen:
                n=int(np.genfromtxt(islice(f,1),delimiter=" ",unpack=False,dtype=np.uint32))
                fcoord=np.vstack((fcoord,np.genfromtxt(islice(f,n),delimiter=" ",unpack=False,dtype=np.float)))
                n_lmnd.append(n)
            fsort=np.argsort(fcoord[:,-1])
            fcoord=fcoord[fsort,:dmn]
            nfnd=sum(n_lmnd)
            nframe=dat_log_rsf[-1]
            for i in range(nframe):
                dat_tmp=np.empty(shape=[0,2],dtype=np.float)
                for f,j in zip(fopen,range(len(fopen))):
                    dat_tmp=np.vstack((dat_tmp,np.loadtxt(islice(f,n_lmnd[j]),unpack=False,dtype=np.float)))
                dat_rsf_tmp=np.vstack((dat_rsf_tmp,dat_tmp[fsort,:]))
                if np.remainder(i+1,24)*np.remainder(i+1,nframe)==0:print "frame " + str(i+1) +"/"+str(nframe)+ " merged"

# Sort quasi-static and waveform by obs 
dmn = dat_seis_tmp.shape[1]
nobs=len(ocoord)
dat_seis = np.empty(shape=[nobs,dmn,dat_seis_tmp.shape[0]/nobs],dtype=np.float)
if alpha: dat_seis_alpha = np.empty([nobs,dmn,dat_seis_alpha_tmp.shape[0]/nobs],dtype=np.float)
dat_qs_sort = np.empty(shape=[nobs,dmn,dat_qs_tmp.shape[0]/nobs],dtype=np.float) 
for i in range(nobs): 
    dat_seis[i,:,:] = np.transpose(dat_seis_tmp[i::nobs,:])
    dat_qs_sort[i,:,:] = np.transpose(dat_qs_tmp[i::nobs,:dmn])
    if alpha: dat_seis_alpha[i,:,:] = np.transpose(dat_seis_alpha_tmp[i::nobs,:])
# Sort fault slip by frame
if slip and rsf==1 and pseudo: # RSF pseudo time
    dat_rsf=np.empty(shape=[nfnd,2,len(dat_rsf_tmp)/nfnd], dtype=np.float)
    for i in range(len(dat_rsf_tmp)/nfnd):
        dat_rsf[:,:,i]=dat_rsf_tmp[i*nfnd:(i+1)*nfnd,:]
# Sort seismic/slip data by event 
dat_seis_sort={}
if alpha: dat_seis_alpha_sort={}
if slip: 
    dat_slip_sort = {} 
    if alpha:
        dat_slip_alpha_sort={}
if len(dat_log_dyn.shape)==0:dat_log_dyn=[dat_log_dyn.item()]
if alpha:
    if len(dat_log_dyn_alpha.shape)==0: dat_log_dyn_alpha=[dat_log_dyn_alpha.item()]
for i in range(len(dat_log_dyn)):
    if i==0:
        start=0
    else:
        start=dat_log_dyn[i-1]
    end=dat_log_dyn[i]
    dat_seis_sort['step '+str(dat_log[i,0])] = dat_seis[:,:,start:end]
    if alpha:
        if i==0:
            start=0
        else:
            start=dat_log_dyn_alpha[i-1]
        end=dat_log_dyn_alpha[i]
        dat_seis_alpha_sort['step '+str(dat_log[i,0])] = dat_seis_alpha[:,:,start:end]
    if slip and dyn:
        if i==0:
            start=0
        else:
            start=dat_log_slip[i-1]
        end=dat_log_slip[i]
        dat_slip_sort['step '+str(dat_log[i,0])] = dat_slip_tmp[:,:,start:end]
        if alpha:
            if i==0:
                start=0
            else:
                start=dat_log_slip_alpha[i-1]
            end=dat_log_slip_alpha[i]
            dat_slip_alpha_sort['step '+str(dat_log[i,0])] = dat_slip_alpha_tmp[:,:,start:end]

# Sort rate state data by quasi-static time step.
if slip and rsf==1 and pseudo: 
    dat_rsf_sort=[]
    for i in range(len(dat_log_rsf)):
        if i==0:
            start=0
        else:
            start=dat_log_rsf[i-1]
        dat_rsf_sort.append(dat_rsf[start:dat_log_rsf[i],:])

# Store sorted data to .mat files
h5file = name_sol+'_fe.h5'

# Save to .h5 file
mdict={'dat_obs_sta': dat_qs_sort,
       'dat_obs_dyn': dat_seis_sort,
       'crd_obs': ocoord,
       'dt': dt,
       'dt_dyn': dt_dyn,
       'dat_log': dat_log,
       'dat_log_dyn': np.array(dat_log_dyn)}
if slip:
    mdict['crd_flt'] = fcoord
    mdict['dat_log_slip'] = dat_log_slip
    mdict['dat_fqs'] = dat_trac_sta
    mdict['dat_slip_sta']=dat_slip_sta
    if dyn:
        mdict['dt_slip'] =  dt_slip
        mdict['dat_slip'] = dat_slip_sort
    if alpha:
        mdict['dt_slip_alpha'] = dt_slip_alpha
        mdict['dat_slip_alpha'] = dat_slip_alpha_sort
    if rsf and pseudo:
        tmp_arr = np.zeros((len(dat_rsf_sort),), dtype=np.object)
        for i in range(len(tmp_arr)):
            tmp_arr[i] = dat_rsf_sort[i]
        dat_rsf_sort = tmp_arr
        mdict['dt_rsf'] = dt_rsf
        mdict['dat_log_rsf'] = dat_log_rsf
        mdict['dat_rsf'] = dat_rsf_sort
if alpha:
    mdict['dt_dyn_alpha'] = dt_alpha
    mdict['dat_obs_alpha'] = dat_seis_alpha_sort
print 'writing to '+ h5file +'...'    
save_dict_to_hdf5(mdict, h5file)
print h5file + ' created'

if clean:
    print 'Cleanup...'
    for f in list(set([log_file,log_dyn_file])|set(files_qs)|set(files_seis)):
        if os.path.isfile(f):
            os.remove(f)
            print 'file '+f+' deleted'
    if slip:
        for f in list(set([fqs_file,log_slip_file])|set(files_slip)):
            if os.path.isfile(f):
                os.remove(f)
                print 'file '+f+' deleted'
        if rsf==1 and pseudo:
            for f in list(set([log_rsf_file])|set(files_rsf)):
                if os.path.isfile(f):
                    os.remove(f)
                    print 'file '+f+' deleted'
        if alpha:
            for f in list(set([log_slip_file_alpha])|set(files_slip_alpha)):
                    if os.path.isfile(f):
                        os.remove(f)
                        print 'file '+f+' deleted'
    if alpha:
        for f in list(set([log_dyn_file_alpha])|set(files_alpha)):
            if os.path.isfile(f):
                os.remove(f)
                print 'file '+f+' deleted'
