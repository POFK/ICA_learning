#!/usr/bin/env python
# coding=utf-8
import h5py 
import numpy as np
import healpy as hp

f1=h5py.File('/home/mtx/data/ICA_sim/21cm_freq200_128.hdf5')
f2=h5py.File('/home/mtx/data/ICA_sim/galaxy_freq200_128.hdf5')
map1=f1['map'].value
map2=f2['map'].value
f1.close()
f2.close()
fr=map1.shape[0]
def fwhm(freq):
    D=40
    return 1.22*1440.*0.21/freq/D
Fwhm=map(fwhm,700)
map0=map1[:,0,:]+map2[:,0,:]
del map1
del map2
for i in range(fr):
    map0[i]=hp.smoothing(map0[i],fwhm=Fwhm[i])
f=h5py.File(nameout,mode='w')
f.create_dataset(name='map',data=map0)
f.close()
##############################################################
1,name2=get_filename()
3=[]
i,j in np.hstack((name1,name2)):
nameout='/home/mtx/data/ICA_sim/data_with_beam/'+i[:5]+j[:3]+j[3:]
i='/home/mtx/data/ICA_sim/'+i
j='/home/mtx/data/ICA_sim/'+j
print i,j,nameout
add_beam(i,j,nameout)
##############################################################
