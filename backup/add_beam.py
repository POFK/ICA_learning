#!/usr/bin/env python
# coding=utf-8
import h5py
import healpy as hp
import numpy as np
def add_beam(name1,name2,nameout):
    f1=h5py.File(name1)
    f2=h5py.File(name2)
    map1=f1['map'].value
    map2=f2['map'].value
    f1.close()
    f2.close()
    fr=map1.shape[0]
    def fwhm(freq):
        D=40
        return 1.22*1440.*0.21/freq/D
    freq=np.linspace(700,800,fr,endpoint=True)
    Fwhm=map(fwhm,freq)
    map0=map1[:,0,:]+map2[:,0,:]
    del map1
    del map2
    for i in range(fr):
        map0[i]=hp.smoothing(map0[i],fwhm=Fwhm[i])
    f=h5py.File(nameout,mode='w')
    f.create_dataset(name='map',data=map0)
    f.close()
def get_filename():
    f1=open('/home/mtx/data/ICA_sim/name_21cm','r')
    f2=open('/home/mtx/data/ICA_sim/name_fg','r')
    a=f1.readlines()
    b=f2.readlines()
    f1.close()
    f2.close()
    c=[]
    d=[]
    for i in range(len(a)):
        c.append(['21cm'+a[i][-18:-1]])
    for j in range(len(b)):
        d.append(['fg'+b[j][-18:-1]])
    return c,d
##################################################################
name1,name2=get_filename()
name3=[]
for i,j in np.hstack((name1,name2)):
    nameout='/home/mtx/data/ICA_sim/data_with_beam/'+i[:5]+j[:3]+j[3:]
    i='/home/mtx/data/ICA_sim/'+i
    j='/home/mtx/data/ICA_sim/'+j
    print i,j,nameout
    add_beam(i,j,nameout)
##################################################################
