#!/usr/bin/env python
# coding=utf-8
import numpy as np

N_bin=5
#v=np.random.randn(100,100,100)
v=np.arange(10**6)
v=v.reshape(100,100,100)
x=np.linspace(0,99,100)
pos_x=x[:,None,None]*np.ones((100,100,100))
pos_y=x[None,:,None]*np.ones((100,100,100))
pos_z=x[None,None,:]*np.ones((100,100,100))
DPM=np.c_[pos_x.reshape(-1),pos_y.reshape(-1),pos_z.reshape(-1)]
binx=np.linspace(0,99.1,N_bin+1)
print 'bin:', binx
### np.histogramdd
n,edges=np.histogramdd(DPM,bins=(binx,binx,binx))
v_bin_np=np.histogramdd(DPM,bins=(binx,binx,binx),weights=v.reshape(-1))[0]
print 'histogram:',v_bin_np
###for
v_bin_for=np.zeros([N_bin,N_bin,N_bin])
n_bin_for=np.zeros([N_bin,N_bin,N_bin])
for i in np.arange(N_bin):
    for j in np.arange(N_bin):
        for k in np.arange(N_bin):
            bool1= (binx[i] <= pos_x) * (pos_x < binx[i+1])
            bool2= (binx[j] <= pos_y) * (pos_y < binx[j+1])
            bool3= (binx[k] <= pos_z) * (pos_z < binx[k+1])
            bool=bool1*bool2*bool3
            v_bin_for[i,j,k]=v[bool].sum()
            n_bin_for[i,j,k]=len(v[bool])
print 'with for:',v_bin_for
print 'compare V:',np.allclose(v_bin_for,v_bin_np)
print 'compare n:',np.allclose(n,n_bin_for)
