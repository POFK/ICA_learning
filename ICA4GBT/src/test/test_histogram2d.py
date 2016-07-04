#!/usr/bin/env python
# coding=utf-8
import numpy as np
import matplotlib.pyplot as plt

N_bin=5
v=np.ones((10,10))
pos_x=np.random.rand(10,10)
pos_y=np.random.rand(10,10)
DPM=np.c_[pos_x.reshape(-1),pos_y.reshape(-1)]
binx=np.linspace(0,1.1,N_bin+1)
print 'bin:', binx
### np.histogramdd
n,edges=np.histogramdd(DPM,bins=(binx,binx))
v_bin_np=np.histogramdd(DPM,bins=(binx,binx),weights=v.reshape(-1))[0]
print 'histogram:',v_bin_np
###for
v_bin_for=np.zeros([N_bin,N_bin])
n_bin_for=np.zeros([N_bin,N_bin])
for i in np.arange(N_bin):
    for j in np.arange(N_bin):
            bool1= (binx[i] <= pos_x) * (pos_x < binx[i+1])
            bool2= (binx[j] <= pos_y) * (pos_y < binx[j+1])
            bool=bool1*bool2
            v_bin_for[i,j]=v[bool].sum()
            n_bin_for[i,j]=len(v[bool])
print 'with for:',v_bin_for
print 'compare V:',np.allclose(v_bin_for,v_bin_np)
print 'compare n:',np.allclose(n,n_bin_for)
### plot
plt.pcolor(binx,binx,v_bin_np.T)
plt.plot(pos_x.reshape(-1),pos_y.reshape(-1),'r.')
plt.colorbar()
plt.xlim([0,1.1])
plt.ylim([0,1.1])
plt.show()

