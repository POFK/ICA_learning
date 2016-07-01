#!/usr/bin/env python
# coding=utf-8
import numpy as np
import matplotlib.pyplot as plt
#==============  read meta  =================================
PATHdataA='/home/ycli/data/gbt/gbt_15hr_41-80_pointcorr/secA_15hr_41-80_pointcorr_clean_map_I_800.npy'
PATHdataB='/home/ycli/data/gbt/gbt_15hr_41-80_pointcorr/secB_15hr_41-80_pointcorr_clean_map_I_800.npy'
from core import algebra
def ReadMeta(data_path):
    '''return  freq ra dec'''
    data = algebra.make_vect(algebra.load(data_path))
    freq = data.get_axis('freq')
    ra = data.get_axis('ra')
    dec = data.get_axis('dec')
    return freq,ra,dec
Freq,Ra,Dec=ReadMeta(PATHdataA)
def plot(data,freq,sigma=False):
    if sigma:
        mean=data[freq].mean()
        std=data[freq].std()
        min=mean-std
        max=mean+std
    else:
        min=None
        max=None
    plt.pcolormesh(Dec,Ra,data[freq],vmin=min,vmax=max)
    plt.colorbar()
    plt.xlim([Dec.min(),Dec.max()])
    plt.ylim([Ra.min(),Ra.max()])
    plt.xlabel('dec')
    plt.ylabel('ra')
    plt.gca().invert_yaxis()


N=6
#cut=0
dataA=np.load(PATHdataA)#[:,cut:-cut,cut:-cut]
Ra=Ra#[cut:-cut]
Dec=Dec#[cut:-cut]
dataB=np.load(PATHdataB)
shape=dataA.shape
mapA=dataA.reshape(shape[0],-1)
mapB=dataB.reshape(shape[0],-1)
C=np.dot(mapA,mapA.T)
U,s,V=np.linalg.svd(C)

S=np.zeros_like(s)
S[:N]=1.
#map_clean=np.dot((1-np.dot(np.dot(U,np.diag(S)),V)),mapA)
Fg=np.dot(np.dot(np.dot(U,np.diag(S)),V),mapA)
map_clean=mapA-Fg
map_clean=map_clean.reshape(shape)
plt.figure('pca',figsize=(24,18))
plt.subplot(212)
plt.title('remove %d mode'%N)
plt.semilogy(s**0.5,'g.-')
plt.semilogy(s[:N]**0.5,'r.')
plt.ylim(0.1,(s**0.5).max())
plt.subplot(231)
plt.title('origin')
plot(dataA,100)
plt.subplot(232)
plt.title('foreground')
plot(Fg.reshape(shape),100)
plt.subplot(233)
plt.title('cleaned')
plot(map_clean,100,sigma=True)

#print map_clean.shape
#map_clean=map_clean.reshape(shape)
#plt.savefig('/home/mtx/ICA_learning/ICA4GBT/pca/data_pca/PCA_remove_%dmode.png'%N)
plt.show()
#np.save('/home/mtx/ICA_learning/ICA4GBT/pca/data_pca/PCA_remove_%dmode'%N,map_clean)
