#!/usr/bin/env python
# coding=utf-8
from IM import IM
from read import ReadMeta
import matplotlib.pyplot as plt
import numpy as np

class PCA(IM):    
    def PcaInit(self,freq,ra,dec,data):
        self.freq=freq
        self.ra=ra
        self.dec=dec
        self.data=data
        self.shape=data.shape
    def RunPca(self):
        print 'Start PCA ...'
        data=self.data
        ra=self.ra
        dec=self.dec
        shape=self.shape
        map=data.reshape(shape[0],-1)
        C=np.dot(map,map.T)
        U,s,V=np.linalg.svd(C)
        S=np.zeros_like(s)
        S[:self.PCn]=1.
        Fg=np.dot(np.dot(np.dot(U,np.diag(S)),V),map)
        data_clean=map-Fg
        self.data_clean=data_clean.reshape(shape)
        self.Fg=Fg.reshape(shape)
        print 'finish PCA ...'
        return self.data_clean
    def Plot(self):
        ra=self.ra
        dec=self.dec
        Ra,Dec=np.meshgrid(dec,ra)
        plt.figure('All frequency map (PCA) ',figsize=(24,12))
        plt.suptitle('plot all frequency and PCn=%d'%self.PCn,x=0.06,y=0.98)
        plt.subplot(1,3,1)
        plt.pcolormesh(Ra,Dec,self.data.mean(axis=0))
        cbar=plt.colorbar()
#       cbar.set_label('mK')
        plt.title('map')
        plt.xlim([dec.min(),dec.max()])
        plt.ylim([ra.min(),ra.max()])
        plt.xlabel('dec')
        plt.ylabel('ra')
        plt.gca().invert_yaxis()

        plt.subplot(1,3,2)
        plt.pcolormesh(Ra,Dec,self.Fg.mean(axis=0))
        cbar=plt.colorbar()
#       cbar.set_label('mK')
        plt.title('removed Fg')
        plt.xlim([dec.min(),dec.max()])
        plt.ylim([ra.min(),ra.max()])
        plt.xlabel('dec')
        plt.ylabel('ra')
        plt.gca().invert_yaxis()

        clean_map_p=self.data_clean.mean(axis=0)
        mean_=clean_map_p.mean()
        std_=clean_map_p.std()
        plt.subplot(1,3,3)
        plt.pcolormesh(Ra,Dec,clean_map_p,vmin=mean_-3*std_,vmax=mean_+3*std_)
        cbar=plt.colorbar()
        cbar.set_label('mK')
        plt.title('cleaned map')
        plt.xlim([dec.min(),dec.max()])
        plt.ylim([ra.min(),ra.max()])
        plt.xlabel('dec')
        plt.ylabel('ra')
        plt.gca().invert_yaxis()

################################################################################
# to test
if (__name__=='__main__'):
    fs=PCA()
    data=np.load(fs.Path1)
    freq,ra,dec=ReadMeta(fs.Path1)
    data=np.load(fs.Path1)[:,10:-10,5:-5]
    ra=ra[10:-10]
    dec=dec[5:-5]
    data*=1000   # converse unit to mK
    fs.PcaInit(freq=freq,ra=ra,dec=dec,data=data)
    cleandata=fs.RunPca()
    fs.Plot()
    plt.show()

