#!/usr/bin/env python
# coding=utf-8
from IM import IM
from read import ReadMeta
from sklearn.decomposition import FastICA
import numpy as np
import matplotlib.pyplot as plt
class ICA(IM):    
    def IcaInit(self,freq,ra,dec,data):
        self.freq=freq
        self.ra=ra
        self.dec=dec
        self.data=data
        self.shape=data.shape
    def RunIca(self):
        print 'Start ICA ...'
        data=self.data
        ICn=self.ICn
        shape=(data.shape[1],data.shape[2])
        S=data.reshape((data.shape[0],-1)).T   #model 2
        ica = FastICA(n_components=self.ICn)
        S_ = ica.fit_transform(S)
        A_ = ica.mixing_
        self.ICs_sum=np.dot(S_,A_.T)
        self.data_clean=(S-self.ICs_sum).T.reshape(self.shape)
        self.ICs_sum=self.ICs_sum.T.reshape(self.shape)
        print 'Finish ICA ...'
        return self.data_clean
    def Plot(self):
        ra=self.ra
        dec=self.dec
        Ra,Dec=np.meshgrid(dec,ra)
        plt.figure('All frequency map',figsize=(24,12))
        plt.suptitle('plot all frequency and ICn=%d'%self.ICn,x=0.06,y=0.98)
        plt.subplot(1,3,1)
        plt.pcolormesh(Ra,Dec,self.data.sum(axis=0))
        cbar=plt.colorbar()
#       cbar.set_label('mK')
        plt.title('map')
        plt.xlim([dec.min(),dec.max()])
        plt.ylim([ra.min(),ra.max()])
        plt.xlabel('dec')
        plt.ylabel('ra')
        plt.gca().invert_yaxis()

        plt.subplot(1,3,2)
        plt.pcolormesh(Ra,Dec,self.ICs_sum.sum(axis=0))
        cbar=plt.colorbar()
#       cbar.set_label('mK')
        plt.title('ICs')
        plt.xlim([dec.min(),dec.max()])
        plt.ylim([ra.min(),ra.max()])
        plt.xlabel('dec')
        plt.ylabel('ra')
        plt.gca().invert_yaxis()

        clean_map_p=self.data_clean.sum(axis=0)
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
    fs=ICA()
    data=np.load(fs.Path1)
    freq,ra,dec=ReadMeta(fs.Path1)
    data=np.load(fs.Path1)[:,10:-10,5:-5]
    ra=ra[10:-10]
    dec=dec[5:-5]
    data*=1000   # converse unit to mK
    fs.IcaInit(freq=freq,ra=ra,dec=dec,data=data)
    cleandata=fs.RunIca()
    fs.Plot()
    plt.show()

