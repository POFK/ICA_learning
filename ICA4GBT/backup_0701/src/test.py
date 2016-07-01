#!/usr/bin/env python
# coding=utf-8
import numpy as np
import matplotlib.pyplot as plt
from PCA import PCA
from ICA import ICA
from getGridPk import getGridPk
from read_par import readPar
from read import ReadMeta

#============================================================
PathPar='Fg.par'# parameter file
par=readPar(PathPar)
print par
Path1=par['Path1']
Path2=par['Path2']
path_data_wigglez='/home/ycli/data/wigglez/gbt_15hr_41-80_pointcorr/reg15data.npy'
#========= Load data ========================================
data_wig=np.load(path_data_wigglez)[:,10:-10,5:-5]
data=np.load(Path1)
freq,ra,dec=ReadMeta(Path1)
data=np.load(Path1)[:,10:-10,5:-5]
ra=ra[10:-10]
dec=dec[5:-5]
data*=1000   # converse unit to mK
#========== run PCA =========================================
pca=PCA(PathPar)
pca.PcaInit(freq=freq,ra=ra,dec=dec,data=data)
data_clean_pca=pca.RunPca()
#pca.Plot()
#plt.savefig('png/meanMapWithPCA.png')

#============================================================
ica=ICA(PathPar)
ica.IcaInit(freq=freq,ra=ra,dec=dec,data=data)
data_clean_ica=ica.RunIca()
#ica.Plot()
#plt.savefig('png/meanMapWithICA.png')

#=============== PCA & ICA ==================================
ica.IcaInit(freq=freq,ra=ra,dec=dec,data=data_clean_pca)
PIca=ica.RunIca()
ica.Plot()
plt.savefig('./PIca.png')
pca.PcaInit(freq=freq,ra=ra,dec=dec,data=data_clean_ica)
IPca=pca.RunPca()
pca.Plot()
plt.savefig('./IPca.png')
#============================================================
getPk=getGridPk(PathPar)
shape=(freq.shape[0],ra.shape[0],dec.shape[0])
ra=ra-ra.mean()
dec=dec-dec.mean()
T_ori=getPk.grid(freq=freq,ra=ra,dec=dec,data=data)
T_ica=getPk.grid(freq=freq,ra=ra,dec=dec,data=data_clean_ica)
T_pca=getPk.grid(freq=freq,ra=ra,dec=dec,data=data_clean_pca)
T_wig=getPk.grid(freq=freq,ra=ra,dec=dec,data=data_wig)
T_IP=getPk.grid(freq=freq,ra=ra,dec=dec,data=IPca)
T_PI=getPk.grid(freq=freq,ra=ra,dec=dec,data=PIca)
#========= plot =============================================
print T_ori.shape,T_ori.dtype
print T_ica.shape,T_ica.dtype

getPk.getPk(T1=T_wig,T2=T_wig)
k,pk_wig,n=getPk.getBin1D()
getPk.getPk(T1=T_ica,T2=T_wig)
k,pk_ica,n=getPk.getBin1D()
getPk.getPk(T1=T_ori,T2=T_wig)
k,pk_ori,n=getPk.getBin1D()
getPk.getPk(T1=T_pca,T2=T_wig)
k,pk_pca,n=getPk.getBin1D()
getPk.getPk(T1=T_IP,T2=T_wig)
k,pk_ip,n=getPk.getBin1D()
getPk.getPk(T1=T_PI,T2=T_wig)
k,pk_pi,n=getPk.getBin1D()

print n
#========================================
plt.close('all')
#plt.loglog(k[3:],pk_wig[3:],'.-',label='wigglez')
plt.loglog(k,np.abs(pk_ori),'.-',label='origin')
plt.loglog(k,np.abs(pk_ica),'.-',label='ica')
plt.loglog(k,np.abs(pk_pca),'.-',label='pca')
plt.loglog(k,np.abs(pk_ip),'.-.',label='ICA > PCA')
plt.loglog(k,np.abs(pk_pi),'.-.',label='PCA > ICA')
#plt.ylim([10**2,10**9])
plt.legend()
plt.show()
#print 'pk_ori',pk_ori
#print 'pk_ica',pk_ica
#print 'pk_pca',pk_pca
#============================================================
