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
#========= Load data ========================================
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
pca.Plot()
plt.savefig('png/meanMapWithPCA.png')

#============================================================
ica=ICA(PathPar)
ica.IcaInit(freq=freq,ra=ra,dec=dec,data=data)
data_clean_ica=ica.RunIca()
ica.Plot()
plt.savefig('png/meanMapWithICA.png')

#============================================================
getPk=getGridPk(PathPar)
shape=(freq.shape[0],ra.shape[0],dec.shape[0])
ra=ra-ra.mean()
dec=dec-dec.mean()
T_ori=getPk.grid(freq=freq,ra=ra,dec=dec,data=data)
T_ica=getPk.grid(freq=freq,ra=ra,dec=dec,data=data_clean_ica)
T_pca=getPk.grid(freq=freq,ra=ra,dec=dec,data=data_clean_pca)
#========= plot =============================================
print T_ori.shape,T_ori.dtype
print T_ica.shape,T_ica.dtype

getPk.getPk(T1=T_ica,T2=T_ica)
k,pk_ica,n=getPk.getBin1D()
getPk.getPk(T1=T_ori,T2=T_ori)
k,pk_ori,n=getPk.getBin1D()
getPk.getPk(T1=T_pca,T2=T_pca)
k,pk_pca,n=getPk.getBin1D()
#========================================
plt.close('all')
plt.loglog(k[3:],pk_ori[3:],'.-',label='origin')
plt.loglog(k[3:],pk_ica[3:],'.-',label='ica')
plt.loglog(k[3:],pk_pca[3:],'.-',label='pca')
#plt.ylim([10**2,10**9])
plt.legend()
#plt.show()
print 'pk_ori',pk_ori
print 'pk_ica',pk_ica
print 'pk_pca',pk_pca
#============================================================
