#!/usr/bin/env python
# coding=utf-8
from read import ReadMeta
from sklearn.decomposition import FastICA
import numpy as np
import matplotlib.pyplot as plt

data_path1 = '/home/ycli/data/gbt/gbt_15hr_41-80_pointcorr/secA_15hr_41-80_pointcorr_clean_map_I_800.npy'
data_path2 = '/home/ycli/data/gbt/gbt_15hr_41-80_pointcorr/secA_15hr_41-80_pointcorr_noise_weight_I_800.npy'

#data_path = '/home/ycli/data/wigglez/gbt_15hr_41-80_pointcorr/reg15data.npy'
freq,ra,dec=ReadMeta(data_path1)
print 'dec',dec.shape
print 'ra',ra.shape
data1=np.load(data_path1)
data2=np.load(data_path2)
data=data1
###################select wedge#############################
#select_ra=5
#select_dec=5
#data=data[:,select_ra:-select_ra,select_dec:-select_dec]
#ra=ra[select_ra:-select_ra]
#dec=dec[select_dec:-select_dec]
################################################################################
N=10
freq_Num=10
#S=data.reshape((data.shape[0],-1))    #model 1
S=data.reshape((data.shape[0],-1)).T   #model 2
print 'S',S.shape
ica = FastICA(n_components=N)
S_ = ica.fit_transform(S)
A_ = ica.mixing_

print 'S_',S_.shape
print 'A_',A_.shape
S_ICs=np.empty_like(data[0])
for i in range(N):
#   ICs=S_[freq_Num,i]*A_[:,i]   #for model 1,  2 is a frequency
    ICs=S_[:,i]*A_[freq_Num,i]   #for model 2,  2 is a frequency
    print ICs.shape
    ICs=ICs.reshape((data.shape[1],data.shape[2]))
    S_ICs+=ICs
    plt.subplot(2,5,i)
    plt.title(str(i))
    plt.imshow(ICs,extent=[dec[0],dec[-1],ra[0],ra[-1]])
#   plt.pcolor(ICs)
    plt.xlabel('dec')
    plt.ylabel('ra')
    plt.colorbar()
#plt.savefig('ICs.png')
#plt.cla()
#plt.clf()
map_cleaned=data[freq_Num]-S_ICs
plt.figure('map_cleaned')
plt.subplot(1,2,1)
plt.title('map_original')
plt.imshow(data[freq_Num],extent=[dec[0],dec[-1],ra[0],ra[-1]])
plt.xlabel('dec')
plt.ylabel('ra')
plt.colorbar()
plt.subplot(1,2,2)
plt.title('map_cleaned')
plt.imshow(map_cleaned-map_cleaned.mean(),extent=[dec[0],dec[-1],ra[0],ra[-1]])
#plt.pcolor(map_cleaned)#,extent=[dec[0],dec[-1],ra[0],ra[-1]])
plt.xlabel('dec')
plt.ylabel('ra')
plt.colorbar()
plt.show()
#plt.savefig('map_cleaned.png')
