#!/usr/bin/env python
# coding=utf-8
from read import ReadMeta
from sklearn.decomposition import FastICA
import numpy as np
import matplotlib.pyplot as plt

N=5

#data_path1 = '/home/ycli/data/gbt/gbt_15hr_41-80_pointcorr/secB_15hr_41-80_pointcorr_clean_map_I_800.npy'
#data_path2 = '/home/ycli/data/gbt/gbt_15hr_41-80_pointcorr/secB_15hr_41-80_pointcorr_noise_weight_I_800.npy'

data_path1 = '/home/ycli/data/gbt/gbt_15hr_41-80_pointcorr/secA_15hr_41-80_pointcorr_clean_map_I_800.npy'
data_path2 = '/home/ycli/data/gbt/gbt_15hr_41-80_pointcorr/secB_15hr_41-80_pointcorr_clean_map_I_800.npy'
data_path3 = '/home/ycli/data/gbt/gbt_15hr_41-80_pointcorr/secC_15hr_41-80_pointcorr_clean_map_I_800.npy'
data_path4 = '/home/ycli/data/gbt/gbt_15hr_41-80_pointcorr/secD_15hr_41-80_pointcorr_clean_map_I_800.npy'

data_wigglez=np.load('/home/ycli/data/wigglez/gbt_15hr_41-80_pointcorr/reg15rand035.npy')
#data_wigglez=np.load('wigglez_gbt_15hr_41-80_pointcorr_reg15data.npy')

#data_path = '/home/ycli/data/wigglez/gbt_15hr_41-80_pointcorr/reg15data.npy'
freq,ra,dec=ReadMeta(data_path1)
print dec
print ra
data1=np.load(data_path1)
data2=np.load(data_path2)
data3=np.load(data_path3)
data4=np.load(data_path4)
#data=data1#*data2
data=(data1+data2+data3+data4)/4.
###################select wedge#############################
#select_ra=15
#select_dec=8
#data=data[:,select_ra:-select_ra,select_dec:-select_dec]
#ra=ra[select_ra:-select_ra]
#dec=dec[select_dec:-select_dec]
#data_wigglez=data_wigglez[:,select_ra:-select_ra,select_dec:-select_dec]
################################################################################
np.save('/home/mtx/ICA_learning/ICA4GBT/output/'+data_path1[-45:-4]+'_cleaned.npy',data)
np.save('/home/mtx/ICA_learning/ICA4GBT/output/'+'wigglez_gbt_15hr_41-80_pointcorr_reg15data.npy',data_wigglez)
#S=data.reshape((data.shape[0],-1))    #model 1
S=data.reshape((data.shape[0],-1)).T   #model 2
print 'S',S.shape
ica = FastICA(n_components=N)
S_ = ica.fit_transform(S)
A_ = ica.mixing_

print 'S_',S_.shape
print 'A_',A_.shape
def ICA_clean(data,N,freq_Num):
    S_ICs=0*np.empty_like(data[-1])
    for i in range(N):
#       ICs=S_[freq_Num,i]*A_[:,i]   #for model 1,  2 is a frequency
        ICs=S_[:,i]*A_[freq_Num,i]   #for model 2,  2 is a frequency
        ICs=ICs.reshape((data.shape[1],data.shape[2]))
        S_ICs+=ICs
    a=data[freq_Num]-S_ICs
    return a-a.mean()
result=[]
for j in np.arange(len(freq)):
    result.append(ICA_clean(data,N,j))
result=np.array(result)
#plt.imshow(ICA_clean(data,N,250),extent=[dec[0],dec[-1],ra[0],ra[-1]])
#plt.colorbar()
#plt.show()
np.save('/home/mtx/ICA_learning/ICA4GBT/output/'+data_path1[-45:-4]+'_cleaned_'+str(N)+'.npy',result)
np.save('/home/mtx/ICA_learning/ICA4GBT/output/'+data_path1[-45:-4]+'_freq.npy',freq)
np.save('/home/mtx/ICA_learning/ICA4GBT/output/'+data_path1[-45:-4]+'_ra.npy',ra)
np.save('/home/mtx/ICA_learning/ICA4GBT/output/'+data_path1[-45:-4]+'_dec.npy',dec)
