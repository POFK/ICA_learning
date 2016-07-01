#!/usr/bin/env python
# coding=utf-8
from parameter import *

np.save('/project/mtx/output/ICA_GBT/freq.npy',freq)
np.save('/project/mtx/output/ICA_GBT/ra.npy',ra)
np.save('/project/mtx/output/ICA_GBT/dec.npy',dec)

f=open('./filename')
filename=f.readlines()
f.close()

data_wig=np.load('/project/mtx/output/ICA_GBT/cleaned_reg15data.npy_N'+str(N)+'.npy')
dataA=np.load('/project/mtx/output/ICA_GBT/cleaned_secA_15hr_41-80_pointcorr_clean_map_I_800.npy_N'+str(N)+'.npy')
dataB=np.load('/project/mtx/output/ICA_GBT/cleaned_secB_15hr_41-80_pointcorr_clean_map_I_800.npy_N'+str(N)+'.npy')
dataC=np.load('/project/mtx/output/ICA_GBT/cleaned_secC_15hr_41-80_pointcorr_clean_map_I_800.npy_N'+str(N)+'.npy')
dataD=np.load('/project/mtx/output/ICA_GBT/cleaned_secD_15hr_41-80_pointcorr_clean_map_I_800.npy_N'+str(N)+'.npy')

################################################################################
def get_data(path,output):
    print path
    data=np.load(path)
    data=data[:,ra_1:ra_2,dec_1:dec_2]
    k,p1,p2,p_cro=ica.Pk(data,dataA,bin=bins)
    return k,p_cro


K_=np.empty(bins)
PK_=np.empty(bins)
for i in filename[5:]:
    i=i[:-1]
    k,p_cro=get_data(i,i.split('/')[-1])
    K_=np.c_[K_,k]
    PK_=np.c_[PK_,p_cro]
K_=K_[:,1:]
PK_=PK_[:,1:]*K_**3/(2*np.pi**2)
np.save('k',K_)
np.save('Pk',PK_)

