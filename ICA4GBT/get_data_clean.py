#!/usr/bin/env python
# coding=utf-8
from parameter import *

np.save('/project/mtx/output/ICA_GBT/freq.npy',freq)
np.save('/project/mtx/output/ICA_GBT/ra.npy',ra)
np.save('/project/mtx/output/ICA_GBT/dec.npy',dec)

#f=open('./filename')
#filename=f.readlines()
#f.close()
#################################################################################
#def get_data(path,output):
#    data=np.load(path)
#    data=data[:,ra_1:ra_2,dec_1:dec_2]
#    map_clean_A=ica.Clean(data,label='a',freq=0,PLOT=0)
#    print 'writing:','/project/mtx/output/ICA_GBT/'+'cleaned_'+output+'_N'+str(N)+'.npy'
#    np.save('/project/mtx/output/ICA_GBT/'+'cleaned_'+output+'_N'+str(N)+'.npy',map_clean_A)


#data_A=ica.LoadData(ica.dataA_path,ica.data_wigglez_path)
#data_B=ica.LoadData(ica.dataB_path,ica.data_wigglez_path)
#data_C=ica.LoadData(ica.dataC_path,ica.data_wigglez_path)
#data_D=ica.LoadData(ica.dataD_path,ica.data_wigglez_path)

data_A=np.load(ica.dataA_path)
data_B=np.load(ica.dataB_path)
data_C=np.load(ica.dataC_path)
data_D=np.load(ica.dataD_path)


#add weight
#weight_A=np.load('/home/ycli/data/gbt/gbt_15hr_41-80_pointcorr/secA_15hr_41-80_pointcorr_noise_weight_I_800.npy')
#weight_B=np.load('/home/ycli/data/gbt/gbt_15hr_41-80_pointcorr/secB_15hr_41-80_pointcorr_noise_weight_I_800.npy')
#weight_C=np.load('/home/ycli/data/gbt/gbt_15hr_41-80_pointcorr/secC_15hr_41-80_pointcorr_noise_weight_I_800.npy')
#weight_D=np.load('/home/ycli/data/gbt/gbt_15hr_41-80_pointcorr/secD_15hr_41-80_pointcorr_noise_weight_I_800.npy')
#data_A=(data_A*weight_A)
#data_B=(data_B*weight_B)
#data_C=(data_C*weight_C)
#data_D=(data_D*weight_D)

data_A=data_A[:,ra_1:ra_2,dec_1:dec_2]
data_B=data_B[:,ra_1:ra_2,dec_1:dec_2]
data_C=data_C[:,ra_1:ra_2,dec_1:dec_2]
data_D=data_D[:,ra_1:ra_2,dec_1:dec_2]

map_clean_A=ica.Clean(data_A,label='a',freq=80,sigma=1,PLOT=0)
map_clean_B=ica.Clean(data_B,label='b',freq=80,sigma=3,PLOT=0)
map_clean_C=ica.Clean(data_C,label='c',freq=80,sigma=5,PLOT=0)
map_clean_D=ica.Clean(data_D,label='d',freq=80,sigma=7,PLOT=0)
#plt.close('all')

np.save('/project/mtx/output/ICA_GBT/'+'cleaned_'+'secA_15hr_41-80_pointcorr_clean_map_I_800.npy'+'_N'+str(N)+'.npy',map_clean_A)
np.save('/project/mtx/output/ICA_GBT/'+'cleaned_'+'secB_15hr_41-80_pointcorr_clean_map_I_800.npy'+'_N'+str(N)+'.npy',map_clean_B)
np.save('/project/mtx/output/ICA_GBT/'+'cleaned_'+'secC_15hr_41-80_pointcorr_clean_map_I_800.npy'+'_N'+str(N)+'.npy',map_clean_C)
np.save('/project/mtx/output/ICA_GBT/'+'cleaned_'+'secD_15hr_41-80_pointcorr_clean_map_I_800.npy'+'_N'+str(N)+'.npy',map_clean_D)





#for i in filename[4:]:
#    i=i[:-1]
#    get_data(i,i.split('/')[-1])
