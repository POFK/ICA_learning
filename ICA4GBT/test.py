#!/usr/bin/env python
# coding=utf-8
from parameter import *
from My_matplotlib_par import *

np.save('/project/mtx/output/ICA_GBT/freq.npy',freq)
np.save('/project/mtx/output/ICA_GBT/ra.npy',ra)
np.save('/project/mtx/output/ICA_GBT/dec.npy',dec)
data_wig=np.load(ica.data_wigglez_path)
################################################################################
data_A=ica.LoadData(ica.dataA_path,ica.data_wigglez_path)
data_B=ica.LoadData(ica.dataB_path,ica.data_wigglez_path)
data_C=ica.LoadData(ica.dataC_path,ica.data_wigglez_path)
data_D=ica.LoadData(ica.dataD_path,ica.data_wigglez_path)

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
data_wig=data_wig[:,ra_1:ra_2,dec_1:dec_2]

map_clean_A=ica.Clean(data_A,label='a',freq=80,sigma=1,PLOT=1,ra=ra,dec=dec)
#plt.savefig('gbt_a_N100_weight.eps')
map_clean_B=ica.Clean(data_A,label='b',freq=80,sigma=3,PLOT=1,ra=ra,dec=dec)
#plt.savefig('gbt_b_N100_weight.eps')
map_clean_C=ica.Clean(data_A,label='c',freq=80,sigma=5,PLOT=1,ra=ra,dec=dec)
#plt.savefig('gbt_c_N100_weight.eps')
map_clean_D=ica.Clean(data_A,label='d',freq=80,sigma=7,PLOT=1,ra=ra,dec=dec)
#plt.savefig('gbt_d_N100_weight.eps')
plt.show()
plt.close('all')

#plt.show()

#np.save('/project/mtx/output/ICA_GBT/SecA_cleaned_'+str(N)+'.npy',map_clean_A)
#np.save('/project/mtx/output/ICA_GBT/SecB_cleaned_'+str(N)+'.npy',map_clean_B)
#np.save('/project/mtx/output/ICA_GBT/SecC_cleaned_'+str(N)+'.npy',map_clean_C)
#np.save('/project/mtx/output/ICA_GBT/SecD_cleaned_'+str(N)+'.npy',map_clean_D)
##################################################################################
##ica.Grid(32,32,128)
#k_0,p1_0,p2_0,p_cro_0=ica.Pk(map_clean_A,data_wig,bin=bins)
#k,p1,p2,p_cro=ica.Pk(map_clean_B,data_wig,bin=bins)
#k_0+=k
#p_cro_0+=p_cro
#k,p1,p2,p_cro=ica.Pk(map_clean_C,data_wig,bin=bins)
#k_0+=k
#p_cro_0+=p_cro
#k,p1,p2,p_cro=ica.Pk(map_clean_D,data_wig,bin=bins)
#k_0+=k
#p_cro_0+=p_cro
#
#std=np.load('std.npy')
#plt.plot(k,p1,'.-.',label='auto_D')
#plt.plot(k,p2,'.-.',label='wig')
##plt.errorbar(k,np.abs(p_cro)*k**3/(2*np.pi**2),yerr=std,fmt='.-',label='cro')
#plt.plot(k,np.abs(p_cro_0)/4.,'o--',label='cro')
#
#k,p1,p2,p_cro=ica.Pk(map_clean_D,map_clean_A,bin=bins)
#plt.plot(k,p_cro,'.-.',label='cro_AD')
#
#print p_cro_0/4
#print std
#a=np.loadtxt('explanatory00_z3_pk.dat')
#a=a[a.shape[0]/3.:,:]
#k=a[:,0]
#pk=a[:,1]
#pk=pk*0.43*0.29*((0.308+(1+0.7)**-3*0.692)/0.37)**-0.5*((1+0.7)/1.8)**0.5*1.48
#plt.plot(k,pk*k**3/(2*np.pi**2),label='delta')
#plt.yscale('log')
#plt.xscale('log')
#plt.legend(loc='upper left')
#plt.show()
##plt.savefig('ps3_weight.eps')
#
#
