#!/usr/bin/env python
# coding=utf-8
from parameter import *
#from My_matplotlib_par import *
f=open('./filename')
filename=f.readlines()
f.close()
################################################################################
data_wig=np.load(ica.data_wigglez_path)
data_wig=data_wig[:,ra_1:ra_2,dec_1:dec_2]
dataA=np.load('/project/mtx/output/ICA_GBT/cleaned_secA_15hr_41-80_pointcorr_clean_map_I_800.npy_N'+str(N)+'.npy')
dataB=np.load('/project/mtx/output/ICA_GBT/cleaned_secB_15hr_41-80_pointcorr_clean_map_I_800.npy_N'+str(N)+'.npy')
dataC=np.load('/project/mtx/output/ICA_GBT/cleaned_secC_15hr_41-80_pointcorr_clean_map_I_800.npy_N'+str(N)+'.npy')
dataD=np.load('/project/mtx/output/ICA_GBT/cleaned_secD_15hr_41-80_pointcorr_clean_map_I_800.npy_N'+str(N)+'.npy')

def get_randPk(data_gbt):
    K_=np.empty(bins)
    PK_=np.empty(bins)
    for i in filename[5:]:
        i=i[:-1]
        output=i.split('/')[-1]
        path=i
        print path
        data=np.load(path)
        data=data[:,ra_1:ra_2,dec_1:dec_2]
        k,p1,p2,p_cro=ica.Pk(data_gbt,data,bin=bins)
        K_=np.c_[K_,k]
        PK_=np.c_[PK_,p_cro]
    K_=K_[:,1:]
    PK_=PK_[:,1:]
#   np.save('/project/mtx/output/ICA_GBT/Pk',PK_)
    return PK_
#================================================================================
#PK_A=get_randPk(dataA)
#PK_B=get_randPk(dataB)
#PK_C=get_randPk(dataC)
#PK_D=get_randPk(dataD)
#mean_A=PK_A.mean(1)
#std_A=PK_A.std(1)
#mean_B=PK_B.mean(1)
#std_B=PK_B.std(1)
#mean_C=PK_C.mean(1)
#std_C=PK_C.std(1)
#mean_D=PK_D.mean(1)
#std_D=PK_D.std(1)
#PK_=np.c_[PK_A,PK_B,PK_C,PK_D]
#mean=PK_.mean(1)
#std=PK_.std(1)
#S_mean=np.c_[mean_A,mean_B,mean_C,mean_D,mean]
#S_std=np.c_[std_A,std_B,std_C,std_D,std]
#np.save('pk_mean.npy',S_mean)
#np.save('pk_std.npy',S_std)

S_mean=np.load('pk_mean.npy')
S_std=np.load('pk_std.npy')
#================================================================================
pk=np.empty(bins)
S_pk=[]
for i in [dataA,dataB,dataC,dataD]:
    k,p1,p2,p_cro=ica.Pk(i,data_wig,bin=bins)
    S_pk.append(p_cro)
    pk+=p_cro
pk/=4.

plt.figure('ps')
plt.plot(k,np.abs(p2),'-',label='Pk_wig')
#plt.errorbar(k,np.abs(S_pk[0]),yerr=S_std[:,0],fmt='.--',label='cro_A_wig')
#plt.errorbar(k,np.abs(S_pk[1]),yerr=S_std[:,1],fmt='.--',label='cro_B_wig')
#plt.errorbar(k,np.abs(S_pk[2]),yerr=S_std[:,2],fmt='.--',label='cro_C_wig')
#plt.errorbar(k,np.abs(S_pk[3]),yerr=S_std[:,3],fmt='.--',label='cro_D_wig')
plt.errorbar(k,np.abs(pk),yerr=S_std[:,4],fmt='o-',label='cro')
ss=np.loadtxt('explanatory00_z3_pk.dat')
k_cla=ss[:,0]
p_cla=ss[:,1]
p_cla*=0.43*0.29*((0.308+(1+0.7)**-3*0.692)/0.37)**-0.5*((1+0.7)/1.8)**0.5*1.48
p_cla*=k_cla**3/(2*np.pi**2)
plt.plot(k_cla[len(k_cla)/3:],p_cla[len(k_cla)/3:],'-',label='T')
plt.yscale('log')
plt.xscale('log')
plt.legend(loc='upper left')

print S_std[:,4]
print pk
plt.show()
#plt.savefig('PS_weight.eps')
plt.close('all')
