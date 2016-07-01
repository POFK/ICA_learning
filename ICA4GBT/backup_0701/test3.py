#!/usr/bin/env python
# coding=utf-8
from parameter import *
bins=20
#data_wig=np.load(ica.data_wigglez_path)
data_wig=np.load(ica.dataA_path)
data_wig=data_wig[:,ra_1:ra_2,dec_1:dec_2]
data_wig=np.random.randn(data_wig.shape[0],data_wig.shape[1],data_wig.shape[2])


ica.Grid(32,64,32)
k,p0,p2,p_cro=ica.Pk(data_wig,data_wig,bin=bins)
plt.plot(k,p2,'.-',label='1')
ica.Grid(64,32,32)
k,p1,p2,p_cro=ica.Pk(data_wig,data_wig,bin=bins)
plt.plot(k,p2,'.-',label='2')
ica.Grid(32,32,64)
k,p1,p2,p_cro=ica.Pk(data_wig,data_wig,bin=bins)
plt.plot(k,p2,'.-',label='3')


ica.Grid(32,32,128)
k,p0,p2,p_cro=ica.Pk(data_wig,data_wig,bin=bins)
plt.plot(k,p_cro,'g.-.',label='32 32 128')
#ica.Grid(32,64,256)
#k,p1,p2,p_cro=ica.Pk(data_wig,data_wig,bin=bins)
#plt.plot(k,p_cro,'gx--',label='32 64 256')
ica.Grid(32,32,256)
k,p1,p2,p_cro=ica.Pk(data_wig,data_wig,bin=bins)
plt.plot(k,p_cro,'bx--',label='32 32 256')
ica.Grid(64,64,256)
k,p1,p2,p_cro=ica.Pk(data_wig,data_wig,bin=bins)
plt.plot(k,p_cro,'rx--',label='64 64 256')

#ica.Grid(32,32,256)
#k,p1,p2,p_cro=ica.Pk(data_wig,data_wig,bin=bins)
#plt.plot(k,p_cro,'b*--',label='32 32 256')
#ica.Grid(64,32,256)
#k,p1,p2,p_cro=ica.Pk(data_wig,data_wig,bin=bins)
#plt.plot(k,p_cro,'g*--',label='64 32 256')
#ica.Grid(16,32,256)
#k,p1,p2,p_cro=ica.Pk(data_wig,data_wig,bin=bins)
#plt.plot(k,p_cro,'r*--',label='16 32 256')





print p0/p1
plt.xscale('log')
plt.yscale('log')
plt.legend(loc='upper left')
plt.show()

