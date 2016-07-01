#!/usr/bin/env python
# coding=utf-8

import numpy as np
from astropy.cosmology import FlatLambdaCDM
import matplotlib.pyplot as plt

h=67.77
ICs_N=50
Nx,Ny,Nz=64,64,256
cosmo = FlatLambdaCDM(H0=h, Om0=0.307)
#data_wigglez=np.load('/home/ycli/data/wigglez/gbt_15hr_41-80_pointcorr/reg15data.npy')
data_wigglez=np.load('wigglez_gbt_15hr_41-80_pointcorr_reg15data.npy')
ra=np.load('secA_15hr_41-80_pointcorr_clean_map_I_800_ra.npy')
dec=np.load('secA_15hr_41-80_pointcorr_clean_map_I_800_dec.npy')
freq=np.load('secA_15hr_41-80_pointcorr_clean_map_I_800_freq.npy')
data0=np.load('secA_15hr_41-80_pointcorr_clean_map_I_800_cleaned.npy')
#data1=np.load('secA_15hr_41-80_pointcorr_clean_map_I_800_cleaned_2.npy')
#data2=np.load('secA_15hr_41-80_pointcorr_clean_map_I_800_cleaned_3.npy')
data3=np.load('secA_15hr_41-80_pointcorr_clean_map_I_800_cleaned_5.npy')
data4=np.load('secA_15hr_41-80_pointcorr_clean_map_I_800_cleaned_10.npy')
data5=np.load('secA_15hr_41-80_pointcorr_clean_map_I_800_cleaned_20.npy')
#data6=np.load('secA_15hr_41-80_pointcorr_clean_map_I_800_cleaned_50.npy')
#data7=np.load('secA_15hr_41-80_pointcorr_clean_map_I_800_cleaned_100.npy')
################################################################################
redshift=1420./(freq/10**6)-1
dc = np.array(cosmo.comoving_distance(redshift))*h/100        #z: Mpc/h
pos_x=-dc[:,None,None]*np.sin(np.pi/180.*dec[None,None,:])*np.cos(np.pi/180.*ra[None,:,None])
pos_y=-dc[:,None,None]*np.sin(np.pi/180.*dec[None,None,:])*np.sin(np.pi/180.*ra[None,:,None])
pos_z=-dc[:,None,None]*np.cos(np.pi/180.*ra[None,:,None])+np.zeros_like(dec[None,None,:])
print pos_x.shape,pos_y.shape,pos_z.shape
print pos_x.max(),pos_x.min()
print pos_y.max(),pos_y.min()
print pos_z.max(),pos_z.min()
DPM=np.c_[pos_x.reshape(-1),pos_y.reshape(-1),pos_z.reshape(-1)]
bin_x=np.linspace(pos_x.min()-10**-5,pos_x.max()+10**-5,Nx+1)
bin_y=np.linspace(pos_y.min()-10**-5,pos_y.max()+10**-5,Ny+1)
bin_z=np.linspace(pos_z.min()-10**-5,pos_z.max()+10**-5,Nz+1)
def f(data,data_wigglez,label):
    data=data*10**-3
    T,edges=np.histogramdd(DPM,bins=(bin_x,bin_y,bin_z),weights=data.reshape(-1))
    T_wigglez=np.histogramdd(DPM,bins=(bin_x,bin_y,bin_z),weights=data_wigglez.reshape(-1))[0]
#   plt.figure('1')
#   plt.imshow(T[1,:,:])
#   plt.figure('2')
    Hx=edges[0][1]-edges[0][0]
    Hy=edges[1][1]-edges[1][0]
    Hz=edges[2][1]-edges[2][0]
    k_x=np.fft.fftfreq(Nx,1./Nx)
    k_y=np.fft.fftfreq(Ny,1./Ny)
    k_z=np.fft.fftfreq(Nz,1./Nz)
    T_k1=np.fft.fftn(T)
    T_k2=np.fft.fftn(T_wigglez)
    Pk_cro=Hx*Hy*Hz/Nx/Ny/Nz*((T_k1.conjugate()*T_k2+T_k2.conjugate()*T_k1)/2).real
    Pk_gbt=Hx*Hy*Hz/Nx/Ny/Nz*(np.abs(T_k1)**2)
    Pk_wig=Hx*Hy*Hz/Nx/Ny/Nz*(np.abs(T_k2)**2)
    k_mag=np.sqrt((2*np.pi/Hx/Nx*k_x[:,None,None])**2+(2*np.pi/Hy/Ny*k_y[None,:,None])**2+(2*np.pi/Hz/Nz*k_z[None,None,:])**2)
    
    #################################################################################
    bin=100
    edges=np.linspace(k_mag.min(),k_mag.max(),bin+1,endpoint=True)
    n=np.histogram(k_mag,edges)[0]
    k_bin=np.histogram(k_mag,edges,weights=k_mag)[0]
    pk_cro=np.histogram(k_mag,edges,weights=Pk_cro)[0]
    pk_gbt=np.histogram(k_mag,edges,weights=Pk_gbt)[0]
    pk_wig=np.histogram(k_mag,edges,weights=Pk_wig)[0]
    plt.figure('Pk')
    plt.semilogy(k_bin/n,pk_cro/n,'.-.',label=label+'_cro')
    plt.semilogy(k_bin/n,pk_gbt/n,'.--',label=label+'_gbt')
    plt.semilogy(k_bin/n,pk_wig/n,'.-',label=label+'_wig')
    plt.legend()
    print pk_cro/np.sqrt(pk_gbt*pk_wig)
    plt.figure('CC')
    plt.plot(k_bin/n,pk_cro/np.sqrt(pk_gbt*pk_wig),'.-.',label=label+'_r')

f(data=data0,data_wigglez=data_wigglez,label='data')
#f(data=data1,data_wigglez=data_wigglez,label='ICs_2')
#f(data=data2,data_wigglez=data_wigglez,label='ICs_3')
f(data=data3,data_wigglez=data_wigglez,label='ICs_5')
f(data=data4,data_wigglez=data_wigglez,label='ICs_10')
f(data=data5,data_wigglez=data_wigglez,label='ICs_20')
#f(data=data6,data_wigglez=data_wigglez,label='ICs_50')
#f(data=data7,data_wigglez=data_wigglez,label='ICs_100')
plt.legend()
plt.show()
