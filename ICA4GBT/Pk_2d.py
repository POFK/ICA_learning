#!/usr/bin/env python
# coding=utf-8
import numpy as np
from astropy.cosmology import FlatLambdaCDM
import matplotlib.pyplot as plt
ra=np.load('secB_15hr_41-80_pointcorr_clean_map_I_800_ra.npy')
dec=np.load('secB_15hr_41-80_pointcorr_clean_map_I_800_dec.npy')
freq=np.load('secB_15hr_41-80_pointcorr_clean_map_I_800_freq.npy')
data=np.load('secB_15hr_41-80_pointcorr_clean_map_I_800_cleaned.npy')
data1=np.load('secB_15hr_41-80_pointcorr_clean_map_I_800_cleaned_2.npy')
data2=np.load('secB_15hr_41-80_pointcorr_clean_map_I_800_cleaned_3.npy')
data3=np.load('secB_15hr_41-80_pointcorr_clean_map_I_800_cleaned_5.npy')
data4=np.load('secB_15hr_41-80_pointcorr_clean_map_I_800_cleaned_10.npy')
data5=np.load('secB_15hr_41-80_pointcorr_clean_map_I_800_cleaned_20.npy')
data6=np.load('secB_15hr_41-80_pointcorr_clean_map_I_800_cleaned_50.npy')
def f(data=data4,ra=ra,dec=dec,freq=freq,label='ICs_10'):
    h=67.77
    cosmo = FlatLambdaCDM(H0=h, Om0=0.307)
    z=1420.4/(freq/10**6)-1
    dc = np.array(cosmo.comoving_distance(z))*h/100        #z: Mpc/h
    N=50
    data=data[N]*10**-3
    freq=freq[N]
    dc=dc[N]
    Hx=(dec[1]-dec[0])/180.*np.pi*dc
    Hy=-(ra[1]-ra[0])/180.*np.pi*dc
    #print Hy
    Lx=Hx*len(dec)
    Ly=Hy*len(ra)
    print np.pi*2/Lx
    freq_x=np.fft.fftfreq(len(dec),1./len(dec))
    freq_y=np.fft.fftfreq(len(ra),1./len(ra))
    deltak=np.fft.fft2(data)
    Pk=np.abs(deltak)**2*(Lx*Ly/len(ra)**2/len(dec)**2)
    k=np.sqrt((2*np.pi/Lx)**2*freq_x[None,:]**2+(2*np.pi/Ly)**2*freq_y[:,None]**2)
    ################################################################################
    bin=10
    edges=np.linspace(k.min(),k.max(),bin+1,endpoint=True)
    n=np.histogram(k,edges)[0]
    print n
    k_bin=np.histogram(k,edges,weights=k)[0]
    print data.shape
    print k.shape
    print Pk.shape
    pk=np.histogram(k,edges,weights=Pk)[0]
    plt.semilogy(k_bin/n,pk/n,'.--',label=label)
f(data=data,label='data')
f(data=data1,label='ICs_2')
f(data=data2,label='ICs_3')
f(data=data3,label='ICs_5')
f(data=data4,label='ICs_10')
f(data=data5,label='ICs_20')
f(data=data6,label='ICs_50')
plt.legend()
plt.show()
