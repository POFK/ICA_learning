#!/usr/bin/env python
# coding=utf-8
from read import ReadMeta
import numpy as np
from astropy.cosmology import FlatLambdaCDM
import matplotlib.pyplot as plt

h=67.77
freq,ra,dec=ReadMeta('/home/ycli/data/wigglez/gbt_15hr_41-80_pointcorr/reg15data.npy')
cosmo = FlatLambdaCDM(H0=h, Om0=0.308)
data_wigglez1=np.load('/home/ycli/data/wigglez/gbt_15hr_41-80_pointcorr/reg15data.npy')
Nx,Ny,Nz=64,64,256
################################################################################
def g(Nx,Ny,Nz):
    redshift=1420./(freq/10**6)-1
    dc = np.array(cosmo.comoving_distance(redshift))*h/100        #z: Mpc/h
    pos_x=dc[:,None,None]*np.sin(np.pi/180.*dec[None,None,:])*np.cos(np.pi/180.*ra[None,:,None])
    pos_y=dc[:,None,None]*np.sin(np.pi/180.*dec[None,None,:])*np.sin(np.pi/180.*ra[None,:,None])
    pos_z=dc[:,None,None]*np.cos(np.pi/180.*ra[None,:,None])+np.zeros_like(dec[None,None,:])
    print pos_x.shape,pos_y.shape,pos_z.shape
    print pos_x.max(),pos_x.min()
    print pos_y.max(),pos_y.min()
    print pos_z.max(),pos_z.min()
    DPM=np.c_[pos_x.reshape(-1),pos_y.reshape(-1),pos_z.reshape(-1)]
    bin_x=np.linspace(pos_x.min()-10**-5,pos_x.max()+10**-5,Nx+1)
    bin_y=np.linspace(pos_y.min()-10**-5,pos_y.max()+10**-5,Ny+1)
    bin_z=np.linspace(pos_z.min()-10**-5,pos_z.max()+10**-5,Nz+1)
    def f(data,data_wigglez,label):
        n,edges=np.histogramdd(DPM,bins=(bin_x,bin_y,bin_z))
        n[n==0]=1
        T=np.histogramdd(DPM,bins=(bin_x,bin_y,bin_z),weights=data.reshape(-1))[0]
        T_wigglez=np.histogramdd(DPM,bins=(bin_x,bin_y,bin_z),weights=data_wigglez.reshape(-1))[0]
    
        T/=n
        T_wigglez/=n
#       T*=Nx*Ny*Nz/T.sum() 
#       T_wigglez*=Nx*Ny*Nz/T_wigglez.sum() 

        Hx=np.abs(edges[0][1]-edges[0][0])
        Hy=np.abs(edges[1][1]-edges[1][0])
        Hz=np.abs(edges[2][1]-edges[2][0])
        k_x=np.fft.fftfreq(Nx,1./Nx)
        k_y=np.fft.fftfreq(Ny,1./Ny)
        k_z=np.fft.fftfreq(Nz,1./Nz)
        T_k1=np.fft.fftn(T)
        T_k2=np.fft.fftn(T_wigglez)
#       window_k = np.sinc(1. / Nx * k_x[:,None,None]) * np.sinc(1. / Ny * k_y[None,:,None]) * np.sinc(1. / Nz * k_z[None,None,:])
#       T_k1/=window_k
#       T_k2/=window_k
        

        Pk_cro=Hx*Hy*Hz/Nx/Ny/Nz*((T_k1.conjugate()*T_k2+T_k2.conjugate()*T_k1)/2).real
        Pk_gbt=Hx*Hy*Hz/Nx/Ny/Nz*(np.abs(T_k1)**2)
        Pk_wig=Hx*Hy*Hz/Nx/Ny/Nz*(np.abs(T_k2)**2)
        k_mag=np.sqrt((2*np.pi/Hx/Nx*k_x[:,None,None])**2+(2*np.pi/Hy/Ny*k_y[None,:,None])**2+(2*np.pi/Hz/Nz*k_z[None,None,:])**2)
        Pk_cro*=(k_mag**3.)/(2*np.pi**2)    # delta(k)^2
        Pk_gbt*=(k_mag**3.)/(2*np.pi**2)    # delta(k)^2
        Pk_wig*=(k_mag**3.)/(2*np.pi**2)    # delta(k)^2
        
        #################################################################################

        bin=50
        edges=np.linspace(k_mag.min(),k_mag.max(),bin+1,endpoint=True)
        n=np.histogram(k_mag,edges)[0]
        k_bin=np.histogram(k_mag,edges,weights=k_mag)[0]
        pk_cro=np.histogram(k_mag,edges,weights=Pk_cro)[0]
        pk_gbt=np.histogram(k_mag,edges,weights=Pk_gbt)[0]
        pk_wig=np.histogram(k_mag,edges,weights=Pk_wig)[0]
        plt.figure('Pk')
        plt.semilogy(k_bin/n,pk_cro/n,'x--',label=label+'_cro')
#       plt.semilogy(k_bin/n,pk_gbt/n,'.-',label=label+'_gbt')
#       plt.semilogy(k_bin/n,pk_wig/n,'o-',label=label+'_wig')
        plt.ylabel('$\Delta ^2(k)$')
        plt.xlabel('$k$')
    
        plt.legend()

        #################################################################################
#       wk=np.histogram(k_mag,edges,weights=window_k)[0]
#       plt.figure('window')
#       plt.plot(k_bin/n,wk/n,'.-',label=str(Nx)+' '+str(Ny)+' '+str(Nz))
#       plt.legend()
        #################################################################################
    f(data=data_wigglez1,data_wigglez=data_wigglez1,label=str(Nx)+' '+str(Ny)+' '+str(Nz))

#g(Nx=32,Ny=32,Nz=32)
#g(Nx=32,Ny=64,Nz=32)
#g(Nx=32,Ny=32,Nz=128)
#g(Nx=64,Ny=64,Nz=64)
#g(Nx=128,Ny=128,Nz=128)
plt.legend()
#plt.xscale('log')
plt.show()
