#!/usr/bin/env python
# coding=utf-8
from read import ReadMeta
import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM
from mpl_toolkits.mplot3d import Axes3D
Nx=128
Ny=64
Nz=32
#==== read data =========================
dataA_path = '/home/ycli/data/gbt/gbt_15hr_41-80_pointcorr/secA_15hr_41-80_pointcorr_clean_map_I_800.npy'
freq,ra,dec=ReadMeta(dataA_path)
data=np.load(dataA_path)
#data=np.load(dataA_path)[:,10:-10,5:-5]
#ra=ra[10:-10]
#dec=dec[5:-5]
data*=1000   # converse mK
#==== griding ===========================
# freq ->> x axis; ra ->> y axis; dec ->> z axis
h=67.77
shape=(freq.shape[0],ra.shape[0],dec.shape[0])
ra=ra-ra.mean()
dec=dec-dec.mean()
cosmo = FlatLambdaCDM(H0=h, Om0=0.307)
redshift=1420./(freq/10**6)-1
dc = np.array(cosmo.comoving_distance(redshift))*h/100        #z: Mpc/h
pos_z=dc[:,None,None]*np.sin(dec[None,None,:]/180.*np.pi)*np.ones(shape)
pos_y=dc[:,None,None]*np.sin(ra[None,:,None]/180.*np.pi)*np.ones(shape)
pos_x=dc[:,None,None]*np.ones(shape)
#==== to test the position ====
#fig = plt.figure()
#ax = Axes3D(fig)
#ax.plot(pos_x.reshape(-1),pos_y.reshape(-1),pos_z.reshape(-1),'.',alpha=0.1)
#ax.set_xlabel('X')
#ax.set_ylabel('y')
#ax.set_zlabel('z')
#plt.show()
#==============================
DPM=np.c_[pos_x.reshape(-1),pos_y.reshape(-1),pos_z.reshape(-1)]
def getPk(Nx=Nx,Ny=Ny,Nz=Nz):
    bin_x=np.linspace(pos_x.min()-10**-5,pos_x.max()+10**-5,Nx+1)
    bin_y=np.linspace(pos_y.min()-10**-5,pos_y.max()+10**-5,Ny+1)
    bin_z=np.linspace(pos_z.min()-10**-5,pos_z.max()+10**-5,Nz+1)
    deltaT=np.empty_like(data)
    for i in np.arange(freq.shape[0]):
        deltaT[i]=data[i]-data[i].mean()
    
    nn=np.histogramdd(DPM,bins=(bin_x,bin_y,bin_z))[0]
    T,edges=np.histogramdd(DPM,bins=(bin_x,bin_y,bin_z),weights=deltaT.reshape(-1))
    Hx=edges[0][1]-edges[0][0]
    Hy=edges[1][1]-edges[1][0]
    Hz=edges[2][1]-edges[2][0]
    print '='*60
    print 'Nx=%d\t\t'%Nx,'Ny=%d\t\t'%Ny,'Nz=%d\t\t'%Nz
    print 'Hx=%f\t'%Hx,'Hx=%f\t'%Hx,'Hx=%f\t'%Hx
    print '1', len(nn[nn!=0])
    print '0', len(nn[nn==0])
    ratio=np.float(len(nn[nn!=0]))/np.float(len(nn[nn==0])+len(nn[nn!=0]))
    print 'ratio:',ratio
    print '='*60
    
    #==== power spectrum estimator ==========
    k_x=np.fft.fftfreq(Nx,1./Nx)
    k_y=np.fft.fftfreq(Ny,1./Ny)
    k_z=np.fft.fftfreq(Nz,1./Nz)
    T_k=np.fft.fftn(T)
    Pk_auto=Hx*Hy*Hz/Nx/Ny/Nz*(np.abs(T_k)**2)
    k_mag=np.sqrt((2*np.pi/Hx/Nx*k_x[:,None,None])**2+(2*np.pi/Hy/Ny*k_y[None,:,None])**2+(2*np.pi/Hz/Nz*k_z[None,None,:])**2)
    Pk_auto*=k_mag**3/(2*np.pi**2)    # delta(k)^2
    #==== log bin 1d ========================
    bin=20
    k_mag[0,0,0]=10**-20
    edges=np.linspace(np.log10(2*10**-2),np.log10(k_mag.max()),bin+1,endpoint=True)
    n=np.histogram(np.log10(k_mag),edges)[0]
    k_bin=np.histogram(np.log10(k_mag),edges,weights=k_mag)[0]
    pk=np.histogram(np.log10(k_mag),edges,weights=Pk_auto)[0]
    n[n==0]=1
#   return k_bin/n,pk/n/ratio
    return k_bin/n,pk/n
#==== plot result =======================
print shape
k_plt,pk_plt=getPk(Nx=128,Ny=64,Nz=32)
plt.loglog(k_plt,pk_plt,'.-',label='256 64 32')
#plt.ylim([10**2,10**9])
plt.legend()
plt.show()
