#!/usr/bin/env python
# coding=utf-8
from IM import IM
from read import ReadMeta
import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM
from read_par import readPar
#from mpl_toolkits.mplot3d import Axes3D
class getGridPk(IM):
    def Pos(self,freq=None,ra=None,dec=None):
        # freq ->> x axis; ra ->> y axis; dec ->> z axis
        self.shape=(freq.shape[0],ra.shape[0],dec.shape[0])
        cosmo = FlatLambdaCDM(H0=self.h, Om0=self.Om0)
        redshift=1420./(freq/10**6)-1
        dc = np.array(cosmo.comoving_distance(redshift))*self.h/100        #z: Mpc/h
        self.pos_z=dc[:,None,None]*np.sin(dec[None,None,:]/180.*np.pi)*np.ones(self.shape)
        self.pos_y=dc[:,None,None]*np.sin(ra[None,:,None]/180.*np.pi)*np.ones(self.shape)
        self.pos_x=dc[:,None,None]*np.ones(self.shape)
        DPM=np.c_[self.pos_x.reshape(-1),self.pos_y.reshape(-1),self.pos_z.reshape(-1)]
        bin_x=np.linspace(self.pos_x.min()-10**-5,self.pos_x.max()+10**-5,self.Nx+1)
        bin_y=np.linspace(self.pos_y.min()-10**-5,self.pos_y.max()+10**-5,self.Ny+1)
        bin_z=np.linspace(self.pos_z.min()-10**-5,self.pos_z.max()+10**-5,self.Nz+1)
        return DPM,bin_x,bin_y,bin_z
    def grid(self,freq=None,ra=None,dec=None):
        print 'gridding ...'
        self.DPM,self.bin_x,self.bin_y,self.bin_z=self.Pos(freq=freq,ra=ra,dec=dec)
        nn,edges=np.histogramdd(self.DPM,bins=(self.bin_x,self.bin_y,self.bin_z))
        self.Hx=edges[0][1]-edges[0][0]
        self.Hy=edges[1][1]-edges[1][0]
        self.Hz=edges[2][1]-edges[2][0]
        print 'Hx=%f\t'%self.Hx,'Hy=%f\t'%self.Hy,'Hz=%f\t'%self.Hz
        print '1', len(nn[nn!=0])
        print '0', len(nn[nn==0])
        ratio=np.float(len(nn[nn!=0]))/np.float(len(nn[nn==0])+len(nn[nn!=0]))
        print 'ratio:',ratio
        print 'finish gridding'
    def grid_21cm(self,data=None):
        deltaT=np.empty_like(data)
        for i in np.arange(freq.shape[0]):
            deltaT[i]=data[i]-data[i].mean()
        nn=np.histogramdd(self.DPM,bins=(self.bin_x,self.bin_y,self.bin_z))[0]
        T,edges=np.histogramdd(self.DPM,bins=(self.bin_x,self.bin_y,self.bin_z),weights=deltaT.reshape(-1))
        return T
    def grid_wig(self,data=None):
        deltax=data/data.mean()
        nn=np.histogramdd(self.DPM,bins=(self.bin_x,self.bin_y,self.bin_z))[0]
        delta,edges=np.histogramdd(self.DPM,bins=(self.bin_x,self.bin_y,self.bin_z),weights=deltax.reshape(-1))
        return delta
    def getPk(self,T1=None,T2=None):
        Hx=self.Hx
        Hy=self.Hy
        Hz=self.Hz
        Nx=self.Nx
        Ny=self.Ny
        Nz=self.Nz
        k_x=np.fft.fftfreq(int(Nx),1./Nx)
        k_y=np.fft.fftfreq(int(Ny),1./Ny)
        k_z=np.fft.fftfreq(int(Nz),1./Nz)
        T_k1=np.fft.fftn(T1)
        T_k2=np.fft.fftn(T2)
        Pk=Hx*Hy*Hz/Nx/Ny/Nz*((T_k1.conjugate()*T_k2+T_k2.conjugate()*T_k1)/2).real
        #for 1D bins
        self.k_mag=np.sqrt((2*np.pi/Hx/Nx*k_x[:,None,None])**2+(2*np.pi/Hy/Ny*k_y[None,:,None])**2+(2*np.pi/Hz/Nz*k_z[None,None,:])**2)
        #for 2D bins
        k_V=np.sqrt(np.zeros_like(2*np.pi/Hx/Nx*k_x[:,None,None])**2+(2*np.pi/Hy/Ny*k_y[None,:,None])**2+(2*np.pi/Hz/Nz*k_z[None,None,:])**2) # k in y-z plane
        k_P=np.sqrt((2*np.pi/Hx/Nx*k_x[:,None,None])**2+np.zeros_like(2*np.pi/Hy/Ny*k_y[None,:,None])**2+np.zeros_like(2*np.pi/Hz/Nz*k_z[None,None,:])**2) # k in x axis
        self.k_2d=np.c_[k_P.reshape(-1),k_V.reshape(-1)]

        self.Pk=Pk*self.k_mag**3/(2*np.pi**2)    # delta(k)^2=k^3/2pi*Pk
        return Pk
    def getBin1D(self):
        bins=self.bins
        self.k_mag[0,0,0]=10**-20
        minK=np.min([2*np.pi/self.Nx/self.Hx,2*np.pi/self.Ny/self.Hy,2*np.pi/self.Nz/self.Hz])  # to find the minimal k
        edges=np.linspace(np.log10(minK),np.log10(self.k_mag.max()+10**-5),bins+1,endpoint=True)
        n=np.histogram(np.log10(self.k_mag),edges)[0]
        k_bin=np.histogram(np.log10(self.k_mag),edges,weights=self.k_mag)[0]
        pk_bin=np.histogram(np.log10(self.k_mag),edges,weights=self.Pk)[0]
        n[n==0]=1
        return k_bin/n,pk_bin/n,n
    def getBin2D(self):
        print 'getting 2d bins:'
        bins=self.bins2d
        Pk2d=self.Pk
        Pk2d[0,0,0]=0.
        minK=np.min([2*np.pi/self.Nx/self.Hx,2*np.pi/self.Ny/self.Hy,2*np.pi/self.Nz/self.Hz])  # to find the minimal k
        edges=np.linspace(np.log10(minK),np.log10(self.k_mag.max()+10**-5),bins+1,endpoint=True)
        edges=10**edges
        edges[0]=0.0
        n=np.histogramdd(self.k_2d,bins=(edges,edges))[0]
        pk_bin=np.histogramdd(self.k_2d,bins=(edges,edges),weights=self.Pk.reshape(-1))[0]
        pk_bin[0,0]=pk_bin[0,0]*n[0,0]/(n[0,0]-1)
        n[0,0]=n[0,0]-1
        print 'edges2d:',edges
        print 'n:',n
        n[n==0]=1
        return edges,pk_bin/n

################################################################################
# to test
if (__name__=='__main__'):
    fs=getGridPk('Fg.par')
    data=np.load(fs.PathWig)
    freq,ra,dec=ReadMeta(fs.Path1)
    #data=np.load(dataA_path)[:,10:-10,5:-5]
    #ra=ra[10:-10]
    #dec=dec[5:-5]
    data*=1000   # converse mK
    shape=(freq.shape[0],ra.shape[0],dec.shape[0])
    ra=ra-ra.mean()
    dec=dec-dec.mean()
    fs.grid(freq=freq,ra=ra,dec=dec)
#   T=fs.grid_21cm(data=data)
    T=fs.grid_wig(data=data)
    fs.getPk(T1=T,T2=T)
    k,pk,n=fs.getBin1D()
    edges2d,pk2d=fs.getBin2D()
    plt.figure('1D')
    plt.loglog(k,pk,'.-',label='256 64 32')
    #plt.ylim([10**2,10**9])
    plt.legend()
    plt.figure('2D')
    plt.pcolor(pk2d.T)
    plt.colorbar()
    plt.xlabel('k_p')
    plt.ylabel('k_v')
    plt.show()
