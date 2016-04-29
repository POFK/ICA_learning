#!/usr/bin/env python
# coding=utf-8
from read import ReadMeta
from sklearn.decomposition import FastICA
import numpy as np
from astropy.cosmology import FlatLambdaCDM
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
class ICA_GBT():    
    def __init__(self,freq,ra,dec,N=5):
        self.N=N  # n_components of ICA
        self.freq,self.ra,self.dec=freq,ra,dec
        self.h=67.77
        self.dataA_path = '/home/ycli/data/gbt/gbt_15hr_41-80_pointcorr/secA_15hr_41-80_pointcorr_clean_map_I_800.npy'
        self.dataB_path = '/home/ycli/data/gbt/gbt_15hr_41-80_pointcorr/secB_15hr_41-80_pointcorr_clean_map_I_800.npy'
        self.dataC_path = '/home/ycli/data/gbt/gbt_15hr_41-80_pointcorr/secC_15hr_41-80_pointcorr_clean_map_I_800.npy'
        self.dataD_path = '/home/ycli/data/gbt/gbt_15hr_41-80_pointcorr/secD_15hr_41-80_pointcorr_clean_map_I_800.npy'
        self.data_wigglez_path='/home/ycli/data/wigglez/gbt_15hr_41-80_pointcorr/reg15data.npy'
    def LoadData(self,PathData,PathWig):
        self.data_wigglez=np.load(PathWig)
        self.data=np.load(PathData)
        return self.data
    def Clean(self,data,label='test',freq=20,sigma=1,PLOT=True):
        shape=(data.shape[1],data.shape[2])
        S=data.reshape((data.shape[0],-1)).T   #model 2
        ica = FastICA(n_components=self.N)
        S_ = ica.fit_transform(S)
        A_ = ica.mixing_
        ICs_sum=np.dot(S_,A_.T)
        data_clean=(S-ICs_sum).T.reshape(data.shape)
        limy,limx=data[freq].shape
        if PLOT== True:
            plt.figure(label)
            plt.subplot(1,3,1)
#           plt.imshow(data[freq],extent=[self.dec[0],self.dec[-1],self.ra[0],self.ra[-1]])
            plt.pcolor(data[freq])
            plt.colorbar()
            plt.title('map')
            plt.xlim([0,limx])
            plt.ylim([0,limy])
            plt.xlabel('dec')
            plt.ylabel('ra')
            plt.subplot(1,3,2)
#           plt.imshow((ICs_sum).T.reshape(data.shape)[freq],extent=[self.dec[0],self.dec[-1],self.ra[0],self.ra[-1]])
            plt.pcolor((ICs_sum).T.reshape(data.shape)[freq])
            plt.colorbar()
            plt.title('ICs')
            plt.xlim([0,limx])
            plt.ylim([0,limy])
            plt.xlabel('dec')
            plt.ylabel('ra')
            plt.subplot(1,3,3)
            mean_=data_clean[freq].mean()
            std_=data_clean[freq].std()
#           plt.imshow(data_clean[freq],vmin=mean_-std_*sigma,vmax=mean_+std_*sigma,extent=[self.dec[0],self.dec[-1],self.ra[0],self.ra[-1]])
            plt.pcolor(data_clean[freq]-mean_,vmin=-std_*sigma,vmax=std_*sigma)
            plt.colorbar()
            plt.title('cleaned map')
            plt.xlim([0,limx])
            plt.ylim([0,limy])
            plt.xlabel('dec')
            plt.ylabel('ra')
        return data_clean
    def Grid(self,Nx=64,Ny=64,Nz=256):
        print 'griding',Nx,Ny,Nz
        h=self.h
        ra=self.ra
        ra=ra-ra.mean()
        dec=self.dec
        dec=dec-dec.mean()
        self.Nx=Nx
        self.Ny=Ny
        self.Nz=Nz
        cosmo = FlatLambdaCDM(H0=h, Om0=0.307)
        redshift=1420./(self.freq/10**6)-1
        dc = np.array(cosmo.comoving_distance(redshift))*h/100        #z: Mpc/h
        pos_x=dc[:,None,None]*np.cos(np.pi/180.*dec[None,None,:])*np.cos(np.pi/180.*ra[None,:,None])
        pos_y=dc[:,None,None]*np.cos(np.pi/180.*dec[None,None,:])*np.sin(np.pi/180.*ra[None,:,None])
        pos_z=dc[:,None,None]*np.sin(np.pi/180.*dec[None,None,:])+np.zeros_like(ra[None,:,None])
#       print pos_x.max(),pos_x.min()
#       print pos_y.max(),pos_y.min()
#       print pos_z.max(),pos_z.min()

        #put redshift direction on z axis. and before this, dec>>z ra>>y redshift>>x
        #after transformation, ra>>x dec>>y redshift>>z
        self.DPM=np.c_[pos_y.reshape(-1),pos_z.reshape(-1),pos_x.reshape(-1)]
        self.bin_x=np.linspace(pos_y.min()-10**-5,pos_y.max()+10**-5,Nx+1)
        self.bin_y=np.linspace(pos_z.min()-10**-5,pos_z.max()+10**-5,Ny+1)
        self.bin_z=np.linspace(pos_x.min()-10**-5,pos_x.max()+10**-5,Nz+1)
    def Pk(self,data1,data2,bin=30):
        DPM=self.DPM
        bin_x=self.bin_x
        bin_y=self.bin_y
        bin_z=self.bin_z
        Nx=self.Nx
        Ny=self.Ny
        Nz=self.Nz
        nn=np.histogramdd(DPM,bins=(bin_x,bin_y,bin_z))[0]
        T1,edges=np.histogramdd(DPM,bins=(bin_x,bin_y,bin_z),weights=data1.reshape(-1))
        T2=np.histogramdd(DPM,bins=(bin_x,bin_y,bin_z),weights=data2.reshape(-1))[0]
        Hx=edges[0][1]-edges[0][0]
        Hy=edges[1][1]-edges[1][0]
        Hz=edges[2][1]-edges[2][0]
        print Hx,Hy,Hz

#       T1*=Nx*Ny*Nz/T1.sum()
#       T2*=Nx*Ny*Nz/T2.sum()

#       nbar=(nn!=0).sum()/(Hx*Nx*Hy*Ny*Hz*Nz)
#       print 'n!=0',nbar*(Hx*Nx*Hy*Ny*Hz*Nz)

        nn[nn==0]=1
        T1/=nn
        T2/=nn

#       plt.figure('0')
#       plt.pcolor(T1[10,:,:])
#       plt.figure('1')
#       plt.pcolor(T1[:,10,:])
#       plt.figure('2')
#       plt.pcolor(T1[:,:,10])
#       plt.show()

#       fig=plt.figure()
#       ss=np.meshgrid(bin_x,bin_y,bin_z)
#       ax=fig.add_subplot(111,projection='3d')
#       ax.scatter3D(DPM[:,0],DPM[:,1],DPM[:,2],s=1,c='r')
#       plt.show()



#       T1/=Hx*Hy*Hz
#       T2/=Hx*Hy*Hz


        k_x=np.fft.fftfreq(Nx,1./Nx)
        k_y=np.fft.fftfreq(Ny,1./Ny)
        k_z=np.fft.fftfreq(Nz,1./Nz)

        T_k1=np.fft.fftn(T1)
        T_k2=np.fft.fftn(T2)
#       self.window_k = np.sinc(1. / Nx * k_x[:,None,None]) * np.sinc(1. / Ny * k_y[None,:,None]) * np.sinc(1. / Nz* k_z[None,None,:])

#       T_k1/=self.window_k
#       T_k2/=self.window_k

        Pk_cro=Hx*Hy*Hz/Nx/Ny/Nz*((T_k1.conjugate()*T_k2+T_k2.conjugate()*T_k1)/2).real
        Pk_1=Hx*Hy*Hz/Nx/Ny/Nz*(np.abs(T_k1)**2)
        Pk_2=Hx*Hy*Hz/Nx/Ny/Nz*(np.abs(T_k2)**2)


        k_mag=np.sqrt((2*np.pi/Hx/Nx*k_x[:,None,None])**2+(2*np.pi/Hy/Ny*k_y[None,:,None])**2+(2*np.pi/Hz/Nz*k_z[None,None,:])**2)
        Pk_cro*=k_mag**3/(2*np.pi**2)    # delta(k)^2
        Pk_1*=k_mag**3/(2*np.pi**2)    # delta(k)^2
        Pk_2*=k_mag**3/(2*np.pi**2)    # delta(k)^2
        
        #################################################################################
#       edges=np.linspace(k_mag.min(),k_mag.max(),bin+1,endpoint=True)
#       n=np.histogram(k_mag,edges)[0]
#       k_bin=np.histogram(k_mag,edges,weights=k_mag)[0]
#       pk_cro=np.histogram(k_mag,edges,weights=Pk_cro)[0]
#       pk_1=np.histogram(k_mag,edges,weights=Pk_1)[0]
#       pk_2=np.histogram(k_mag,edges,weights=Pk_2)[0]
##############log bin ##########################################################
        k_mag[0,0,0]=10**-20
        edges=np.linspace(np.log10(2*10**-2),np.log10(k_mag.max()),bin+1,endpoint=True)
        n=np.histogram(np.log10(k_mag),edges)[0]
        k_bin=np.histogram(np.log10(k_mag),edges,weights=k_mag)[0]
        pk_cro=np.histogram(np.log10(k_mag),edges,weights=Pk_cro)[0]
        pk_1=np.histogram(np.log10(k_mag),edges,weights=Pk_1)[0]
        pk_2=np.histogram(np.log10(k_mag),edges,weights=Pk_2)[0]
        n[n==0]=1
        return k_bin/n,pk_1/n,pk_2/n,pk_cro/n
