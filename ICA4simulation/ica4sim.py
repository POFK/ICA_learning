#!/usr/bin/env python
# coding=utf-8
import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
from sklearn.decomposition import FastICA
from subprocess import call
import h5py


class ICA():

    @classmethod
    def ReadMap(self, filename=''):
        self.f = h5py.File(filename)
        self.map = self.f['map'].value
        self.f.close()
#       for self.i in range(len(self.map)):
#           self.fig = plt.figure(filename)
#           hp.mollview(self.map[self.i][0])
#           plt.savefig(filename + '_freq_%d.png' % self.i)
#           plt.cla()
        return self.map

    @classmethod
    def ReadData(self, filename=''):
        self.f = h5py.File(filename)
        self.map = self.f['map'].value
        self.f.close()
        return self.map

    @classmethod
    def rebuild(self, N, A_, S_, freq_num=0):
        self.s = 0
        for i in np.arange(N):
            self.s = self.s + S_[:, i] * A_[freq_num, i]
        hp.mollview(self.s,min=self.s.mean()-self.s.std(),max=self.s.mean()+self.s.std(), title='rebuild_freq_%d' % freq_num)
        plt.savefig('N_%d_rebuid_freq_%d.eps' % (N, freq_num))
        plt.cla()
        return self.s

    @classmethod
    def Residuals(self, a, b):
        '''calculation residuals of a and b, with sum |a-b|'''
        return np.abs(a - b).sum()/a.shape[0]

    @classmethod
    def GetComponent(self, N):
        plt.clf()
        for i in range(N):
            mean=S_[:,i].mean()
            std=S_[:,i].std()
            hp.mollview(S_[:, i],min=mean-std,max=mean+std,title='component $%d$' % (i + 1))
            plt.savefig('N_%d_result%i.eps' % (N, i))
            plt.cla()

    @classmethod
    def RESULT(self, freq, *args):  # N=N,pol=pol,A_=A_,S=S_,a=map1,b=map3):
        '''args: (N=N,pol=pol,A_=A_,S_=S_,a=map1,b=map3)'''
        (N, pol, A_, S_, a, b) = args
        plt.clf()
        i = freq
        s = ICA.rebuild(N=N, A_=A_, S_=S_, freq_num=i)

        plt_data=a[i, pol] + b[i, pol]
        plt_mean=plt_data.mean()
        plt_std=plt_data.std()
        hp.mollview(plt_data,min=plt_mean-plt_std,max=plt_mean+plt_std,title='simulation_all_freq%d' % i)
        plt.savefig('N_%d_simulation_freq%d.eps' % (N, i))
        plt.cla()

        result_21cm = a[i, pol] + b[i, pol] - s
        result_21cm = result_21cm - result_21cm.mean()

        def PlotResult():
            plt_data=result_21cm
            plt_mean=plt_data.mean()
            plt_std=plt_data.std()
            hp.mollview(plt_data,min=plt_mean-plt_std,max=plt_mean+plt_std, title='rebuild_21cm_freq%d' % i)
            plt.savefig('N_%d_rebuild_21cm_%d.eps' % (N, i))
            plt.cla()
            # plot rebuilt 21cm map, which has subtract its average value
            plt_data=map1[i,pol]-result_21cm
            plt_mean=plt_data.mean()
            plt_std=plt_data.std()
            hp.mollview(plt_data,min=plt_mean-plt_std,max=plt_mean+plt_std, title='residuals_21cm_freq%d' %i)
            plt.savefig('N_%d_residuals_21cm%d.eps' % (N, i))

            plt_data=map1[i,pol]
            plt_mean=plt_data.mean()
            plt_std=plt_data.std()
            hp.mollview(plt_data,min=plt_mean-plt_std,max=plt_mean+plt_std,title='sim_21cm_freq%d' % i)
            plt.savefig('N_%d_sim_21cm_%d.eps' % (N, i))
            plt.cla()

            plt_data=map3[i,pol]
            plt_mean=plt_data.mean()
            plt_std=plt_data.std()
            hp.mollview(plt_data,min=plt_mean-plt_std,max=plt_mean+plt_std,title='sim_syn_%d' % i)
            plt.savefig('N_%d_sim_syn_%d.eps' % (N, i))
            plt.clf()

        PlotResult()
        print ICA.Residuals(map1[i, pol], result_21cm)
        return ICA.Residuals(map1[i, pol], result_21cm)

#########################read file##################################

map1 = ICA.ReadMap(
    filename='/home/mtx/data/ICA_sim/21cm_freq200_128.hdf5')
#map3 = ICA.ReadMap(filename='/home/mtx/data/ICA_sim/galaxy_freq200_128.hdf5')
map3 = 1*ICA.ReadMap(filename='/home/mtx/data/ICA_sim/fg_freq200_128.hdf5')

Freq_num = map1.shape[0]
##########################load data##################################
# hp.mollview(map3[0][0])
pol = 0
S = []
for i in range(Freq_num):
    S.append(map1[i, pol] + map3[i,pol])
S = np.array(S).T
##########################ICA components#############################
N = 4  # components number
#####################################################################
ica = FastICA(n_components=N)
S_ = ica.fit_transform(S)
A_ = ica.mixing_
#################################plot ###############################
call('rm *.eps', shell=True)


def get_simpic(N):
    for i in range(N):
        hp.mollview(map1[i, pol])
        plt.title('$signal-21cm-%i$' % i)
        plt.savefig('N_%d_signal_21cm_%i.eps' % (N, i))
        plt.cla()
    for j in range(N):
        hp.mollview(map3[j, pol])
        plt.title('$signal-galaxy-%i$' % j)
        plt.savefig('N_%d_signal_galaxy_%i.eps' % (N, j))
        plt.cla()

plt.close('all')

ICA.GetComponent(N=N)
ICA.RESULT(100,N, pol, A_, S_, map1, map3)

#call('./write_tex.pyc', shell=True)
#call('cp result.pdf 0.000001.pdf', shell=True)
#call('cp result.pdf result_component_%d.pdf'%N,shell=True)
