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
    def rebuild(self, N, A_, S_, freq_num=100):
        self.s = 0
        for i in range(N):
            self.s = self.s + S_[:, i] * A_[freq_num, i]
        hp.mollview(self.s, title='rebuild_freq_%d' % freq_num)
        plt.savefig('N_%d_rebuid_freq_%d.eps' % (N, freq_num))
        plt.cla()
        return self.s


map1 = ICA.ReadMap(
    filename='/home/mtx/data/ICA_sim/21cm_freq300_128.hdf5')
# map3=ICA.ReadMap(filename='/home/fileF/projects/python/test/ICA/data/fg_freq20.hdf5')
map3 = ICA.ReadMap(
    filename='/home/mtx/data/ICA_sim/galaxy_freq300_128.hdf5')
# map3=ICA.ReadMap(filename='/home/mtianxiang/desktop/proj/ICA/data/point_freq10.hdf5')
Freq_num = map1.shape[0]
##########################load data##################################
# hp.mollview(map3[0][0])
pol = 0
S = []
for i in range(Freq_num):
    S.append(map1[i][pol] + map3[i][pol])
S = np.array(S).T


##########################ICA components#############################
N = 7  # components number


ica = FastICA(n_components=N)
S_ = ica.fit_transform(S)
A_ = ica.mixing_

call('rm *.eps', shell=True)
########plot ###############################


def get_simpic(N):
    for i in range(N):
        hp.mollview(map1[i][pol])
        plt.title('$signal-21cm-%i$' % i)
        plt.savefig('N_%d_signal_21cm_%i.eps' % (N, i))
        plt.cla()
    for j in range(N):
        hp.mollview(map3[j][pol])
        plt.title('$signal-galaxy-%i$' % j)
        plt.savefig('N_%d_signal_galaxy_%i.eps' % (N, j))
        plt.cla()


def get_component():
    plt.clf()
    for i in range(N):
        #        hp.mollview(S_[:, i], title='$%d$' %
        #                    (i + 1), sub=((N / 3 + bool(N % 3)) * 100 + 30 + i + 1))
        #    plt.savefig('component.eps')
        #    plt.cla()
        hp.mollview(S_[:, i], title='component $%d$' % (i + 1))
        plt.savefig('N_%d_result%i.eps' % (N, i))
        plt.cla()


# get_simpic(N)
get_component()
for i in range(100, 101):
    s = ICA.rebuild(N=N, A_=A_, S_=S_, freq_num=i)
    hp.mollview(map1[i][pol] + map3[i][pol], title='simulation_all_freq%d' % i)
    plt.savefig('N_%d_simulation_freq%d.eps' % (N, i))
    plt.cla()
    a=(map1[i][pol] +map3[i][pol] -s)
    a=a-a.mean()
    hp.mollview(a,
        title='rebuild_21cm_freq%d' %
        i)
    plt.savefig('N_%d_rebuild_21cm_%d.eps' % (N, i))
    plt.cla()
    hp.mollview((map1[i][pol]-a), title='residuals_21cm_freq%d' %i)
    plt.savefig('N_%d_residuals_21cm%d.eps' % (N, i))
    hp.mollview(map1[i][pol], title='sim_21cm_freq%d' % i)
#    print a.mean(),map1[i][pol].mean()
    plt.savefig('N_%d_sim_21cm_%d.eps' % (N, i))
    plt.cla()
    hp.mollview(map3[i][pol], title='sim_galactic_syn%d' % i)
    plt.savefig('N_%d_sim_syn_%d.eps' % (N, i))
    plt.clf()
#print A_
#print A_.shape
#call('./write_tex.pyc', shell=True)
#print A_[0,:]
#print S_[:,1]
