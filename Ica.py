#!/usr/bin/env python
# coding=utf-8

import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
import h5py


class ICA():
    @classmethod
    def Filename(self,freq,pixel):
        self.freq=freq
        self.pixel=pixel
        self.filename_21cm='/home/mtx/data/ICA_sim/21cm_freq%d_%d.hdf5'%(freq,pixel)
        self.filename_fg='/home/mtx/data/ICA_sim/data_with_beam/21cm_fg_freq%d_%d.hdf5'%(freq,pixel)
        print '#'*80+ '\n#Open:\n\t'+self.filename_21cm+'\n\t'+self.filename_fg+'\n#\n'+'#'*80+'\n\n'
        return (self.filename_21cm,self.filename_fg)

    @classmethod
    def ReadMap(self, filename=''):
        self.f = h5py.File(filename)
        self.map = self.f['map'].value
        self.f.close()
        return self.map

    @classmethod
    def Pre(self,data):
        '''centering and whitening'''
        data=data-data.mean()
        data=data/data.var()
        return data

    @classmethod
    def rebuild(self, N, A_, S_, freq_num=0, Plot=True):
        self.s = 0
        for i in range(N):
            self.s = self.s + S_[:, i] * A_[freq_num, i]
        if Plot==True:
            hp.mollview(self.s, title='rebuild_freq_%d' % freq_num)
            plt.savefig('N_%d_rebuid_freq_%d.eps' % (N, freq_num))
            plt.cla()
        return self.s

    @classmethod
    def Residuals(self, a, b):
        '''calculation residuals of a and b, with sum |a-b|'''
        return np.abs(a - b).sum()/a.shape[0]

    @classmethod
    def GetComponent(self, N, S_):
        plt.clf()
        for i in range(N):
            hp.mollview(S_[:, i], title='component $%d$' % (i + 1))
            plt.savefig('N_%d_component%i.eps' % (N, i))
            plt.cla()

    @classmethod
    def RESULT(self, freq, *args):  # N=N,pol=pol,A_=A_,S=S_,a=map1,b=map3,c=map_withbeam,Plot=True):
        '''args: (N=N,pol=pol,A_=A_,S_=S_,a=map1,b=map2,c=map3,Plot=True,withbeam=True)'''
        (N, pol, A_, S_, a, b, Plot, withbeam) = args
        i = freq
        s = ICA.rebuild(N=N, A_=A_, S_=S_, freq_num=i,Plot=False)
        if Plot==True:
            hp.mollview(a[i, pol] + b[i, pol], title='simulation_all_freq%d' % i)
            plt.savefig('N_%d_simulation_freq%d.eps' % (N, i))
            plt.cla()
        if withbeam:
            result_21cm = b[i] - s
        else :
            result_21cm = a[i,pol]+b[i,pol]-s

        result_21cm = result_21cm - result_21cm.mean()#To average rebuilt 21cm signal 

        def PlotResult():
            hp.mollview(result_21cm, title='rebuild_21cm_freq%d' % i)
            plt.savefig('N_%d_rebuild_21cm_%d.eps' % (N, i))
            plt.cla()
            # plot rebuilt 21cm map, which has subtract its average value
            hp.mollview(
                (a[
                    i,
                    pol] -
                    result_21cm),
                title='residuals_21cm_freq%d' %
                i)
            plt.savefig('N_%d_residuals_21cm%d.eps' % (N, i))

            hp.mollview(a[i, pol], title='sim_21cm_freq%d' % i)
            plt.savefig('N_%d_sim_21cm_%d.eps' % (N, i))
            plt.cla()
            hp.mollview(b[i, pol], title='sim_foreground_%d' % i)
            plt.savefig('N_%d_sim_fg_%d.eps' % (N, i))
            plt.clf()

        if Plot==True:
            PlotResult()
        print ICA.Residuals(a[i, pol], result_21cm)
        return ICA.Residuals(a[i, pol], result_21cm)
