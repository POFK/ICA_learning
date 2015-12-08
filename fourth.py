#!/usr/bin/env python
# coding=utf-8
import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
from sklearn.decomposition import FastICA
from subprocess import call
import Ica as tt
#####################################################
'''For data with beam.'''

def Get_Result(filename1,filename2,N=4,withbeam=True,powerspectrum=False,noise=False):
    '''filename1: 21cm
    filename2 :foreground or 21cm+fg+beam
    '''

#   call('rm *.eps', shell=True)

    Plot = False
    pol = 0
#   N = 4
#   withbeam=True
    ############read data#####################
    map1 = tt.ICA.ReadMap(filename1)
    map2 = tt.ICA.ReadMap(filename2)
    if withbeam==True: 
        map = map2
        Freq_num = map.shape[0]
        S = map.T
    else:
        map=map1[:,pol]+map2[:,pol]
        Freq_num = map.shape[0]
        S=map.T
#####################################################
    if noise==True:
        print map[0]
	map=map+10*np.random.normal(loc=1,scale=2,size=map.shape)
        print '#'*80
        print map[0]
    ############FastICA#######################
    ica = FastICA(
        n_components=N,
        algorithm='parallel',
        whiten=True,
        fun='logcosh',
        fun_args=None,
        max_iter=200,
        tol=0.0001,
        w_init=None,
        random_state=None)
    S_ = ica.fit_transform(S)
    A_ = ica.mixing_
    ##########################################
#    tt.ICA.GetComponent(N, S_)
#    re=tt.ICA.rebuild(N, A_, S_, 0, Plot=False)
#    residuals=tt.ICA.RESULT(0,N,pol,A_,S_,map1,map2,map3,Plot)
    res = []
#    return (map1[100,pol],map1[100,pol]+map2[100,pol]-re)
    ##########################################

    for i in range(Freq_num):
        res.append(tt.ICA.RESULT(i, N, pol, A_, S_, map1, map2, Plot, withbeam))
        print res[-1]
        plt.close('all')
    res = np.array(res)
    resx = np.linspace(700, 800, Freq_num, endpoint=True)
#   plt.plot(resx,res,label='freq_%d pixel_%d'%(F,P))
#   plt.show()
    def GetPs(freq):
	name='21cm_Syn_35beam_noise'
#       plt.close('all')
        re=tt.ICA.rebuild(N,A_,S_,freq,Plot)
        rebuild_21cm=map[freq,:]-re
        rebuild_21cm=rebuild_21cm-rebuild_21cm.mean()
        cl1=hp.anafast(map1[freq,pol])
        cl2=hp.anafast(rebuild_21cm)
        l=np.arange(len(cl1))
        plt.figure('APS')
        plt.semilogy(l,l*(l+1)*cl1,'r-',label='simulation')
        plt.semilogy(l,l*(l+1)*cl2,'b-',label='result')
        plt.xlabel('$l$')
        plt.ylabel('$l(l+1)C_l$')
	plt.title(name)
	plt.legend()
   	plt.savefig('./result/APS/'+name+'.eps')
    if powerspectrum==True:
        GetPs(0)    
    return np.c_[resx,res]
def add_sim_21cm_in_freq():
    map1 = tt.ICA.ReadMap('../data/ICA_sim/21cm_freq200_128.hdf5')
    x = np.linspace(700, 800, 200, endpoint=True)
    return (x,map1[:,0,0])
#a=Get_Result('../data/ICA_sim/21cm_freq200_128.hdf5','../data/ICA_sim/data_with_beam/D2021cm_fg_freq200_128.hdf5',withbeam=True,powerspectrum=True)
#a2=Get_Result('../data/ICA_sim/21cm_freq200_128.hdf5','../data/ICA_sim/data_with_beam/D4021cm_fg_freq200_128.hdf5',withbeam=True,powerspectrum=True)
#a3=Get_Result('../data/ICA_sim/21cm_freq300_128.hdf5','../data/ICA_sim/data_with_beam/D4021cm_fg_freq300_128.hdf5',withbeam=True,powerspectrum=False)
#b=Get_Result('../data/ICA_sim/21cm_freq200_128.hdf5','../data/ICA_sim/fg_freq200_128.hdf5',withbeam=True,powerspectrum=True,noise=True)
b2=Get_Result('../data/ICA_sim/21cm_freq200_128.hdf5','../data/ICA_sim/data_with_beam/_35arcmin_21cm_fg_freq200_128.hdf5',withbeam=True,powerspectrum=True,noise=True)
#b3=Get_Result('../data/ICA_sim/21cm_freq200_128.hdf5','../data/ICA_sim/data_with_beam/D4021cm_fg_freq200_128.hdf5',withbeam=True,powerspectrum=False)
####################################################################################
#plt.figure('200_128',figsize=(18,10))
#sim=add_sim_21cm_in_freq()
#plt.plot(sim[0],sim[1],label='sim')
####################################################################################
#plt.plot(a2[:,0],a2[:,1],label='freq_%d pixel_%d_withbeam'%(100,128))
#plt.plot(b3[:,0],b3[:,1],label='freq_%d pixel_%d_withbeam'%(200,128))
#plt.plot(a3[:,0],a3[:,1],label='freq_%d pixel_%d_withbeam'%(300,128))
#plt.plot(b[:,0],b[:,1],label='syn+ps+beam+noise')
#plt.plot(b2[:,0],b2[:,1],label='syn+ps+beam')
#plt.plot(b3[:,0],b3[:,1],label='syn+ps+beam')
#plt.xlabel('Frequency(MHz)')
#plt.ylabel('Residuals')
#plt.legend(loc='upper left')
#plt.show()
#plt.savefig('./result/compare_noise2.eps')
##########################################
#call('./write_tex.pyc', shell=True)
