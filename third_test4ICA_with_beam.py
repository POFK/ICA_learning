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

def Get_Result(filename1,filename2,withbeam=True):
    '''filename1: 21cm
    filename2 :foreground or 21cm+fg+beam
    '''

    call('rm *.eps', shell=True)

    Plot = False
    pol = 0
    N = 4
#   withbeam=True
    ############read data#####################
    map1 = tt.ICA.ReadMap(filename1)
    map2 = tt.ICA.ReadMap(filename2)
    if withbeam==True: 
        map = map2
        Freq_num = map.shape[0]
        S = map.T
        del map
    else:
        map=map1[:,pol]+map2[:,pol]
        Freq_num = map.shape[0]
        S=map.T
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
    return np.c_[resx,res]
plt.figure(1)
a=Get_Result('../data/ICA_sim/21cm_freq200_128.hdf5','../data/ICA_sim/data_with_beam/D2021cm_fg_freq200_128.hdf5',withbeam=True)
b=Get_Result('../data/ICA_sim/21cm_freq200_128.hdf5','../data/ICA_sim/fg_freq200_128.hdf5',withbeam=False)
a2=Get_Result('../data/ICA_sim/21cm_freq100_128.hdf5','../data/ICA_sim/data_with_beam/D2021cm_fg_freq100_128.hdf5',withbeam=True)
b2=Get_Result('../data/ICA_sim/21cm_freq100_128.hdf5','../data/ICA_sim/fg_freq100_128.hdf5',withbeam=False)
a3=Get_Result('../data/ICA_sim/21cm_freq200_256.hdf5','../data/ICA_sim/data_with_beam/D2021cm_fg_freq200_256.hdf5',withbeam=True)
b3=Get_Result('../data/ICA_sim/21cm_freq200_256.hdf5','../data/ICA_sim/fg_freq200_256.hdf5',withbeam=False)
plt.plot(a[:,0],a[:,1],label='freq_%d pixel_%d_withbeam'%(200,128))
plt.plot(b[:,0],b[:,1],label='freq_%d pixel_%d'%(200,128))
plt.plot(a2[:,0],a2[:,1],label='freq_%d pixel_%d_withbeam'%(100,128))
plt.plot(b2[:,0],b2[:,1],label='freq_%d pixel_%d'%(100,128))
plt.plot(a3[:,0],a3[:,1],label='freq_%d pixel_%d_withbeam'%(200,256))
plt.plot(b3[:,0],b3[:,1],label='freq_%d pixel_%d'%(200,256))

plt.legend()
plt.show()
##########################################
#call('./write_tex.pyc', shell=True)
