#!/usr/bin/env python
# coding=utf-8
import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
from sklearn.decomposition import FastICA
from subprocess import call
import Ica as tt


def Get_Result(F, P):
    call('rm *.eps', shell=True)

    Plot = False
    pol = 0
    N = 4
    ############read data#####################
#   map1, map2 = tt.ICA.Filename(F, P)  # freq,pixel
    map1 = tt.ICA.ReadMap(map1)
    map2 = tt.ICA.ReadMap(map2)
    map = map1 + map2
    Freq_num = map.shape[0]
    # map=map-map.mean()
    S = map[:, pol].T
    del map
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
#    residuals=tt.ICA.RESULT(0,N,pol,A_,S_,map1,map2,Plot)
    res = []
#    return (map1[100,pol],map1[100,pol]+map2[100,pol]-re)
    ##########################################

    for i in range(Freq_num):
        res.append(tt.ICA.RESULT(i, N, pol, A_, S_, map1, map2, Plot))
        print res[-1]
        plt.close('all')
    res = np.array(res)
    resx = np.linspace(700, 800, Freq_num, endpoint=True)
#   plt.plot(resx,res,label='freq_%d pixel_%d'%(F,P))
#   plt.show()
    return np.c_[resx,res]
plt.figure(1)
a=Get_Result(200,128)
b=Get_Result(300,128)
plt.plot(a[:,0],a[:,1],label='freq_%d pixel_%d'%(200,128))
plt.plot(b[:,0],b[:,1],label='freq_%d pixel_%d'%(300,128))
plt.legend()
plt.show()
##########################################
#call('./write_tex.pyc', shell=True)
