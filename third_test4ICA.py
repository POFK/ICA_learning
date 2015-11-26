#!/usr/bin/env python
# coding=utf-8
import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import FastICA
from subprocess import call
import Ica as tt

call('rm *.eps', shell=True)

Plot=False
pol=0
N=5
############read data#####################
map1,map2=tt.ICA.Filename(100,128) #freq,pixel
map1=tt.ICA.ReadMap(map1)
map2=tt.ICA.ReadMap(map2)
map=map1+map2
Freq_num = map.shape[0]
#map=map-map.mean()
S=map[:,pol].T
del map
############FastICA#######################
ica=FastICA(n_components=N, algorithm='parallel', whiten=True, fun='logcosh', fun_args=None, max_iter=200, tol=0.0001, w_init=None, random_state=None)
S_=ica.fit_transform(S)
A_ = ica.mixing_
##########################################
tt.ICA.GetComponent(N,S_)
tt.ICA.rebuild(N,A_,S_,0,Plot=True)
residuals=tt.ICA.RESULT(0,N,pol,A_,S_,map1,map2,False)
##########################################
from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

res=[]
for i in range(rank*Freq_num/4,(rank+1)*Freq_num/4):
    res.append(tt.ICA.RESULT(i, N, pol, A_, S_, map1, map2,False))
    plt.close('all')

if rank!=0:
    comm.send(res,dest=0)
elif rank==0:
    a1=comm.recv(source=1)
    a2=comm.recv(source=2)
    a3=comm.recv(source=3)
    res=res+a1+a2+a3
    res=np.array(res)

    f=np.linspace(700,800,Freq_num,endpoint=True)
    p=np.c_[f,res]
#   np.savetxt('test_pixel_%d_Of_freq200'%Freq_num,p)

resx=np.linspace(700,800,Freq_num,endpoint=True)
plt.plot(resx,res)
plt.show()

###########################################
#del map1
#del map2
###########################################
#call('./write_tex.pyc', shell=True)
