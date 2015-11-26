#!/usr/bin/env python
# coding=utf-8
from second_test4ICA import *
from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

res=[]
for i in range(rank*Freq_num/4,(rank+1)*Freq_num/4):
    res.append(ICA.RESULT(i, N, pol, A_, S_, map1, map3))
    plt.close('all')

if rank!=0:
    comm.send(res,dest=0)
elif rank==0:
    a1=comm.recv(source=1)
    a2=comm.recv(source=2)
    a3=comm.recv(source=3)
#   a4=comm.recv(source=4)
#   a5=comm.recv(source=5)
#   a6=comm.recv(source=6)
#   a7=comm.recv(source=7)
    res=res+a1+a2+a3#+a4+a5+a6+a7
    res=np.array(res)

    f=np.linspace(700,800,Freq_num,endpoint=True)
    p=np.c_[f,res]
    np.savetxt('Freq_%d_withFG'%Freq_num,p)

#   plt.clf()
#   plt.plot(f,res)
#   plt.show()

