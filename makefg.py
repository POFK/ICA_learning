#!/usr/bin/env python
# coding=utf-8
import numpy as np
import matplotlib.pyplot as plt
import scipy.special as special
import mpi4py.MPI as MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
theta = np.linspace(0, 30/180.*np.pi, 256,endpoint=False)+30/180.*np.pi/256/2
phi = np.linspace(0.,30/180.*np.pi, 256,endpoint=False)+30/180.*np.pi/256/2
X,Y=np.meshgrid(theta,phi)
#l=range(1,128)
l=range(1,2)
def delta(theta,phi):
    s=0
    def fg(l, v, v2, A=700., vf=130. * 10**6, e=4.0):
        return A * (1000. / l)**2.4 * (vf**2. / (v * v2))**(2 * 2.8) * np.exp(-(np.log(v / v2))**2 / (2 * e**2))
    for i in l:
        m=range(-i,i+1)
        Cl=fg(i,830.78,829.22)
        for j in m:
            Ylm=abs(special.sph_harm(n=i,m=j,theta=theta,phi=phi))
            s=s+Ylm**2*Cl
    return s
I=[]
for i in range(64*rank,(rank+1)*64):
#for i in range(len(theta)):
    S=[]
    for j in range(len(phi)):
        S.append(delta(theta[i],phi[j]))
    I.append(S)
if rank==0:
    I1=comm.recv(source=1)
    I2=comm.recv(source=2)
    I3=comm.recv(source=3)
    out =[]
    for i in I:
        out.append(i)
    for i in I1:
        out.append(i)
    for i in I2:
        out.append(i)
    for i in I3:
        out.append(i)
    print len(out),len(out[0])
    plt.pcolor(X,Y,np.array(out))
    plt.colorbar()
    plt.xlim(np.min(theta),np.max(theta))
    plt.xlim(np.min(phi),np.max(phi))
    #plt.savefig('./fg.png')
    plt.show()
    #plt.cla()
    np.savetxt('data_fg',out)
    
if rank==1:
    comm.send(I,dest=0)
    print 'rank1 ok'
if rank ==2:
    comm.send(I,dest=0)
    print 'rank2 ok'
if rank==3:
    comm.send(I,dest=0)
    print 'rank3 ok'



#print out

