#!/usr/bin/env python
# coding=utf-8
from read_par import readPar
class IM():
    def __init__(self,parFile='Fg.par'):
        print 'init ...'
        par=readPar(parFile)
        self.Path1=par['Path1']
        self.Path2=par['Path2']
        self.Nx=par['Nx']
        self.Ny=par['Ny']
        self.Nz=par['Nz']
        self.h =par['h']
        self.Om0=par['Om0']
        self.bins=par['bins']
        self.ICn=par['ICn']
        self.PCn=par['PCn']
        print 'Nx=%d\t\t'%self.Nx,'Ny=%d\t\t'%self.Ny,'Nz=%d\t\t'%self.Nz
        print 'Number of ICs :\t\t%d'%self.ICn
        print 'init finished ...'


