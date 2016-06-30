#!/usr/bin/env python
# coding=utf-8
def readPar(path):
    f=open(path,'r')
    data=f.readlines()
    f.close()
    P=[]
    for i in data:
        if '#' in i:
            continue
        elif i=='\n':
            continue
        else:
            P.append(i.split('\n')[0].split('='))
    par={}
    for i in P:
        if 'Path' in i[0]:
            i[1]=str(i[1])
        else:
            i[1]=float(i[1])
        par[i[0]]=i[1]
    print '='*80
    print 'Read par ...'
    print par
    return par
