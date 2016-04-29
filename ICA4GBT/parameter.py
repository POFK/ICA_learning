#!/usr/bin/env python
# coding=utf-8
from ICA import *
N=40
#ra_1,ra_2,dec_1,dec_2=16,-16,9,-9
#ra_1,ra_2,dec_1,dec_2=10,-10,9,-9
ra_1,ra_2,dec_1,dec_2=0,-1,0,-1
freq,ra,dec=ReadMeta('/home/ycli/data/wigglez/gbt_15hr_41-80_pointcorr/reg15data.npy')
#-------------------
ra=ra[ra_1:ra_2]
dec=dec[dec_1:dec_2]
#-------------------

print 'ICs number:',N
ica=ICA_GBT(freq,ra,dec,N=N)
ica.Grid(32,32,128)
#ica.Grid(64,64,128)
bins=30
