#!/usr/bin/env python
# coding=utf-8
data_path = '/home/ycli/data/gbt/gbt_15hr_41-80_pointcorr/secA_15hr_41-80_pointcorr_clean_map_I_800.npy'
from core import algebra
def ReadMeta(data_path):
    '''return  freq ra dec'''
    data = algebra.make_vect(algebra.load(data_path))
    freq = data.get_axis('freq')
    ra = data.get_axis('ra')
    dec = data.get_axis('dec')
    return freq,ra,dec

