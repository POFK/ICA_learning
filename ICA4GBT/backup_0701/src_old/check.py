#!/usr/bin/env python
# coding=utf-8
import numpy as np
import matplotlib.pyplot as plt
data=np.load('secB_15hr_41-80_pointcorr_clean_map_I_800_cleaned_2.npy')
plt.imshow(data[8])
plt.colorbar()
plt.show()
