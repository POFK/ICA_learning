#!/usr/bin/env python
# coding=utf-8
import matplotlib.pyplot as plt
import numpy as np

a=np.random.randn(3,5)
plt.pcolor(a)
plt.colorbar()
plt.show()
