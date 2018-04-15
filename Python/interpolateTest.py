# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 09:53:40 2016

@author: peterweilnboeck
"""

import matplotlib.pyplot as plt
from numpy import genfromtxt
from numpy import transpose
from scipy.interpolate import interp1d

my_data = genfromtxt('EOS-numbers1.csv', delimiter=',')

my_data1 = transpose(my_data)
rho = my_data1[0]
p = my_data1[1]

rho[0]
p[0]
rhoFromP = interp1d(p, rho)
pFromRho = interp1d(rho, p, kind='cubic')

plt.plot(p, rhoFromP(p))
plt.show()

print my_data1
print "Rho: ",rho
print "p: ",p