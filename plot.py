#!/usr/bin/python

from __future__ import print_function
try:    
    from future_builtins import map, zip
except ImportError:
    pass

import matplotlib.pyplot as plt

X, Y, Z = [], [], []
#for line1, line2 in zip(open("data/array_coord.dat"), open("Q.dat")):
#    x, y = map(float, line1.split())
#    z = float(line2)
#    X.append(x)
#    Y.append(y)
#    Z.append(z)
    
for line in open("a.dat"):
    x, y, z = map(float, line.split())
    X.append(x)
    Y.append(y)
    Z.append(z)
    
    
plt.hexbin(X, Y, Z)
plt.show()
