#!/usr/bin/python

from __future__ import print_function
try:    
    from future_builtins import map, zip
except ImportError:
    pass

import matplotlib.pyplot as plt


X, Y, Z = [], [], []
fn = "data/500.current.dat"
for line in open(fn):
    try:
        x, y, z = map(float, line.split()[:3])
        X.append(x)
        Y.append(y)
        Z.append(z)
    except ValueError:
        pass
    
    
plt.hexbin(X, Y, Z)
plt.show()
