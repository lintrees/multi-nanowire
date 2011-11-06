#!/usr/bin/python

from __future__ import print_function
try:    
    from future_builtins import map, zip
except ImportError:
    pass

from math import hypot

fnin = "barrier.debug"
fnout = "barrier.dat"

origin = True

curve = []
for line in open(fnin):
    x, y, U = map(float, line.split())
    if origin:
        x0, y0 = x, y
        origin = False
    d = hypot(x-x0, y-y0)
    curve.append((d, -U))

if __name__ == "__main__":
    with open(fnout, 'w') as fout:
        for d, U in curve:
            print("{:.15e}  {:12.8f}".format(d, U), file=fout)
            
del fnin, fnout, origin
