#!/usr/bin/python

from __future__ import print_function
try:    
    from future_builtins import map, zip
except ImportError:
    pass

from sys import argv
from math import pi, sin, cos, sqrt
from random import uniform, triangular, normalvariate

fn_array_coord = "data/array_coord.dat"

N = int(argv[1])
R = float(argv[2])

n = sqrt(N)+2
d = R/n
sigma = .1

print(d, n)

coords = []
for i in range(1, int(n)+1):
    dtheta = 2*pi/i
    r0 = i*d
    for j in range(i):
        r =  normalvariate(r0, d*sigma)
        theta = normalvariate(dtheta*j, dtheta*sigma)
        coords.append((r*sin(theta), r*cos(theta)))
    
with open(fn_array_coord, 'w') as fout:
    for x, y in coords:
        print("{:20.12e}   {:20.12e}".format(x, y), file=fout)
