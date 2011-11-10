#!/usr/bin/python

from __future__ import print_function
try:    
    from future_builtins import map, zip
except ImportError:
    pass

from sys import argv
from math import pi, sin, cos, sqrt, hypot
from random import uniform, triangular, normalvariate

fn_array_coord = "data/array_coord.dat"

R = float(argv[1])
min_d = float(argv[2])
N = int(argv[3])

coords = []

def consistent(x, y):
    return all(hypot(x-xi, y-yi) > min_d for xi, yi in coords)
    
while len(coords) < N:
    theta, r = uniform(0, 2*pi), triangular(0, R, R)
    x, y = r*sin(theta), r*cos(theta)
    if consistent(x, y):
        coords.append((x, y))
    
#for i in range(1, int(n)+1):
#    dtheta = 2*pi/i
#    r0 = i*d
#    theta0 = uniform(0, 2*pi)
#    for j in range(i):
#        r =  normalvariate(r0, d*sigma)
#        theta = normalvariate(dtheta*j, dtheta*sigma) + theta0
#        coords.append((r*sin(theta), r*cos(theta)))
    
with open(fn_array_coord, 'w') as fout:
    for x, y in coords:
        print("{0:20.12e}   {1:20.12e}".format(x, y), file=fout)
