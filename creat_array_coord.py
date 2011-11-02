#!/usr/bin/python

from __future__ import print_function
try:    
    from future_builtins import map, zip
except ImportError:
    pass

from sys import argv
from math import pi, sin, cos
from random import uniform, triangular

fn_array_coord = "array_coord.dat"

n = int(argv[1])
maxr = float(argv[2])

coords = []
for i in range(n):
    r, theta = triangular(0, maxr, maxr), uniform(0, 2*pi)
    coords.append((r*sin(theta), r*cos(theta)))
    
with open(fn_array_coord, 'w') as fout:
    for x, y in coords:
        print("{:20.12e}   {:20.12e}".format(x, y), file =fout)
