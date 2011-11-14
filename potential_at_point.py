#!/usr/bin/python

from __future__ import print_function
try:    
    from future_builtins import map, zip
except ImportError:
    pass

from sys import argv
from operator import attrgetter

from scipy.interpolate import LinearNDInterpolator
from numpy import array


datadir = "data/"

if len(argv) == 1:
    fn_potential = "single.potential.out"
    fn_triangle = "single.triangle.out"
elif len(argv) == 2 and argv[1] == "background":
    fn_potential = "single.background.potential.out"
    fn_triangle = "single.background.triangle.out"
else:
    raise

x0, y0 = map(float, input().split()[:2])


class Point(object): pass

allPoints = []
for line in open(datadir+fn_potential):
    p = Point()
    p.x, p.y, p.phi, p.doc, p.po = map(float, line.split())
    allPoints.append(p)
    
triangles = []
for line in open(datadir+fn_triangle):
    triNo = map(lambda no: int(no)-1, line.split()[:3])
    tri = frozenset(map(allPoints.__getitem__, triNo))
    assert len(tri) == 3    
    triangles.append(tri)

On_same_side = lambda x, y, x3, y3, x1, y1, x2, y2: (x-x1)*(y1-y2) == (y-y1)*(x1-x2) or ((x-x1)*(y1-y2) < (y-y1)*(x1-x2)) == ((x3-x1)*(y1-y2) < (y3-y1)*(x1-x2))
def Is_in_tri(x, y, tri):
    x1, x2, x3 = map(attrgetter('x'), tri)
    y1, y2, y3 = map(attrgetter('y'), tri)
    return (On_same_side(x, y, x3, y3, x1, y1, x2, y2)
        and On_same_side(x, y, x1, y1, x2, y2, x3, y3)
        and On_same_side(x, y, x2, y2, x3, y3, x1, y1))
        
for tri in triangles:
    if Is_in_tri(x0, y0, tri):
        x = tuple(map(attrgetter('x'), tri))
        y = tuple(map(attrgetter('y'), tri))
        phi = tuple(map(attrgetter('phi'), tri))
        m = array((x, y))
        m = m.transpose()        
        print(LinearNDInterpolator(m, phi)((x0, y0)))
        break
        

