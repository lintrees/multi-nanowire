#!/usr/bin/python

from __future__ import print_function
try:    
    from future_builtins import map, zip
except ImportError:
    pass

from scipy.constants import epsilon_0, pi
from scipy.optimize import newton_krylov, broyden1, broyden2

from modify_potential_curve import curve

coeff = 1/(4*pi*epsilon_0)

def f(X):
    R, Q, U0 = X
    fi = lambda xi: coeff*Q/(xi+R)+U0
    print(R, Q, U0)
    return sum(map(lambda xy: (fi(xy[0])-xy[1])**2, curve))
    
guess = [1e-8, -1e-16, 5]
sol = broyden2(f, guess, verbose=1)

print(sol)


