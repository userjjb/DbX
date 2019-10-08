import numpy as np
import scipy as sp

import matplotlib.pyplot as plt
import numpy.polynomial.legendre as leg

n = 16

def f(x,y):
    return np.sin(2*np.pi*x) * np.sin(2*np.pi*y)

def curve(a, b, c, x):
    return a + b*x + c*x**2

Qx,Qw = leg.leggauss(n)
nx,ny = np.meshgrid(Qx,Qx)

Le = leg.legvander(Qx, n-1)
#La = sp.interpolate.lagrange(n)

Q,R = np.linalg.qr(Le)