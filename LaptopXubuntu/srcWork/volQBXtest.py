import numpy as np

order = 5
nelems_x = 10
nelems_y = 10
B = [-1,1,-1,1] #left,right,bottom,top

def potential(x,y,s,K):
    return np.exp(-(x**2/(s*K) + y**2/s))

def density(x,y,s,K): #Analytically determined from potential
    return (1/(2*np.pi)) * potential(x,y,s,K) * ((4*(K**2*y**2+x**2)-(2*s*K)*(K+1)) / (s*K)**2)

xx = np.linspace(-1,1,300)
xx,yy = np.meshgrid(xx,xx)

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure(figsize=(9,5))
ax = fig.gca(projection='3d')
ax.plot_surface(xx,yy,potential(xx,yy,np.sqrt(0.001),1))

Qx,Qw = np.polynomial.legendre.leggauss(order+1)
vx,vy = np.meshgrid(np.linspace(B[0],B[1],nelems_x+1),np.linspace(B[2],B[3],nelems_y+1))
hx = vx[1:,1:] - vx[:-1,:-1]
hy = vy[1:,1:] - vy[:-1,:-1]
eQx,eQy = np.meshgrid((Qx+1)/2 , (Qx+1)/2)

na = np.newaxis
nx = eQx*hx.ravel()[:,na,na] + vx[:-1,:-1].ravel()[:,na,na]
ny = eQy*hy.ravel()[:,na,na] + vy[:-1,:-1].ravel()[:,na,na]