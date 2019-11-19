#from mpl_toolkits.mplot3d import Axes3D
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.scatter(nx.ravel(), ny.ravel(), nz.ravel(), c=Qwp.ravel(), cmap=plt.hot())
#plt.show()

import numpy as np
import matplotlib.pyplot as plt

na = np.newaxis
pi = np.pi
sin = np.sin

d = 0.01 # expansion center distance
quad_order = 11
h = 1/4
domain = np.array([-1,1]) # square domain for simplicity
domain_width = domain[1]-domain[0]

Qx, Qw = np.polynomial.legendre.leggauss(quad_order+1)
Qxm = (h/domain_width) * (Qx+1) # mapped quad points run from 0 to h

nx, ny, nz = np.meshgrid(Qxm, Qxm, Qxm) # nodes of the tensor-product element
Qwp = (h/domain_width)**3 * Qw * Qw[:, na] * Qw[:, na, na] # quad weights for 3D element

elem_origins = np.arange(domain[0], domain[1], h)
ex, ey, ez = np.meshgrid(elem_origins, elem_origins, elem_origins)
# construct the per element nodes, ravel goes in the order z, x, then y within 
#  each element and across elements in the domain
# raveled nodes are easier to handle when we do the dot product with quad weights
# nex etc. indexing is nex[element number, intra-element node number]
nex = nx.ravel()[na, :] + ex.ravel()[:, na]
ney = ny.ravel()[na, :] + ey.ravel()[:, na]
nez = nz.ravel()[na, :] + ez.ravel()[:, na]

#def f(x, y, z):
#    return sin(pi*x) * sin(pi*y) * sin(pi*z)
#
#def exact_u(x, y, z):
#    return sin(pi*x) * sin(pi*y) * sin(pi*z) / (-3*pi**2)

def f(x, y, z):
    return (-3*np.exp(-x**2/32 - y**2/32 - z**2/32))/16 + (np.exp(-x**2/32 - y**2/32 - z**2/32)*x**2)/256 + (np.exp(-x**2/32 - y**2/32 - z**2/32)*y**2)/256 + (np.exp(-x**2/32 - y**2/32 - z**2/32)*z**2)/256

def exact_u(x, y, z):
    return np.exp(-(x**2/32 + y**2/32 + z**2/32))

# naming goes as K(order de-sing kernel)S(order Taylor expansion of kernel)
def K2S3(d, r):
    return (8*d**6 + 8*d**4*r**2 + 7*d**2*r**4 + 2*r**6)/(2.*(d**2 + r**2)**3.5)

def K4S5(d, r):
    k = d**2 + r**2
    return (3*d**4)/(2.*k**2.5) + (3*d**2 + 2*r**2)/(2.*k**1.5) + (3*d**2*(2*d**4 - 3*d**2*r**2))/(4.*k**3.5) + (d**4*(6*d**4 - 23*d**2*r**2 + 6*r**4))/(4.*k**4.5) + (3*d**6*(8*d**6 - 88*d**4*r**2 + 115*d**2*r**4 - 20*r**6))/(16.*k**6.5) + (3*d**4*(8*d**6 - 56*d**4*r**2 + 39*d**2*r**4 - 2*r**6))/(16.*k**5.5)

density = f(nex, ney, nez)
tgt = np.array([0.5,0.5,0.5]) # example test target
r = np.sqrt((nex - tgt[0])**2 + (ney - tgt[1])**2 + (nez - tgt[2])**2)

kernel_order2 = K2S3(d, r)
kernel_order4 = K4S5(d, r)

Qwpr = Qwp.ravel()
p2 = np.sum( (kernel_order2 * density) @ Qwpr)
p4 = np.sum( (kernel_order4 * density) @ Qwpr)

print(p2)
print(p4)
print(exact_u(tgt[0], tgt[1], tgt[2]))