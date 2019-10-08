import numpy as np
import matplotlib.pyplot as plt
import cseAnalytical
from timeit import default_timer as timer

na = np.newaxis
pi = np.pi
sin = np.sin

def f(x, y, z):
    return 1

# naming goes as K(order de-sing kernel)S(order Taylor expansion of kernel)
def K2S3(d, r):
    return (8*d**6 + 8*d**4*r**2 + 7*d**2*r**4 + 2*r**6)/(2.*(d**2 + r**2)**3.5)

def K2S3E(d, r):
    k = d**2 + r**2
    return d**2/k**1.5 + 1/np.sqrt(k) + (d**2*(2*d**2 - r**2))/(2.*k**2.5) - (d**3*(-2*d**3 + 3*d*r**2))/(2.*k**3.5)

def K2S6(d, r):
    k = d**2 + r**2
    return d**2/k**1.5 + 1/np.sqrt(k) + (d**4*(2*d**2 - 3*r**2))/(2.*k**3.5) + (d**2*(2*d**2 - r**2))/(2.*k**2.5) + (d**4*(8*d**4 - 24*d**2*r**2 + 3*r**4))/(8.*k**4.5) + (d**6*(8*d**4 - 40*d**2*r**2 + 15*r**4))/(8.*k**5.5) + (d**6*(16*d**6 - 120*d**4*r**2 + 90*d**2*r**4 - 5*r**6))/(16.*k**6.5)

def K4S4(d, r):
    k = d**2 + r**2
    return (3*d**4)/(2.*k**2.5) + (3*d**2 + 2*r**2)/(2.*k**1.5) + (3*d**2*(2*d**4 - 3*d**2*r**2))/(4.*k**3.5) - (d**3*(-6*d**5 + 23*d**3*r**2 - 6*d*r**4))/(4.*k**4.5) + (3*d**4*(8*d**6 - 56*d**4*r**2 + 39*d**2*r**4 - 2*r**6))/(16.*k**5.5)

def K4S5(d, r):
    k = d**2 + r**2
    return (3*d**4)/(2.*k**2.5) + (3*d**2 + 2*r**2)/(2.*k**1.5) + (3*d**2*(2*d**4 - 3*d**2*r**2))/(4.*k**3.5) + (d**4*(6*d**4 - 23*d**2*r**2 + 6*r**4))/(4.*k**4.5) + (3*d**6*(8*d**6 - 88*d**4*r**2 + 115*d**2*r**4 - 20*r**6))/(16.*k**6.5) + (3*d**4*(8*d**6 - 56*d**4*r**2 + 39*d**2*r**4 - 2*r**6))/(16.*k**5.5)

def K4S6(d, r):
    k = d**2 + r**2
    return (3*d**4)/(2.*k**2.5) + (3*d**2 + 2*r**2)/(2.*k**1.5) + (3*d**2*(2*d**4 - 3*d**2*r**2))/(4.*k**3.5) + (d**4*(6*d**4 - 23*d**2*r**2 + 6*r**4))/(4.*k**4.5) + (3*d**6*(8*d**6 - 88*d**4*r**2 + 115*d**2*r**4 - 20*r**6))/(16.*k**6.5) + (3*d**4*(8*d**6 - 56*d**4*r**2 + 39*d**2*r**4 - 2*r**6))/(16.*k**5.5) + (d**6*(48*d**8 - 760*d**6*r**2 + 1590*d**4*r**4 - 585*d**2*r**6 + 20*r**8))/(32.*k**7.5)

def K6S7(d, r):
    k = d**2 + r**2
    return (15*d**6)/(8.*k**3.5) + (15*d**2*(2*d**6 - 5*d**4*r**2))/(16.*k**4.5) + (15*d**4 + 20*d**2*r**2 + 8*r**4)/(8.*k**2.5) + (5*d**3*(6*d**7 - 37*d**5*r**2 + 20*d**3*r**4))/(16.*k**5.5) + (15*d**4*(8*d**8 - 88*d**6*r**2 + 115*d**4*r**4 - 20*d**2*r**6))/(64.*k**6.5) + (3*d**5*(40*d**9 - 680*d**7*r**2 + 1563*d**5*r**4 - 680*d**3*r**6 + 40*d*r**8))/(64.*k**7.5) + (5*d**6*(48*d**10 - 1160*d**8*r**2 + 4078*d**6*r**4 - 3195*d**4*r**6 + 520*d**2*r**8 - 8*r**10))/(128.*k**8.5) + (15*d**7*(16*d**11 - 520*d**9*r**2 + 2578*d**7*r**4 - 3143*d**5*r**6 + 980*d**3*r**8 - 56*d*r**10))/(128.*k**9.5)
# d=0.15*h, quad_order = 50 + 5*i
# np.polyfit(np.log10([1/2,1/4,1/6,1/8]),np.log10([1.88E-5,3.38E-7,7.5608705812231233e-09,8.69E-11]),1)

def map_bounds(tgt):
    return np.array([-0.5-tgt[0], 0.5-tgt[0], -0.5-tgt[1], 0.5-tgt[1], -0.5-tgt[2], 0.5-tgt[2]])


domain = np.array([-0.5, 0.5]) # square domain for simplicity

test_tgts = np.array([[0,0,0],[0.1,0,0],[0,0.3,0],[-0.4,-0.4,0],[0.1,0.25,-0.38]])
n_refine = 4
error2 = np.zeros((n_refine, test_tgts.shape[0]))
error4 = np.zeros((n_refine, test_tgts.shape[0]))

q = []
start = timer()
for i in range(n_refine):
    print(i)
    quad_order = 35 + 5*i
    h = 1/3
    d = 0.15*h # expansion center distance
    
    Qx, Qw = np.polynomial.legendre.leggauss(quad_order+1)
    Qxm = h * (Qx+1)/2 # mapped quad points in interval [0, h]
    
    nx, ny, nz = np.meshgrid(Qxm, Qxm, Qxm) # nodes of the tensor-product element
    Qwp = (h/2)**3 * Qw * Qw[:, na] * Qw[:, na, na] # quad weights for 3D element
    
    elem_origins = np.arange(domain[0], domain[1], h)
    ex, ey, ez = np.meshgrid(elem_origins, elem_origins, elem_origins)
    # construct the per element nodes, ravel goes in the order z, x, then y within 
    #  each element and across elements in the domain
    # raveled nodes are easier to handle when we do the dot product with quad weights
    # nex etc. indexing is nex[element number, intra-element node number]
    nex = nx.ravel()[na, :] + ex.ravel()[:, na]
    ney = ny.ravel()[na, :] + ey.ravel()[:, na]
    nez = nz.ravel()[na, :] + ez.ravel()[:, na]
    
    density = f(nex, ney, nez)
    
    for j, tgt in enumerate(test_tgts):
        r = np.sqrt((nex - tgt[0])**2 + (ney - tgt[1])**2 + (nez - tgt[2])**2)
        Qwpr = Qwp.ravel()
        exact = cseAnalytical.exact_u(map_bounds(tgt))
        
        #kernel_order2 = K2S6(d, r)
        #p2 = np.sum( (kernel_order2 * density) @ Qwpr)
        #error2[i, j] = np.abs((p2-exact)/exact)
        kernel_order4 = K6S7(d, r)
        p4 = np.sum( (kernel_order4 * density) @ Qwpr)
        error4[i, j] = np.abs((p4-exact)/exact)
        # pr = np.sum( (density/r) @ Qwpr)
        # errorR[i, j] = np.abs((pr-exact)/exact)
    last = quad_order
    q.append(quad_order)

end = timer()
print(end - start)

plt.figure(figsize=(7.8,5.2))
#plt.semilogy(q, np.max(error2, axis=1), label = 'K2S3, h: ' + np.str(h) + ', d: ' + np.str(d/h))
plt.semilogy(q, np.max(error4, axis=1), label = 'K4S5, h: ' + np.str(h) + ', d: ' + np.str(d/h))
plt.text(last*0.9, np.max(error4[-1])*1.1, '%.2E' % np.max(error4[-1]))
plt.legend()
plt.tight_layout()