import numpy as np
import matplotlib.pyplot as plt
import cseAnalytical

na = np.newaxis
pi = np.pi
sin = np.sin

def f(x, y, z):
    return 1

# naming goes as K(order de-sing kernel)S(order Taylor expansion of kernel)
def K2S3(d, r):
    return (8*d**6 + 8*d**4*r**2 + 7*d**2*r**4 + 2*r**6)/(2.*(d**2 + r**2)**3.5)

def K4S5(d, r):
    k = d**2 + r**2
    return (3*d**4)/(2.*k**2.5) + (3*d**2 + 2*r**2)/(2.*k**1.5) + (3*d**2*(2*d**4 - 3*d**2*r**2))/(4.*k**3.5) + (d**4*(6*d**4 - 23*d**2*r**2 + 6*r**4))/(4.*k**4.5) + (3*d**6*(8*d**6 - 88*d**4*r**2 + 115*d**2*r**4 - 20*r**6))/(16.*k**6.5) + (3*d**4*(8*d**6 - 56*d**4*r**2 + 39*d**2*r**4 - 2*r**6))/(16.*k**5.5)

def map_bounds(tgt):
    return np.array([-0.5-tgt[0], 0.5-tgt[0], -0.5-tgt[1], 0.5-tgt[1], -0.5-tgt[2], 0.5-tgt[2]])


domain = np.array([-0.5, 0.5]) # square domain for simplicity

test_tgts = np.array([[0,0,0],[0.1,0,0],[0,0.3,0],[-0.4,-0.4,0],[0.1,0.25,-0.38]])
n_refine = 9
error2 = np.zeros((n_refine, tgts.shape[0]))
error4 = np.zeros((n_refine, tgts.shape[0]))

for i in range(n_refine):
    print(i)
    quad_order = 6 + 2**i
    h = 1
    d = 0.001*h # expansion center distance
    
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
        
        kernel_order2 = K2S3(d, r)
        kernel_order4 = K4S5(d, r)
        
        Qwpr = Qwp.ravel()
        p2 = np.sum( (kernel_order2 * density) @ Qwpr)
        p4 = np.sum( (kernel_order4 * density) @ Qwpr)
        
        exact = cseAnalytical.exact_u(map_bounds(tgt))
        
        error2[i, j] = np.abs((p2-exact)/exact)
        error4[i, j] = np.abs((p4-exact)/exact)

plt.loglog(6 + 2**np.arange(n_refine), np.max(error2, axis=1), label = '2')
plt.loglog(6 + 2**np.arange(n_refine), np.max(error4, axis=1), label = '4')
plt.legend()