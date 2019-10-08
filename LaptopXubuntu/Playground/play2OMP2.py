import numpy as np
import matplotlib.pyplot as plt
import cseAnalytical
from timeit import default_timer as timer
import kernels2OMP
import kernels2

na = np.newaxis

# d=0.15*h, quad_order = 5 + 5*i
# np.polyfit(np.log10([1/2, 1/3, 1/4, 1/6, 1/8]),np.log10([1.88E-5, 1.84E-6, 3.38E-7, 7.56e-09, 8.69E-11]),1)
# d=0.2*h, quad_order = 5 + 3*i
# np.polyfit(np.log10([1/2, 1/3, 1/4, 1/6, 1/8]),np.log10([1.18E-4, 3.23E-6, 1.84E-6, 1.28e-07, 7.76E-9]),1)
# 0.03 1.07E-10 290, 0.05 7.76E-9 140, 0.075 3.38E-7 80, 0.1 1.85E-6 55, 0.15 1.88E-5 35, 0.2 1.18E-4 25
# 0.026 4.95E-10, 0.028 1.14E-10, 0.029 9.13E-11 366s
# d=0.015 <600 1.15E-10 423s
def map_bounds(tgt):
    return np.array([-0.5-tgt[0], 0.5-tgt[0], -0.5-tgt[1], 0.5-tgt[1], -0.5-tgt[2], 0.5-tgt[2]])


domain = np.array([-0.5, 0.5]) # square domain for simplicity

test_tgts = np.array([[0.1,0,0],[0,0.3,0],[-0.4,-0.4,0],[0.1,0.25,-0.38]])
n_refine = 1
error = np.zeros((n_refine, test_tgts.shape[0]))

q = []
start = timer()
#@profile
#def dummy():
for i in range(n_refine):
    print('\n'+np.str(i)+' ', end='')
    quad_order = 30 + 20*i
    h = 1/2
    d = 0.025*h # expansion center distance
    
    Qx, Qw = np.polynomial.legendre.leggauss(quad_order+1)
    Qxm = h * (Qx+1)/2 # mapped quad points in interval [0, h]
    
    nodeX = np.array(np.meshgrid(Qxm, Qxm, Qxm)).reshape(3, -1) # nodes of the tensor-product element
    Qwp = (h/2)**3 * Qw * Qw[:, na] * Qw[:, na, na] # quad weights for 3D element
    
    elem_origins = np.arange(domain[0], domain[1], h)
    elemX = np.array(np.meshgrid(elem_origins, elem_origins, elem_origins)).reshape(3, -1)
    # construct the per element nodes, ravel goes in the order z, x, then y within 
    #  each element and across elements in the domain
    # raveled nodes are easier to handle when we do the dot product with quad weights
    # nex etc. indexing is nex[element number, intra-element node number]
#    nex = nx.ravel()[na, :] + ex.ravel()[:, na]
#    ney = ny.ravel()[na, :] + ey.ravel()[:, na]
#    nez = nz.ravel()[na, :] + ez.ravel()[:, na]
    
    density = 1 # future implementations with non-trivial density will need a way to compute this in the F90 module
    Qwp = Qwp.ravel()
    print('\n' + np.str(timer() - start))
    for j, tgt in enumerate(test_tgts):
        print(j,end='')
        exact = cseAnalytical.exact_u(map_bounds(tgt))
        kern = kernels2OMP.k6s7(8, d, nodeX, elemX, tgt, density, Qwp)
        #kern = kernels2.k6s7(d, nex, ney, nez, tgt)
        pot = np.sum( (kern * density) @ Qwp)
        error[i, j] = np.abs((pot-exact)/exact)

    q.append(quad_order)

print('\n' + np.str(timer() - start))
print(np.max(error[-1,:]))

#dummy()

if n_refine>1:
    plt.figure(figsize=(7.8,5.2))
    plt.semilogy(q, np.max(error, axis=1), label = 'K4S5, h: ' + np.str(h) + ', d: ' + np.str(d/h))
    plt.text(q[-1]*0.9, np.max(error[-1])*1.1, '%.2E' % np.max(error[-1]))
    plt.legend()
    plt.tight_layout()
