import numpy as np
import matplotlib.pyplot as plt

error = np.array([
                 [1.19E-5, 1.39E-4, 2.53E-4, 3.61E-4, 4.10E-4, 1.99E-3, 4.42E-3], 
                 [1.96E-7, 1.25E-5, 3.53E-5, 7.83E-5, 3.65E-4, 3.51E-4, 5.61E-4], 
                 [1.84E-7, 5.56E-7, 3.78E-6, 1.19E-5, 1.35E-4, 3.61E-4, 3.78E-4], 
                 [3.20E-8, 1.86E-7, 2.54E-7, 1.90E-7, 1.27E-5, 7.70E-5, 2.14E-4], 
                 [6.61E-9, 5.55E-8, 1.13E-7, 1.86E-7, 5.38E-7, 1.19E-5, 5.35E-5]])

max_quad = np.array([
                 [37, 20, 17, 15, 10, 6, 5], 
                 [50, 25, 19, 17, 10, 8, 7], 
                 [48, 33, 23, 20, 10, 7, 7], 
                 [50, 33, 26, 25, 11, 7, 7], 
                 [50, 26, 25, 19, 11, 9, 6]])

dd = np.array([0.1, 0.15, 0.175, 0.2, 0.3, 0.4, 0.5])
hh = np.array([1/2, 1/3, 1/4, 1/6, 1/8])

plt.loglog(dd, error.T)
plt.figure(2)
plt.loglog(hh, error)

dofs = np.tile(1/hh,(7,1)).T**3 * (max_quad+1)**3

# Second order kernel: error_d is with the dd2 refinement for h=1/4
 # error_h is for h-refinement with d = 0.1
dd2 = [0.1, 0.05, 0.025, 0.0125, 0.00625, 0.003125, 0.0015625]
error_d = [  6.22525439e-03,   1.64214309e-03,   3.89685824e-04,
         9.46775673e-05,   1.71911336e-05,   1.15075688e-05,
         1.39282609e-05]
error_h = [1.64E-3, 7E-4, 3.96E-4, 1.71E-4, 9.5E-5]

# Lazy convg fit
def lzp(x, y):
    return np.polyfit(np.log10(x),np.log10(y),1)