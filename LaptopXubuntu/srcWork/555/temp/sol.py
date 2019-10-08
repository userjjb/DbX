import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def eval(x,y):
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y,x)
    theta += (2 * np.pi)*(theta < 0.0)
    return r**(2.0/3) * np.sin(2 * theta / 3.0)

xx = np.linspace(-1,1,50)
XX,YY = np.meshgrid(xx,xx)
#plt.plot(XX,YY,eval(XX,YY),projection='3d')
fig = plt.figure()
ax = fig.gca(projection='3d')

surf = ax.plot_surface(XX, YY, eval(XX,YY), linewidth=0)



