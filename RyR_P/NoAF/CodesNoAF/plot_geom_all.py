import numpy as np
import matplotlib.pyplot as plt
from sys import platform

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

# Use latex figures
plt.rc('text', usetex = True)
plt.rc('font', **{'family': 'serif', 'serif': ['Times New Roman']})

dx = 0.1

fig = plt.figure(figsize=(20,6))

posTT = np.loadtxt('../data/posTT-LCC.data').T
posLCC = np.loadtxt('../data/posLCC.data').T
plt.scatter(dx*posLCC[1,:], dx*posLCC[0,:], s=20, marker='x', linewidths=3.5, c='green', label='LCC')
#plt.scatter(dx*posTT[1,:], dx*posTT[0,:], s=10, marker='o', c='black', label='TT')
plt.xlabel('X [$\mu$m]')
plt.ylabel('Y [$\mu$m]')
#print(np.mean(posTT[0,:]))

#x1, x2, y1, y2 = plt.axis()
#plt.axis([-0.5, 100.5, 0, 15])

#plt.axes().set_aspect('equal', 'datalim')
#plt.legend(loc=0)

if platform == 'darwin': # mac
    fig.set_tight_layout(True)
else:
    fig.tight_layout()

plt.show()

