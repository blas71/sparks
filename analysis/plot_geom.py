import numpy as np
import matplotlib.pyplot as plt
from sys import platform

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

# Use latex figures
plt.rc('text', usetex = True)
plt.rc('font', **{'family': 'serif', 'serif': ['Times New Roman']})

dx = 1
zoom = 0
save = 0

wd = 'RyR_P/12x80/noAF/RyRP_noAF1/'
posRyR = np.loadtxt('../' + wd + 'data/posRyR.data').T
posRyR_sarco = np.loadtxt('../' + wd + 'data/posRyR-sarco.data').T
posTT = np.loadtxt('../' + wd + 'data/posTT-LCC.data').T
posAT = np.loadtxt('../' + wd + 'data/posAT-LCC.data').T
# posLCC = np.loadtxt('../' + wd + 'data/posLCC.data').T
# posNCX = np.loadtxt('../' + wd + 'data/posNCX.data').T
# posSR = np.loadtxt('../' + wd + 'data/posSR.data').T
# posTnC = np.loadtxt('../' + wd + 'data/posTnC.data').T

fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(111)

#ax.scatter(dx*posSR[1,:], dx*posSR[0,:], s=10, marker='o', c='orange', edgecolors=None, alpha=0.4, label='SR')
#ax.scatter(dx*posNCX[1,:], dx*posNCX[0,:], s=20, marker='o', c='orange', alpha=0.4, label='TT')
#ax.scatter(dx*posLCC[1,:], dx*posLCC[0,:], s=100, marker='x', linewidths=3.5, c='green', label='LCC')
ax.scatter(dx*posRyR[1,:], dx*posRyR[0,:], s=1, marker='s', label='RyR')
ax.scatter(dx*posRyR_sarco[1,:], dx*posRyR_sarco[0,:], s=3, marker='s', c='r', alpha=0.3, label='RyR membrane')
ax.scatter(dx*posTT[1,:], dx*posTT[0,:], s=5, marker='o', c='black', label='T-tubules')
ax.scatter(dx*posAT[1,:], dx*posAT[0,:], s=5, marker='o', c='black', label='A-tubules')

plt.xlabel('X [$\mu$m]')
plt.ylabel('Y [$\mu$m]')

#x1, x2, y1, y2 = plt.axis()
#plt.axis([-0.5, 100.5, 0, 15])

plt.axes().set_aspect('equal', 'datalim')
plt.legend(loc=0)

if save:
    plt.savefig('../tubules2/plots/TT_geom.pdf')

plt.show()

