import numpy as np
import matplotlib.pyplot as plt
import imp
import RyR_cluster
import os

RyR_cluster = imp.reload(RyR_cluster)
clustering = RyR_cluster.clustering

# Use latex figures
plt.rc('text', usetex = True)
plt.rc('font', **{'family': 'serif', 'serif': ['Times New Roman']})

dx = 1
save = 0

Nx = 800
Ny = 120

wd = ['NoAF', 'AF']
nconf = {wd[0]: [1, 2, 3, 4, 5],
         wd[1]: [1, 2, 3, 4, 5]}

gg = 'RyR_P'

r_ring = 10
AREA = []
i = 0
while 2*i*r_ring < Ny:
    rr1 = 2*i*r_ring
    rr2 = 2*(i+1)*r_ring
    Aout = (Nx-rr1)*(Ny-rr1)
    Ain = (Nx-rr2)*(Ny-rr2)
    AREA.append(Aout-Ain)
    i += 1

AREA = 0.1**2 * np.array(AREA)
CLUSTERS = {}
RYRS = {}
i = 0
for wdi in wd:
    CLUSTERS[wdi] = np.zeros(len(AREA))
    RYRS[wdi] = np.zeros(len(AREA))
    for nf in nconf[wdi]:
        fol = '../' + gg + '/' + wdi + '/' + wdi + str(nf) + '/'
        posRyR = np.loadtxt(fol + 'data/posRyR.data').T
        posRyR_sarco = np.loadtxt(fol + 'data/posRyR-sarco.data').T
        posRyR = np.hstack((posRyR,posRyR_sarco))
        nRyR = posRyR.shape[1]

        nRyRsPixel = np.loadtxt(fol + 'data/nRyRpixel.data').T
        nRyRsPixel = np.hstack((nRyRsPixel, 9*np.ones(posRyR_sarco.shape[1])))

        cluster_size, min_dist, Vcluster, Vdict, VRyRsPixel = clustering(0.15, 0.1, posRyR.T, nRyRsPixel)
        nClusters = len(list(Vcluster.keys()))
        
        for j in Vcluster:
            center = np.mean(Vcluster[j], axis=0)
            d1 = center[0]
            d2 = Ny-center[0]
            d3 = center[1]
            d4 = Nx-center[1]
            d = np.min([d1, d2, d3, d4])
            ring = int(d/r_ring)
            if ring==6:
                ring=5
            CLUSTERS[wdi][ring] += 1
            for k in Vcluster[j]:
                ryrs_k = VRyRsPixel[k][1]
                d1 = k[0]
                d2 = Ny-k[0]
                d3 = k[1]
                d4 = Nx-k[1]
                d = np.min([d1, d2, d3, d4])
                ring = int(d/r_ring)
                if ring==6:
                    ring=5
                RYRS[wdi][ring] += ryrs_k
                 
fig = plt.figure(figsize=(5,5))
ax = plt.subplot(221)
for wdi in wd:
    ax.plot(CLUSTERS[wdi], '--o', label=wdi)
ax.set_title('Number of clusters per layer')
ax.legend()

ax = plt.subplot(223)
i = 0
for wdi in wd:
    ax.plot(CLUSTERS[wdi]/AREA[i], '--o', label=wdi)
    i += 1
ax.set_title('Number of clusters per layer / Area layer')
ax.legend()

ax = plt.subplot(222)
for wdi in wd:
    ax.plot(RYRS[wdi], '--o', label=wdi)
ax.set_title('Number of RyRs per layer')
ax.legend()

ax = plt.subplot(224)
i = 0
for wdi in wd:
    ax.plot(RYRS[wdi]/AREA[i], '--o', label=wdi)
    i += 1
ax.set_title('Number of RyRs per layer / Area layer')
ax.legend()

fig.tight_layout()
plt.show()

