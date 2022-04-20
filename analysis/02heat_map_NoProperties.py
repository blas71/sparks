import numpy as np
import matplotlib.pyplot as plt
import imp
import RyR_cluster
import scipy as sp
import scipy.ndimage
from scipy.interpolate import interp2d
import os

RyR_cluster = imp.reload(RyR_cluster)
clustering = RyR_cluster.clustering

# Use latex figures
plt.rc('text', usetex = True)
plt.rc('font', **{'family': 'serif', 'serif': ['Times New Roman']})

dx = 1
save = 0

nryrs_th = 0

Nx = 800
Ny = 120

gg = 'RyR_P'

wd = ['NoAF', 'AF']
nconf = {wd[0]: [1, 2, 3],
         wd[1]: [1, 2, 3]}
vmax = [2, 5]
tdes = 5000
folder_name = '../plots/' + gg + '/'

cmap = 'coolwarm'
dt = 0.012

colors = ['darkgray', 'royalblue']
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

tmax = {}

for foli in wd:
    tmax[foli] = 0

EVENTS = {}

for foli in wd:
    EVENTS[foli] = {}
    for i in nconf[foli]:
        ff = '../' + gg + '/' + foli + '/' + foli + str(i)
        file_follow = open(ff + '/data/follow_RyR_1.data', 'r')
        lines_ryrs = file_follow.readlines()
        file_follow.close()
        EVENTS[foli][i] = np.array(lines_ryrs)

fig2 = plt.figure(figsize=(5,4))
ax2 = plt.subplot(111)
ww = 0.4
ss = -1

i = 0
for wdi in wd:
    MAT = np.zeros((Ny,Nx))
    X = []
    Y = []
    for nf in nconf[wdi]:
        fol = '../' + gg + '/' + wdi + '/' + wdi + str(nf) + '/'
        cj, ryrk, t2, t1 = np.loadtxt(fol + 'data/timeRyR-open.data').T
        posRyR = np.loadtxt(fol + 'data/posRyR.data').T
        posRyR_sarco = np.loadtxt(fol + 'data/posRyR-sarco.data').T
        posRyR = np.hstack((posRyR,posRyR_sarco))
        nRyR = posRyR.shape[1]

        nRyRsPixel = np.loadtxt(fol + 'data/nRyRpixel.data').T
        nRyRsPixel = np.hstack((nRyRsPixel, 9*np.ones(posRyR_sarco.shape[1])))

        cj = cj.astype(int) - 1
        ryrk = ryrk.astype(int)
        t2 = t2.astype(int)
        t1 = t1.astype(int)
        posRyR = posRyR.astype(int)
        
        idx_t1 = np.argsort(t1)
        t1 = t1[idx_t1]
        t2 = t2[idx_t1]
        cj = cj[idx_t1]
        
        idx2 = t1>int(tdes/dt)
        t1 = t1[idx2]
        t2 = t2[idx2]
        cj = cj[idx2]
        

        cluster_size, min_dist, Vcluster, Vdict, VRyRsPixel = clustering(0.15, 0.1, posRyR.T, nRyRsPixel)
        nClusters = len(list(Vcluster.keys()))
        color = np.zeros(nClusters)

        Nlines = len(EVENTS[wdi][nf])

        for l in range(Nlines):
            line = EVENTS[wdi][nf][l]
            line = line.split(' \n')[0]
            line = line.split(', ')
            line = np.array([int(integer) for integer in line])
            cli = []
            x1 = 0
            x2 = 0
            for li in line:
                x1 += posRyR[0,cj[li]] / len(line)
                x2 += posRyR[1,cj[li]] / len(line)
                pos = (posRyR[0,cj[li]], posRyR[1,cj[li]])
                cli.append(Vdict[pos])
            X.append(x1)
            Y.append(x2)
            cli = list(set(cli))
            for clj in cli:
                color[clj-1] += 1
        
        tmax[wdi] += np.max(t2)*dt/1000 - tdes/1000

        cond = color>0
        maxi = max(color)
        for l in range(nClusters):
            cli = color[l]
            if cli>0:
                x, y = zip(*Vcluster[l+1])
                MAT[x,y] += cli*np.ones(len(x))
    
    MAT = sp.ndimage.filters.gaussian_filter(MAT, 1, mode='constant')
    #MAT[MAT==0] = None
    dist = []
    dlim = 0
    for l in range(len(X)):
        d1 = X[l]
        d2 = Ny-X[l]
        d3 = Y[l]
        d4 = Nx-Y[l]
        d = np.min([d1, d2, d3, d4])
        if d>=dlim:
            dist.append(int(d/r_ring))

    bins = np.arange(len(AREA))
    a1, a2 = np.histogram(dist, bins=len(AREA))
    a1[0] = a1[0] * AREA[0] / ((Nx-2*dlim)*(Ny-2*dlim) - (Nx-2*r_ring)*(Ny-2*r_ring)) / 0.1**2
    ff = a1/AREA/tmax[wdi]*1000
    norm = np.sum(ff)
    #ff /= norm
    ax2.bar(bins+ss*ww*0.5, ff, width=ww, color=colors[i])
    ax2.set_xlabel('distance to membrane [$\mu$m]')
    ax2.set_ylabel('sparks/s/1000$\mu$m$^2$')

    fig = plt.figure(figsize=(15,3.5))
    ax = fig.add_subplot(111)
    bb = ax.imshow(MAT, cmap=cmap, origin="lower", vmin=0, vmax=vmax[i], aspect='equal')
    fig.colorbar(bb, ax=ax)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    ax.set_xticks([])
    ax.set_yticks([])
    
    ax.set_aspect('equal', 'datalim')
    fig.set_facecolor("white")

    fig.tight_layout()

    if save:
        fig.savefig(folder_name + '02heat_map_sum_' + wdi + '.pdf')

    i += 1
    ss += 2

ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

xlabel = []
for k in range(len(AREA)):
    xlabel.append(str(int(k*0.1*r_ring)) + '-' + str(int((k+1)*0.1*r_ring)))
ax2.set_xticks(np.arange(len(AREA)))
ax2.set_xticklabels(xlabel)

fig2.tight_layout()

if save:
    fig2.savefig(folder_name + '02frequency_layers.pdf')

plt.show()
