import numpy as np
import matplotlib.pyplot as plt
import imp
import RyR_cluster
import matplotlib.cm as cm
RyR_cluster = imp.reload(RyR_cluster)
clustering = RyR_cluster.clustering

def fun(Nx,Ny,px,py,Roip):
    ROI = []
    for i in range(2*Roip+1):
        for j in range(2*Roip+1):
            key = np.array([i+px,j+py])
            if np.all(key<np.array([Nx,Ny])) and np.all(key>np.array([0,0])):
                ROI.append((i+px,j+py))
    return ROI

# Use latex figures
plt.rc('text', usetex = True)
plt.rc('font', **{'family': 'serif', 'serif': ['Times New Roman']})

save = 1

Nt = 36
Nsave = 1000
dt = 0.0034
Nx = 800
Ny = 120
Roip = 12
tdes = -1 
Nprev = 4
dx = 0.1

wd = ['noAF', 'AF']
colors = ['grey', 'blue']
nf = len(wd)

figs = []
for i in range(nf):
    figs.append(plt.figure(figsize=(8,2)))

ii = 0
for wdi in wd:
    ff = '../RyR_P/12x80/' + wdi + '/RyRP_' + wdi + '1'
    cj, ryrk, t2, t1 = np.loadtxt(ff + '/data/timeRyR-open.data').T
    posRyR = np.loadtxt(ff + '/data/posRyR.data').T
    posRyR_sarco = np.loadtxt(ff + '/data/posRyR-sarco.data').T
    tf, Ca, Casr, CaTnC, CaCM, CabSR, CSQ = np.loadtxt(ff + '/data/calcium.data').T
    posRyR = np.hstack((posRyR,posRyR_sarco))
    nRyR = posRyR.shape[1]

    cj = cj.astype(int) - 1
    ryrk = ryrk.astype(int)
    t2 = t2.astype(int)
    t1 = t1.astype(int)
    posRyR = posRyR.astype(int)
    
    Vcluster, Vdict = clustering(0.15, 0.1, posRyR)
    Vsparks = {}

    idx_t1 = np.argsort(t1)
    t1 = t1[idx_t1]
    t2 = t2[idx_t1]
    cj = cj[idx_t1]
    
    idx2 = t1>int(tdes/dt)
    t1 = t1[idx2]
    t2 = t2[idx2]
    cj = cj[idx2]

    sparks = []
    count = 0

    SparksNotCounted = list(np.arange(len(cj)))

    while SparksNotCounted:
    
        k = SparksNotCounted[0]
    
        px = posRyR[1,cj[k]]
        py = posRyR[0,cj[k]]
        ROI = fun(Nx,Ny,px,py,Roip)
        nROI = len(ROI)
        
        sparks.append((px,py))

        # RyRs in the open state
        corner_ld = np.array((py-Roip,px-Roip))
        corner_ru = np.array((py+Roip,px+Roip))
        
        ToDel = []
        for i in range(len(SparksNotCounted)):
            si = SparksNotCounted[i]
            keyi = posRyR[:,cj[si]]
            tt1 = t1[si]
            if np.all(keyi>corner_ld) and np.all(keyi<corner_ru) and t1[k] <= tt1 <= t1[k]+Nt*Nsave:
                cli = Vdict[tuple(posRyR[:,cj[si]])]
                if cli in Vsparks:
                    Vsparks[cli] += 1
                else:
                    Vsparks[cli] = 1
                ToDel.append(i)
    
        if ToDel:
            NToDel = len(ToDel)
            for i in range(len(ToDel)):
                iDel = ToDel[NToDel-i-1]
                del SparksNotCounted[iDel]
    
    sparks = np.array(sparks).T

    ax = figs[ii].add_subplot(111)
    for cli in Vsparks:
        for clj in Vcluster[cli]:
            c = cm.coolwarm(float(Vsparks[cli])/7,1)
            ax.scatter(dx*clj[1], dx*clj[0], s=1, marker='s', c=c)

    #ax.scatter(dx*sparks[0,:], dx*sparks[1,:], s=10, marker='s', c=colors[ii])

    ax.set_aspect('equal', 'datalim')
    ax.set_xlim([0,80])
    ax.set_ylim([0,12])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ii += 1

for i in range(nf):
    figs[i].tight_layout()

if save:
    for i in range(nf):
        figs[i].savefig('../plots/12x80/spatial_' + wd[i] + '.pdf')

plt.show()
