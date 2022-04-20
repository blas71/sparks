import numpy as np
import matplotlib.pyplot as plt
import imp
import RyR_cluster
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

Nx = 1000
Ny = 150
dx = 0.1
Roip = 12
dt = 0.0034
tdes = 2000 # ms

Nt_RyR = int(70/dt)
nconf = 5
wd = np.arange(nconf) + 1

FREQ = []
TMAX = []

for nci in wd:
    print(nci)
    #ff = '../CTRL/CTRL' + str(nci)
    ff = '../RyR_P/AF/RyRP_AF' + str(nci)
    cj, ryrk, t2, t1 = np.loadtxt(ff + '/data/timeRyR-open.data').T
    posRyR = np.loadtxt(ff + '/data/posRyR.data').T
    posRyR_sarco = np.loadtxt(ff + '/data/posRyR-sarco.data').T
    posRyR_all = np.hstack((posRyR,posRyR_sarco))
    nRyR = posRyR.shape[1]

    t, ca, casr, cab1, cab2, cab3, cab4 = np.loadtxt(ff + '/data/calcium.data').T
    TMAX.append(t[-1]-tdes)

    cj = cj.astype(int) - 1
    ryrk = ryrk.astype(int)
    t2 = t2.astype(int)
    t1 = t1.astype(int)
    posRyR = posRyR.astype(int)
    posRyR_sarco = posRyR_sarco.astype(int)
    posRyR_all = posRyR_all.astype(int)
    
    idx_t1 = np.argsort(t1)
    t1 = t1[idx_t1]
    t2 = t2[idx_t1]
    cj = cj[idx_t1]
    
    idx2 = t1>int(tdes/dt)
    t1 = t1[idx2]
    t2 = t2[idx2]
    cj = cj[idx2]
    
    SparksNotCounted = list(np.arange(len(cj)))
    nsparks = 0

    ii = 0
    while SparksNotCounted:

        k = SparksNotCounted[0]
        nsparks += 1

        px = posRyR_all[1,cj[k]]
        py = posRyR_all[0,cj[k]]
        ROI = fun(Nx,Ny,px,py,Roip)
        nROI = len(ROI)

        # RyRs in the open state
        corner_ld = np.array((py-Roip,px-Roip))
        corner_ru = np.array((py+Roip,px+Roip))
        
        ToDel = []
        for i in range(len(SparksNotCounted)):
            si = SparksNotCounted[i]
            keyi = posRyR_all[:,cj[si]]
            tt1 = t1[si]
            if np.all(keyi>corner_ld) and np.all(keyi<corner_ru) and t1[k] <= tt1 <= t1[k]+Nt_RyR:
                ToDel.append(i)
    
        if ToDel:
            NToDel = len(ToDel)
            for i in range(len(ToDel)):
                iDel = ToDel[NToDel-i-1]
                del SparksNotCounted[iDel]
        
        ii += 1

    FREQ.append(nsparks)

TMAX = np.array(TMAX) / 1000
FREQ = np.array(FREQ)
print('frequency:', np.sum(FREQ)/np.sum(TMAX), 'sparks/s')
