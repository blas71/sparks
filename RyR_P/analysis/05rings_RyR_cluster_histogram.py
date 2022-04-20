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

save = 1

def get_ring(px, py):
    p = np.array([px,py])
    dist = np.min(np.sqrt(np.sum((membrane - p)**2, axis=1)))
    return int(dist/r_ring)

wdv = ['noAF']
colors = ['dodgerblue', 'darkorange']

Nsave = 1500
Nx = 1000
Ny = 150
dx = 0.1
Roip = 12
tdes = -1
dt = 0.0034
Nt_RyR = 17*Nsave
nconf = 5 
tmax = 17 * nconf

r_ring = 7
membrane = []
for i in range(Nx):
    membrane.append((i,0))
for i in range(Nx):
    membrane.append((i,Ny-1))
for i in range(Ny):
    membrane.append((0,i))
for i in range(Ny):
    membrane.append((Ny-1,i))
membrane = np.array(membrane)

i = 0
AREA = []
while 2*i*r_ring < Ny:
    rr1 = 2*i*r_ring
    rr2 = 2*(i+1)*r_ring
    Aout = (Nx-rr1)*(Ny-rr1)
    Ain = (Nx-rr2)*(Ny-rr2)
    AREA.append(Aout-Ain)
    i += 1

AREA = dx**2 * np.array(AREA)
SPARKS = {}

fig = plt.figure(figsize=(8,3.5))
ax = plt.subplot(111)
ww = 0.3
ss = 0
jj = 0

for wd in wdv:
    for nci in range(nconf):
        nci += 1
        print(nci)
        ff = '../' + wd + '/RyR_phospho_' + wd + str(nci)
        cj, ryrk, t2, t1 = np.loadtxt(ff + '/data/timeRyR-open.data').T
        posRyR = np.loadtxt(ff + '/data/posRyR.data').T
        posRyR_sarco = np.loadtxt(ff + '/data/posRyR-sarco.data').T
        posRyR_all = np.hstack((posRyR,posRyR_sarco))
    
        cj = cj.astype(int) - 1
        ryrk = ryrk.astype(int)
        t2 = t2.astype(int)
        t1 = t1.astype(int)
        posRyR = posRyR.astype(int)
        posRyR_sarco = posRyR_sarco.astype(int)
        posRyR_all = posRyR_all.astype(int)
        
        nRyR = posRyR.shape[1]
        nRyR_sarco = posRyR_sarco.shape[1]
    
        idx_t1 = np.argsort(t1)
        t1 = t1[idx_t1]
        t2 = t2[idx_t1]
        cj = cj[idx_t1]
        
        idx2 = t1>int(tdes/dt)
        t1 = t1[idx2]
        t2 = t2[idx2]
        cj = cj[idx2]
        
        SparksNotCounted = list(np.arange(len(cj)))
    
        ii = 0
        while SparksNotCounted:
    
            k = SparksNotCounted[0]
            nRyRs_spark = 0
    
            px = posRyR_all[1,cj[k]]
            py = posRyR_all[0,cj[k]]
            ROI = fun(Nx,Ny,px,py,Roip)
            nROI = len(ROI)
    
            ring = get_ring(px, py)
            
            # RyRs in the open state
            corner_ld = np.array((py-Roip,px-Roip))
            corner_ru = np.array((py+Roip,px+Roip))
            
            ToDel = []
            for i in range(len(SparksNotCounted)):
                si = SparksNotCounted[i]
                keyi = posRyR_all[:,cj[si]]
                tt1 = t1[si]
                if np.all(keyi>corner_ld) and np.all(keyi<corner_ru) and t1[k] <= tt1 <= t1[k]+Nt_RyR:
                    nRyRs_spark += 1
                    ToDel.append(i)
        
            if ToDel:
                NToDel = len(ToDel)
                for i in range(len(ToDel)):
                    iDel = ToDel[NToDel-i-1]
                    del SparksNotCounted[iDel]
            
            if ring in SPARKS:
                SPARKS[ring].append(nRyRs_spark)
            else:
                SPARKS[ring] = [nRyRs_spark]
            ii += 1

    for k in SPARKS:
        nRyRs_spark = len(SPARKS[k])
        ax.bar(k+ss*ww*0.5, nRyRs_spark/AREA[k]/tmax, width=ww, color=colors[jj])
    
    ss += 2
    jj += 1

# ax.set_yscale('log')
# ax.set_ylim([0.3, 1200])
ax.set_yticks([0, 0.1, 0.2]) 
ax.set_title('spark density [events/s/$\mu$m$^2$]') 
# bb = nbins[:-1][aa[0]>0]
# ax.set_xticks(bb[::2])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

xlabel = []
for k in range(len(AREA)):
    xlabel.append(str(np.round(k*dx*r_ring,1)) + '-' + str(np.round((k+1)*dx*r_ring,1)))

ax.set_xticks(np.arange(len(AREA)))
ax.set_xticklabels(xlabel)
ax.set_xlabel('Distance to membrane [$\mu$m]')

fig.tight_layout()
if save:
    plt.savefig('../plots/Fig05_2_ctrl.pdf')

# fig = plt.figure(figsize=(7,9))
# 
# nbins = np.arange(31).astype(int)+1
# 
# for k in range(len(AREA)):
#     
#     ax = plt.subplot(len(AREA),1,k+1)
# 
#     nRyRs_spark = SPARKS[k]
#     aa = np.histogram(nRyRs_spark, bins=nbins)
#     ax.hist(nRyRs_spark, bins=nbins, align='left', rwidth=0.8)
# 
#     ax.set_yscale('log')
#     # ax.set_ylim([0.3, 1200])
#     ax.set_yticks([10, 1000]) 
#     # bb = nbins[:-1][aa[0]>0]
#     # ax.set_xticks(bb[::2])
#     ax.spines['top'].set_visible(False)
#     ax.spines['right'].set_visible(False)
#     xlabel = str(np.round(k*dx*r_ring,1)) + '-' + str(np.round((k+1)*dx*r_ring,1))
#     ax.set_ylabel(xlabel, rotation=0, labelpad=20)
# 
# fig.tight_layout()
# if save:
#     plt.savefig('../plots/Fig05_3.pdf')

plt.show()
