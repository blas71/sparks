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

save = 0
#wdv = ['NCX', 'LCC']
wdv = ['CTRL']

fig = plt.figure(figsize=(4,3.2))
ax = plt.subplot(111)

nmax = 17
nbins = np.arange(nmax+1).astype(int)+1
ww = 0.3   
ss = -1

color = ['skyblue', 'limegreen']
kk = 0

Nsave = 1500
Nx = 1000
Ny = 150
dx = 0.1
Roip = 12
tdes = -1
dt = 0.0034
Nt_RyR = 17*Nsave
nconf = 5 

for wd in wdv:
    SPARKS = []
    TMAX = []

    for nci in range(nconf):
        nci += 1
        print(nci)
        #direc = '../' + wd + '/' + wd + '_' + str(nci) 
        direc = '../' + wd + '/' + wd + str(nci) 
        cj, ryrk, t2, t1 = np.loadtxt(direc + '/data/timeRyR-open.data').T
        posRyR = np.loadtxt(direc + '/data/posRyR.data').T
        posRyR_sarco = np.loadtxt(direc + '/data/posRyR-sarco.data').T
        posRyR = np.hstack((posRyR,posRyR_sarco))
        nRyR = posRyR.shape[1]

        # t, ca, casr, cab1, cab2, cab3, cab4 = np.loadtxt(direc + '/data/calcium.data').T
        # TMAX.append(t[-1]-tdes)
    
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
        
        SparksNotCounted = list(np.arange(len(cj)))
    
        ii = 0
        while SparksNotCounted:
    
            k = SparksNotCounted[0]
            nRyRs_spark = 0
    
            px = posRyR[1,cj[k]]
            py = posRyR[0,cj[k]]
            ROI = fun(Nx,Ny,px,py,Roip)
            nROI = len(ROI)
            
            # RyRs in the open state
            corner_ld = np.array((py-Roip,px-Roip))
            corner_ru = np.array((py+Roip,px+Roip))
            
            ToDel = []
            for i in range(len(SparksNotCounted)):
                si = SparksNotCounted[i]
                keyi = posRyR[:,cj[si]]
                tt1 = t1[si]
                if np.all(keyi>corner_ld) and np.all(keyi<corner_ru) and t1[k] <= tt1 <= t1[k]+Nt_RyR:
                    nRyRs_spark += 1
                    ToDel.append(i)
        
            if ToDel:
                NToDel = len(ToDel)
                for i in range(len(ToDel)):
                    iDel = ToDel[NToDel-i-1]
                    del SparksNotCounted[iDel]
            
            SPARKS.append(nRyRs_spark)
    
            ii += 1
     
    hist_sparks, edges = np.histogram(SPARKS, bins=nbins)
    ax.bar(nbins[:-1] + ss*ww*0.5, hist_sparks, width=ww, color=color[kk])
    ss += 2
    kk += 1

ax.set_xticks([1, 5, 10, 15])
ax.set_yscale('log')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_title('open RyR2')

fig.tight_layout()

if save:
    plt.savefig('../plots_part3/Fig02.pdf')
plt.show()
