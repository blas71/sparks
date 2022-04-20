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
Nx = 1000
Ny = 150
Roip = 12
tdes = -1 
Nprev = 4

wd = ['noAF', 'AF']
nf = len(wd)

figs = []
for i in range(nf):
    figs.append(plt.figure(figsize=(4,3)))

ii = 0
for wdi in wd:
    ff = '../RyR_P/' + wdi + '/FILMS/RyRP_' + wdi + '1'
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
    
    idx_t1 = np.argsort(t1)
    t1 = t1[idx_t1]
    t2 = t2[idx_t1]
    cj = cj[idx_t1]
    
    idx2 = t1>int(tdes/dt)
    t1 = t1[idx2]
    t2 = t2[idx2]
    cj = cj[idx2]

    t = []
    Ca_all = []
    count = 0
    cci = 0

    SparksNotCounted = list(np.arange(len(cj)))
    
    while SparksNotCounted:
    
        k = SparksNotCounted[0]
    
        px = posRyR[1,cj[k]]
        py = posRyR[0,cj[k]]
        ROI = fun(Nx,Ny,px,py,Roip)
        nROI = len(ROI)

        ti = int(t1[k]/Nsave)

        if ii==1:
            print(count)
            if np.mod(count,7)==0:
                if dt*Nsave*(ti+Nt)<tf[-1] and dt*Nsave*(ti-Nprev)>0:
                    ROI = fun(Nx,Ny,px,py,Roip)
                    nROI = len(ROI)
                    Cai = []
                    for i in range(Nt):
                        ca = np.loadtxt(ff + '/codes/fort.' + str(1002+ti+i-Nprev)).T
                        ca = np.reshape(ca, (Nx,Ny))
                        ca_mean = 0
                        for j in ROI:
                            ca_mean += ca[j] / nROI
                        Cai.append(ca_mean)
                    t.append((np.arange(Nt)-Nprev) * dt * Nsave + t1[k]*dt)
                    Ca_all.append(Cai)
                    cci += 1
        else:
            if dt*Nsave*(ti+Nt)<tf[-1] and dt*Nsave*(ti-Nprev)>0:
                ROI = fun(Nx,Ny,px,py,Roip)
                nROI = len(ROI)
                Cai = []
                for i in range(Nt):
                    ca = np.loadtxt(ff + '/codes/fort.' + str(1002+ti+i-Nprev)).T
                    ca = np.reshape(ca, (Nx,Ny))
                    ca_mean = 0
                    for j in ROI:
                        ca_mean += ca[j] / nROI
                    Cai.append(ca_mean)
                t.append((np.arange(Nt)-Nprev) * dt * Nsave + t1[k]*dt)
                Ca_all.append(Cai)
                cci += 1
        count += 1
    
        # RyRs in the open state
        corner_ld = np.array((py-Roip,px-Roip))
        corner_ru = np.array((py+Roip,px+Roip))
        
        ToDel = []
        for i in range(len(SparksNotCounted)):
            si = SparksNotCounted[i]
            keyi = posRyR[:,cj[si]]
            tt1 = t1[si]
            if np.all(keyi>corner_ld) and np.all(keyi<corner_ru) and t1[k] <= tt1 <= t1[k]+Nt*Nsave:
                ToDel.append(i)
    
        if ToDel:
            NToDel = len(ToDel)
            for i in range(len(ToDel)):
                iDel = ToDel[NToDel-i-1]
                del SparksNotCounted[iDel]

    ax = figs[ii].add_subplot(111)
    for i in range(cci):
        if ii==0:
            c = cm.Greys(float(i)/cci,1)
        else:
            c = cm.Blues(float(i)/cci,1)
        ax.plot(t[i], np.array(Ca_all[i]) + 0.1*i, color=c)

    ax.set_xlabel('t [ms]')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # ax.spines['left'].set_visible(False)
    # ax.spines['bottom'].set_visible(False)
    
    ii += 1

for i in range(nf):
    figs[i].tight_layout()

if save:
    for i in range(nf):
        figs[i].savefig('../plots/Ca_traces_' + wd[i] + '.pdf')

plt.show()
