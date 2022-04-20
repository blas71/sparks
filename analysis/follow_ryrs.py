import numpy as np
import matplotlib.pyplot as plt

# Use latex figures
plt.rc('text', usetex = True)
plt.rc('font', **{'family': 'serif', 'serif': ['Times New Roman']})

Nx = 800
Ny = 120
Rroi = 0.7 # microns
dx = 0.1
Roip2 = int(np.round(Rroi / dx))**2
dt = 0.012
tdes = 8000 #ms
Nt = int(10/dt)

#Type = 'AF'
#dist_folder = ['RyR_P']
Type = 'CTRL'
dist_folder = ['.']

nconf = np.arange(6) + 9

for nci in nconf:
    Vsparks = {}
    print(nci)
    for foli in dist_folder:
        Vsparks[foli] = {}
    
        #wd1 = foli + '/CONTINUATION/' + Type + str(nci)
        wd1 = foli + '/' + Type + '/' + Type + str(nci)
    
        cj, ryrk, t2, t1 = np.loadtxt('../' + wd1 + '/data/timeRyR-open.data').T
        posRyR = np.loadtxt('../' + wd1 + '/data/posRyR.data').T
        posRyR_sarco = np.loadtxt('../' + wd1 + '/data/posRyR-sarco.data').T
        
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
        
        SparksNotCounted = list(np.arange(len(cj)))
    
        while SparksNotCounted:
            k = SparksNotCounted[0]
            Vsparks[foli][k] = []
            contour = [k]
            ToDel = [0]
            while contour:
                sj = contour[0]
                key_sj = posRyR[:,cj[sj]]
                Nsj = t1[sj]
                for i in range(len(SparksNotCounted)):
                    si = SparksNotCounted[i]
                    keyi = posRyR[:,cj[si]]
                    Ni = t1[si]
                    if Ni>Nsj and Ni<Nsj+Nt:
                        dist2 = np.sum((key_sj-keyi)**2)
                        if dist2<Roip2:
                            Vsparks[foli][k].append(si)
                            contour.append(si)
                            ToDel.append(i)
                del contour[0]
                contour = list(set(contour))
        
            ToDel = list(set(ToDel))
            ToDel.sort()
            if ToDel:
                NToDel = len(ToDel)
                for i in range(len(ToDel)):
                    iDel = ToDel[NToDel-i-1]
                    del SparksNotCounted[iDel]
        
            Vsparks[foli][k] = list(set(Vsparks[foli][k]))
            Vsparks[foli][k].sort()

        ff = open('../' + wd1 + '/data/follow_RyR_1.data', 'w')
        for k in Vsparks[foli]:
            line = str(k)
            for sj in Vsparks[foli][k]:
                line += ", " + str(sj)
            ff.write(line + " \n")
        ff.close()
