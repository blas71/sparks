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
fraction = False
folder = '_3_20'
amp_th = 0.1 # 0, 0.015, 0.05, 0.1
FDHM_th = 20

def get_ring(px, py):
    p = np.array([px,py])
    dist = np.min(np.sqrt(np.sum((membrane - p)**2, axis=1)))
    return int(dist/r_ring)

wd2 = 'RyRP_'
wdv3 = ['noAF', 'AF']

colors = ['grey', 'dodgerblue']

Nsave = 1000
Nx = 800
Ny = 120
dx = 0.1
Roip = 12
tdes = 2500
dt = 0.0034
Nt = 22
nconf = [6, 6]
Nprev = 3 

r_ring = 10
membrane = []
for i in range(Nx):
    membrane.append((i,0))
for i in range(Nx):
    membrane.append((i,Ny-1))
for i in range(Ny):
    membrane.append((0,i))
for i in range(Ny):
    membrane.append((Nx-1,i))
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

fig = plt.figure(figsize=(4,2))
ax = plt.subplot(111)
ww = 0.4
ss = -1
jj = 0

AAA = []

for wd3 in wdv3:
    tmax = (17. - tdes/1000.) * nconf[jj]
    SPARKS = {}
    AAA.append([])
    for nci in range(nconf[jj]):
        nci += 1
        print(nci)
        ff = '../RyR_P/12x80/' + wd3 + '/' + wd2 + wd3 + str(nci)
        ca_traces = np.loadtxt(ff + '/data/spark_traces.data')
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
        #ca_traces = ca_traces[idx2,:] 
        
        SparksNotCounted = list(np.arange(len(cj)))
    
        ii = 0
        while SparksNotCounted and ii<ca_traces.shape[0]:
    
            k = SparksNotCounted[0]
            Ni = int(t1[k]/Nsave)-Nprev
            nRyRs_spark = 0
    
            px = posRyR_all[1,cj[k]]
            py = posRyR_all[0,cj[k]]
            ROI = fun(Nx,Ny,px,py,Roip)
            nROI = len(ROI)

            # Calcium trace in that ROI
            Cai = ca_traces[ii,:]
    
            ring = get_ring(px, py)

            # absloute amplitude and amplitude wrt to the baseline
            ind_max = np.argmax(Cai[Nprev:])+Nprev
            amp = Cai[ind_max]
            bline = np.quantile(Cai[:ind_max],0.02)
            amp_rel = amp-bline
            AAA[jj].append(amp_rel)

            # FDHM
            max2 = (amp-bline)/2 + bline
            FDHM = np.sum(Cai>max2)*dt*Nsave
            
            # RyRs in the open state

            corner_ld = np.array((py-Roip,px-Roip))
            corner_ru = np.array((py+Roip,px+Roip))
            
            ToDel = []
            for i in range(len(SparksNotCounted)):
                si = SparksNotCounted[i]
                keyi = posRyR_all[:,cj[si]]
                tt1 = t1[si]
                #if np.all(keyi>corner_ld) and np.all(keyi<corner_ru) and t1[k] <= tt1 <= t1[k]+Nt_RyR:
                if np.all(keyi>corner_ld) and np.all(keyi<corner_ru) and Ni <= tt1/Nsave <= Ni+Nt:
                    nRyRs_spark += 1
                    ToDel.append(i)
        
            if ToDel:
                NToDel = len(ToDel)
                for i in range(len(ToDel)):
                    iDel = ToDel[NToDel-i-1]
                    del SparksNotCounted[iDel]
            
            if amp_rel>amp_th and FDHM>FDHM_th:
                if ring in SPARKS:
                    SPARKS[ring].append(nRyRs_spark)
                else:
                    SPARKS[ring] = [nRyRs_spark]
            ii += 1

    if fraction:
        if jj==0:
            nRyRs_spark = len(SPARKS[5])
            NORM = nRyRs_spark/AREA[5]
    else:
        NORM=tmax
    
    save_array = np.zeros(6)
    for k in SPARKS:
        nRyRs_spark = len(SPARKS[k])
        #ax.bar(k+ss*ww*0.5, nRyRs_spark/AREA[k]/NORM, width=ww, color=colors[jj])
        save_array[k] = nRyRs_spark/AREA[k]/NORM
    if fraction:
        np.savetxt('../RyR_P/12x80/data/Filter' + folder + '/fraction_' + wd3 + '.data', save_array, fmt='%1.3f')
    else:
        np.savetxt('../RyR_P/12x80/data/Filter' + folder + '/frequency_' + wd3 + '.data', save_array, fmt='%1.5f')

    if jj==0:
        SPARKS_noAF = SPARKS
    ss += 2
    jj += 1

# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)
# 
# xlabel = []
# for k in range(len(AREA)):
#     xlabel.append(str(int(k*dx*r_ring)) + '-' + str(int((k+1)*dx*r_ring)))
# 
# ax.set_xticks(np.arange(len(AREA)))
# ax.set_xticklabels(xlabel)
# ax.set_xlabel('Distance to membrane [$\mu$m]')
# if fraction:
#     ax.set_ylabel('Fraction') 
#     ax.set_ylim([0, 20])
#     ax.set_yticks([0, 10, 20])
# else:
#     #ax.set_ylim([0, 12])
#     ax.set_ylim([0, 6])
#     ax.set_yticks([0, 3, 6])
#     ax.set_ylabel('spark density \n [events/s/1000$\mu$m$^2$]') 
# 
# fig.tight_layout()
# if save:
#     if fraction:
#         plt.savefig('../../../Dropbox/PFM_Marchena/2D/plots_phospho_csq/panel_CSQ_RyR/12x80/Filter' + folder + '/subfigB/Fig05_CSQ_RYR_fraction.pdf')
#     else:
#         plt.savefig('../../../Dropbox/PFM_Marchena/2D/plots_phospho_csq/panel_CSQ_RyR/12x80/NoFilter/subfigB/Fig05_CSQ_RYR_3.pdf')

#plt.show()
