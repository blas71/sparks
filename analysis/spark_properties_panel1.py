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
folder = '_3_20'
amp_th = 0.1 # 0, 0.015, 0.05, 0.1
FDHM_th = 20

Nsave = 1000
Nx = 800
Ny = 120
dx = 0.1
Roip = 12
tdes = 2500
dt = 0.0034
Nt = 22
nconf = [4, 4]
Nprev = 3
fol = ['noAF', 'AF']
AREA_total = Nx*Ny*dx**2
Ni_max = 6001

SPARKS = {}
ll = ['AMP', 'Time to peak', 'FDHM', 'Frequency']
titles = ['Amplitude', 'Time to peak', 'FDHM', 'Frequency']
ylabel = ['$\Delta c / c_o$', 't$_{peak}$ [ms]', 'FDHM [ms]', 'f$_{sparks}$ [events/s/1000$\mu$m$^2$]']

for li in ll:
    SPARKS[li] = []

jj = 0
for foli in fol:
    wd = np.arange(nconf[jj]) + 1
    for li in ll:
        SPARKS[li].append([])
    for nci in wd:
        ff = '../12x80/' + foli + '/RyRP_CSQ_' + foli + str(nci)
        #for aaa in range(1):
        ca_traces = np.loadtxt(ff + '/data/spark_traces.data')
        cj, ryrk, t2, t1 = np.loadtxt(ff + '/data/timeRyR-open.data').T
        posRyR = np.loadtxt(ff + '/data/posRyR.data').T
        posRyR_sarco = np.loadtxt(ff + '/data/posRyR-sarco.data').T
        posRyR = np.hstack((posRyR,posRyR_sarco))
        nRyR = posRyR.shape[1]
        
        cj = cj.astype(int) - 1
        ryrk = ryrk.astype(int)
        t2 = t2.astype(int)
        t1 = t1.astype(int)
        posRyR = posRyR.astype(int)
        
        Vcluster, Vdict = clustering(0.15, 0.1, posRyR)
        
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
        nsparks = 0
        while SparksNotCounted and ii<ca_traces.shape[0]:

            k = SparksNotCounted[0]
            Ni = int(t1[k]/Nsave)-Nprev

            if Ni+Nt+1002<Ni_max:
                t = dt * Nsave * (Ni+np.arange(Nt))
                px = posRyR[1,cj[k]]
                py = posRyR[0,cj[k]]
                ROI = fun(Nx,Ny,px,py,Roip)
                nROI = len(ROI)
                
                # Calcium trace in that ROI
                Cai = ca_traces[ii,:]
        
                ## PROPERTIES ##
                # absloute amplitude and amplitude wrt to the baseline
                ind_max = np.argmax(Cai[Nprev:])+Nprev
                amp = Cai[ind_max]
                bline = np.quantile(Cai[:ind_max],0.02)
                amp_rel = amp-bline

                # FDHM
                max2 = (amp-bline)/2 + bline
                FDHM = np.sum(Cai>max2)*dt*Nsave

                if amp_rel>amp_th and FDHM>FDHM_th:
                    SPARKS['AMP'][jj].append(amp_rel)

                    nsparks += 1
                    
                    # time to peak
                    SPARKS['Time to peak'][jj].append((ind_max-Nprev)*Nsave*dt)
                    
                    # FDHM
                    max2 = (amp-bline)/2 + bline
                    SPARKS['FDHM'][jj].append(FDHM)
                
                # RyRs in the open state
                corner_ld = np.array((py-Roip,px-Roip))
                corner_ru = np.array((py+Roip,px+Roip))
                
                cl_in_fire = []
                nCaRUi = []
                nRyRsEventi = []
                ToDel = []
                for i in range(len(SparksNotCounted)):
                    si = SparksNotCounted[i]
                    keyi = posRyR[:,cj[si]]
                    tt1 = t1[si]
                    if np.all(keyi>corner_ld) and np.all(keyi<corner_ru) and Ni <= tt1/Nsave <= Ni+Nt:
                        tt2 = t2[si]
                        cl_in_fire.append([tt1,tt2])
                        nRyRsEventi.append(tuple(posRyR[:,cj[si]]))
                        nCaRUi.append(Vdict[tuple(posRyR[:,cj[si]])])
                        ToDel.append(i)
        
                if ToDel:
                    NToDel = len(ToDel)
                    for i in range(len(ToDel)):
                        iDel = ToDel[NToDel-i-1]
                        del SparksNotCounted[iDel]
                ii += 1
            else:
                break

        SPARKS['Frequency'][jj].append(nsparks)

    jj += 1

colors = ['grey', 'dodgerblue']

# fig = plt.figure(figsize=(3,2))
# ax = plt.subplot(111)

PROPERTIES = np.zeros((4,2))
ERROR = np.zeros((4,2))
for i in range(len(fol)):
    tmax = 17. * nconf[i]
    PROPERTIES[0,i] = np.mean(SPARKS['Frequency'][i]) / tmax * 1000 / AREA_total
    #ax.bar(i, x, color=colors[i])

# ax.set_title(titles[-1])
# ax.set_xticks([])
# ax.set_ylabel(ylabel[-1])
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)
# 
# fig.tight_layout()
# 
# if save:
#     plt.savefig('../../../Dropbox/PFM_Marchena/2D/plots_phospho_csq/panel_CSQ_RyR/12x80/Filter' + folder + '/subfigA/sparks-frequency.pdf')

# fig = plt.figure(figsize=(4,2))
# 
# axs = []
# for i in range(len(ll[:-1])):
#     axi = plt.subplot(1,3,i+1)
#     axs.append(axi)

for i in range(len(fol)):
    j = 0
    for li in ll[:-1]:
        x = np.mean(SPARKS[li][i])
        xerr = np.std(SPARKS[li][i])
        PROPERTIES[j+1,i] = x
        ERROR[j+1,i] = xerr
        #axs[j].bar(i, x, color=colors[i])
        #axs[j].errorbar(i, x, yerr=np.array([[0],[xerr]]), fmt='-', capsize=2, lw=0.5, markeredgewidth=0.5, color='k')
        j += 1


np.savetxt('../../../Dropbox/PFM_Marchena/2D/plots_phospho_csq/panel_CSQ_RyR/12x80/Filter' + folder + '/subfigA/spark_properties.data', PROPERTIES, fmt='%1.4f')
np.savetxt('../../../Dropbox/PFM_Marchena/2D/plots_phospho_csq/panel_CSQ_RyR/12x80/Filter' + folder + '/subfigA/spark_properties_err.data', ERROR, fmt='%1.4f')

# i = 0
# y1 = [2, 30, 30]
# y2 = [4, 60, 60]
# ylim = [4.5, 60, 75]
# for axi in axs:
#     axi.set_title(titles[i])
#     axi.set_xticks([])
#     axi.set_ylim([0, ylim[i]])
#     axi.set_yticks([0, y1[i], y2[i]])
#     axi.set_ylabel(ylabel[i])
#     axi.spines['top'].set_visible(False)
#     axi.spines['right'].set_visible(False)
#     i += 1
# 
# fig.tight_layout()
# 
# if save:
#     plt.savefig('../../../Dropbox/PFM_Marchena/2D/plots_phospho_csq/panel_CSQ_RyR/12x80/Filter' + folder + '/subfigA/sparks-properties.pdf')
# 
# plt.show()

