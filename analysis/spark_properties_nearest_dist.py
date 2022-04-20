import numpy as np
import matplotlib.pyplot as plt

# Use latex figures
plt.rc('text', usetex = True)
plt.rc('font', **{'family': 'serif', 'serif': ['Times New Roman']})

save = 0

amp_th = 0
FDHM_th = 10
FWHM_th = 0.5
bl_th = 0
nryrs_th = 0

amp_th_max = 10
FDHM_th_max = 120
FWHM_th_max = 300
time_to_peak_th_max = 60
bl_th_max = 0.3

gg = 'CSQ'
fol = ['NoAF', 'AF']
label = ['no AF', 'AF']
nconf = {fol[0]: [1, 2, 3, 4, 5],
         fol[1]: [1, 2, 3, 4]}

dx = 0.1
Ny = 120
Nx = 800

PROPERTIES = {}

for foli in fol:
    PROPERTIES[foli] = {}
    for i in nconf[foli]:
        ff = '../' + gg + '/' + foli + '/' + foli + str(i)
        PROPERTIES_i = np.loadtxt(ff + '/data/spark_properties_1.data')
        PROPERTIES[foli][i] = PROPERTIES_i

for foli in fol:
    for i in nconf[foli]:
        x = PROPERTIES[foli][i]
        cond = (x[0,:]>FWHM_th) & (x[1,:]>amp_th) & (x[2,:]>FDHM_th) & (x[0,:]<FWHM_th_max) & (x[1,:]<amp_th_max) & (x[2,:]<FDHM_th_max) & (x[3,:]<time_to_peak_th_max) & (x[5,:]<bl_th_max) & (x[3,:]>bl_th) & (x[4,:]>nryrs_th)
        PROPERTIES[foli][i] = PROPERTIES[foli][i][:,cond]


ll = ['position']
titles = ['Distance to nearest spark']
ylabel = ['Distance [$\mu$m]']

count = {}
dist_min = {}
bw = 10
nbars = int(Ny/2/bw)
kk = 0

for foli in fol:
    count[foli] = np.zeros(nbars)
    dist_min[foli] = np.zeros(nbars)
    for nf in nconf[foli]:
        print(foli, nf)
        positions = np.vstack((PROPERTIES[foli][nf][7], PROPERTIES[foli][nf][6]))
        positions = positions.T
        nsparks = positions.shape[0]
        for i in range(len(positions)):
            dist = np.zeros(len(positions)-1)
            posi = positions[i,:]
            px, py = positions[i,:]
        
            # distance to membrane
            pxx = (px, Nx-1-px)
            pyy = (py, Ny-1-py)
            dist_mem = np.min((np.min(pxx),np.min(pyy)))
            index = int(dist_mem/bw)
        
            jj = 0
            posi = positions[i,:]
            cond = np.ones(nsparks).astype(bool)
            cond[i] = False
            dist = np.sqrt(np.sum((posi - positions[cond,:])**2, axis=1))
            dist = dist[dist>2]
            dist_min[foli][index] += dx*min(dist)
            count[foli][index] += 1
        kk += 1

fig = plt.figure(figsize=(4,3))
ax = plt.subplot(111)

colors = ['gainsboro', 'dodgerblue']

i = 0
ss = -1
for foli in fol:
    #ax.bar(np.arange(nbars)+0.4*ss*0.5, dist_min[foli]/count[foli], width=0.4, color=colors[i])
    ax.plot(np.arange(nbars), dist_min[foli]/count[foli], color=colors[i])
    i += 1
    ss += 2

xlabel = []
for k in range(nbars):
    xlabel.append(str(int(k*dx*bw)) + '-' + str(int((k+1)*dx*bw)))

#ax.set_ylim([0, 0.2])
#ax.set_yticks([0, 1, 2])
ax.set_xticks(np.arange(nbars))
ax.set_xticklabels(xlabel)
ax.set_xlabel('Distance to membrane [$\mu$m]')
ax.set_ylabel('Distance to nearest spark [$\mu$m]')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

fig.tight_layout()

if save:
    plt.savefig('../' + gg + '/plots/spark-nearest-dist.pdf')

plt.show()

