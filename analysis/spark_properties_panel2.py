import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import os

# Use latex figures
plt.rc('text', usetex = True)
plt.rc('font', **{'family': 'serif', 'serif': ['Times New Roman']})

save = 1

amp_th = 0
FDHM_th = 10.
FWHM_th = 0.
bl_th = 0
nryrs_th = 0

amp_th_max = 10
FDHM_th_max = 120
FWHM_th_max = 300
time_to_peak_th_max = 60
bl_th_max = 0.3

# gg = 'RyR_P'
# fol = ['NoAF', 'AF']
# label = ['no AF', 'AF']
# nconf = {fol[0]: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17],
#          fol[1]: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]}
gg = '.'
fol = ['CTRL']
label = ['CTRL']
nconf = {fol[0]: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]}

tmax = {}
for foli in fol:
    tmax[foli] = 28. * len(nconf[foli])

PROPERTIES = {}

for foli in fol:
    for i in nconf[foli]:
        ff = '../' + gg + '/' + foli + '/' + foli + str(i)
        PROPERTIES_i = np.loadtxt(ff + '/data/spark_properties_1.data')
        print(foli, PROPERTIES_i.shape)
        if foli not in PROPERTIES:
            PROPERTIES[foli] = PROPERTIES_i
        else:
            PROPERTIES[foli] = np.hstack((PROPERTIES[foli], PROPERTIES_i))

for foli in fol:
    x = PROPERTIES[foli]
    cond = (x[0,:]>FWHM_th) & (x[1,:]>amp_th) & (x[2,:]>FDHM_th) & (x[0,:]<FWHM_th_max) & (x[1,:]<amp_th_max) & (x[2,:]<FDHM_th_max) & (x[3,:]<time_to_peak_th_max) & (x[5,:]<bl_th_max) & (x[3,:]>bl_th) & (x[4,:]>nryrs_th)
    PROPERTIES[foli] = PROPERTIES[foli][:,cond]

idx = [0, 1, 2, 3, 5]
nbins = []
number = [25, 25, 15, 7, 20]
for i in range(len(idx)):
    mins = []
    maxs = []
    for foli in fol:
        mins.append(min(PROPERTIES[foli][idx[i],:]))
        maxs.append(max(PROPERTIES[foli][idx[i],:]))
    mins = min(mins)
    maxs = max(maxs)
    nbins.append(np.linspace(mins, maxs, number[i]))

titles = ['FWHM', 'Amplitude', 'FDHM', 'Time to peak', 'Baseline']
ylabel = ['FWHM [$\mu$m]', '$\Delta c /c_o$ [$\mu$M]', 'FDHM [ms]', 't$_{peak}$ [ms]', 'baseline [$\mu$M]']

color = ['darkgray', 'royalblue']
aligns = ['left', 'right']
fig = plt.figure(figsize=(6,4))
for i in range(5):
    X = []
    ax = plt.subplot(3,2,i+1)
    for foli in fol:
        X.append(PROPERTIES[foli][idx[i],:])
    ax.hist(X, bins=nbins[i], color=color[:len(fol)], align='left', rwidth=0.8, edgecolor='black', linewidth=0.7)

    if i==4:
        ax.set_yscale('log')
        ax.set_yticks([1, 10, 100, 1000])
    ax.set_title(titles[i])
    ax.set_xlabel(ylabel[i])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

fig.tight_layout()

########

fig4 = plt.figure(figsize=(7,5))
fig5 = plt.figure(figsize=(7,5))
fig6 = plt.figure(figsize=(6,5))

for i in range(5):
    ax1 = fig4.add_subplot(3,2,i+1)
    ax2 = fig5.add_subplot(3,2,i+1)
    ax3 = fig6.add_subplot(3,2,i+1)
    jj = -1
    for j in range(len(fol)):
        x = PROPERTIES[fol[j]][idx[i],:]
        nryrs = PROPERTIES[fol[j]][4,:]
        cond = nryrs<45
        ax1.scatter(nryrs+0.15*jj, x, color=color[j], s=10, edgecolor='black', linewidth=0.1, alpha=0.5)
        ax2.scatter(nryrs[cond]+0.15*jj, x[cond], color=color[j], s=10, edgecolor='black', linewidth=0.1, alpha=0.5)
        nryrs_set = list(set(nryrs))
        for xset_i in nryrs_set:
            yvalue = x[nryrs==xset_i]
            ax3.scatter(xset_i+0.15*jj, np.mean(yvalue), s=20, color=color[j])
        jj += 2

    #ax3.set_xlim([0,26])
    for ax in [ax1, ax2, ax3]:
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.set_title(titles[i])
        ax.set_xlabel('$N_{RyRs}$')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

fig4.tight_layout()
fig5.tight_layout()
fig6.tight_layout()

########

ww = 0.4
fig2 = plt.figure(figsize=(4,4))
for i in range(4):
    ax = fig2.add_subplot(2,2,i+1)
    k = 0
    for foli in fol:
        x = np.mean(PROPERTIES[foli][i])
        xerr = np.std(PROPERTIES[foli][i])
        ax.bar(k, x, yerr=xerr, color=color[k], width=ww, error_kw=dict(lw=0.5, capsize=4, capthick=0.5))
        k += 1
    ax.set_title(titles[i])
    ax.set_xlim([-1,2])
    ax.set_xticks([])
    ax.set_ylabel(ylabel[i])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

fig2.tight_layout()

########

fig3 = plt.figure(figsize=(4,3))
ax = plt.subplot(111)

X = []
for foli in fol:
    x = PROPERTIES[foli][4]
    X.append(x)

nbins = np.arange(max(X[-1])+2)+1
nbins = nbins[::2]

ww = 0.4 * (nbins[1]-nbins[0])
i = -1
j = 0
for foli in fol:
    hist, ssss = np.histogram(X[j], bins=nbins)
    #bins, hist = ax.hist(X, bins=nbins, color=color[:len(fol)], align='left', rwidth=0.8, edgecolor='black', linewidth=0.7, label=label)
    ax.bar(nbins[:-1]-i*ww*0.5, hist/tmax[foli], color=color[j], width=ww, edgecolor='black', linewidth=0.7, label=label[j])
    i += 2
    j += 1

ax.set_xticks(nbins[:-1])
#ax.set_xticklabels((nbins[::3]+1).astype(int))
ax.set_yscale('log')
ax.set_xlabel('O$_{RyRs}$')
ax.set_ylabel('events/s')
ax.legend()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

fig3.tight_layout()

for foli in fol:
    print(foli + '- frequency:', PROPERTIES[foli].shape[1]/tmax[foli])

if save:
    # folder = '../plots/' + gg + '/amp' + str(amp_th) + '_FDHM' + str(FDHM_th) + '_nRyRs' + str(nryrs_th)
    # print(folder)
    # os.mkdir(folder)
    # filefreq = open(folder + '/freq.txt', 'w')
    # for foli in fol:
    #     filefreq.write(foli + '- frequency:' + str(PROPERTIES[foli].shape[1]/tmax[foli]) + ' \n')
    # filefreq.close()
    folder = '../CTRL/plots'
    #folder = '../' + gg + '/plots'
    fig.savefig(folder + '/09properties_dist.pdf')
    fig2.savefig(folder + '/09properties_mean.pdf')
    fig3.savefig(folder + '/09properties_RyRs.pdf')
    fig4.savefig(folder + '/09properties_RyRs_scatter.pdf')
    fig5.savefig(folder + '/09properties_RyRs_scatter_zoom.pdf')
    fig6.savefig(folder + '/09properties_RyRs_scatter_mean.pdf')

plt.show()
