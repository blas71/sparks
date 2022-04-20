import numpy as np
import matplotlib.pyplot as plt
from sys import platform

# Use latex figures
plt.rc('text', usetex = True)
plt.rc('font', **{'family': 'serif', 'serif': ['Times New Roman']})

save = 1

gg = 'CSQ'
wd2v = ['AF', 'NoAF']
wd3v = np.arange(1)+2

fig = plt.figure(figsize=(6,4))
ax1 = fig.add_subplot(2,1,1)
ax2 = fig.add_subplot(2,1,2)
axs = [ax1, ax2]

for wd2 in wd2v:
    for wd3 in wd3v:
        ff = '../' + gg + '/' + wd2 + '/' + wd2 + str(wd3)
        t, Ca, Casr, CaTnC, CaCM, CabSR, CaSrTot, CabRhod = np.loadtxt(ff + '/data/calcium.data').T
        t /= 1000
        cond = t>0
        ax1.plot(t[cond], Ca[cond], label=wd2+str(wd3), alpha=0.6)
        ax2.plot(t[cond], Casr[cond])

ax1.axvline(x=8, c='k')
ax2.axvline(x=8, c='k')

#ax1.set_ylim([0.07, 0.15])
#ax2.set_ylim([750, 950])

titles = ['Free calcium [$\mu M$]', 'SR calcium', 'currents']
i = 0
ax1.legend()
for ax in axs:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.locator_params('x', nbins=6, tight='True')
    ax.locator_params('y', nbins=6, tight='True')
    ax.set_xlabel('time [s]')
    ax.set_ylabel('calcium [$\mu M$]')
    ax.set_title(wd2v[0])
    i += 1
    
if platform == 'darwin': # mac
    fig.set_tight_layout(True)
else:
    fig.tight_layout()

if save:
    plt.savefig('../' + gg + '/plots/07concentrations.pdf')

plt.show()
