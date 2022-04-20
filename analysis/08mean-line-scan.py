import numpy as np
import matplotlib.pyplot as plt

# Use latex figures
plt.rc('text', usetex = True)
plt.rc('font', **{'family': 'serif', 'serif': ['Times New Roman']})

Mi = 1000
Mj = 150
dx = 0.1

dt = 0.0034 * 1500 / 1000
Nt = 480
#Ni = 4297
Ni = 10002

save = 0

dt2 = 0.0034 * 100 / 1000
Ni2 = (Ni-1002) * int(dt/dt2)
Nt2 = int(dt*Nt/dt2)
t2 = np.arange(Nt2)*dt2

# wd = ['05', '25', '65']
# Li = [0.5, 2.5, 6.5]
wd1 = ['with_kCaSR', 'without_kCaSR']
#wd1 = ['SR_diff']
nconf = 2
nf = len(wd1)

k = 0

fig1 = plt.figure(1, figsize=(7,2))
# fig2 = plt.figure(2, figsize=(7,2))
# fig3 = plt.figure(3, figsize=(9,2))

for j in range(nf):
    ca_scan = np.zeros((Nt,Mj))
    for i in range(Nt):
        ca = np.loadtxt('../RyR_P/AF/' + wd1[j] + '/RyRP_AF1/codes/fort.' + str(i+Ni)).T
        ca = ca.reshape((Mi, Mj))
        ca_scan[i,:] = np.sum(ca[400:800,:], axis=0)/400
    
    ax1 = fig1.add_subplot(1,nf,j+1)
    surf = ax1.imshow(ca_scan.T, origin="lower", vmin=300, vmax=700, aspect='auto', cmap=plt.get_cmap('jet'), extent=[0, dt*Nt, 0, Mj*dx])
    #ax1.set_title(wd1[j])
    #if j==1:
    #    fig1.colorbar(surf, ticks=[0, 0.05], aspect=10)
    #ax1.set_xticks([0, 0.2, 0.4])
    ax1.set_xlabel('t [s]')
    ax1.set_ylabel(r'Y [$\mu$m]')

    # ax2 = fig2.add_subplot(1,nf,j+1)
    # SS_cy, SS_sr, CT_cy, CT_sr = np.loadtxt('L0' + wd[j] + '/data/calcium-SS-CT.data').T
    # ax2.plot(t2, SS_cy[Ni2:Ni2+Nt2], c='r', lw=1)
    # ax2.plot(t2, CT_cy[Ni2:Ni2+Nt2], c='k', lw=1)
    # ax2.set_xticks([0, 0.2, 0.4])
    # ax2.set_yticks([0, 0.5, 1])
    # ax2.set_xlabel('t [s]')
    # ax2.set_ylabel(r'$c_i$ [$\mu$M]')
    # x1, x2, y1, y2 = ax2.axis()
    # ax2.axis([x1, t2[-1], 0, 1.4])

    # # Hide the right and top spines
    # ax2.spines['right'].set_visible(False)
    # ax2.spines['top'].set_visible(False)
    # 
    # # Only show ticks on the left and bottom spines
    # ax2.yaxis.set_ticks_position('left')
    # ax2.xaxis.set_ticks_position('bottom')

# fig1.tight_layout()
# fig2.tight_layout()
# fig3.tight_layout()

if save:
    fig1.savefig('L005/RyR2_5/plots/08line-scan-long.pdf')
    fig2.savefig('L005/RyR2_5/plots/08line-scan-traces.pdf')
    fig3.savefig('L005/RyR2_5/plots/08line-scan-TT.pdf')

plt.show()
