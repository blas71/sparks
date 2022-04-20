import numpy as np
import matplotlib.pyplot as plt

# Use latex figures
plt.rc('text', usetex = True)
plt.rc('font', **{'family': 'serif', 'serif': ['Times New Roman']})

wd2v = ['AF', 'noAF']
    
fig = plt.figure(figsize=(6,4))
save = 1

for wd2 in wd2v:
    ff = '../CSQ_RyRp/WithPotential/' + wd2 + '/RyRP_CSQ_' + wd2 + '2'
    
    posRyR = np.loadtxt(ff + '/data/posRyR.data').T
    posRyR_sarco = np.loadtxt(ff + '/data/posRyR-sarco.data').T
    posRyR = np.hstack((posRyR,posRyR_sarco))
    nRyR = 9*posRyR.shape[1]
    
    cj, kk, t2, t1 = np.loadtxt(ff + '/data/timeRyR-open.data').T
    
    t1 = (t1/100).astype(int)
    t2 = (t2/100).astype(int)
    
    Nt = max(t2)+1
    t = np.arange(Nt)*100*0.0034/1000
    
    OPEN = np.zeros(Nt)
    
    for i in range(len(t1)):
        OPEN[t1[i]:t2[i]] += 1
    
    cond = (t>12)
    plt.plot(t[cond], OPEN[cond]/nRyR*100, label=wd2, alpha=0.7)

plt.legend()
plt.xlabel('t [ms]')
plt.ylabel('Fraction $O_{RyR}$ ($\%$)')
plt.yticks([0, 0.25, 0.5])

if save:
    plt.savefig('../plots/12x80/OPEN_RYR_2Hz.pdf')
plt.show()
