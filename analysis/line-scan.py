import numpy as np
import matplotlib.pyplot as plt

gg = 'RyR_P'
fol = ['NoAF', 'AF']
Nt = 300
Nx = 800
Ny = 120
Ni = 3000

fig = plt.figure(figsize=(5,5))
ax = plt.subplot(111)

for Type in fol:
    wd1 = '../' + gg + '/' + Type + '/' + Type + '1'
    path_cab = wd1 + '/Codes' + Type + '/Cabc'
    line_scan = np.zeros(Ny)
    for i in range(Nt):
        cab = np.loadtxt(path_cab + str(i+Ni).zfill(4) + '.data').T
        cab = np.reshape(cab, (Nx,Ny))
        line_scan += np.mean(cab[500:600,:], axis=0)/Nt
    ax.plot(line_scan, label=Type)

ax.legend()
plt.show()

    
