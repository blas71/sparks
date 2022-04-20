import numpy as np
import matplotlib.pyplot as plt
import random

Nx = 15
Ny = 100

Vnn = {}
for i in range(Nx-2):
    i += 1
    for j in range(Ny-2):
        j += 1
        Vnn[(i,j)] = [(i+1,j), (i-1,j), (i,j+1), (i,j-1)]

ii = [1, Nx-2]
k = 0
for i in [0,Nx-1]:
    for j in range(Ny-2):
        j += 1
        Vnn[(i,j)] = [(ii[k],j), (i,j+1), (i,j-1)]
    k += 1

jj = [1, Ny-2]
k = 0
for j in [0,Ny-1]:
    for i in range(Nx-2):
        i += 1
        Vnn[(i,j)] = [(i-1,j), (i+1,j), (i, jj[k])]
    k += 1

Vnn[(0,0)] = [(1,0), (0,1)]
Vnn[(Nx-1,0)] = [(Nx-2,0), (Nx-1,1)]
Vnn[(0,Ny-1)] = [(1,Ny-1), (0,Ny-2)]
Vnn[(Nx-1,Ny-1)] = [(Nx-2,Ny-1), (Nx-1,Ny-2)]

latt = np.zeros((Nx,Ny))

nCaRU = 10
nCLmax = 20 # maximum number of clusters
ci = 0 

CLidx = []

for i in range(nCaRU):
    CLidx.append([])
    origin = tuple((np.array([Nx,Ny]) * np.random.rand(2)).astype(int))
    CLidx[i].append(origin)
    latt[origin] = 1
    contour = Vnn[origin]
    nCL = 1

    while contour and nCL<nCLmax+1:
        # choose randomly an element:
        new_cl = random.choice(contour)
        # add to list and matrix
        latt[new_cl] = 1
        CLidx[i].append(new_cl)
        # add NNs to contour
        contour = contour + Vnn[new_cl]
        # remove exisiting clusters
        contour.remove(new_cl)
        # remove duplicated NNs
        contour = list(set(contour))

        nCL += 1

for i in Vnn:
    latt[i] += 1

plt.imshow(latt, interpolation=None, origin='lower')
plt.colorbar()
plt.show()
