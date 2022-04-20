import numpy as np
import random 

def SR_strips(Lx, Ly, Tx, Ty, dx, lambd):

    Nx = int(Lx / dx)
    Ny = int(Ly / dx)

    Ti = int(Tx / dx)
    Tj = int(Ty / dx)
    
    pos_strips = []

    TT_mat = np.zeros((Nx,Ny))
    NCX_TT = []

    i = 1
    while Tj*i<Ny:
        Tji = Tj*i

        Ldown = np.random.exponential(lambd)
        Lup = np.random.exponential(lambd)
        while Lup + Ldown > Lx:
            Ldown = np.random.exponential(lambd)
            Lup = np.random.exponential(lambd)

        Ldown = int(Ldown/dx)
        Lup = int(Lup/dx)

        for k in [-1,0,1]:
            if Tji+k>-1 and Tji+k<Ny:
                for j in range(Nx-2):
                    pos_strips.append((j+1,Tji+k))
                    if j+1<Ldown or j+1>Nx-Lup:
                        if j+1>1 and j+1<Nx-2:
                            NCX_TT.append((j+1,Tji+k))
                            TT_mat[(j+1,Tji+k)] += 1
        i += 1

    latt = np.zeros((Nx,Ny))
    for i in pos_strips:
        latt[i] = 1

    pos_TnC = []
    for i in range(Nx):
        for j in range(Ny):
            if latt[i,j]==0:
                pos_TnC.append((i,j))

    return pos_strips, pos_TnC, NCX_TT, TT_mat

