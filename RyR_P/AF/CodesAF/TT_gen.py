import numpy as np
import random 

def tubules_gen(Lx, Ly, dx, pos_zplanes, pos_NCX, lambd, densityTT):
    """
    pos_zplanes: list of tuples with the z-planes nodes
    lambd: parameter in the exponential distribution
    densityTT: density of T-tubules
    """

    Nx = int(Lx/dx)
    Ny = int(Ly/dx)
    
    # total number of TT
    nTT = densityTT * len(pos_zplanes)

    # create matrix with available positions
    lattice = np.zeros((Nx, Ny))
    
    # 1 for the z-planes
    for i in pos_zplanes:
        lattice[i] = 1
    
    nTTi = 0
    ttubules = []
    while nTTi<nTT:
        # choose vertical length L/2
        L2 = np.random.exponential(lambd)

        L2 = int(L2/dx)
        L = L2*2 + 1

        # grow t-tubule up & down L/2
        tti = np.zeros((2,L)).astype(int)
        for i in range(L):
            tti[0,i] = i
        tti[0,:] -= L2

        # choose origin
        pos_o = random.choice(pos_zplanes)
        for i in range(L):
            tti[:,i] += np.array(pos_o)
            # PBC 
            tti[0,i] -= Nx * int(tti[0,i] / Nx)

        #Vcluster[jj] -= Mx * (Vcluster[jj] / Mx).astype(int)

        # check if this space is empty
        empty = L
        for i in range(L):
            empty -= lattice[tuple(tti[:,i])]

        if empty==0:
            for i in range(L):
                # add new T-tubules
                ttubules.append(tuple(tti[:,i]))
                pos_NCX.append(tuple(tti[:,i]))
                # change the available matrix
                lattice[tuple(tti[:,i])] = 0
                lattice[(tti[0,i], tti[1,i]+1)] = 0
                lattice[(tti[0,i], tti[1,i]-1)] = 0
                lattice[(tti[0,i], tti[1,i]+2)] = 0
                lattice[(tti[0,i], tti[1,i]-2)] = 0
                
            nTTi += L

    return np.array(ttubules), pos_NCX


