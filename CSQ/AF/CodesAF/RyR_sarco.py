import numpy as np

def sarco_grid(Lx, Ly, Tx, Ty, dx):
    """"
    Lx, Ly: internal spatial lengths (microns)
    dx: spatial grid (microns)
    """

    RyR = np.array([(0,0), (0,1), (1,0), (1,1)])
    nsubCaRU = RyR.shape[0]

    Nx = int(Lx / dx)
    Ny = int(Ly / dx)
    
    Ti = int(Tx / dx)
    Tj = int(Ty / dx)

    nCaRUx = np.ceil(Lx/Tx).astype(int) # number of CaRUs on the cell membrane on the short axis
    nCaRUy = np.ceil(Ly/Ty).astype(int) - 1 # number of CaRUs on the cell membrane on the long axis

    lattice = np.zeros((Nx, Ny))
    posRR = []
    posLCC = []

    # CaRUs - long axis 
    for i in [0,Nx-2]:
        for j in range(nCaRUy):
            v = np.array([i, (j+1)*Tj])
            if i==0:
                posLCC.append((i, (j+1)*Tj))
            elif i==Nx-2:
                posLCC.append((i+1, (j+1)*Tj))
            for k in range(nsubCaRU):
                posRRi = tuple(v + RyR[k,:])
                lattice[posRRi] = 1
                posRR.append(posRRi) 

    # CaRUs - short axis
    for j in [0,Ny-2]:
        for i in range(nCaRUx):
            v = np.array([i*Ti, j])
            if j==0:
                posLCC.append((i*Ti, j))
            elif j==Ny-2:
                posLCC.append((i*Ti, j+1))
            for k in range(nsubCaRU):
                posRRi = tuple(v + RyR[k,:])
                lattice[posRRi] = 1
                posRR.append(posRRi) 

    # upper corners RyR
    v = np.array([Nx-2, 0])
    for k in range(nsubCaRU):
        posRRi = tuple(v + RyR[k,:])
        lattice[posRRi] = 1
        posRR.append(posRRi) 
    v = np.array([Nx-2, Ny-2])
    for k in range(nsubCaRU):
        posRRi = tuple(v + RyR[k,:])
        lattice[posRRi] = 1
        posRR.append(posRRi) 

    # upper corners LCC
    posLCC.append((Nx-1,0))
    posLCC.append((Nx-1,Ny-1))

    return np.array(posRR), np.array(posLCC), lattice


