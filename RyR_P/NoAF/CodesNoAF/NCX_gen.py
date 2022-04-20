import numpy as np

def NCX_fun(Lx, Ly, dx):

    Nx = int(Lx / dx)
    Ny = int(Ly / dx)

    pos_NCX = []
    
    latt = np.ones((Nx,Ny))
    latt_int = np.ones((Nx-2,Ny-2))

    for i in range(Nx-2):
        for j in range(Ny-2):
            latt[i+1,j+1] = latt[i+1,j+1] - latt_int[i,j]

    for i in range(Nx):
        for j in range(Ny):
            if latt[i,j]==1:
                pos_NCX.append((i,j))
    
    return pos_NCX
