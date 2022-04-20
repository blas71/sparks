import numpy as np
import random
import bisect
from scipy.optimize import fsolve
from scipy.special import erf

def RyR_gen(Lx, Ly, Ty, sx, sy, mu, sigma, lambd, dx, nRyR, TT_mat):
    """
    Lx: box size in the X axis (microns)
    Ly: box size in the Y axis (microns)
    dx: spatial grid (microns)
    """

    Nx = int(Lx / dx)
    Ny = int(Ly / dx)

    Vnn = Vnn_fun(Nx, Ny)

    posRR = []
    LCC_TT = [] 
    latt = np.zeros((Nx,Ny))
    nRyRsPixel = np.zeros((Nx,Ny))

    nRyRsCluster = []
    while np.sum(nRyRsCluster)<4508*9:
        U = np.random.rand()
        x_soli = fsolve(CDF, 0.5, args=(mu, sigma, lambd, U))[0]
        if x_soli>0:
            nRyRsCluster.append(int(np.floor(x_soli))+1)

    print(np.mean(nRyRsCluster))
    nCaRUy = int(np.floor(Ly / Ty))
    nCaRUx = len(nRyRsCluster) / nCaRUy
    Ti = int(Nx / nCaRUx) + 1
    Tj = int(Ty / dx)

    nCaRUx = int(nCaRUx)

    zi = np.random.normal(0, sx, nCaRUx*nCaRUy)
    zj = np.random.normal(0, sy, nCaRUx*nCaRUy)

    k = 0
    for j in range(nCaRUy):
        for i in range(nCaRUx):
            posi = (1+i)*Ti + np.rint(zi[k]).astype(int)
            posj = (1+j)*Tj + np.rint(zj[k]).astype(int)
            if (posi<Nx-1) & (posj<Ny-1) & (posi>0) & (posj>0) :
                if TT_mat[(posi,posj)]==1:
                    LCC_TT.append((posi,posj))
                v = (posi, posj)
                latt, nRyRsPixel = cluster(Nx, Ny, nRyR, nRyRsCluster[k], Vnn, v, latt, nRyRsPixel)
            k += 1

    nRyRpixel_vec = []
    for i in range(Nx):
        for j in range(Ny):
            if latt[i,j] == 1:
                posRR.append((i,j))
                nRyRpixel_vec.append(nRyRsPixel[i,j])

    return posRR, nRyRsPixel, nRyRpixel_vec, LCC_TT

def CDF(x, mu, sigma, lambd, U):
    return 0.25*(1 + erf((x-mu)/(np.sqrt(2)*sigma))) + 0.5*(1 - np.exp(-lambd * x)) - U

def cluster(Nx, Ny, nRyR, nRyRsCluster, Vnn, origin, latt, nRyRsPixel):
    nCL = 0
    CLidx = [origin]
    contour = Vnn[origin]
    latt[origin] = 1
    nCL = 1 # current number of pixels

    res = np.mod(nRyRsCluster,nRyR)
    
    if res==0:
        nCLmax = int(nRyRsCluster / nRyR) # maximum number of pixels
        nRyRsPixel[origin] = nRyR
    else:
        nCLmax = int(np.floor(nRyRsCluster / nRyR)) + 1 # number of pixels
        nRyRsPixel[origin] = res

    while contour and nCL<nCLmax:
        check = 0
        for ci in contour:
            check += latt[ci]
        if check == len(contour):
            break
        # choose randomly an element of the contour
        new_cl = random.choice(contour)
        if latt[new_cl] == 0:
            # add to list of pixels of the cluster and to the matrix
            CLidx.append(new_cl)
            latt[new_cl] = 1
            nRyRsPixel[new_cl] += nRyR
            # add NNs to contour
            contour = contour + Vnn[new_cl]
            # remove the current point
            contour.remove(new_cl)
            # remove duplicated NNs
            contour = list(set(contour))

            nCL += 1

    return latt, nRyRsPixel

def Vnn_fun(Nx, Ny):

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

    return Vnn

def AT_gen(L_, std_, Lx, Ly, Tx, dx, NAT_, NCX_TT, posRR, lattice, nRyRpixel_vec):
    L = int(L_ / dx)
    std = int(std_ / dx)
    Nx = int(Lx / dx)
    Ny = int(Ly / dx)
    NAT = int(NAT_ / dx)
    Ti = int(Tx / dx)
    
    Lxv = np.arange(Nx)
    Lyv = np.arange(Ny)

    RyR = np.array([(0,0), (0,1), (1,0), (1,1)])
    nsubCaRU = RyR.shape[0]

    AT_latt = np.zeros((Nx,Ny))
    NATi = 0
    posAT = []

    while NATi<NAT:
        zi = int(np.random.normal(L, std))
        if zi>0:
            px = random.choice(Lxv)
            py = random.choice(Lyv)
            if py+zi<Ny-3 and py+zi>2 and px>2 and px<Nx-3:
                if np.sum(AT_latt[px,py:py+zi])==0:
                    AT_latt[px,py:py+zi] = 1
                    NATi += zi
                    for l in range(zi):
                        for k in [-1,0,1]:
                            if px+k>0 and px+k<Nx-1:
                                NCX_TT.append((px+k,py+l))
                        if np.mod(l,Ti)==0:
                            posAT.append((px, py+l))
                            v = np.array([px, py+l])
                            for ll in range(nsubCaRU):
                                posRRi = tuple(v + RyR[ll,:])
                                lattice[posRRi] = 1
                                posRR.append(posRRi)
                                nRyRpixel_vec.append(9)
    
    return posAT, NCX_TT, np.array(posRR), lattice, nRyRpixel_vec
