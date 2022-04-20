import numpy as np
import ast

def rotation(Lx, Ly, dx, gap_m, width_m, Vcluster):
    """"
    Given a set of clusterized points, the code applies a random rotation with PBC and then maps the result to lattice.
    Lx: box size in the X axis (microns)
    Ly: box size in the Y axis (microns)
    dx: spatial grid (microns)
    gap_m: gap between sarcolemma and internal space (microns)
    width_m: membrane width (microns)
    Vcluster: labelling of each cluster
    """

    Mx = round(Lx / dx)
    My = round(Ly / dx)
    
    width = int(width_m / dx)
    gap = int(gap_m / dx)

    lattice = np.zeros((Mx,Mx))

    for ii in Vcluster:
        Vcluster[ii] = np.array(Vcluster[ii])
    positions = []

    for jj in Vcluster:
        # get positions of the j-th cluster
        pos_j = Vcluster[jj]
        # center of mass
        centerj = np.mean(pos_j, axis=0)
        # relative positions to CM
        pos_j = pos_j - centerj
        # pick a random angle in (0, 2pi)
        phi = 2 * np.pi * np.random.rand()
        # apply the matrix rotation
        for ii in range(pos_j.shape[0]):
            pos_j[ii,:] = np.array([pos_j[ii,0]*np.cos(phi) - pos_j[ii,1]*np.sin(phi), \
                    pos_j[ii,0]*np.sin(phi) + pos_j[ii,1]*np.cos(phi)])
        # save the new rotated cluster in the dictionary
        Vcluster[jj] = pos_j + centerj
        # map to the discrete lattice
        Vcluster[jj] = Vcluster[jj].astype(int)
        # # PBC 
        Vcluster[jj] -= Mx * (Vcluster[jj] / Mx).astype(int)

    # transform to a matrix form
    for jj in Vcluster:
        for ii in range(Vcluster[jj].shape[0]):
            lattice[tuple(Vcluster[jj][ii])] += 1

    # save duplicated points only one time
    for i in range(Mx):
        for j in range(Mx):
            if lattice[i,j] > 0:
                positions.append((i+width+gap,j+width+gap))

    return positions





