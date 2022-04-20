import numpy as np

def clustering(TH, dx, ind, nRyRsPixel):
    """
    given a set of points the clustering algorithm returns the same system packed into clusters
    TH: (threshold) maximum distance between two RyRs in the same cluster (microns)
    dx: spatial grid (microns)
    ind: points in the lattice
    """
    
    ind = np.array(ind)
    ind = ind.astype(int)
    nRR = ind.shape[0]
    Vhash = []
    for i in range(nRR):
        Vhash.append(tuple(ind[i,:]))
    
    Vcluster = {}
    Vdict = {}
    VRyRsPixel = {}
    count = 0
    for i in Vhash:
        Vdict[i] = 0
        VRyRsPixel[i] = [0, nRyRsPixel[count]]
        count += 1

    ci = 0
    positions = []   
    # get cluster value
    values = np.fromiter(iter(Vdict.values()), dtype=int)
    Vhash = list(Vdict.keys())
    V = np.asarray(Vhash)
    maxd = TH / dx
    
    # check if there is any node without cluster
    while np.any(values == 0):
        # initial node where we initiate the construction of a cluster
        # init_point = Vhash[np.where(values == 0)[0][0]]
        for j in Vdict:
            cj = Vdict[j]
            if cj==0:
                init_point = j
                break
        # add this initial point to the contour
        contour = [init_point] 
        ci += 1
        # assign the cluster value ci
        Vdict[init_point] = ci
        VRyRsPixel[init_point][0] = ci
        # Vcluster: cluster (i) --> all elements
        Vcluster[ci] = []
        Vcluster[ci].append(init_point)
        ## ALGORITHM ##
        # check if there's any element on the contour
        while contour:
            contour2 = [] # In this list we will save the new contour generated
            for elem_c in contour:
                # calculate the distance to the other nodes
                dist = np.sqrt(np.sum((elem_c - V)**2, axis=1))
                # get those inside a radius maxd
                elem_new = V[(dist<maxd) & (dist>0.5)]
                elem_new = tuple(map(tuple, elem_new))
                # if these elements do not belong to the cluster, then we include them
                for elem_new_i in elem_new:
                    if Vdict[elem_new_i] == 0:
                        Vdict[elem_new_i] = ci
                        VRyRsPixel[elem_new_i][0] = ci
                        Vcluster[ci].append(elem_new_i)
                        contour2.append(elem_new_i)
            contour = contour2
        # get cluster value
        values = np.fromiter(iter(Vdict.values()), dtype=int)

    ## PROPERTIES ##
    cluster_size = np.zeros(ci)
    cluster_center = np.zeros((ci,2))
    for ii in VRyRsPixel:
        cli = VRyRsPixel[ii][0]-1
        cluster_size[cli] += VRyRsPixel[ii][1]
        cluster_center[cli,:] += np.array(ii) * VRyRsPixel[ii][1]
    
    for ii in range(ci):
        cluster_center[ii,:] /= cluster_size[ii]
    
    # Nearest neighbour distance between cluster centers
    min_dist = np.zeros(ci)
    mask = np.ones(ci).astype(bool)
    for i in range(ci):
        cli = cluster_center[i,:]
        mask[i] = False
        dist = np.sqrt(np.sum((cli-cluster_center[mask,:])**2, axis=1))
        min_dist[i] = np.min(dist)
        mask[i] = True

    return cluster_size, min_dist, Vcluster, Vdict, VRyRsPixel
    # return Vcluster, Vdict
