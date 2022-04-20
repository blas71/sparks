import numpy as np

def clustering(TH, dx, ind):
    
    nRR = ind.shape[1]
    Vhash = []
    for i in range(nRR):
        Vhash.append(tuple(ind[:,i]))
    
    Vcluster = {}
    Vdict = {}
    for i in Vhash:
        Vdict[i] = 0
     
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
        init_point = Vhash[np.where(values == 0)[0][0]]
        # add this initial point to the contour
        contour = [init_point] 
        ci += 1
        # assign the cluster value ci
        Vdict[init_point] = ci
        # cluster --> all elements
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
                elem_new = V[dist<maxd]
                elem_new = tuple(map(tuple, elem_new))
                # if these elements do not belong to the cluster, then we include them
                for elem_new_i in elem_new:
                    if Vdict[elem_new_i] == 0:
                        Vdict[elem_new_i] = ci
                        Vcluster[ci].append(elem_new_i)
                        contour2.append(elem_new_i)
            contour = contour2
        # get cluster value
        values = np.fromiter(iter(Vdict.values()), dtype=int)
    
    return Vcluster, Vdict
