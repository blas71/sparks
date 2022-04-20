import random
import numpy as np

def LCC_fun(Lx, Ly, Tx, Ty, dx):
    """
    Lx, Ly: spatial dimensions
    """
    
    # generation loop 
    nLCC = 0
    count = 0
    VLCC = {}
    pos_LCC = []
    while nLCC < nLCCtotal:
        # choose randomly a cluster
        cli = random.choice(accepted)
        # choose randomly between 1 and 5 points
        nLCCi = int((1 + 5*np.random.random()))
        # choose the corresponding points
        VLCC[count] = random.sample(Vsarco[cli], nLCCi)
        # remove current cluster
        accepted.remove(cli)
        # update number of LCC
        nLCC += nLCCi
        count += 1

    for ci in VLCC:
        for lcci in VLCC[ci]:
            pos_LCC.append(lcci)

    return VLCC, np.array(pos_LCC)
