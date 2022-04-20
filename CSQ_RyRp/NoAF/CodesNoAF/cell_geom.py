import RyR_generation_exp 
import RyR_cluster
import RyR_rotation
import RyR_plot
import RyR_sarco
import LCC_gen
import NCX_gen
import strips_gen
import TT_gen
import imp
import numpy as np
import random
import matplotlib.pyplot as plt

RyR_generation_exp = imp.reload(RyR_generation_exp)
RyR_gen = RyR_generation_exp.RyR_gen
AT_gen = RyR_generation_exp.AT_gen

RyR_plot = imp.reload(RyR_plot)
plot_scatter = RyR_plot.plot_scatter

RyR_sarco = imp.reload(RyR_sarco)
sarco_grid = RyR_sarco.sarco_grid

NCX_gen = imp.reload(NCX_gen)
NCX_fun = NCX_gen.NCX_fun

strips_gen = imp.reload(strips_gen)
SR_strips = strips_gen.SR_strips

RyR_cluster = imp.reload(RyR_cluster)
clustering = RyR_cluster.clustering

save = 1

#### PARAMETERS ####
# RyR distribution parameters
Ty = 1.6 # microns
Tx = 0.5 # microns
sx = 0.4 # dispersion - microns
sy = 0.4 # dispersion - microns
mu = 60
sigma = 25
lambd = 1

# internal region
Lx = 12 # microns
Ly = 80 # microns
dx = 0.1 # microns
nRyR = 9

# T-tubules
lambd2 = 3 # mean value 1/lambd

# A-tubules
BranchLength = 1.7 # microns
BranchSTD = 1 # microns
ATdensity = 0.06 # microns / microns**2
NAT = int(ATdensity * Lx* Ly)

#### SARCOLEMMA REGION ####
# External RyRs and LCCs
pos_RR_sarco, pos_LCC, lattice = sarco_grid(Lx, Ly, Tx, Ty/2, dx)

# Exchanger
pos_NCX = NCX_fun(Lx, Ly, dx)
 
#### INTERNAL REGIONAL ####
# generate SR strips and T-tubules
pos_zplanes, pos_TnC, NCX_TT, TT_mat = SR_strips(Lx, Ly, Tx, Ty, dx, lambd2)

# generate a random distribution of RyR 
pos_RR_int, nRyRsPixel, nRyRpixel_vec, LCC_TT = RyR_gen(Lx, Ly, Ty, sx, sy, mu, sigma, lambd, dx, nRyR, TT_mat)

# generate Axial Tubules
posAT, NCX_TT, pos_RR_int, lattice, nRyRpixel_vec = AT_gen(BranchLength, BranchSTD, Lx, Ly, Tx, dx, NAT, NCX_TT, pos_RR_int, lattice, nRyRpixel_vec)

#print(pos_RR_int)
pos_strips = pos_zplanes
# add cell membrane
for i in pos_NCX:
    pos_strips.append(i)

# add TT to NCX
for i in NCX_TT:
    pos_NCX.append(i)

pos_RR = np.vstack((pos_RR_int, pos_RR_sarco))
pos_RR = list(map(tuple, pos_RR))
pos_background = []
for i in pos_NCX:
    if i not in pos_RR:
        pos_background.append(i)

#pos_background = np.array(pos_background)
#pos_RR = np.array(pos_RR)
#plt.scatter(pos_background[:,1], pos_background[:,0], s=4)
#plt.scatter(pos_RR[:,1], pos_RR[:,0], s=0.5)
#plt.axes().set_aspect('equal', 'datalim')
#plt.show()
#pos_RR = np.concatenate((pos_RR_int, pos_RR_sarco), axis=0)
#cl_i, min_dist_i, Vcluster, Vdict = clustering(0.15, dx, pos_RR_int) 

if save:
    # save coordinates
    np.savetxt('../data/posRyR.data', pos_RR_int, fmt='%i')
    np.savetxt('../data/posRyR-sarco.data', pos_RR_sarco, fmt='%i')
    np.savetxt('../data/nRyRpixel.data', nRyRpixel_vec, fmt='%i')
    np.savetxt('../data/posLCC.data', pos_LCC, fmt='%i')
    np.savetxt('../data/posTT-LCC.data', LCC_TT, fmt='%i')
    np.savetxt('../data/posAT-LCC.data', posAT, fmt='%i')
    np.savetxt('../data/posNCX.data', pos_NCX, fmt='%i')
    np.savetxt('../data/posSR.data', pos_strips, fmt='%i')
    np.savetxt('../data/posTnC.data', pos_TnC, fmt='%i')
    np.savetxt('../data/posBackground.data', pos_background, fmt='%i')
    
    # save parameters
    f = open('../data/parameters.data', 'w')
    f.write(str(Lx) + '\n')
    f.write(str(Ly) + '\n')
    f.write(str(dx) + '\n')
    f.write(str(len(pos_RR_int)) + '\n')
    f.write(str(len(pos_RR_sarco)) + '\n') 
    f.write(str(len(pos_LCC)) + '\n') 
    f.write(str(len(LCC_TT)) + '\n') 
    f.write(str(len(posAT)) + '\n') 
    f.write(str(len(pos_NCX)) + '\n') 
    f.write(str(len(pos_strips)) + '\n') 
    f.write(str(len(pos_TnC)) + '\n') 
    f.write(str(len(pos_background)) + '\n') 
    f.write(' ' + '\n')
    f.write('-------------' + '\n')
    f.write(' ' + '\n')
    f.write('Lx (microns)' + '\n')
    f.write('Ly (microns)' + '\n')
    f.write('dx (microns)' + '\n')
    f.write('Number of internal RyR' + '\n')
    f.write('Number of RyR in the sarcolemma' + '\n')
    f.write('Number of LCC' + '\n')
    f.write('Number of LCC T-tubules' + '\n')
    f.write('Number of NCX' + '\n')
    f.write('Number of SR grid' + '\n')
    f.write('Number of TnC grid' + '\n')
    f.write('Number of background grid' + '\n')
    f.close()

# plot
#plot_mcs(0, Vcluster, Vsarco)
#plot_scatter(0, dx, pos_RR_int, pos_RR_sarco, pos_LCC, pos_NCX, pos_strips, dx*np.mean(min_dist_i))
