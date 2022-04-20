import numpy as np
import matplotlib.pyplot as plt

# Use latex figures
plt.rc('text', usetex = True)
plt.rc('font', **{'family': 'serif', 'serif': ['Times New Roman']})


def plot_scatter(save, dx, pos_RR_int, pos_RR_sarco, pos_LCC, pos_NCX, pos_strips, mean):
    """"
    Given three dictionaries s.t. the keys are the cluster sizes and the values 
    are points in the lattice, the code plot these distributions 
    """
    
    pos_RR_int = dx * pos_RR_int
    pos_RR_sarco = dx * pos_RR_sarco
    pos_LCC = dx * pos_LCC
    pos_NCX = dx * np.array(pos_NCX)
    pos_strips = dx * np.array(pos_strips)

    fig = plt.figure(figsize=(12,7))

    #plt.scatter(pos_NCX[:,1], pos_NCX[:,0], s=24, marker='o', alpha=0.4, color='blue', label='NCX')
    plt.scatter(pos_RR_sarco[:,1], pos_RR_sarco[:,0], s=12, marker='^', color='orange', label='RyR mem.')
    plt.scatter(pos_RR_int[:,1], pos_RR_int[:,0], s=12, marker='o', color='blue', label='RyR int')
    #plt.scatter(pos_LCC[:,1], pos_LCC[:,0], s=60, marker='x', color='black', label='LCC')
    plt.scatter(pos_strips[:,1], pos_strips[:,0], s=6, marker='o', color='black', alpha=0.3, label='z-lines')

    plt.legend(loc=0)
    plt.xlabel(r'X [$\mu m$]')
    plt.ylabel(r'Y [$\mu m$]')
    plt.axes().set_aspect('equal', adjustable='box')
    plt.title('center to center CaRU distance: ' + str(round(mean,2)) + '$\mu$m')
    fig.set_tight_layout(True)
    
    if save:
        plt.savefig('../plots/cell-geometry.pdf')
    
    plt.show()

def histo(Vcluster):
    Vhist = {}
    count = 0
    mean_size = 0
    for ci in Vcluster:
        size = len(Vcluster[ci])
        mean_size += size
        count += 1
        if size in Vhist:
            Vhist[size] += 1
        else:
            Vhist[size] = 1

    mean_size /= count
    ci = np.array(list(Vhist.keys()))
    dist = np.array(list(Vhist.values()))

    return ci, dist, mean_size


def plot_mcs(save, Vcluster, Vsarco):
    """"
    Given two cluster distrutions the code plots two histogram with the cluster sizes
    Vcluster: keys=label, values=coordinates
    save=1 save the figure
    """
    
    ci, dist, mean = histo(Vcluster)
    ci_sarco, dist_sarco, mean_sarco = histo(Vsarco)

    print('plotting...')

    fig = plt.figure(figsize=(4,4))
    plt.semilogy(ci, dist, 'o')
    plt.xlabel('Number of RyRs')
    plt.title('Internal space - mean cluster size: ' + str(np.round(mean, 3)))
    fig.set_tight_layout(True)

    if save:
        plt.savefig('../plots/cluster-dist.pdf')

    fig = plt.figure(figsize=(4,4))
    plt.semilogy(ci_sarco, dist_sarco, 'o')
    plt.xlabel('Number of RyRs')
    plt.title('Sarcolemma - mean cluster size: ' + str(np.round(mean_sarco, 3)))
    fig.set_tight_layout(True)

    if save:
        plt.savefig('../plots/cluster-dist-sarco.pdf')

    plt.show()
