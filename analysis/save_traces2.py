import numpy as np
import matplotlib.pyplot as plt
import os

def fun(Nx,Ny,px,py,Roip):
    ROI = []
    for i in range(2*Roip+1):
        for j in range(2*Roip+1):
            key = np.array([i+px-Roip,j+py-Roip])
            if np.all(key<np.array([Nx,Ny])) and np.all(key>np.array([0,0])):
                ROI.append((i+px-Roip,j+py-Roip))
    return ROI

ff = '1'
gg = '.'
Type = 'CTRL'
wd1 = '../' + gg + '/' + Type + '/' + Type + '14'
#path_calc = wd1 + '/Codes' + Type + '/Calc'
#path_cab = wd1 + '/Codes' + Type + '/Cabc'
path_calc = wd1 + '/codes/Calc'
path_cab = wd1 + '/codes/Cabc'

Nsave = 833
Nsave_ryr = 28
factor = int(Nsave/Nsave_ryr)+1
Nx = 800
Ny = 120
dx = 0.1
dt = 0.012
tdes = 8000 #ms
Ni_max = 2801
Nprev = 2

file_ryrs = open(wd1 + '/data/follow_RyR_' + ff + '.data', 'r')
lines_ryrs = file_ryrs.readlines()
Nlines = len(lines_ryrs)

cj, ryrk, t2, t1 = np.loadtxt(wd1 + '/data/timeRyR-open.data').T
posRyR = np.loadtxt(wd1 + '/data/posRyR.data').T
posRyR_sarco = np.loadtxt(wd1 + '/data/posRyR-sarco.data').T
posRyR = np.hstack((posRyR,posRyR_sarco))
nRyR = posRyR.shape[1]

cj = cj.astype(int) - 1
ryrk = ryrk.astype(int)
t2 = t2.astype(int)
t1 = t1.astype(int)
posRyR = posRyR.astype(int)

idx_t1 = np.argsort(t1)
t1 = t1[idx_t1]
t2 = t2[idx_t1]
cj = cj[idx_t1]

idx2 = t1>int(np.round(tdes/dt))
t1 = t1[idx2]
t2 = t2[idx2]
cj = cj[idx2]

t1 -= int(np.round(tdes/dt))
t2 -= int(np.round(tdes/dt))

Ca_all = []
Cab_all = []
OPENS = []
CLOSE = []
PROPERTIES = [[], [], [], [], [], [], [], []]

for l in range(Nlines):
    line = lines_ryrs[l]
    line = line.split(' \n')[0]
    line = line.split(', ')
    line = np.array([int(integer) for integer in line])

    t1_min = []
    t2_max = []
    topen = []
    posk = []
    for ryrs_k in line:
        t1_min.append(t1[ryrs_k])
        t2_max.append(t2[ryrs_k])
        posk.append(posRyR[:,cj[ryrs_k]])
        topen.append(t2[ryrs_k] - t1[ryrs_k])

    t1_min = np.min(t1_min)
    t2_max = np.max(t2_max)
    posk = np.array(posk)
    topen = np.array(topen)

    center = np.average(posk, axis=0, weights=topen)
    center = np.round(center).astype(int)
    dist_max = dx*np.sqrt(np.max(np.sum((center-posk)**2,axis=1)))

    Rroi = max(dist_max+0.5, 1.9)
    Rroip = int(Rroi/dx)
    Ni = int(t1_min/Nsave)+1 - Nprev
    if Ni<=0:
        Ni = 1
    Nif = int(t2_max/Nsave) 
    Nif += min(30, Ni_max-Nif)
    Nt = Nif - Ni

    if Nif<Ni_max:
        # center of ROI
        PROPERTIES[6].append(center[0])
        PROPERTIES[7].append(center[1])

        open_trace = np.zeros(Nt*factor)
        close_trace = np.zeros(Nt*factor)
        for ryrs_k in line:
            t1i = int((t1[ryrs_k] - t1_min)/Nsave_ryr) + 1
            t2i = int((t2[ryrs_k] - t1_min)/Nsave_ryr) + 1
            open_trace[t1i] += 1
            close_trace[t2i] -= 1
        open_trace = open_trace.astype(int)
        close_trace = close_trace.astype(int)
        OPENS.append(list(open_trace))
        CLOSE.append(list(close_trace))

        ROI = fun(Nx,Ny,center[1],center[0],Rroip)
        nROI = len(ROI)

        cab_mean = np.zeros((nROI,Nt))
        ca_mean = np.zeros((nROI,Nt))
        # Calcium trace in that ROI
        for i in range(Nt):
            cab = np.loadtxt(path_cab + str(i+Ni).zfill(4) + '.data').T
            ca = np.loadtxt(path_calc + str(i+Ni).zfill(4) + '.data').T
            cab = np.reshape(cab, (Nx,Ny))
            ca = np.reshape(ca, (Nx,Ny))
            jj = 0
            for j in ROI:
                cab_mean[jj,i] = cab[j]
                ca_mean[jj,i] = ca[j]
                jj += 1

        Caib = np.mean(cab_mean, axis=0)
        Cai = np.mean(ca_mean, axis=0)
        cab_max = cab_mean.flatten()
        cab_max.sort()
        cabmax = np.mean(cab_max[-5:])
        FWHM = []
        for i in range(Nt):
            FWHM.append(2 * np.sqrt(np.sum(cab_mean[:,i]>cabmax/2) / np.pi))

        Cab_all.append(list(Caib))
        Ca_all.append(list(Cai))
        PROPERTIES[0].append(dx * np.max(FWHM))

        # relative amplitude
        ind_max = np.argmax(Caib)
        if ind_max==0:
            ind_max = 2
        amp = Caib[ind_max]
        bline = np.quantile(Caib[:ind_max],0.02)
        amp_rel = (amp-bline) / bline
        PROPERTIES[1].append(amp_rel)
        PROPERTIES[5].append(bline)

        # FDHM
        max2 = (amp-bline)/2 + bline
        FDHM = np.sum(Caib>max2)*dt*Nsave
        PROPERTIES[2].append(FDHM)

        # time to peak
        t_to_peak = (ind_max-Nprev)*Nsave*dt
        PROPERTIES[3].append(t_to_peak)

        # Number of opened RyRs
        PROPERTIES[4].append(len(topen))


file_traces = open(wd1 + '/data/spark_traces_' + ff + '.data', 'w')
for trace in Ca_all:
    file_traces.write(" ".join(map(str, trace)))
    file_traces.write("\n")

file_traces.close()

file_traces = open(wd1 + '/data/spark_traces_Rhod_' + ff + '.data', 'w')
for trace in Cab_all:
    file_traces.write(" ".join(map(str, trace)))
    file_traces.write("\n")

file_traces.close()

file_opens = open(wd1 + '/data/open_ryr_traces_' + ff + '.data', 'w')
for trace in OPENS:
    file_opens.write(" ".join(map(str, trace)))
    file_opens.write("\n")

file_opens.close()

file_close = open(wd1 + '/data/close_ryr_traces_' + ff + '.data', 'w')
for trace in CLOSE:
    file_close.write(" ".join(map(str, trace)))
    file_close.write("\n")

file_close.close()

np.savetxt(wd1 + '/data/spark_properties_' + ff + '.data', PROPERTIES)
