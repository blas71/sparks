import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.special import erf

mu = 75
sigma = 20
lambd = 1/10

def CDF(x,U):
    return 0.5 * (0.5 * (1 + erf((x-mu)/(np.sqrt(2)*sigma))) + 1 - np.exp(-lambd * x)) - U

# x = np.linspace(0, 2*mu, 200)
# dx = x[1]-x[0]
# 
# f = np.gradient(CDF(x,0.1), dx)
# print(np.sum(x*f*dx))
# plt.semilogy(x, f)

Uvec = []
x_sol = []
nRyRs = 9

while np.sum(x_sol)<84375:
    U = np.random.rand()
    Uvec.append(U)
    x_soli = fsolve(CDF, 0.5, args=U)[0]
    if x_soli>0:
        x_sol.append(x_soli)

x_sol = np.array(x_sol).astype(int)+1
x = np.linspace(np.min(x_sol), np.max(x_sol), 200)

nCL = (x_sol / nRyRs).astype(int) + 1

print('RyRs: ', np.sum(x_sol))
print('RyRs/cluster: ', int(np.mean(x_sol)))
print('clusters: ', len(nCL))
print('CaRU(x): ', len(nCL) / 100 * 1.6)

plt.figure()
plt.subplot(2,1,1)
plt.plot(x, CDF(x,Uvec[0])+Uvec[0])

plt.subplot(2,1,2)
plt.hist(x_sol, bins=35)
plt.show()
