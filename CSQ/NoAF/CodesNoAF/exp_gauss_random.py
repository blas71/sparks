import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf
from scipy.optimize import fsolve

mu = 100
sigma = 20
lambd = 1/10

def CDF(x, U):
    return 0.5 * (0.5 * (1 + erf((x-mu)/(np.sqrt(2)*sigma))) + 1 - np.exp(-lambd * x)) - U

nRyR = []
for i in range(1800):
    U = np.random.rand()
    x_soli = fsolve(CDF, 0.5, args=(U))[0]
    nRyR.append(int(x_soli)+1)

plt.hist(nRyR, bins=50)

plt.show()

# beta = 5 # mean value = beta
# mu = 50
# sigma = 20
# r = []
# rexp = []
# rnorm = []
# 
# N = 80000
# 
# for i in range(N):
#     rr = np.random.rand()
#     if rr<0.5:
#         r.append(np.random.exponential(beta))
#     else:
#         r.append(np.random.normal(mu, sigma))
# 
# for i in range(N):
#     rexp.append(np.random.exponential(beta))
#     rnorm.append(np.random.normal(mu, sigma))
# 
# rhist = np.histogram(r, bins=100)
# rhist_exp = np.histogram(rexp, bins=100)
# rhist_norm = np.histogram(rnorm, bins=100)
# 
# plt.plot(rhist[0])
# plt.plot(rhist_exp[0])
# plt.plot(rhist_norm[0])
# plt.show()
