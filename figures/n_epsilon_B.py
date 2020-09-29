'''n_epsilon_B plot.'''
from sys import exit
import numpy as np
from matplotlib import pyplot as plt

from helpers import *

EXPLO_EXTENT = 20
# Parameter ranges:
N_ = np.linspace(BEST_N - 1., BEST_N + 1., EXPLO_EXTENT)
EPSILONB_ = np.linspace(BEST_EB - 1., BEST_EB + 1., EXPLO_EXTENT)


X, Y = np.meshgrid(np.array([10 ** n for n in N_]),
                   np.array([10 ** eb for eb in EPSILONB_]))
Z = np.zeros_like(X)
for i, n in enumerate(N_):
    for j, eb in enumerate(EPSILONB_):
        Z[i,j] = log(fit_err(cocoon_dic([n, BEST_E, eb, BEST_GB_MIN,
                                         BEST_GB_MAX])))
        print(Z[i, j])
        print(n, eb)

plt.figure(9)
plt.xlabel("$n$"  + " " + r"($\rm{cm}^{-3}$)")
plt.ylabel("$\epsilon_b$")
plt.xscale("log")
plt.yscale("log")
plt.contourf(X, Y, Z, 400, cmap='jet')
cbar = plt.colorbar()
cbar.set_label(r'$\log \chi^2$')
plt.show()
exit(0)
plt.savefig("plots/n_epsilon_B.pdf", bbox_inches = 'tight')
