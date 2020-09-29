'''umn plot.'''

import numpy as np
from matplotlib import pyplot as plt

from helpers import *


EXPLO_EXTENT = 15
# Parameter ranges:
N_ = np.linspace(BEST_N - 2., BEST_N + 2., EXPLO_EXTENT) #[1.1, 1.15, 1.2, 1.25, 1.3]]
GBMIN_ = np.linspace(10 ** BEST_GB_MIN - 1.5, 10 ** BEST_GB_MIN + 1.2, EXPLO_EXTENT)#[1. + x for x in [0.25, 0.3, 0.35]]#[10. ** x for x in np.linspace(-2., 0.3, 9)]

X, Y = np.meshgrid(np.array([10 ** n for n in N_]), np.array([um for um in GBMIN_]))
Z = np.zeros_like(X)
for i, n in enumerate(N_):
    for j, um in enumerate(GBMIN_):
        Z[i,j] = log(fit_err(cocoon_dic([n, BEST_E, BEST_EB, log(um), BEST_GB_MAX])))
        print(Z[i, j])
        print(n, um)

plt.figure(9)
plt.xlabel("$n$"  + " " + r"($\rm{cm}^{-3}$)")
plt.ylabel(r"$u_m$")
plt.xscale("log")
plt.contourf(X, Y, Z, 100, cmap='jet')
#plt.scatter(x, y, c=e, cmap='jet',vmin=e.min(), vmax=e.max(), label = r"$\log \chi^2$")
#plt.legend()
c = plt.colorbar()
c.set_label(r'$\log \chi^2$')
plt.savefig("plots/umn.pdf", bbox_inches = 'tight')
