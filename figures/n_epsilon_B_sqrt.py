'''n_epsilon_B_sqrt plot.'''

import numpy as np
from matplotlib import pyplot as plt

from helpers import *

EXPLO_EXTENT = 20
# Parameter ranges:s
N_ = np.linspace(-4, -2, EXPLO_EXTENT) #[1.1, 1.15, 1.2, 1.25, 1.3]]
SQRTNEPSILONB_ = np.linspace(-4, -2, EXPLO_EXTENT)#[x * 1.e-2 for x in [0.3, 0.35, 0.4, 0.45]]#[10. ** x for x in np.linspace(-2., -2., 5)


X, Y = np.meshgrid(np.array([10 ** n for n in N_]), np.array([10 ** sneb for sneb in SQRTNEPSILONB_]))
Z = np.zeros_like(X)
for i, n in enumerate(N_):
    for j, sneb in enumerate(SQRTNEPSILONB_):
        Z[i,j] = log(fit_err(cocoon_dic([n, BEST_E, 2 * sneb - n, BEST_GB_MIN, BEST_GB_MAX])))
        print(Z[i, j])
        print(n, sneb)

plt.figure(9)
plt.xlabel("$n$"  + " " + r"($\rm{cm}^{-3}$)")
plt.ylabel("$\sqrt{n \epsilon_b}$")
plt.xscale("log")
plt.yscale("log")
plt.contourf(X, Y, Z, 400, cmap='jet')
cbar = plt.colorbar()
cbar.set_label(r'$\log \chi^2$')
plt.savefig("plots/n_epsilon_B_sqrt.pdf", bbox_inches = 'tight')
