'''umnbis.'''

import numpy as np
from matplotlib import pyplot as plt

from helpers import *

EXPLO_EXTENT = 15
# Parameter ranges:
N_ = np.linspace(-4, -2, EXPLO_EXTENT) #[1.1, 1.15, 1.2, 1.25, 1.3]]
NGBMINA_ = np.linspace(-6, -2, EXPLO_EXTENT)#[1. + x for x in [0.25, 0.3, 0.35]]#[10. ** x for x in np.linspace(-2., 0.3, 9)]

X, Y = np.meshgrid(np.array([10 ** n for n in N_]), np.array([10 ** numa for numa in NGBMINA_]))
Z = np.zeros_like(X)
for i, n in enumerate(N_):
    for j, numa in enumerate(NGBMINA_):
        if 0.2 * n - 0.2 * numa < 0.48:
            Z[i,j] = log(fit_err(cocoon_dic([n, BEST_E, BEST_EB, 0.2 * n - 0.2 * numa, BEST_GB_MAX])))
        else:
            Z[i,j] = 1.
        print(n, 0.2 * n - 0.2 * numa, Z[i, j])

plt.figure(9)
plt.xlabel("$n$"  + " " + r"($\rm{cm}^{-3}$)")
plt.ylabel(r"$n u_m^{-\alpha}$")
plt.xscale("log")
plt.yscale("log")
plt.contourf(X, Y, Z, 100, cmap='jet')
#plt.scatter(x, y, c=e, cmap='jet',vmin=e.min(), vmax=e.max(), label = r"$\log \chi^2$")
#plt.legend()
c = plt.colorbar()
c.set_label(r'$\log \chi^2$')
plt.savefig("plots/umnbis.pdf", bbox_inches = 'tight')
