#!/bin/env python3
"""prob plot."""

import numpy as np
from matplotlib import pyplot as plt
from numpy import sqrt

f = lambda r: r/(sqrt(r ** 2 - 2.12 ** 2) * sqrt(9**2 - 2.12**2)) if (r > 2.12 and r < 9.) else 0.


x = np.linspace(0, 10., 123338)
print("Plotting distribution.")
plt.figure(4)
plt.xlabel('$R$ (kpc)')
plt.ylim(-0.1, 1.1)
plt.ylabel(r'$f_R(R$ | $R_{\rm{proj}} = 2.12$ kpc$)$ (1/kpc)')
plt.plot(x, np.array([f(a) for a in x]) )
plt.text(2.12, 0., r"$R_{\rm{proj}}$ = 2.12 kpc", fontsize=13, withdash=True, dashdirection=1,
     rotation=45, dashlength=55, color='brown')
plt.text(9, 0., r"$R_G = 9$ kpc", fontsize=13, withdash=True, dashdirection=0,
     rotation=-35, dashlength=35, color = 'brown')
plt.grid(True, which='major')

plt.savefig("plots/prob.pdf", bbox_inches='tight')
