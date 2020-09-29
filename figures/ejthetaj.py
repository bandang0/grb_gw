'''ejthetaj plot.'''

import numpy as np
from matplotlib import pyplot as plt

from helpers import *

E_j01 = log(5.4e49)
E_j1 = log(5.3e48)
E_j10 = log(4.82e47)
E_j50 = log(5.36e46)
E_j75 = log(1.8e46)

print(E_j01)

theta_max = 0.25 * 10 ** 0.799 * 180 / np.pi
theta_gw = 0.25 * 10 ** 0.291 * 180 / np.pi
f = lambda a, b, c, t: (c - b * log(np.pi * t / (180. * 0.25))/a) + 52
x = np.array([10 ** a for a in np.linspace(0.5, 1.95, 25)])

plt.figure(9)
plt.text(min(x), E_j01, "$E_j(f_\gamma = 0.1\%)$", fontsize=13, withdash=True, dashdirection=1,
     rotation=0, dashlength=15, color='brown')
plt.text(min(x), E_j1, "$E_j(f_\gamma = 1\%)$", fontsize=13, withdash=True, dashdirection=1,
     rotation=0, dashlength=15, color='brown')
plt.text(min(x), E_j10, "$E_j(f_\gamma = 10\%)$", fontsize=13, withdash=True, dashdirection=1,
     rotation=0, dashlength=15, color='brown')
plt.text(min(x), E_j50, "$E_j(f_\gamma = 50\%)$", fontsize=13, withdash=True, dashdirection=1,
     rotation=0, dashlength=15, color='brown')

plt.plot(x, f(0.736, 2., 2.10, x), "-", color='blue', label=r"$\theta_{\rm{obs}} = \theta_{\rm{obs}}^{GW}$")
plt.plot(x, f(1.43, 2., 0.857, x),  "-", color='blue')

plt.plot(x, f(0.736, 2., 0.489, x), "-", color='red', label=r"$\theta_{\rm{obs}} = \theta_{\rm{obs}}^{GW}/2$")
plt.plot(x, f(1.43, 2., 0.71, x),  "-", color='red')

plt.plot(x, f(0.736, 2., -1.12, x), "-", color='green', label=r"$\theta_{\rm{obs}} = \theta_{\rm{obs}}^{GW}/4$")
plt.plot(x, f(1.43, 2., 0.561, x),  "-", color='green')

plt.fill_between(x, [min(f(1.43, 2., 0.561, a),f(0.736, 2., -1.12, a)) for a in x] , 100 * f(1.43, 2., 0.561, x), color = 'green', alpha = '0.6')#facecolor = "none", hatch = "|", edgecolor = 'green')

plt.fill_between(x, [min(f(0.736, 2., 0.489, a),f(1.43, 2., 0.71, a)) for a in x] , 100 * f(1.43, 2., 0.561, x), facecolor = "none", hatch = "/", edgecolor = 'red', linewidth = 0.0)

plt.fill_between(x, [min(f(0.736, 2., 2.10, a),f(1.43, 2., 0.857, a)) for a in x] , 100 * f(1.43, 2., 0.561, x), facecolor = "none", hatch = ".", edgecolor = 'blue', linewidth = 0.0)

plt.xlabel(r"$\theta_j (\rm{deg})$")
plt.ylabel(r"$\log(E_j /\rm{erg})$")
#plt.xscale('log')
#plt.yscale('log')
plt.legend()
plt.xlim(min(x), max(x))
plt.ylim(E_j50 - 0.1, 56)

plt.savefig("plots/ejthetaj.pdf", bbox_inches = 'tight')
