#!/bin/env python3
"""n_vary plot"""

from matplotlib import pyplot as plt

from helpers import *

print("Plotting light curve.")
plt.figure(4)
plt.errorbar(OBS_TIMES, OBS_FLUXES, yerr = OBS_ERRORS, color = "grey", fmt="+")
plt.xscale('log')
plt.yscale('log')
plt.xlim([10, 400.])
plt.ylim([1.e-3, 4.e-1])
plt.xlabel('Observer Time (days)')
plt.ylabel(r'$S_{\nu} \times (D/40Mpc)^{-2}$ (mJy)')
#plt.title('Light Curve (at 3 GHz)')
plt.grid(True, which='major')
lab = { '-2':r'$n = \bar{n}$' + ' / 4',
        '-1':r'$n = \bar{n}$' + ' / 2',
        '0': r'$n = \bar{n}$',
        '1': r'$n = \bar{n} \times 2$',
        '2': r'$n = \bar{n} \times 4$'}
for s in ['-2', '-1', '0', '1', '2']:
    t_obs, F_nu_obs = data_list("data/n{}.out".format(s))
    t_obs_days = t_obs / day
    F_nu_obs_mJy = F_nu_obs * 1.e26
    plt.plot(t_obs_days, F_nu_obs_mJy, label=lab[s])
plt.legend()
plt.savefig("plots/n_vary.pdf", bbox_inches='tight')
