#!/bin/env python3
"""jet_theta plot."""

from matplotlib import pyplot as plt

from helpers import *

print("Plotting jet_theta.")
plt.figure(4)
plt.errorbar(OBS_TIMES, OBS_FLUXES, yerr = OBS_ERRORS, color = "grey", fmt="+")
plt.xscale('log')
plt.yscale('log')
plt.xlim([0.005, 400.])
plt.ylim([2.e-3, 8])
plt.xlabel('Observer Time (days)')
plt.ylabel(r'$S_{\nu} \times (D/40Mpc)^{-2}$ (mJy)')
#plt.title('Light Curve (at 3 GHz)')
plt.grid(True, which='major')
lab = { '20':r'$\theta_{\rm{obs}} = 20$ deg (hidden)',
        '15':r'$\theta_{\rm{obs}} = 15$ deg',
        '10': r'$\theta_{\rm{obs}} = 10$ deg',
        '5': r'$\theta_{\rm{obs}} = 5$ deg',
        '0': r'$\theta_{\rm{obs}} = 0$ deg'}
for s in ['20', '15', '10', '5', '0']:
    t_obs, F_nu_obs = data_list("data/j{}.out".format(s))
    t_obs_days = t_obs / day
    F_nu_obs_mJy = F_nu_obs * 1.e26
    plt.plot(t_obs_days, F_nu_obs_mJy, "-", label=lab[s])
plt.legend()
plt.savefig("plots/jet_theta.pdf", bbox_inches='tight')
