#!/bin/env python3
"""Just plot the data from files listed on CL on the same plot."""

from sys import argv
from matplotlib import pyplot as plt

from helpers import *

print("Plotting light curve.")
plt.figure(4)
plt.errorbar(OBS_TIMES_d, OBS_FLUXES_mJy, yerr = OBS_ERRORS_mJy, color = "grey", fmt="+")
plt.xscale('log')
plt.yscale('log')
#plt.xlim([10, 1000.])
#plt.ylim([1.e-3, 4.e-1])
plt.xlabel('Observer Time (days)')
plt.ylabel(r'$S_{\nu} \times (D/40Mpc)^{-2}$ (mJy)')
#plt.title('Light Curve (at 3 GHz)')
plt.grid(True, which='major')
for s in range(int((len(argv) - 1) / 2)):
    t_obs, F_nu_obs = data_list("data/{}.out".format(argv[2 * s + 1]))
    t_obs_days = t_obs / day
    F_nu_obs_mJy = F_nu_obs * 1.e26
    plt.plot(t_obs_days, F_nu_obs_mJy, label=argv[2*s + 2])
plt.legend()
plt.savefig("plots/over_plot.pdf", bbox_inches='tight')
