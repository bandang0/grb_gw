#!/bin/env python3
"""remnants plot."""

from matplotlib import pyplot as plt

from helpers import *

print("Plotting light curve.")
plt.figure(4)
#plt.errorbar(OBS_TIMES, OBS_FLUXES, yerr = OBS_ERRORS, color = "grey", fmt="+")
plt.xscale('log')
plt.yscale('log')
plt.xlim([7.3e-5, 2.9e2])
plt.ylim([6.e-5, 1.9e10])
plt.xlabel('Observer Time (days)')
plt.ylabel(r'$S_{0.3-10\rm{ keV}} \times (D/40Mpc)^{-2}$ ($\mu$Jy)')
#plt.title('Light Curve (at 3 GHz)')
plt.grid(True, which='major')
plt.savefig("plots/remants.pdf", bbox_inches='tight')
