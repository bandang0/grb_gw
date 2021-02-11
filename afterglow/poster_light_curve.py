#!/bin/env python3
"""Light curve plot for FS poster."""

from matplotlib import pyplot as plt

from helpers import *

print("Plotting light_curve plot.")
plt.figure(4)
#plt.errorbar(OBS_TIMES, OBS_FLUXES, yerr = OBS_ERRORS, color = "grey", fmt="+")
plt.xscale('log')
plt.yscale('log')
plt.xlim([10, 400.])
#plt.ylim([2.e-7, 3])
plt.xlabel('Observer Time (days)')
plt.ylabel(r'$S_{\nu} \times (D/40Mpc)^{-2} (\mu Jy)$')
#plt.title('Light Curve (at 3 GHz)')
plt.grid(True, which='major')
plt.errorbar(OBS_TIMES_d, OBS_FLUXES_muJy, yerr = OBS_ERRORS_muJy, fmt = "+", color = "grey")

# configuration box
conf_string = ("Exterior medium:\n$n$ = {:.0e} $cm^{{-3}}$\n$\epsilon_B$ = {:.0e}\n\nJet\n$E_0$ = {:.0e} erg\n"+r"$\theta_v$ = {} deg"+"\n"+r"$\theta_j$ = {} deg").format(1.e-3, 5.e-3, 2.e52, 20, 5)
conf_string = r"Exterior medium:"+"\n"+r"$n = 10^{-3} cm^{{-3}}$"+"\n"+r"$\epsilon_B = 5 \times 10^{-4}$"

plt.text(15, 400, conf_string, fontsize = 12, weight='bold')

# light curves
add_plot(plt, "best_cocoon", linestyle= "--", color = "black", label = "Cocoon (increasing phase)")
add_plot(plt, "best_jet_decrease", color = "blue", label = "Best-fit jet (decreasing phase)")
add_plot(plt, "best_jet_decrease_n", color = "orange", label = r"$n = \bar{n} \times 10$")
add_plot(plt, "best_jet_decrease_e0", color = "red", label = r"$E_0 = \bar{E_0} \times 10$")
add_plot(plt, "best_jet_decrease_theta_v", color = "green", label = r"$\theta_v = \bar{\theta_v} / 2$")

plt.ylim([10, 1000])
plt.xlim([10, 400])
plt.legend()
plt.title("3 GHz light curve")
plt.savefig("plots/poster_light_curve.pdf", bbox_inches='tight')
