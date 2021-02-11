from sys import exit
import matplotlib.pyplot as plt
from matplotlib import colors
import pandas as pd
from numpy import sin, cos, sqrt
from numpy import log10 as log
import numpy as np
from helpers import *
from sys import argv

PLOT_DIR = f"plots/flares"
plt.style.use('presentation')

# Fixed parameters
Eiso = 1.e52
tj = 0.1
tv = 0.12
a = 0.
b = -1.5
Ep = 10 * keV
nup = Ep / hPlanck
nuobs = 5.0 * keV / hPlanck
Gamma = float(argv[1])
beta = sqrt(1 - 1 / Gamma ** 2)

# LOS material parameters
r_los = 5 / Gamma
re_los = 1.e14
te_los = re_los / (beta * cLight)
tobs_los = te_los - re_los / cLight
tm_los = te_los - re_los / cLight
tM_los = te_los - re_los * cos(r_los) / cLight

# Constrained tobs and dtobs
tobsf = float(argv[2])
dtobsf = float(argv[3])

# Plot LOS emission
FILE_NAME = f"new-single_constrained_patch-{Gamma}-{tobsf}-{dtobsf}"
plt.figure()
T = 10 ** np.linspace(log(tm_los), log(tM_los), 2000)
F_los = np.array([simple_hle_patch(t, nuobs, 0, r_los, te_los, re_los, Gamma, Eiso, nup, a, b)
      for t in T])
plt.plot(T - tobs_los, 12 * F_los, label = "LOS and plateau ", color="black")

T = 10 ** np.linspace(log(tM_los - tobs_los), 4, 2000)
Fp = np.array([toy_plateau(t)
      for t in T])
plt.plot(T, Fp / 10, color="black")
# Plot the off-axis patches
for i, r in enumerate([5 / Gamma]):
    for chi in [tv - tj + r, tv, tv + tj - r]:
        # Calculate t_delay and re to have the right tobsf
        re = cLight * dtobsf / (cos(chi - r) - cos(chi + r))
        t_delay = tobsf + tobs_los - re * (1 / beta - cos(chi - r)) / cLight

        # Report re and t_delay
        print(f"r = {r * Gamma} / Gamma, chi = {chi / tv:4.3f} * tv, re = {re:3.2e}, tdelay = {t_delay:3.2f}")

        te = t_delay + re / (beta * cLight)
        tm = te - re * cos(max(0., chi - r)) / cLight
        tM = te - re * cos(min(Pi, chi + r)) / cLight

        T = 10 ** np.linspace(log(tm), log(tM), 2000)
        F = np.array([simple_hle_patch(t, nuobs, chi, r, te, re, Gamma, Eiso, nup, a, b)
              for t in T])

        plt.plot(T - tobs_los, F, linewidth=0.8, linestyle = ls_l[i],
                 color=cmap(chi / (tv + tj)),
                 label=r"($\delta \theta + \psi / 2) / \theta_v = $" + f"{chi / tv:4.2f}")

plt.ylabel(f"Spectral luminosity ({hPlanck * nuobs/keV} keV)")
plt.xlabel("Time since last LOS pulse (s)")
plt.xscale('log')
plt.yscale('log')
plt.xlim(1)
plt.ylim(max(F_los) * 1.e-9, max(F_los))
plt.legend(loc='lower left')
#plt.title(r"$\Gamma$ = " + f"{Gamma}, " + r"$t_{obs}^F$ = " + f"{tobsf}, " + r"$\Delta t_{obs}^F$ = " + f"{dtobsf}")
plt.savefig(f"{PLOT_DIR}/{FILE_NAME}.pdf", bbox_inches='tight')
