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

Eiso = 1.e52
tj = 0.1
a = 0.
b = -1.5
Gamma = float(argv[1])
t_delay = float(argv[2])
beta = sqrt(1 - 1 / Gamma ** 2)

r_los = 1 / Gamma
re_los = 1.e14
te_los = re_los / (beta * cLight)
tobs_los = te_los - re_los / cLight
tm_los = te_los - re_los / cLight
tM_los = te_los - re_los * cos(r_los) / cLight
FILE_NAME = f"single_delayed_patch-{Gamma}-{t_delay}"

nuobs = 5.0 * keV / hPlanck
nug = 300. * keV / hPlanck
Ep = 10 * keV
nup = Ep / hPlanck

r = 5 / Gamma
plt.figure()
for tv in [tj / 2, tj, 2 * tj, 5 *  tj]:
    for re in [re_los / 2, re_los, 2 * re_los]:
        ls = '-' if (re == re_los / 2) else (':' if (re == re_los) else "-.")
        te = t_delay + re / (beta * cLight)
        tm = te - re * cos(max(0., tv - r)) / cLight
        tM = te - re * cos(min(Pi, tv + r)) / cLight

        T = 10 ** np.linspace(log(tm), log(tM), 2000)
        F = np.array([simple_hle_patch(t, nuobs, tv, r, te, re, Gamma, Eiso, nup, a, b)
              for t in T])
        plt.plot(T - tobs_los, F, linestyle = ls,
                       color= cmap(0.5 * tv / tj), label=f"{tv/tj},{re/re_los}")
T = 10 ** np.linspace(log(tm_los), log(tM_los), 2000)
F = np.array([simple_hle_patch(t, nuobs, 0, r_los, te_los, re_los, Gamma, Eiso, nup, a, b)
      for t in T])

plt.plot(T - tobs_los, F, label = "LOS", color="black")
plt.ylabel(f"Spectral luminosity ({hPlanck * nuobs/keV} keV)")
plt.xlabel("Time since last LOS pulse (s)")
plt.xscale('log')
plt.yscale('log')
plt.xlim((tM_los - tm_los) / 2)
plt.ylim(max(F) * 1.e-8, max(F))
plt.legend()
plt.savefig(f"{PLOT_DIR}/{FILE_NAME}.pdf", bbox_inches='tight')
