from sys import exit
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.colors import Normalize
import pandas as pd
from numpy import sin, cos, sqrt
from numpy import log10 as log
from scipy.interpolate import interp1d
from astropy.io import fits
import numpy as np
from helpers import *
from sys import argv

PLOT_DIR = f"plots/flares"
plt.style.use('presentation')

# Fixed parameters
Eiso = 1.e52
tj = 0.1
tv = 0.13
a = 0.
b = -1.5
Ep = 1. * keV
nup = Ep / hPlanck
nuobs = 5.0 * keV / hPlanck

# Plateau data
with fits.open('data/plateau_03.fits') as hdul:
    data = hdul[1].data
plateau = interp1d(data['tobs'], data['dtheta=0.03'], bounds_error=False, fill_value=0.)

# Flares to plot
letters =['A', 'B', 'C','D', 'E',  'F', 'G',  'I',  'J', 'K', 'L']
taus =   [1,   3,   0.8, 0.1, 1,   .7,  .2,   1,    4,    0.2, 9]
tejs =   [60,  60,  60,  60,  15,  20,  0,    5,    0,    60,  80]
gammas = [400, 200, 200, 200, 200, 400, 400,  400,  400,  400, 400]
dts =    [.03, .08, .08, .08, .03, .08, .08,  0.03, 0.03, 0.03, 0.03]

print("Plotting flares over afterglow")
plt.figure()
tobs_l = np.logspace(-1.5, 4.5, 400)
afterglow = interp1d(tobs_l, [plateau(u) + toy_esd(u) for u in tobs_l], bounds_error=False, fill_value = 0.)
plt.plot(tobs_l, afterglow(tobs_l), color="black")
plt.plot(tobs_l, toy_esd(tobs_l), color="grey", linestyle="--", linewidth=0.8)
plt.plot(tobs_l, plateau(tobs_l), color="grey", linestyle="--",  linewidth=0.8)

for i in range(len(letters)):
    letter = letters[i]
    Gamma = gammas[i]
    tau = taus[i]
    t_delay = tejs[i]
    chi = dts[i]
    beta = sqrt(1 - 1 / Gamma ** 2)
    re = 2 * cLight * tau * Gamma ** 2
    dia = 3 / Gamma
    r = dia / 2
    S = (1 - beta * cos(chi))/(1 - beta)
    L_on = BAT_c * Eiso / (2 * Gamma * tau * nup)
    te = t_delay + re / (beta * cLight)
    tm = te - re * cos(max(0., chi - r)) / cLight
    tM = te - re * cos(min(Pi, chi + r)) / cLight
    eta_t = (tM - tm) / (tm)
    peak_patch = XRT_c * peak_patch_luminosity(nuobs, chi, r, te, re, Gamma, Eiso, nup, a, b)
    patch_afterglow = log(peak_patch / afterglow(tm))
    print(f"{letter} & {Gamma} & {chi / tv:3.2f} & {tau} & {t_delay:3.2f} & {S:3.2f} & {re:3.2e} & {eta_t:3.2f} & -- & {L_on:3.2e} \\\\")

    T = 10 ** np.linspace(log(tm), log(tM), 2000)
    L_p = XRT_c * np.array([simple_hle_patch(t, nuobs, chi, r, te, re, Gamma, Eiso, nup, a, b)
          for t in T])
    plt.plot(T, L_p, linestyle=ls_l[i % 4], label=f"{letter}")

plt.ylabel(f"Luminosity [0.3-30] keV (erg/s)")
plt.xlabel("Time (s)")
plt.xscale('log')
plt.yscale('log')
plt.ylim([4.e44, 1.e51])
plt.xlim(20, 1.e4)
plt.legend(loc='upper right')
plt.savefig(f"{PLOT_DIR}/flares.pdf", bbox_inches='tight')

print("Plotting activity conditions")
tej_l = np.linspace(-10, 100, 200)
tau_l = np.logspace(-3, 2, 200)

for g in [200, 400]:
    r = 1.5 / g
    beta = sqrt(1 - 1 / g ** 2)
    for chi in [tv - tj, tv - tj / 2]:
        FILE_NAME = f"activity-{g}-{chi/tv:.2f}"
        dcos = cos(chi - r) - cos(chi + r)
        ccm = cos(chi - r)
        fig, ax = plt.subplots()

        # Off axis luminosity
        L = log(np.array([[XRT_c * peak_patch_luminosity(nuobs, chi, r, tej + 2 * tau * g ** 2 / beta, 2 * cLight * tau * g ** 2, g, Eiso, nup, a, b)/afterglow(tej + 2 * tau * g ** 2 * (1 / beta - ccm)) for tau in tau_l] for tej in tej_l]))
        CS = ax.contour(tau_l, tej_l, L, 9, colors="k")
        ax.clabel(CS, inline=1, inline_spacing=0, fmt='%1.1f', rightside_up=0, fontsize="smaller", use_clabeltext=1)
        #plt.pcolor(tau_l, tej_l, L, cmap="Greys", norm=Normalize(-3, 3))
        #plt.colorbar(label="$\log L_{off}$ / AG")

        # Arrival time off
        Y1 = [30 - 2 * tau * g**2 * (1 - beta * ccm) / beta for tau in tau_l]
        Y2 = [1000 - 2 * tau * g**2 * (1 - beta * ccm) / beta for tau in tau_l]
        plt.plot(tau_l, Y1,
                 color = "blue", linestyle=":")
        plt.plot(tau_l, [100 - 2 * tau * g**2 * (1 - beta * ccm) / beta for tau in tau_l],
                 color = "blue", linestyle="--")
        plt.plot(tau_l, [300 - 2 * tau * g**2 * (1 - beta * ccm) / beta for tau in tau_l],
                 label=r"$t_{off}$ = 30, 100, 300, 500, 1000s", color = "blue")
        plt.plot(tau_l, [500 - 2 * tau * g**2 * (1 - beta * ccm) / beta for tau in tau_l],
                 color = "blue", linestyle="--")
        plt.plot(tau_l, Y2,
                  color = "blue", linestyle=":")
        plt.fill_between(tau_l, Y1, Y2, color="blue", alpha=0.2)

        # Aspect ratio
        Y1 = [(1/0.03) * 2 * tau * g **2 * dcos - 2 * tau * g ** 2  * (1 - beta * ccm) / beta for tau in tau_l]
        Y2 = [(1/1) * 2 * tau * g **2 * dcos - 2 * tau * g ** 2  * (1 - beta * ccm) / beta for tau in tau_l]
        plt.plot(tau_l, Y1 ,
                 linestyle=":", color = "red")
        plt.plot(tau_l, [(1/0.1) * 2 * tau * g **2 * dcos - 2 * tau * g ** 2  * (1 - beta * ccm) / beta for tau in tau_l],
                 linestyle="--", color = "red")
        plt.plot(tau_l, [(1/0.3) * 2 * tau * g **2 * dcos - 2 * tau * g ** 2  * (1 - beta * ccm) / beta for tau in tau_l],
                 label=r"$\eta$ = .03, .1, .3, .5, 1", color = "red")
        plt.plot(tau_l, [(1/0.5) * 2 * tau * g **2 * dcos - 2 * tau * g ** 2  * (1 - beta * ccm) / beta for tau in tau_l],
                 linestyle="--", color = "red")
        plt.plot(tau_l, Y2,
                 linestyle=":", color = "red")
        plt.fill_between(tau_l, Y1, Y2, color="red", alpha=0.2)

        # On axis luminosity
        plt.vlines(BAT_c * Eiso / (2 * g * nup * 1.e49), ymin = -10, ymax = 100, color="green", linestyle=":")
        plt.vlines(BAT_c * Eiso / (2 * g * nup * 1.e50), ymin = -10, ymax = 100, color="green", linestyle="--")
        plt.vlines(BAT_c * Eiso / (2 * g * nup * 1.e51), ymin = -10, ymax = 100, color="green", label=r"$\log L_{on, BAT}$ = 49, 50, 51, 52, 53")
        plt.vlines(BAT_c * Eiso / (2 * g * nup * 1.e52), ymin = -10, ymax = 100, color="green", linestyle="--")
        plt.vlines(BAT_c * Eiso / (2 * g * nup * 1.e53), ymin = -10, ymax = 100, color="green", linestyle=":")
        plt.axvspan(BAT_c * Eiso / (2 * g * nup * 1.e49), BAT_c * Eiso / (2 * g * nup * 1.e53), color="green", alpha=0.2)

        # Flare points
        for i in range(len(letters)):
            if gammas[i] == g and dts[i] == chi:
                plt.plot([taus[i]], [tejs[i]], marker='o', color="black")
                plt.text(taus[i], tejs[i] + 5, letters[i])
        plt.xscale("log")
        plt.ylim(-10, 100)
        plt.xlim(1.e-3, 100)
        plt.ylabel(r"$t_{ej}$ [s]")
        plt.xlabel(r"$\tau~(= \Delta t_{on})$ [s]")
        plt.title(r"$\Gamma$ = " + f"{g}, "+ r"$\delta \theta / \theta_v$ = s" + f"{chi/tv:.2f}")
        if g >= 1000:
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='xx-small')
        plt.savefig(f"{PLOT_DIR}/{FILE_NAME}.pdf", bbox_inches='tight')
