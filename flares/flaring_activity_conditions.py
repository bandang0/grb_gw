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
import sys

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

PLOT_DIR = f"plots/flares"

# Fixed parameters
Eiso = 1.e54
tj = 0.1
tv = 0.13
a = -0.5
b = -2.0
Ep = 1. * keV
nup = Ep / hPlanck
nuobs = 5.0 * keV / hPlanck
nug = 75. * keV / hPlanck

if "wind" in argv:
    case = "wind"
elif "ism" in argv:
    case = "ism"
else:
    print("Must specify case...")
    exit(0)

afterglow_df = pd.read_csv(f"data/afterglow_{case}.csv", sep=" ", names = ['tobs', 'fag'])
afterglow_df = afterglow_df.sort_values("tobs")
if case == "wind":
    log_afterglow = interp1d(afterglow_df.tobs, afterglow_df.fag, fill_value='extrapolate')
    afterglow_spectral = lambda t: 10 ** log_afterglow(log(t)) / XRT_c
elif case == "ism":
    # Must place at density of 1.e-3
    time_f = -1.
    flux_f = 3 * 0.38
    afterglow_l = 10 ** afterglow_df.fag
    toy_esd_l = XRT_c * toy_esd(10 ** afterglow_df.tobs) * 0.6
    only_plateau_l = log(afterglow_l - toy_esd_l)
    nans = np.count_nonzero(np.isnan(only_plateau_l))
    print(nans)
    log_plateau = interp1d(afterglow_df.tobs[nans:] - time_f, only_plateau_l[nans:] - flux_f, fill_value='extrapolate')
    afterglow_spectral = lambda t: (toy_esd(t) + 10 ** log_plateau(log(t))/XRT_c)

# Flares to plot
letters = ['A', 'B', 'C', 'D', 'E', 'F', 'I', 'J', 'K']
taus = {'A': 1, 'B': 4, 'C' : 0.5, 'D': 0.5, 'E': 1, 'F': 9, 'G': 1, 'H': 1, 'I':  4, 'J': 2, 'K': 9}
tejs = {'A': 50, 'B': 60, 'C': 60, 'D': 60, 'E': 50, 'F': 80, 'G': 20, 'H': 5, 'I': 0, 'J':  60, 'K': 80}
gammas = {'A': 100, 'B': 200, 'C': 200, 'D': 200, 'E': 200, 'F': 400, 'G':200, 'H': 400, 'I': 400, 'J': 400, 'K': 400}
dts =  {'A':.03, 'B': .08, 'C': .08, 'D': .08, 'E': .03, 'F': .03, 'G': .03, 'H': 0.03, 'I': 0.03, 'J': 0.03, 'K': 0.03}
sizes = {'A':9, 'B': 6, 'C': 3, 'D': 9, 'E': 3, 'F': 9, 'G': 3, 'H': 3, 'I': 3, 'J': 6, 'K': 3}

print("Plotting flares over afterglow")
plt.figure()
tobs_l = np.logspace(0, 6, 1000)

plt.plot(tobs_l, XRT_c * afterglow_spectral(tobs_l), color="black")
for i, letter in enumerate(letters):
    Gamma = gammas[letter]
    tau = taus[letter]
    t_delay = tejs[letter]
    dia = sizes[letter] / Gamma
    r = dia / 2
    chi = dts[letter] + r
    beta = sqrt(1 - 1 / Gamma ** 2)
    re = 2 * cLight * tau * Gamma ** 2
    S = (1 - beta * cos(chi))/(1 - beta)
    te = t_delay + re / (beta * cLight)
    tm = te - re * cos(max(0., chi - r)) / cLight
    tM = te - re * cos(min(Pi, chi + r)) / cLight
    T_patch = 10 ** np.linspace(log(tm), log(tM), 2000)
    L_patch = np.array([simple_hle_patch_band(t, nuobs, chi, r, te, re, Gamma, Eiso, nup, a, b, XRT_0, XRT_1) for t in T_patch])
    t_above = [T_patch[j] for j in range(len(T_patch)) if L_patch[j] > (XRT_c * afterglow_spectral(T_patch[j]))]
    eta_t = (tM - tm) / tm
    eta_p = 4 * r / chi
    if len(t_above) == 0:
        eta_o = 0.
    else:
        eta_o = (max(t_above) - min(t_above)) / max(t_above)
    L_on = peak_patch_luminosity_band(nug, 0, r, te, re, Gamma, Eiso, nup, a, b, BAT_0, BAT_1)
    peak_patch = peak_patch_luminosity_spectral(nuobs, chi, r, te, re, Gamma, Eiso, nup, a, b)
    patch_afterglow = log(peak_patch / afterglow_spectral(tm))

    plt.plot(T_patch, L_patch, linestyle=ls_l[i % 4], color=colors_l[i], label=f"{letter}")

    print(f"{letter} & {Gamma} & {sizes[letter]} & {chi / tv:3.2f} & {tau:3.2f} & {t_delay:6.2f} & {S:7.2f} & {re:3.2e} & {eta_t:3.2f} & {eta_p:3.2f} & {eta_o:3.2f} & -- & {L_on:3.2e} \\\\")

plt.ylabel(f"Luminosity [0.3-30] keV (erg/s)")
plt.xlabel("Time (s)")
plt.xscale('log')
plt.yscale('log')
plt.ylim([4.e43, 1.e51])
plt.xlim(30, 10**4.5)
plt.legend(loc='upper right')
plt.savefig(f"{PLOT_DIR}/{case}_flares.pdf", bbox_inches='tight')

#plot patches
plt.figure()
fig = plt.gcf()
ax = fig.gca()
plt.gca().set_aspect('equal', adjustable='box')
ax.add_artist(plt.Circle((0, 0), tj, edgecolor='black', fill=False))

for i, letter in enumerate(letters):
    ax.add_artist(plt.Circle((dts[letter] + 0.5 * sizes[letter] / gammas[letter] - tv, 0), 0.5 * sizes[letter] / gammas[letter], linewidth=2, linestyle=ls_l[i % 4], label=f"{letter}", fill=False, edgecolor=colors_l[i ]))
plt.scatter([-tv], [0], marker = '+', color='red')
plt.xlim([-1.1 * tv, tv])
plt.ylim([-tv, tv])
plt.savefig(f"{PLOT_DIR}/positions.pdf")

print("Plotting activity conditions")
tej_l = np.linspace(-10, 100, 200)
tau_l = np.logspace(-3, 2, 200)

for g in [100, 200, 400]:
    r = 1.5 / g
    beta = sqrt(1 - 1 / g ** 2)
    for chi in [tv - tj, tv - tj / 2]:
        FILE_NAME = f"activity-{g}-{chi/tv:.2f}"
        dcos = cos(chi - r) - cos(chi + r)
        ccm = cos(chi - r)
        fig, ax = plt.subplots()

        # Off axis luminosity
        L = log(np.array([[peak_patch_luminosity_spectral(nuobs, chi, r, tej + 2 * tau * g ** 2 / beta, 2 * cLight * tau * g ** 2, g, Eiso, nup, a, b)/afterglow_spectral(tej + 2 * tau * g ** 2 * (1 / beta - ccm)) for tau in tau_l] for tej in tej_l]))
        CS = ax.contour(tau_l, tej_l, L, 7, colors="k")
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
                 label=r"$t_{\rm flare}$ = 30, 100, 300, 500, 1000 [s]", color = "blue")
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
                 label=r"$\eta$ = 0.03, 0.1, 0.3, 0.5, 1", color = "red")
        plt.plot(tau_l, [(1/0.5) * 2 * tau * g **2 * dcos - 2 * tau * g ** 2  * (1 - beta * ccm) / beta for tau in tau_l],
                 linestyle="--", color = "red")
        plt.plot(tau_l, Y2,
                 linestyle=":", color = "red")
        plt.fill_between(tau_l, Y1, Y2, color="red", alpha=0.2)

        # On axis luminosity
        plt.vlines(BAT_c * Eiso / (2 * g * nup * 1.e50), ymin = -10, ymax = 100, color="green", linestyle=":")
        plt.vlines(BAT_c * Eiso / (2 * g * nup * 1.e51), ymin = -10, ymax = 100, color="green", linestyle="--")
        plt.vlines(BAT_c * Eiso / (2 * g * nup * 1.e52), ymin = -10, ymax = 100, color="green", label=r"$L_{\rm on, BAT} = 10^{50}, 10^{51}, 10^{52}, 10^{53}, 10^{54}$ [erg/s]")
        plt.vlines(BAT_c * Eiso / (2 * g * nup * 1.e53), ymin = -10, ymax = 100, color="green", linestyle="--")
        plt.vlines(BAT_c * Eiso / (2 * g * nup * 1.e54), ymin = -10, ymax = 100, color="green", linestyle=":")
        plt.axvspan(BAT_c * Eiso / (2 * g * nup * 1.e50), BAT_c * Eiso / (2 * g * nup * 1.e54), color="green", alpha=0.2)

        # Flare points
        for letter in letters:
            if gammas[letter] == g and dts[letter] == chi:
                plt.plot([taus[letter]], [tejs[letter]], marker='o', color="black")
                plt.text(taus[letter], tejs[letter] + 5, letter)
        plt.xscale("log")
        plt.ylim(-10, 100)
        plt.xlim(1.e-3, 100)
        if g == 100 or chi/tv > 0.6:
            plt.ylabel(r"$t_{ej}$ [s]")
        if g == 400 or chi/tv > 0.6:
            plt.xlabel(r"$\tau~(= \Delta t_{on})$ [s]")
        plt.title(r"$\Gamma$ = " + f"{g}, "+ r"$\delta \theta / \theta_v$ = " + f"{chi/tv:.2f}")
        if g == 400 and chi/tv > 0.6:
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.savefig(f"{PLOT_DIR}/{case}_{FILE_NAME}.pdf", bbox_inches='tight')
