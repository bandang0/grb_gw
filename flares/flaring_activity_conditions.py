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
Eiso = 1.e53
Eiso_los = 1.2e52
tj = 0.1
tv = 0.13
a = -0.1
b = -1.2
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

#afterglow_df = pd.read_csv(f"data/afterglow_{case}.csv", sep=" ", names = ['tobs', 'fag'])
#afterglow_df = afterglow_df.sort_values("tobs")
plateau_wind = lambda x: 1.8e45
plateau_ism = lambda x: 6.6e45
# if case == "wind":
#     log_afterglow = interp1d(afterglow_df.tobs, afterglow_df.fag, fill_value='extrapolate')
#     afterglow_spectral = lambda t: 10 ** log_afterglow(log(t)) / XRT_c
#     afterglow_band = lambda t: (1.8e45 + 1.e48 * (t / 100)**(-3))
#     plateau_level = 1.8e45
# elif case == "ism":
#     # Must place at density of 1.e-3
#     time_f = -1.
#     flux_f = 3 * 0.38
#     afterglow_l = 10 ** afterglow_df.fag
#     toy_esd_l = XRT_c * toy_esd(10 ** afterglow_df.tobs) * 0.6
#     only_plateau_l = log(afterglow_l - toy_esd_l)
#     nans = np.count_nonzero(np.isnan(only_plateau_l))
#     log_plateau = interp1d(afterglow_df.tobs[nans:] - time_f, only_plateau_l[nans:] - flux_f, fill_value='extrapolate')
#     afterglow_spectral = lambda t: (toy_esd(t) + 10 ** log_plateau(log(t))/XRT_c)
#     afterglow_band = lambda t: (6.6e45 + 1.e48 * (t / 100)**(-3))
#     plateau_level = 6.6e45

# Flares to plot
letters = ['A', 'C', 'D', 'G', 'H', 'I', 'K']
taus = {'A': 1, 'B': 1, 'C' : 30, 'D': 10., 'E': 2, 'F': 9, 'G': 10, 'H': 3., 'I': 1, 'J': 10, 'K': 10}
tejs = {'A': 50, 'B': 50, 'C': 60, 'D': 0., 'E': 50, 'F': 80, 'G': 50, 'H': 80, 'I': 40, 'J': 60, 'K': 80}
gammas = {'A': 100, 'B': 50, 'C': 200, 'D': 100, 'E': 100, 'F': 200, 'G': 100, 'H': 100, 'I': 200, 'J': 250, 'K': 200}
tis =  {'A': .03, 'B': .08, 'C': .03, 'D': .03, 'E': .03, 'F': .03, 'G': .03, 'H': 0.03, 'I': 0.03, 'J': 0.02, 'K': 0.03}
sizes = {'A': 9, 'B': 6, 'C': 6, 'D': 3, 'E': 5, 'F': 9, 'G': 3, 'H': 3, 'I': 3, 'J': 3, 'K': 3}

print("Plotting flares over afterglow")
tej_los = 30
tau_los = 2
g_los = 50
r_los = 1
chi_los = 0
beta_los = sqrt(1 - 1 / g_los ** 2)
re_los = 2 * cLight * g_los ** 2 * tau_los
te_los = tej_los + re_los / (cLight * beta_los)
tm_los = te_los - re_los * cos(max(0., - r_los)) / cLight
tM_los = te_los - re_los * cos(min(Pi, r_los)) / cLight
print(f"{tm_los}, {tM_los}, {re_los:.2e}")
tobs_l = np.logspace(0, 6, 1000)
L_los = np.array([simple_hle_patch_band_zero(t, nuobs, chi_los, r_los, te_los, re_los, g_los, Eiso_los, nup, a, b, XRT_0, XRT_1) for t in tobs_l])
L_wind = np.array([plateau_wind(tobs) for tobs in tobs_l])
L_ism = np.array([plateau_ism(tobs) for tobs in tobs_l])
plt.figure()
plt.plot(tobs_l, L_wind , color="black", linestyle=":")
plt.plot(tobs_l, L_ism, color="black", linestyle=":")
plt.plot(tobs_l, L_los, color="black", linestyle=":")
plt.plot(tobs_l, L_los + L_wind, color="black", linestyle="-")
plt.plot(tobs_l, L_los + L_ism, color="black", linestyle="-")
for i, letter in enumerate(letters):
    Gamma = gammas[letter]
    tau = taus[letter]
    t_delay = tejs[letter]
    dia = sizes[letter] / Gamma
    r = dia / 2
    chi = tis[letter] + r
    beta = sqrt(1 - 1 / Gamma ** 2)
    re = 2 * cLight * tau * Gamma ** 2
    S = (1 - beta * cos(chi))/(1 - beta)
    te = t_delay + re / (beta * cLight)
    tm = te - re * cos(max(0., chi - r)) / cLight
    tM = te - re * cos(min(Pi, chi + r)) / cLight
    T_patch = 10 ** np.linspace(log(tm), log(tM), 2000)
    L_patch = np.array([simple_hle_patch_band(t, nuobs, chi, r, te, re, Gamma, Eiso, nup, a, b, XRT_0, XRT_1) for t in T_patch])
    t_above_wind = [T_patch[j] for j in range(len(T_patch)) if L_patch[j] >  (plateau_wind(T_patch[j]) + simple_hle_patch_band_zero(T_patch[j], nuobs, chi_los, r_los, te_los, re_los, g_los, Eiso_los, nup, a, b, XRT_0, XRT_1))]
    t_above_ism = [T_patch[j] for j in range(len(T_patch)) if L_patch[j] >  (plateau_ism(T_patch[j]) + simple_hle_patch_band_zero(T_patch[j], nuobs, chi_los, r_los, te_los, re_los, g_los, Eiso_los, nup, a, b, XRT_0, XRT_1))]
    eta_t = (tM - tm) / tm
    eta_p = 4 * r / tis[letter]
    if len(t_above_wind) == 0:
        eta_o_wind = 0.
    else:
        eta_o_wind = (max(t_above_wind) - min(t_above_wind)) / max(t_above_wind)

    if len(t_above_ism) == 0:
        eta_o_ism = 0.
    else:
        eta_o_ism = (max(t_above_ism) - min(t_above_ism)) / max(t_above_ism)
    L_on = peak_patch_luminosity_band(nug, 0, r, te, re, Gamma, Eiso, nup, a, b, BAT_0, BAT_1)
    #peak_patch = peak_patch_luminosity_band(nuobs, chi, r, te, re, Gamma, Eiso, nup, a, b, XRT_0, XRT_1)
    #patch_afterglow = log(peak_patch / afterglow_band(tm))

    plt.plot(T_patch, L_patch, linestyle=ls_l[i % 4], color=colors_l[i % 10], label=f"{letter}")

    print(f"{letter} & {Gamma:3d} & {sizes[letter]:3d} & {tis[letter] / tv:3.2f} & {tau:5.1f} & {int(t_delay):3d} & {int(S):4d} & {re:3.1e} & {eta_t:3.2f} & {eta_p:3.2f} & {eta_o_wind:3.2f} & {eta_o_ism:3.2f} & {L_on:3.1e} \\\\")

plt.ylabel(f"Luminosity [0.3-30] keV (erg/s)")
plt.xlabel("Time (s)")
plt.xscale('log')
plt.yscale('log')
plt.ylim([4.e44, 1.e49])
plt.xlim([50, 7.e3])
plt.legend(loc='upper right')
plt.savefig(f"{PLOT_DIR}/both_flares.pdf", bbox_inches='tight')

#plot patches
plt.figure()
fig = plt.gcf()
ax = fig.gca()
plt.gca().set_aspect('equal', adjustable='box')
ax.add_artist(plt.Circle((0, 0), tj, edgecolor='black', fill=False))

for i, letter in enumerate(letters):
    ax.add_artist(plt.Circle((tis[letter] + 0.5 * sizes[letter] / gammas[letter] - tv, 0), 0.5 * sizes[letter] / gammas[letter], linewidth=2, linestyle=ls_l[i % 4], label=f"{letter}", fill=False, edgecolor=colors_l[i % 10]))
plt.scatter([-tv], [0], marker = '+', color='red')
plt.scatter([0], [0], marker = '.', color='black')
plt.xlim([-1.1 * tv, tv])
plt.ylim([-tv, tv])
plt.savefig(f"{PLOT_DIR}/positions.pdf")

print("Plotting activity conditions")
tej_l = np.linspace(-10, 100, 200)
tau_l = np.logspace(-3, 2, 200)

for case in ['wind', 'ism']:
    if case == 'wind':
        afterglow_band = lambda x: (plateau_wind(x) + 1.e47 * (x/100) ** (-3))
    if case == 'ism':
        afterglow_band = lambda x: (plateau_ism(x) + 1.e47 * (x/100) ** (-3))
    for g in [50, 100, 200]:
        r = 4 / g
        beta = sqrt(1 - 1 / g ** 2)
        for ti in [tv - tj, tv - tj / 2]:
            print(f"{g}, {ti}")
            FILE_NAME = f"activity-{g}-{ti:.2f}"
            chi = ti + r
            dcos = cos(chi - r) - cos(chi + r)
            ccm = cos(chi - r)
            fig, ax = plt.subplots()

            # Off axis luminosity
            L = log(np.array([[peak_patch_luminosity_band(nuobs, chi, r, tej + 2 * tau * g ** 2 / beta, 2 * cLight * tau * g ** 2, g, Eiso, nup, a, b, XRT_0, XRT_1)/afterglow_band(tej + 2 * tau * g ** 2 * (1 / beta - ccm)) for tau in tau_l] for tej in tej_l]))
            CS = ax.contour(tau_l, tej_l, L, 13, colors="red")
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
            # Y1 = [(1/0.1) * 2 * tau * g **2 * dcos - 2 * tau * g ** 2  * (1 - beta * ccm) / beta for tau in tau_l]
            # Y2 = [(1/2) * 2 * tau * g **2 * dcos - 2 * tau * g ** 2  * (1 - beta * ccm) / beta for tau in tau_l]
            # plt.plot(tau_l, Y1 ,
            #          linestyle=":", color = "red")
            # plt.plot(tau_l, [(1/0.3) * 2 * tau * g **2 * dcos - 2 * tau * g ** 2  * (1 - beta * ccm) / beta for tau in tau_l],
            #          linestyle="--", color = "red")
            # plt.plot(tau_l, [(1/0.5) * 2 * tau * g **2 * dcos - 2 * tau * g ** 2  * (1 - beta * ccm) / beta for tau in tau_l],
            #          label=r"$\Delta t / t$ = 0.1, 0.3, 0.5, 1., 2.", color = "red")
            # plt.plot(tau_l, [(1/1.) * 2 * tau * g **2 * dcos - 2 * tau * g ** 2  * (1 - beta * ccm) / beta for tau in tau_l],
            #          linestyle="--", color = "red")
            # plt.plot(tau_l, Y2,
            #          linestyle=":", color = "red")
            # plt.fill_between(tau_l, Y1, Y2, color="red", alpha=0.2)

            # On axis luminosity
            plt.vlines(BAT_c * Eiso / (2 * g * nup * 1.e50), ymin = -10, ymax = 100, color="green", linestyle=":")
            plt.vlines(BAT_c * Eiso / (2 * g * nup * 1.e51), ymin = -10, ymax = 100, color="green", linestyle="--")
            plt.vlines(BAT_c * Eiso / (2 * g * nup * 1.e52), ymin = -10, ymax = 100, color="green", label=r"$L_{\rm on, BAT} = 10^{50}, 10^{51}, 10^{52}, 10^{53}, 10^{54}$ [erg/s]")
            plt.vlines(BAT_c * Eiso / (2 * g * nup * 1.e53), ymin = -10, ymax = 100, color="green", linestyle="--")
            plt.vlines(BAT_c * Eiso / (2 * g * nup * 1.e54), ymin = -10, ymax = 100, color="green", linestyle=":")
            plt.axvspan(BAT_c * Eiso / (2 * g * nup * 1.e50), BAT_c * Eiso / (2 * g * nup * 1.e54), color="green", alpha=0.2)

            # Flare points
            for letter in letters:
                if gammas[letter] == g and tis[letter] == ti:
                    plt.plot([taus[letter]], [tejs[letter]], marker='o', color="black")
                    plt.text(taus[letter], tejs[letter] + 5, letter)
            plt.xscale("log")
            plt.ylim(-10, 100)
            plt.xlim(1.e-3, 100)
            if g == 100:
                plt.ylabel(r"$t_{ej}$ [s]")
            if case == 'ism':
                plt.xlabel(r"$\tau~(= \Delta t_{on})$ [s]")
            if case == "wind":
                plt.title(r"$\Gamma$ = " + f"{g}, "+ r"$\theta_i$ = " + f"{ti:.2f}" + " rad")
            if g == 50 and ti == 0.03:
                plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.savefig(f"{PLOT_DIR}/{case}_{FILE_NAME}.pdf", bbox_inches='tight')
