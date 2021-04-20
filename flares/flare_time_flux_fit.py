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
case = argv[1]

tau_l = np.logspace(-1.5, 2, 200)

if case == "060204B":
    LF = 4.28e49
    TF = 317.5
    T90 = 139.4
if case == "060413":
    LF = 3.20e48
    TF = 640.1
    T90 = 147.7
if case == "070704":
    LF = 3.68e+49
    TF = 313.2
    T90 = 384.2
if case == "081008":
    LF = 8.33e+49
    TF = 301.7
    T90 = 185.5
if case == "081102":
    LF = 3.08e+48
    TF = 959.5
    T90 = 62.7
if case == "100619A":
    LF = 1.83e+49
    TF = 962.0
    T90 = 97.5
if case == "110407A":
    LF =  5.06e+48
    TF = 443.0
    T90 = 145.5
if case == "120728A":
    LF = 7.13e+47
    TF = 521.7
    T90 = 22.0

def Eiso_fit(tf, lf, tau, tej):
    return tau * lf * ((tf - tej)/tau) ** 3 / 0.14

def S_factor(lf, tau, e):
    return (tau * lf / (e * 0.14)) ** (-1/3)

Emin = 1.e60
Emax = 1.50
fig, ax = plt.subplots()
for tej in filter(lambda x: x < TF, [0, T90 / 2, T90, 2 * T90, 3 * T90, 5 * T90]):
    Eiso_l = [Eiso_fit(TF, LF, tau, tej) for tau in tau_l]
    Emin = min(Emin, min(Eiso_l))
    Emax = max(Emax, max(Eiso_l))
    plt.plot(tau_l, Eiso_l, label = r"$t_{\rm ej} /T_{90}$ = " + f"{tej / T90:.1f}")
print(Emin, Emax)
E_range = np.logspace(log(Emin), log(Emax), 200)
S = [[S_factor(LF, tau, e) for tau in tau_l] for e in E_range]
CS = ax.contour(tau_l, E_range, S, levels = [1, 5, 10, 20, 30, 50, 100, 200, 300, 500, 1000], colors="black", label=r"$\mathcal{S}$")
ax.clabel(CS, inline=1, inline_spacing=0, fmt='%1.1f', rightside_up=0, fontsize="smaller", use_clabeltext=1)
plt.text(0.5, 10 ** (log(Emin) + 0.90 * log(Emax / Emin)), r"$L^{\rm flare}$" + f" = {LF} erg/s,\n" + r"$t^{\rm flare}$" + f" = {TF} s")
plt.ylim([Emin, Emax])
plt.xlabel(r"$\tau$ [s]")
plt.ylabel(r"$E_{\rm iso}$ [erg]")
plt.xscale('log')
plt.yscale('log')
plt.legend(loc="upper left")
plt.title(f"GRB{case} " + r"($T_{90}$" +f" = {T90} s)")
plt.savefig(f"{PLOT_DIR}/{case}_2dfit.pdf", bbox_inches='tight')
