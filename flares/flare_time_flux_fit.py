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

if case == "050406":
    LF = 1.e48
    TF = 218
    T90 = 22
if case == "050713A":
    LF = 5.e49
    TF = 100
    T90 = 125
if case == "050713Abis":
    LF = 1.e49
    TF = 158
    T90 = 125
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
if case == "A":
    LF = 3.e48
    TF = 100
    T90 = 30
if case == 'average':
    LF = 3.e49
    TF = 50
    T90 = 30
if case == 'average2':
    LF = 3.e48
    TF = 50
    T90 = 30
if case == 'median':
    LF = 1.8e49
    TF = 85
    T90 = 40
if case == 'early':
    LF = 2.8e49
    TF = 74
    T90 = 40
fgeo = 0.5
def Eiso_fit(tf, lf, tau, tej):
    return tau * lf * ((tf - tej)/tau) ** 3 / fgeo

def S_factor(lf, tau, e):
    return (tau * lf / (e * fgeo)) ** (-1/3)

tau_l = np.logspace(log(0.09), log(20), 200)
Emin = 1.e60
Emax = 1.50
fig, ax = plt.subplots()
for tej in filter(lambda x: x < TF, [0, T90 / 2, T90, 1.5 * T90, 2 * T90, 3 * T90, 5 * T90]):
    Eiso_l = [Eiso_fit(TF, LF, tau, tej) for tau in tau_l]
    Emin = min(Emin, min(Eiso_l))
    Emax = max(Emax, max(Eiso_l))
    plt.plot(tau_l, Eiso_l, label = r"$t_{\rm ej} /T_{90}$ = " + f"{tej / T90:.1f}")

tau_p = np.logspace(-1, log(5))
E_p = [max(10**(3) * LF * tau/fgeo,  Eiso_fit(TF, LF, tau, 2 * T90)) for tau in tau_p]
E_sup = [3.e53 for t in tau_p]
plt.fill_between(tau_p, E_p, np.maximum(E_sup, E_p), color="purple", alpha=0.5)
E_range = np.logspace(log(Emin), log(Emax), 200)
S = [[S_factor(LF, tau, e) for tau in tau_l] for e in E_range]
CS = ax.contour(tau_l, E_range, S, levels = [1, 2, 5, 10, 20, 30, 50], colors="black", label=r"$\mathcal{S}$")
ax.clabel(CS, inline=1, inline_spacing=0, fmt='%1.1f', rightside_up=0, fontsize="smaller", use_clabeltext=1)
plt.text(1, 10 ** (55), r"$L^{\rm flare}$" + f" = {LF} erg/s,\n" + r"$t^{\rm flare}$" + f" = {TF} s", bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
plt.ylim([1.e51, 5.e55])
plt.xlabel(r"$\tau$ [s]")
plt.ylabel(r"$E_{\rm iso}$ [erg]")
plt.xscale('log')
plt.yscale('log')
plt.legend(loc="upper left")
plt.title(f"Median XRT flare " + r"($T_{90}$" +f" = {T90} s)")
plt.savefig(f"{PLOT_DIR}/{case}_2dfit.pdf", bbox_inches='tight')
