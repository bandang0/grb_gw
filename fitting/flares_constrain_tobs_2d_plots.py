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
Ep = 0.1 * keV
nup = Ep / hPlanck
nuobs = 5.0 * keV / hPlanck

def t_delay(tobs, dtobs, chi, r, g):
    beta = sqrt(1 - 1 / g ** 2)
    ret = tobs + re_los * (1 / beta - 1) / cLight - dtobs * (1 / beta - cos(chi - r)) / (cos(chi - r) - cos(chi + r))
    return max(-100, min(100, ret))

def Re(dtobs, chi, r):
    return cLight * dtobs / (cos(chi - r) - cos(chi + r))

def t_ar(t_ej, chi, r, re, beta):
    return t_ej + re / (beta * cLight) - re * cos(max(0., chi - r)) / cLight

def dtobs(chi, r, re, beta):
    return re * (cos(max(0, chi-r) - cos(chi + r)))/cLight

def delta_t_t(t_ej, chi, r, re, beta):
    t = t_ar(t_ej, chi, r, re, beta)
    if t < 1.:
        return 0
    else:
        return dtobs(chi, r, re, beta)/t

tej_l = np.linspace(-100, 200, 200)
re_l = 10 ** np.linspace(13, 16, 200)

for g in [250, 125]:
    r = 0.5 / g
    beta = sqrt(1 - 1 / g ** 2)
    for chi in [tv - tj + r, tv - tj / 2, tv]:
        FILE_NAME_1 = f"dtt-{g}-{chi/tv:.2f}"
        plt.figure()
        T = log(np.array([[max(0,t_ar(t_ej, chi, r, re, beta)) for t_ej in tej_l] for re in re_l]))
        plt.contourf(tej_l, re_l, T, cmap=cmap, levels=20, vmin=0, vmax=3.5)
        #plt.colorbar(label=r"$t_{ar}$ (s)", extend="both")

        U = log(np.array([[XRT_c * peak_patch_luminosity(nuobs, chi, r, t_ej + re/(beta * cLight), re, g, Eiso, nup, a, b) for t_ej in tej_l] for re in re_l]))
        plt.contour(tej_l, re_l, U, cmap='jet', levels=10, vmin=42, vmax=53)
        #plt.colorbar(label="Peak Luminosity [0.3-30] keV")

        V = log(np.array([[dtobs(chi, r, re, beta) for t_ej in tej_l] for re in re_l]))
        plt.contour(tej_l, re_l, V, cmap='Greys', levels=10, vmin = 2, vmax = 6)
        #plt.colorbar(label="$\Delta t_{obs}$ (s)")

        plt.ylabel(r"$R_e$ (cm)")
        plt.xlabel(r"$t_{ej}$ (s)")
        plt.yscale('log')
        plt.title(r"$\Gamma$ = " + f"{g}, " + r"$\delta \theta / \theta_v$ = " + f"{chi/tv:.2f},  "+ r"$r \Gamma$ = " + f"{r * g:.2f}" )
        plt.savefig(f"{PLOT_DIR}/{FILE_NAME_1}.pdf", bbox_inches='tight')

plt.figure()
T = np.array([[ 0 if t_ej > 0 else 3.5 for t_ej in tej_l] for re in re_l])
plt.contourf(tej_l, re_l, T, cmap=cmap, levels=20, vmin=0, vmax=3.5)
plt.colorbar(label=r"$\log t_{ar}$ (s)", extend="both")

U = np.array([[42 if t_ej > 0  else 53 for t_ej in tej_l] for re in re_l])
plt.contour(tej_l, re_l, U, cmap='jet', levels=10, vmin=42, vmax=53)
plt.colorbar(label=r"$\log L_{p}^{XRT}$ (erg/s)")

V = np.array([[2 if t_ej > 0 else 6 for t_ej in tej_l] for re in re_l])
plt.contour(tej_l, re_l, V, cmap='Greys', levels=10, vmin = 2, vmax = 6)
plt.colorbar(label="$\log \Delta t_{obs}$ (s)")

plt.ylabel(r"$R_e$ (cm)")
plt.xlabel(r"$t_{ej}$ (s)")
plt.yscale('log')
plt.savefig(f"{PLOT_DIR}/axis.pdf", bbox_inches='tight')
