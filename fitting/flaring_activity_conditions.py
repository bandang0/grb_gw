from sys import exit
import matplotlib.pyplot as plt
from matplotlib import colors
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
Ep = 0.1 * keV
nup = Ep / hPlanck
nuobs = 5.0 * keV / hPlanck

# Plateau data
with fits.open('data/plateau_03.fits') as hdul:
    data = hdul[1].data
print(min(data['tobs']))
plateau = interp1d(data['tobs'], data['dtheta=0.03'], bounds_error=False, fill_value=0.)

tej_l = np.linspace(0, 100, 200)
tau_l = 10 ** np.linspace(-3, 2, 200)

for g in [100, 200, 400]:
    theta_p = 1.5 / g
    beta = sqrt(1 - 1 / g ** 2)
    for dtheta in [tv - tj, tv - tj / 2, tv]:
        FILE_NAME = f"activity-{g}-{dtheta/tv:.2f}"
        dcos = cos(dtheta - theta_p) - cos(dtheta + theta_p)
        ccm = cos(dtheta - theta_p)
        plt.figure()

        # Arrival time off
        plt.plot(tau_l, [100 - 2 * tau * g**2 * (1 - beta * cos(dtheta - theta_p)) / beta for tau in tau_l],
                 color = "blue", linestyle=":")
        plt.plot(tau_l, [300 - 2 * tau * g**2 * (1 - beta * cos(dtheta - theta_p)) / beta for tau in tau_l],
                 label=r"$t_{off}$ = 100, 300, 1000s", color = "blue")
        plt.plot(tau_l, [1000 - 2 * tau * g**2 * (1 - beta * cos(dtheta - theta_p)) / beta for tau in tau_l],
                  color = "blue", linestyle=":")

        # Arrival time on
        plt.plot(tau_l, [0 - tau for tau in tau_l],
                 linestyle=":", color = "grey")
        plt.plot(tau_l, [50 - tau for tau in tau_l],
                 label=r"$t_{on}$ = 0, 50, 100s", color = "grey")
        plt.plot(tau_l, [100 - tau for tau in tau_l],
                  linestyle=":", color = "grey")

        # Aspect ratio
        plt.plot(tau_l, [(1/0.1) * 2 * tau * g **2 * dcos - 2 * tau * g ** 2  * (1 - beta * ccm) / beta for tau in tau_l],
                 linestyle=":", color = "orange")
        plt.plot(tau_l, [(1/0.3) * 2 * tau * g **2 * dcos - 2 * tau * g ** 2  * (1 - beta * ccm) / beta for tau in tau_l],
                 label=r"$\eta$ = .1, .3, 1", color = "orange")
        plt.plot(tau_l, [(1/1) * 2 * tau * g **2 * dcos - 2 * tau * g ** 2  * (1 - beta * ccm) / beta for tau in tau_l],
                 linestyle=":", color = "orange")

        # On axis luminosity
        plt.vlines(BAT_c * Eiso / (2 * g * nup * 1.e51), ymin = -10, ymax = 100, color="green", linestyle=":")
        plt.vlines(BAT_c * Eiso / (2 * g * nup * 1.e52), ymin = -10, ymax = 100, color="green", label=r"$L_{on, BAT}$ = 51, 52, 53 erg/s")
        plt.vlines(BAT_c * Eiso / (2 * g * nup * 1.e53), ymin = -10, ymax = 100, color="green", linestyle=":")

        # Off axis luminosity
        L = log(np.array([[XRT_c * peak_patch_luminosity(nuobs, dtheta, theta_p, tej + 2 * tau * g ** 2, 2 * cLight * tau * g ** 2, g, Eiso, nup, a, b)/(toy_esd(tej + 2 * tau * g ** 2 * (1 - beta * ccm) / beta) + plateau(tej + 2 * tau * g ** 2 * (1 - beta * ccm) / beta)) for tau in tau_l] for tej in tej_l]))
        plt.contour(tej_l, tau_l, np.clip(L, -2, 2), cmap=cmap, vmin = -2, vmax = 2, levels=3)
        plt.colorbar(label="$\log L_{off}$ / (ESD + pl.)")

        plt.xscale("log")
        plt.ylim(-10, 100)
        plt.xlim(1.e-3, 100)
        plt.ylabel(r"$t_{ej}$ [s]")
        plt.xlabel(r"$\tau (= \Delta t_{on})$ [s]")
        plt.title(r"$\Gamma =$" + f"{g}, "+ r"$\delta \theta / \theta_v = $" + f"{dtheta/tv:.2f}")
        if g == 100:
            plt.legend()
        plt.savefig(f"{PLOT_DIR}/{FILE_NAME}.pdf", bbox_inches='tight')


exit(0)
T = log(np.array([[max(0,t_ar(t_ej, chi, r, re, beta)) for t_ej in tej_l] for re in re_l]))

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
