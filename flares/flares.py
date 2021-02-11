from sys import exit
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.cm import rainbow
import pandas as pd
from numpy import sin, cos, sqrt
from numpy import log10 as log
import numpy as np
from helpers import *
from sys import argv
from numpy import random
import corner
import emcee

rin = float(argv[1])
tv = float(argv[2])
Gamma0 = float(argv[3])
DGamma = float(argv[4])
Dtej = float(argv[5])
tj = 0.1

cci = pd.read_csv("data/circle_packings.data", sep=" ", names=['r', 'n'])
r = min(cci['r'], key=lambda x:abs(x-rin/tj))
N_patches = list(cci.loc[cci['r'] == r]['n'])[0]
r = r * tj
PLOT_DIR = f"plots/flares"
CIRCLES_FILE = f"data/cci_coords/cci{N_patches}.txt"
FILE_NAME = f"{r}-{tv}-{Gamma0}-{DGamma}-{Dtej}"

print(f"r: {r}, tv: {tv}, Gamma0: {Gamma0}, DGamma: {DGamma}, Dtej: {Dtej}, N_p: {N_patches}")

plt.style.use('presentation')

np.random.seed(52)

# physics
tvar = 1.
Eiso = 1.e50
a = 0
b = -1.5
Ep = 1.00000 * keV
nup = Ep / hPlanck
nuobs = 5.0 * keV / hPlanck
nug = 300. * keV / hPlanck

df = pd.read_csv(CIRCLES_FILE, sep=" ",
        names=["i", "x", "y"])

df['theta'] = tj * sqrt(df['x'] ** 2 + df['y'] ** 2)
df['phi'] = np.arctan2(df['y'], df['x'])
df['chic'] = np.arccos(cos(tv) * cos(df['theta']) + sin(tv) * sin(df['theta']) * cos(df['phi']))
df['chim'] = np.maximum(df['chic'] - r, 0.)
df['chiM'] = np.minimum(df['chic'] + r, Pi)
df['Gamma'] = Gamma0 + DGamma * random.uniform(0, 1, df.shape[0])
df['beta'] = np.sqrt(1 - 1/df['Gamma']**2)
df['Re'] = 2 * df['Gamma'] ** 2 * cLight * tvar
df['t0IS'] = df['Re'] / (df['beta'] * cLight)
df['Tej'] = Dtej * random.uniform(0, 1, df.shape[0])
df['Te'] = df['Re'] / (df['beta'] * cLight) + df['Tej']
df['tm'] = df['Te'] - df['Re'] * cos(df['chim']) / cLight
df['tM'] = df['Te'] - df['Re'] * cos(df['chiM']) / cLight
chicmin = min(df['chic'])
chicmax = max(df['chic'])
tmmin = min(df['tm'])
tMmax = max(df['tM'])

# plt.figure()
# for chic in [0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.45]:
#     alpha = 0.05
#     tm = t0IS - RIS * cos(max(chic - alpha, 0)) / cLight
#     tM = t0IS - RIS * cos(min(chic + alpha, Pi)) / cLight
#     x = 10 ** np.linspace(log(tm), log(tM), 600)
#     y = np.array([Fnu(chic, alpha, t0IS, RIS, nuobs, a) for a in x])
#     plt.plot(x, y, color = cmap((chic-0)/(0.45 - 0)), label=r"$\chi_c = $" + f"{chic}")
#     alpha = 0.005
#     tm = t0IS - RIS * cos(max(chic - alpha, 0)) / cLight
#     tM = t0IS - RIS * cos(min(chic + alpha, Pi)) / cLight
#     x = 10 ** np.linspace(log(tm), log(tM), 600)
#     y = np.array([Fnu(chic, alpha, t0IS, RIS, nuobs, a) for a in x])
#     plt.plot(x, y, "--", color = cmap((chic - 0)/(0.45 - 0)))
# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel(r"$t_{\rm obs}$ [s]")
# plt.ylabel(r"$L^{\rm obs}$ [erg/s/Hz]")
# plt.legend()
# plt.savefig(f"{PLOT_DIR}/pulses{Gamma}.png", bbox_inches='tight')

fig, axs = plt.subplots(2, 2, figsize=(15,15))
sd = axs[0,0].scatter(df['x'], df['y'], c=(df['Tej']), s=490 * (r / 0.01) ** 2,
            cmap=cmap, label="Ejection time [s]")
fig.colorbar(sd, ax=axs[0,0])
axs[0,0].set_aspect('equal', adjustable='box')
axs[0,0].axis('off')
axs[0,0].legend()

sd = axs[0,1].scatter(df['x'], df['y'], c=(df['tm']), s=490 * (r / 0.01) ** 2,
            cmap=cmap, label=r"$t_m$ [s]")
axs[0,1].plot(tv/tj, 0, "+r")
axs[0,1].set_aspect('equal', adjustable='box')
fig.colorbar(sd, ax=axs[0,1])
axs[0,1].axis('off')
axs[0,1].legend()

for i, p in df.iterrows():
    te = p['Te']
    chic = p['chic']
    re = p['Re']
    tm = p['tm']
    tM = p['tM']
    g = p['Gamma']
    b = p['beta']
    x = np.linspace(tm, tM, 300)
    y = np.array([simple_hle_patch(to, nuobs, chic, r, te, re, g, Eiso, nup, a, b) for to in x])
    axs[1,0].plot(x, (10 - 0.3) * keV * y / hPlanck, color = cmap((chic - chicmin)/(chicmax - chicmin)))
axs[1,0].set_xscale('log')
axs[1,0].set_yscale('log')
axs[1,0].set_xlabel(r"$t_{\rm obs}$ [s]")
axs[1,0].set_ylabel(r"$L^{\rm iso}$ [erg/s]")

N = 400
T = 10 ** np.linspace(log(tmmin), log(tMmax), N)
Fx = np.zeros(N)
Fg = np.zeros(N)
for i, p in df.iterrows():
    te = p['Te']
    chic = p['chic']
    re = p['Re']
    tm = p['tm']
    tM = p['tM']
    g = p['Gamma']
    b = p['beta']
    y = np.array([(simple_hle_patch(to, nuobs, chic, r, te, re, g, Eiso, nup, a, b) if ((to < tM) and (tm < to)) else 0.) for to in T])
    Fx = Fx + y
axs[1,1].plot(T, (10 - 0.3) * keV * Fx / hPlanck, label=f"patch = {r/(0.5/Gamma0):3.2f} x 1/" + r"$\Gamma$, $\Gamma$ = " + f"{Gamma0}, " + r"$\theta_v/\theta_j$ = " + f"{tv/tj}" + r", $\Delta t_{ej}$ = " + f"{Dtej} s")
axs[1,1].set_xscale('log')
axs[1,1].set_yscale('log')
axs[1,1].set_xlabel(r"$t_{\rm obs}$ [s]")
axs[1,1].set_ylabel(r"$L^{\rm iso}$ [erg/s]")
axs[1,1].legend()
plt.savefig(f"{PLOT_DIR}/{FILE_NAME}_total.png", bbox_inches='tight')

# plt.figure()
# x = np.linspace(-1, 1, 50)
# thvl = [0.1, 0.15, 0.2, 0.25, 0.3]
# for thv in thvl:
#     y = [2 * sin(thv + a * tj) * 100  ** 2 * max(0.5/100,
#            sin(thv + a * tj)) * np.arccos((cos(tj) -
#            cos(thv) * cos(thv + a * tj))/(sin(thv) * sin(thv + a * tj)))
#            for a in x]
#     plt.plot(x, y, color=cmap((thv-0.1)/(0.3-0.1)), label=r"$\theta_v$" + f" = {thv}")
# plt.xlabel(r"$(\theta_v - \chi_0)/\theta_j$")
# plt.ylabel(r"$\Delta^{\rm min} t_{ej}$ [s]")
# plt.legend()
# plt.savefig(f"{PLOT_DIR}/success100.png", bbox_inches='tight')
#
# plt.figure()
# for thv in thvl:
#     y = [2 * sin(thv + a * tj) * 400  ** 2 * max(0.5/400,
#            sin(thv + a * tj)) * np.arccos((cos(tj) -
#            cos(thv) * cos(thv + a * tj))/(sin(thv) * sin(thv + a * tj)))
#            for a in x]
#     plt.plot(x, y, "--", color=cmap((thv-0.1)/(0.3-0.1)), label=r"$\theta_v$" + f" = {thv}")
# plt.xlabel(r"$(\theta_v - \chi_0)/\theta_j$")
# plt.ylabel(r"$\Delta^{\rm min} t_{ej}$ [s]")
# plt.savefig(f"{PLOT_DIR}/success400.png", bbox_inches='tight')
