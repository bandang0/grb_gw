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

tj = 0.1
tv = 0.12
Np = 1395
Gamma_core = 400
Gamma_LOS = Gamma_core / 4
r = 1 / Gamma_core
T90 = 20

PLOT_DIR = f"plots/flares"
CIRCLES_FILE = f"data/cci_coords/cci{Np}.txt"
FILE_NAME = f"flares_with_plateaus_{Np}-{T90}-{r * Gamma_core}"

print(f"r: {r}, tv: {tv}, Gamma_core: {Gamma_core}, T90: {T90}, N_p: {Np}")
plt.style.use('presentation')
np.random.seed(52)

# physics
tvar = 1.
Eiso_gamma = 1.e52
a = 0
b = -1.5
Ep = 10.0 * keV
nup = Ep / hPlanck
nuobs = 10.0 * keV / hPlanck
nug = 300. * keV / hPlanck

df = pd.read_csv(CIRCLES_FILE, sep=" ",
        names=["i", "x", "y"])

#LOS material
R_los = 2 * Gamma_LOS**2 * cLight * tvar
beta_los = sqrt(1 - 1 / Gamma_LOS**2)
tm_los = R_los / (beta_los * cLight) - R_los / cLight

# Patches
df['theta'] = tj * sqrt(df['x'] ** 2 + df['y'] ** 2)
df['phi'] = np.arctan2(df['y'], df['x'])
df['chic'] = np.arccos(cos(tv) * cos(df['theta']) + sin(tv) * sin(df['theta']) * cos(df['phi']))
df['chim'] = np.maximum(df['chic'] - r, 0.)
df['chiM'] = np.minimum(df['chic'] + r, Pi)
df['Gamma'] = Gamma_core
df['beta'] = np.sqrt(1 - 1/df['Gamma']**2)
df['Re'] = 2 * df['Gamma'] ** 2 * cLight * tvar
df['Tej'] = T90 * random.uniform(-0.5, 0.5, df.shape[0])
df['Te'] = df['Tej'] + df['Re'] / (df['beta'] * cLight)
df['tm'] = df['Te'] - df['Re'] * cos(df['chim']) / cLight
df['tM'] = df['Te'] - df['Re'] * cos(df['chiM']) / cLight
df['deltat'] = df['tM'] - df['tm']
df['fmax'] = 0.
ssize = 70
print(tm_los)
for j in range(Np):
    df.at[j, 'fmax'] = simple_hle_patch(10**((log(df.at[j, 'tM']) + log(df.at[j, 'tm']))/2),
             nuobs, df.at[j, 'chic'], r, df.at[j,'Te'], df.at[j,'Re'], df.at[j,'Gamma'], Eiso_gamma, nup, a, b)

fig, axs = plt.subplots(2, 2, figsize=(15,15))
sd = axs[0,0].scatter(df['tm'], df['deltat'], c=(log(df['fmax'])), s=ssize,
            cmap=cmap, label="Maximum spectral flux @ 10keV [erg/s/Hz]")
fig.colorbar(sd, ax=axs[0,0])
axs[0,0].set_aspect('equal', adjustable='box')
axs[0,0].axis('on')
axs[0,0].set_xscale('log')
axs[0,0].set_yscale('log')
axs[0,0].set_xlabel(r"$t_m$ [s]")
axs[0,0].set_ylabel(r"$\Delta t$ [s]")
axs[0,0].set_xlim([1, 10000])
axs[0,0].set_ylim([1, 10000])
axs[0,0].legend()

sd = axs[0,1].scatter(df['x'], df['y'], c=(df['tm']), s=ssize,
            cmap=cmap, label=r"$t_m$ [s]")
axs[0,1].plot(tv/tj, 0, "+r")
axs[0,1].set_aspect('equal', adjustable='box')
fig.colorbar(sd, ax=axs[0,1])
axs[0,1].axis('off')
axs[0,1].legend()

x = 10 ** np.linspace(2, 6, 300)
axs[1,0].plot(x, [toy_plateau(u) for u in x])
axs[1,0].set_xscale('log')
axs[1,0].set_yscale('log')
axs[1,0].set_xlabel(r"$t_{\rm obs}$ [s]")
axs[1,0].set_ylabel(r"Plateau emission @ 10keV [erg/s/Hz]")


sd = axs[1,1].scatter(df['x'], df['y'], c=(df['Tej']), s=ssize,
            cmap=cmap, label=r"$t_{ej}$ [s]")
axs[1,1].plot(tv/tj, 0, "+r")
axs[1,1].set_aspect('equal', adjustable='box')
fig.colorbar(sd, ax=axs[1,1])
axs[1,1].axis('off')
axs[1,1].legend()
plt.savefig(f"{PLOT_DIR}/{FILE_NAME}.png", bbox_inches='tight')
