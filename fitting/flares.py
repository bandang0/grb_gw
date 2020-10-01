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

r = float(argv[1])
tv = float(argv[2])
Gamma = float(argv[3])
DGamma = float(argv[4])
Dtej = float(argv[5])
tj = 0.1
N_patches = int(10 ** (-2.05273509 * log(r/tj) - 0.14390698))

PLOT_DIR = f"plots/flares"
CIRCLES_FILE = f"data/cci_coords/cci{N_patches}.txt"
FILE_NAME = f"{r}-{tv}-{Gamma}-{DGamma}-{Dtej}"

print(f"r: {r}, tv: {tv}, Gamma: {Gamma}, DGamma: {DGamma}, Dtej: {Dtej}, N_p: {N_patches}")
plt.style.use('presentation')
BINS = 70

cdict = {'red':   ((0.0,  0.0, 0.0),
                   (1.,  1.0, 1.0),
                   (1.0,  1., 1.0)),

         'green': ((0.0,  0., 0.0),
                   (0.,  .0, .0),
                   (1.0,  1.00, 1.0)),

         'blue':  ((0.0,  0.72, 0.0),
                   (0.5,  0.0, 0.0),
                   (1.0,  0.11, 1.0))}

cmap = colors.LinearSegmentedColormap.from_list('custom',rainbow(np.linspace(0, 1, 10)))
np.random.seed(52)

# physics
tvar = 1.
Eiso = 1.e50
a = 0.
b = 1.5
Ep = 1. * keV
nup = Ep / hPlanck
nuobs = 5.0 * keV / hPlanck
nug = 300. * keV / hPlanck
s = lambda nu: ((nu / nup) ** (a if (nu < nup) else (-b)))/(nup * (1/(a+1) + 1/(b-1)))
def Fnu(chic, r, te, re, nuobs, t, Gamma, beta):
    return ((Eiso / 4 * Pi) * (re / cLight) * deltaphi(chic, r, np.arccos((te - t) * cLight / re))
        * s(nuobs * Gamma * (1 - beta * cLight * (te - t) / re))
        / (Gamma * (1 - beta * cLight * (te - t) / re)) ** 2)


df = pd.read_csv(CIRCLES_FILE,
        names=["i", "x", "y"])


df['theta'] = tj * sqrt(df['x'] ** 2 + df['y'] ** 2)
df['phi'] = np.arctan2(df['y'], df['x'])
df['chic'] = np.arccos(cos(tv) * cos(df['theta']) + sin(tv) * sin(df['theta']) * cos(df['phi']))
df['chim'] = np.maximum(df['chic'] - r, 0.)
df['chiM'] = np.minimum(df['chic'] + r, Pi)
df['Gamma'] = Gamma + DGamma * random.uniform(0, 1, df.shape[0])
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
print(tmmin, tMmax)

plt.figure()
x = np.linspace(-1, 1, 50)
thvl = [0.1, 0.15, 0.2, 0.25, 0.3]
for thv in thvl:
    y = [2 * sin(thv + a * tj) * 100  ** 2 * max(0.5/100,
           sin(thv + a * tj)) * np.arccos((cos(tj) -
           cos(thv) * cos(thv + a * tj))/(sin(thv) * sin(thv + a * tj)))
           for a in x]
    plt.plot(x, y, color=cmap((thv-0.1)/(0.3-0.1)), label=r"$\theta_v$" + f" = {thv}")
plt.xlabel(r"$(\theta_v - \chi_0)/\theta_j$")
plt.ylabel(r"$\Delta^{\rm min} t_{ej}$ [s]")
plt.legend()
plt.savefig(f"{PLOT_DIR}/success100.png", bbox_inches='tight')

plt.figure()
for thv in thvl:
    y = [2 * sin(thv + a * tj) * 400  ** 2 * max(0.5/400,
           sin(thv + a * tj)) * np.arccos((cos(tj) -
           cos(thv) * cos(thv + a * tj))/(sin(thv) * sin(thv + a * tj)))
           for a in x]
    plt.plot(x, y, "--", color=cmap((thv-0.1)/(0.3-0.1)), label=r"$\theta_v$" + f" = {thv}")
plt.xlabel(r"$(\theta_v - \chi_0)/\theta_j$")
plt.ylabel(r"$\Delta^{\rm min} t_{ej}$ [s]")
plt.savefig(f"{PLOT_DIR}/success400.png", bbox_inches='tight')

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


plt.figure()
plt.scatter(df['x'], df['y'], c=(df['Tej']), s=100 * (r / 0.01) ** 2,
            cmap=cmap, label="Ejection time [s]")
plt.colorbar()
plt.axis('off')
plt.legend()
plt.savefig(f"{PLOT_DIR}/{FILE_NAME}_ejtime.png", bbox_inches='tight')

plt.figure()
plt.scatter(df['x'], df['y'], c=(df['tm']), s=100 * (r / 0.01) ** 2,
            cmap=cmap, label=r"$t_m$ [s]")
plt.plot(tv/tj, 0, "+r")
plt.gca().set_aspect('equal', adjustable='box')
plt.colorbar()
plt.axis('off')
plt.legend()
plt.savefig(f"{PLOT_DIR}/{FILE_NAME}_tm.png", bbox_inches='tight')

plt.figure()
for i, p in df.iterrows():
    te = p['Te']
    chic = p['chic']
    re = p['Re']
    tm = p['tm']
    tM = p['tM']
    g = p['Gamma']
    b = p['beta']
    x = np.linspace(tm, tM, 300)
    y = np.array([Fnu(chic, r, te, re, nuobs, a, g, b) for a in x])
    plt.plot(x, (10 - 0.3) * keV * y / hPlanck, color = cmap((chic - chicmin)/(chicmax - chicmin)))
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$t_{\rm obs}$ [s]")
plt.ylabel(r"$L^{\rm iso}$ [erg/s]")
#plt.ylim([1.e9, 1.e20])
plt.savefig(f"{PLOT_DIR}/{FILE_NAME}_individual.png", bbox_inches='tight')

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
    y = np.array([(Fnu(chic, r, te, re, nuobs, a, g, b) if ((a < tM) and (tm < a)) else 0.) for a in T])
    Fx = Fx + y
plt.figure()
plt.plot(T, (10 - 0.3) * keV * Fx / hPlanck)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$t_{\rm obs}$ [s]")
plt.ylabel(r"$L^{\rm iso}$ [erg/s]")
#plt.ylim([1.e10, 1.e20])
plt.savefig(f"{PLOT_DIR}/{FILE_NAME}_summed.png", bbox_inches='tight')
