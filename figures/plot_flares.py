'''Make plots for populations.'''

from sys import exit
from os.path import basename, splitext
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d
from numpy import sin, cos
import numpy as np
from helpers import *
from sys import argv

DATA_FILE = argv[1]
plt.style.use('presentation')
NAME = f"{basename(splitext(argv[1])[0])}"
PLOT_DIR = f"plots/flares"
BINS = 100

def simplateau(t):
    if t < 10**2:
        return 1.e48 * (t / 10**2)**(-3)
    elif t < 10 ** (3.5):
        return 1.e48
    else:
        return 1.e48 * (t / 10 ** (3.5))**(-1)

Dl = 1000. * Mpc
# Read event data and create pandas data frame.
print(f"Collecting data from {DATA_FILE}")
df = pd.read_csv(DATA_FILE, sep=" ",
        names=['t', 'fx', 'fg'])
df['F'] = 1.e-25 * df['fx'] * (30 - 0.3) * keV / (hPlanck)

#dfp = pd.read_csv("data/plateau_011.data", sep=" ",
#dfp['fp'] = dfp['fp'] * 1.e26
#plateau = interp1d(log(dfp['t']), log(dfp['fp']))

T = df['t']
plateau_flux = [simplateau(x) for x in T]


plt.figure()
plt.plot(df['t'], df['F'], '+-', label="flare")
plt.plot(df['t'], plateau_flux, '+-', label="plateau")
plt.plot(T, plateau_flux  + df['F'], '+-', label="both")
#plt.xlim([10, 1.e4])
#plt.ylim([1.e16, 1.e18])
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$t_{obs}$ [s]")
plt.ylabel(r"$L_{0.3-30{\rm keV}, iso}$ [erg/s/cm2]")
plt.legend()
plt.savefig(f"{PLOT_DIR}/{NAME}.pdf", bbox_inches = 'tight')
