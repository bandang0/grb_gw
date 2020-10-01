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


DATA_FILE = "data/circle_packings.data"
PLOT_DIR = f"plots/flares"

plt.style.use('presentation')

df = pd.read_csv(DATA_FILE, sep=" ",
        names=["r", "n"])

fit = np.polyfit(log(df['r']), log(df['n']), 1)
print(fit)
plt.figure()
plt.plot(df['r'], df['n'], label="data")
plt.plot(df['r'], 10**(fit[0] * log(df['r']) + fit[1]), label="fit")
plt.xlabel(r"$r$")
plt.ylabel(r"$N$")
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.savefig(f"{PLOT_DIR}/packings.png", bbox_inches='tight')
