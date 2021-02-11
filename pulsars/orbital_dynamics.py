from sys import exit
import matplotlib.pyplot as plt
import pandas as pd
from numpy import sin, cos
import numpy as np
from scipy.interpolate import interp1d
from helpers import *
from sys import argv

DATA_FILE = argv[1]
PLOT_DIR = "plots/orbital_dynamics"
plt.style.use('palmerio')
BINS = 50

df = pd.read_csv(DATA_FILE, sep = " ", names =['t', 'a', 'i'])

plt.figure()
plt.plot(df['t'], df['a'])
plt.xscale('log')
#plt.yscale('log')
plt.xlabel('t / yr')
plt.ylabel('a / au')
plt.show()
