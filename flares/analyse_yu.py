from os import system
from sys import exit
import pickle
import numpy as np
import matplotlib.pyplot as plt
import emcee
import sys
from sys import argv
import pandas as pd
from helpers import *
import corner
from multiprocessing import Pool
from scipy.stats import kstest, norm
from astropy.cosmology import FlatLambdaCDM

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

df = pd.read_csv(f"data/yuflares.data", sep=" ", names=['n', 'z', 'ts', 'ts_', 'tp', 'tp_', 'te', 'te_', 'fp', 'fp_', 'sf', 'sf_', 'lp', 'lp_', 'e', 'e_'])

df['tpz'] = df['tp']/(1 + df['z'])
print(f"Median arrival time: {df.tpz.median()} s")
print(f"Median peak XRT luminosity: {df.lp.median() * 1.e48} erg/s")

print(f"Mean arrival time: {df.tpz.mean()} s")
print(f"Mean peak XRT luminosity: {df.lp.mean() * 1.e48} erg/s")

df_early = df.loc[df['tp'] < 1000]

print(f"Median arrival time: {df_early.tpz.median()} s (early flares)")
print(f"Median peak XRT luminosity: {df_early.lp.median() * 1.e48} erg/s (early flares)")

print(f"Mean arrival time: {df_early.tpz.mean()} s (early flares)")
print(f"Mean peak XRT luminosity: {df_early.lp.mean() * 1.e48} erg/s (early flares)")

plt.figure()
plt.hist(log(df.tpz), label="Arrival time (log)", bins=20)
plt.xlabel("log(TP)")
plt.show()
