'''Likelihood-free inference for BNSM population.'''

from sys import exit
from os.path import basename, splitext
import matplotlib.pyplot as plt
import pandas as pd
from numpy import sin, cos
import numpy as np
from helpers import *
from sys import argv

DATA_FILE = argv[1]
PLOT_DIR = f"plots/{basename(splitext(argv[1])[0])}"
BINS = 100 

# Constants
GWH = 500. * 1.58
AH = GWH / 1.58
vlas = 15.
skas = 3.
rad10s = 10.
ngvlas = 1.
g0 = 100.

# Measurements
data = {'TP': [164],
        'FP': [120],
        'D': [41],
        'BA': [4.1]}

# Errors
errors = {'TP': [7],
        'FP': [9],
        'D': [3],
        'BA': [0.5]}

# Horizons and sensitivities
AH = [210]
S = [10]

# Choices
espilon = 10.
N = 100

# Read event data and create pandas data frame.
print(f"Collecting data from {DATA_FILE}")
df = pd.read_csv(DATA_FILE, sep=" ",
        names=['ti', 'a1', 'a2', 'es', 'ln0', 'leb0', 'd', 'tv', 'pt', 'pf',
        'ba'])
NTOT = len(df)

for i in range(len(AH)):
    tmp_df = df.loc[df['pf'] > S[i] \
            && 1. + 6 * cos(df['tv']) ** 2 + cos(df['tv']) ** 4\
            > 8 * df['d'] ** 2 / AH[i] ** 2]
    data_real = dict([(tag, np.random.normal(data[tag][i], errors[tag][i], N))\
            for tag in data])
    for j in range(N):
        loc_df = tmp_df.loc[(tmp_df['d'] - data_real['D'][j]) ** 2 \
                        +   (tmp_df['pt'] - data_real['TP'][j]) ** 2\
                        +   (tmp_df['pf'] - data_real['FP'][j]) ** 2\
                        +   (tmp_df['ba'] - data_real['BA'][j]) ** 2 > epsilon]

    
# Add other columns to frame
df['lpt'] = log(df['pt'])
df['ln'] = log(df['n'])
df['lni'] = log(df['ni'])
df['lnc'] = log(df['nc'])
df['lpf'] = log(df['pf'])
df['le0'] = log(df['e0'])
df['be'] = df['ba'] / (sin(df['tv']) + cos(df['tv']) * df['ba'])
df['gp'] = 1. / np.sqrt(1 - df['be'] ** 2)
df['ld'] = log(df['d'])
df['ig'] = 1. + 6 * cos(df['tv']) ** 2 + cos(df['tv']) ** 4 \
        > 8 * df['d'] ** 2 / AH ** 2
df['iv'] = (df['pf'] > vlas)
df['is'] = (df['pf'] > skas)
df['ir'] = (df['pf'] > rad10s)
df['in'] = (df['pf'] > ngvlas)
df['igv'] = (df['ig'] & df['iv'])
df['tod'] = (3 * df['e0'] / (4 * Pi * df['n'] * mProton * g0 ** 2 \
        * cLight ** 2)) ** (1. / 3.) / (2. * g0 * cLight ** 2) / day
df['ltod'] = log(df['tod'])
df['lpt-ltod'] = df['lpt'] - df['ltod']
df['aas'] = 86400 * 1000 * 3600 * 180 * cLight * df['ba'] \
        / (Pi * 1.e6 * 3.086e18 * df['d'])
df['laas'] = log(df['aas'])
df['lna'] = log(3.e9 * (1 - df['be'] * cos(df['tv'])) / (1 - df['be']))
df['dln'] = df['lna'] - df['lni']

# Collect event data according to detections in GW, AG, etc.
print("Creating derived dataframes.")
gw_vla = df.loc[df['igv']]
vla = df.loc[df['iv']]
gw = df.loc[df['ig']]
ska = df.loc[df['is']]
rad = df.loc[df['ir']]

print(f"GW+VLA: {len(gw_vla)}, VLA: {len(vla)}, GW: {len(gw)}, All: {len(df)}")

# Fractions as function of horizon
print("Functions of horizon.")
n_gw_vla = list()
n_d = list()
theta_v_bar = list()
theta_v_std = list()
gwh_l = np.linspace(10, GWH, 20)
for gwh in gwh_l:
    ah = gwh / 1.58
    tmp_df = vla.loc[1. + 6 * cos(df['tv']) ** 2 + cos(df['tv']) ** 4\
            > 8 * df['d'] ** 2 / ah ** 2]
    theta_v_bar.append(tmp_df['tv'].mean())
    theta_v_std.append(np.std(tmp_df['tv']))
    n_gw_vla.append(len(tmp_df))
    n_d.append(len(df.loc[df['d'] < gwh]))

n_gw_vla = np.array(n_gw_vla)
n_d = np.array(n_d)

print("-- theta_obs")
plt.figure()
plt.errorbar(gwh_l, theta_v_bar, yerr=theta_v_std)
plt.xlabel("Horizon/Mpc")
plt.ylabel("mean theta obs")
plt.title(DATA_FILE + " " + "tv")
plt.savefig(f"{PLOT_DIR}/tv_horizon.png")

print("-- Number")
plt.figure()
plt.plot(log(gwh_l), log(n_gw_vla))
plt.xlabel("log(Horizon/Mpc)")
plt.ylabel("log(#(GW+VLA))")
plt.title(DATA_FILE + " " + "numbers")
plt.savefig(f"{PLOT_DIR}/numbers_horizon.png")

print("-- Fraction")
plt.figure()
plt.plot(gwh_l, n_gw_vla/n_d)
plt.xlabel("Horizon/Mpc")
plt.ylabel("#(GW+VLA)/#tot [= alpha * #(GW+VLA)/#GW]")
plt.title(DATA_FILE + " " + "fractions")
plt.savefig(F"{PLOT_DIR}/fractions_horizon.png")

# Collect fractions
print("Calculating fractions")
c = np.array([float(i) + 1. for i in range(len(df['d']))])
gw_vlan = np.cumsum(df['igv'])
vlan = np.cumsum(df['iv'])
gwn = np.cumsum(df['ig'])

# Titles for plots
titles = {'d': 'Distance/Mpc',
        'ln': 'log(n/cm-3)',
        'n': 'n/cm-3',
        'le0': 'log(E0/erg)',
        'e0': 'E0/erg',
        'gp': 'peak gamma',
        'tj': 'jet angle/rad',
        'tv': 'viewing angle/rad',
        'dln': 'log(nuobs/a) - log(ni)',
        'lpt': 'log(time of peak/day)',
        'dlpt': 'pl time - time',
        'dlpf': 'pl flux - flux',
        'lppt': 'log(pl time of peak/day',
        'lppf': 'log(pl peak flux/muJy',
        'lpf': 'log(peak flux/muJy',
        'ba': 'apparent beta at peak',
        'laas': 'log(apparent angular motion/ (mas/day))',
        'aas': 'proper motion (mas/d)',
        'ltod': 'log(obs dec time/day)',
        'lni': 'log(injection frequency/Hz)',
        'lnc': 'log(synchrotron cut frequency/Hz)',
        'lna': 'log(doppler-shifted obs freq / Hz)'}

print("Fractions.")
plt.figure()
plt.plot(df['d'], gw_vlan / c, label='GW+VLA')
plt.plot(df['d'], gwn / c, label='GW')
plt.plot(df['d'], vlan / c, label='VLA')
plt.xlabel('D/Mpc')
plt.ylabel('fraction of all events < D detected')
plt.legend()
plt.savefig(f"{PLOT_DIR}/fractions.png")

for i, lab in enumerate(['d','laas', 'ln', 'gp', 'le0', 'ltod', 'tj', 'tv',
        'lpt', 'lpf', 'dln', 'ba', 'aas',
        'lni', 'lnc', 'lna']):
    print(f"{lab}")
    plt.figure()
    plt.hist([df[lab], gw[lab], vla[lab], gw_vla[lab]],
            label=['All', 'GW', 'VLA', 'GW+VLA'],
            histtype='step', bins=BINS, density=True)
    plt.xlabel(f"{lab}")
    plt.ylabel(f"dN/d{lab}")
    plt.title(DATA_FILE+ " " + titles[lab])
    plt.legend()
    plt.savefig(f"{PLOT_DIR}/d{lab}.png", bbox_inches='tight')
    plt.close()

print("ld")
plt.figure()
plt.hist([df['ld'], gw['ld'], vla['ld'], gw_vla['ld']],
        label=['All', 'GW', 'VLA', 'GW+VLA'],
        histtype='step', bins=BINS, cumulative=True,
        log=True)
plt.xlabel("log(D/Mpc)")
plt.ylabel("N(<D)")
plt.legend()
plt.savefig(f"{PLOT_DIR}/ld.png", bbox_inches='tight')

print("dltoddga")
plt.figure()
plt.axes = gw_vla.plot.hexbin('lpt-ltod', 'gp', gridsize=BINS, cmap='magma')
plt.xlabel("log(pt/tod)")
plt.ylabel("gp")
plt.title(DATA_FILE+ " " + "GW+VLA")
plt.savefig(f"{PLOT_DIR}/dlpt-ltoddgp.png", bbox_inches='tight')

print("dlndgp")
plt.figure()
plt.axes = gw_vla.plot.hexbin('ln', 'gp', gridsize=BINS, cmap='magma')
plt.xlabel('ln')
plt.ylabel('gp')
plt.title(DATA_FILE+ " " + "GW + VLA")
plt.savefig(f"{PLOT_DIR}/dlndgp.png", bbox_inches='tight')

print("dtvdgp")
plt.figure()
plt.axes = gw_vla.plot.hexbin('tv', 'gp', gridsize=BINS, cmap='magma')
plt.xlabel('tv')
plt.ylabel('gp')
plt.title(DATA_FILE+ " " + "GW + VLA")
plt.savefig(f"{PLOT_DIR}/dtvdgp.png", bbox_inches='tight')

print("dtvdd")
plt.figure()
plt.axes = gw_vla.plot.hexbin('d', 'tv', gridsize=BINS,cmap='magma')
plt.title(DATA_FILE+ " " + "GW + VLA")
plt.savefig(f"{PLOT_DIR}/dtvdd_cp.png", bbox_inches='tight')
plt.close('all')

plt.figure()
plt.axes = gw.plot.hexbin('d', 'tv', gridsize=BINS,cmap='magma')
plt.title(DATA_FILE+ " " + 'GW')
plt.savefig(f"{PLOT_DIR}/dtvdd_gw.png", bbox_inches='tight')
plt.close('all')

print("dtvdlpf")
plt.figure()
plt.axes = gw.plot.hexbin('lpf', 'tv', gridsize=BINS,cmap='magma')
plt.title(DATA_FILE+ " " + 'GW')
plt.savefig(f"{PLOT_DIR}/dlpfdtv_gw.png", bbox_inches='tight')
plt.close('all')

print("dddlpf")
plt.figure()
plt.axes = gw.plot.hexbin('d', 'lpf', gridsize=BINS,cmap='magma')
plt.title(DATA_FILE+ " " + 'GW')
plt.savefig(f"{PLOT_DIR}/dddlpf_gw.png", bbox_inches='tight')
plt.close('all')

print("dlptldpf")
plt.figure()
plt.axes = gw.plot.hexbin('lpt', 'lpf', gridsize=BINS,cmap='magma')
plt.title(DATA_FILE+ " " + 'GW')
plt.savefig(f"{PLOT_DIR}/dlptdlpf_gw.png", bbox_inches='tight')
plt.close('all')

print("dlndpf")
plt.figure()
plt.axes = gw.plot.hexbin('ln', 'lpf', gridsize=BINS,cmap='magma')
plt.title(DATA_FILE+ " " + 'GW')
plt.savefig(f"{PLOT_DIR}/dlndlpf_gw.png", bbox_inches='tight')
plt.close('all')

print("dtvdaas")
plt.figure()
plt.axes = gw_vla.loc[gw_vla['aas'] < 0.02].plot.hexbin('tv', 'aas', gridsize=BINS,cmap='jet')
plt.title(DATA_FILE+ " " + "GW+VLA")
plt.savefig(f"{PLOT_DIR}/dtvdaas_cp.png", bbox_inches='tight')
plt.close('all')

print("dddass")
plt.figure()
plt.axes = gw_vla.loc[gw_vla['aas'] < 0.02].plot.hexbin('d', 'aas', gridsize=BINS,cmap='jet')
plt.title(DATA_FILE+ " " + "GW+VLA")
plt.savefig(f"{PLOT_DIR}/dddaas_cp.png", bbox_inches='tight')
plt.close('all')
print(max(gw_vla['be']), min(gw_vla['be']))

print("dgpdlpt")
plt.figure()
plt.axes = gw_vla.plot.hexbin('gp', 'lpt', gridsize=BINS,cmap='jet')
plt.title(DATA_FILE+ " " + "GW+VLA")
plt.savefig(f"{PLOT_DIR}/dgpdlpt_cp.png", bbox_inches='tight')
plt.close('all')

print("dtvdlpt")
plt.figure()
plt.axes = gw_vla.plot.hexbin('tv', 'lpt', gridsize=BINS,cmap='jet')
plt.title(DATA_FILE+ " " + "GW+VLA")
plt.savefig(f"{PLOT_DIR}/dtvdlpt_cp.png", bbox_inches='tight')
plt.close('all')

gw_vla['t'] = gw_vla[['tj', 'tv']].max(axis=1)
print("dgpdt")
plt.figure()
plt.axes = gw_vla.plot.hexbin('gp', 't', gridsize=BINS,cmap='jet')
plt.title(DATA_FILE+ " " + "GW+VLA")
plt.savefig(f"{PLOT_DIR}/dgpdt_cp.png", bbox_inches='tight')
plt.close('all')


