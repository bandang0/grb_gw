'''Make plots for populations.'''

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
GWH =143. * 1.58
AH = float(argv[2])
rs = float(argv[3])
g0 = 100.

# Read event data and create pandas data frame.
print(f"Collecting data from {DATA_FILE}")
df = pd.read_csv(DATA_FILE, sep=" ",
        names=['d', 'n', 'e0', 'eb',
        'tv', 'tj', 'pt', 'pf', 'na', 'ni', '1', 'b', 'c'])
df = df.sort_values(by=['d'], ascending=[True])
NTOT = len(df)

# Add other columns to frame
df['ig'] = (1. + 6 * cos(df['tv']) ** 2 + cos(df['tv']) ** 4 \
        > 8 * df['d'] ** 2 / AH ** 2)
df['iv'] = (df['pf'] > rs)
df['igv'] = (df['ig'] & df['iv'])

# Collect event data according to detections in GW, AG, etc.
print("Creating derived dataframes.")
gw_radio = df.loc[df['igv']]
radio = df.loc[df['iv']]
gw = df.loc[df['ig']]

print(f"AH: {AH}, rs: {rs}")
print(f"Joint/GW: {float(len(gw_radio))/len(gw)}, Radio/Joint: {float(len(radio))/len(gw_radio)}, GW/Tot: {float(len(gw))/len(df)}, All: {len(df)}")
exit(0)

# Titles for plots
titles = {'d': 'Distance/Mpc',
        'ln': 'log(n/cm-3)',
        'n': 'n/cm-3',
        'le0': 'log(E0/erg)',
        'e0': 'E0/erg',
        'g': 'peak gamma',
        'tj': 'jet angle/rad',
        'tv': 'viewing angle/rad',
        'lpt': 'log(time of peak/day)',
        'lpf': 'log(peak flux/muJy',
        'pf' : 'peak flux / muJy',
        'ba': 'apparent beta at peak',
        'laas': 'log(apparent angular motion/ (mas/day))',
        'aas': 'proper motion (mas/d)',
        'ltod': 'log(obs dec time/day)',
        'lni': 'log(injection frequency/Hz)',
        'lnc': 'log(synchrotron cut frequency/Hz)',
        'lna': 'log(doppler-shifted obs freq / Hz)'}

for i, lab in enumerate(['d','lpf', 'ln', 'g', 'le0', 'tj', 'tv',
        'lpt', 'lpf', 'lna', 'lni', 'lnc']):
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

exit(0)
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

print("Fractions.")
plt.figure()
plt.plot(df['d'], gw_vlan / c, label='GW+VLA')
plt.plot(df['d'], gwn / c, label='GW')
plt.plot(df['d'], vlan / c, label='VLA')
plt.xlabel('D/Mpc')
plt.ylabel('fraction of all events < D detected')
plt.legend()
plt.savefig(f"{PLOT_DIR}/fractions.png")

for i, lab in enumerate(['d','pf', 'laas', 'ln', 'gp', 'le0', 'ltod', 'tj', 'tv',
        'lpt', 'lpf', 'dln', 'ba', 'aas',
        'lni', 'lnc', 'lna']):
    print(f"{lab}")
    plt.figure()
    plt.hist([df[lab], gw[lab], vla[lab], gw_vla[lab]],
            label=['All', 'GW', 'VLA', 'GW+VLA'],
            histtype='step', bins=BINS, density=True)
    plt.xlabel(f"{lab}")
    plt.ylabel(f"dN/d{lab}")
    if lab == 'gp':
        plt.xlim([0, 11])
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

print("dlndd")
plt.figure()
plt.axes = gw_vla.plot.hexbin('ln', 'd', gridsize=BINS, cmap='magma')
plt.xlabel('ln')
plt.ylabel('d')
plt.title(DATA_FILE+ " " + "GW + VLA")
plt.savefig(f"{PLOT_DIR}/dlndd.png", bbox_inches='tight')


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
