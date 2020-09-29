from sys import exit
import matplotlib.pyplot as plt
import pandas as pd
from numpy import sin, cos, log10
import numpy as np
from helpers import *
from sys import argv
import corner
import pickle
import emcee
from itertools import product

dLmax = argv[1]
energy = argv[2]

DATA_FILE = f"data/ecm_{dLmax}_{energy}.data"
PLOT_DIR = f"plots/hubble_{dLmax}_{energy}"
OUTPUT = "output"
plt.style.use('presentation')
BINS = 300

# Thresholds
MAGLIM_ZTF = 21.
MAGLIM_LSST = 24.5
VLAS = 15. * microJy
TJ = 0.1
HSAS = 2 * mas
adec = -2.2
ainc = 0.8
dlogtm = 9.e-2

df = pd.read_csv(DATA_FILE, sep=" ",
        names=["d", "n", "eiso", "eb", "tv",
                "pt", "pf", "ppm", "dtvlbi",
                "kng", "knr"])


#df['np'] = (1/adec - 1/ainc) * log(VLAS / df['pf'])/dlogtm
df['le'] = log10(df['eiso'])
df['ci'] = cos(df['tv'])
df['dtheta'] = df['dtvlbi'] * df['ppm']
df['kn'] = ((df['kng'] < MAGLIM_ZTF) | (df['knr'] < MAGLIM_ZTF))
df['age'] = (df['pf'] > 10. * VLAS)
df['vlbi'] = (df['dtheta'] > HSAS)
df['vlbi0'] = (df['ppm'] > 3. * muas/day)


df1 = df.loc[df['kn']]
df2 = df.loc[df['kn'] & df['age']]
df3 = df.loc[df['kn'] & df['age'] & df['vlbi']]
dfp = df2.loc[df2['ppm'] > 0.]
dfd = df2.loc[df2['dtheta'] > 0.]

print(f"Data read from: {DATA_FILE}")

print(f"Total events in sample: {len(df['d'])}")
print(f"Total events in first distance bin: {len(df['d'])/BINS}")
print(f"Total events in first distance-cosi bin: {len(df['d'])/BINS**2}")
print(f"Average log(Eiso) in sample: {df['le'].mean()}")
print(f"Largest distance for all:   {max(df['d']) / Mpc}")
print(f"Largest distance for KN:    {max(df1['d']) / Mpc}")
print(f"Largest distance for AG(e): {max(df2['d']) / Mpc}")
print(f"Largest distance for VLBI:  {max(df3['d']) / Mpc}")

print(f"Distance histogram")
plt.figure()
plt.hist([df1['d'] / Mpc, df2['d'] / Mpc,
    df3['d'] / Mpc],
        label=['L1', 'L2', 'L3'],
        histtype='step', bins=BINS)
plt.xlabel(r"$D$")
plt.legend()
plt.savefig(f"{PLOT_DIR}/d.pdf", bbox_inches='tight')

print(f"cos iota histogram")
plt.figure()
plt.hist([df1['ci'], df2['ci'],
    df3['ci']],
        label=['L1', 'L2', 'L3'],
        histtype='step', bins=BINS)
plt.xlabel(r"$\cos \iota$")
plt.legend()
plt.savefig(f"{PLOT_DIR}/ci.pdf", bbox_inches='tight')

print(f"Magnitude histogram")
plt.figure()
plt.hist([df['knr'], df['kng']],
        label=['r', 'g'],
        histtype='step', bins=BINS)
plt.xlabel(r"$r, g$")
plt.legend()
plt.savefig(f"{PLOT_DIR}/m.pdf", bbox_inches='tight')

print(f"PPM histogram")
plt.figure()
plt.hist([log(dfp['ppm'] / (muas/day))],
        label=['PPM / muas/day'],
        histtype='step', bins=BINS)
plt.xlabel(r"$\log$ PPM / (muas/day)")
plt.legend()
plt.savefig(f"{PLOT_DIR}/ppm.pdf", bbox_inches='tight')

print(f"Dtheta histogram")
plt.figure()
plt.hist([log(dfd['dtheta'] / (mas))],
        label=['dtheta / mas'],
        histtype='step', bins=BINS)
plt.xlabel(r"$\log \Delta \theta$ / mas")
plt.legend()
plt.savefig(f"{PLOT_DIR}/dtheta.pdf", bbox_inches='tight')


print(f"Proper motion vs. total displacement scatter plot (log)")
plt.figure()
plt.scatter(dfd['ppm'] / (muas/day), dfd['dtvlbi'] * dfd['ppm'] / mas)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"PPM / (muas/day)")
plt.ylabel(r"$\Delta \theta$ / mas")
plt.savefig(f"{PLOT_DIR}/s.png", bbox_inches='tight')

print(f"Total displacement vs. cos(i) scatter plot")
plt.figure()
plt.scatter(dfd['dtvlbi'] * dfd['ppm'] / mas, dfd['ci'])
plt.xlabel(r"$\Delta \theta$ / mas")
plt.ylabel(r"$\cos \iota$")
plt.xscale('log')
plt.savefig(f"{PLOT_DIR}/t.png", bbox_inches='tight')

print(f"Total displacement vs. D scatter plot")
plt.figure()
plt.scatter(dfd['dtvlbi'] * dfd['ppm'] / mas, dfd['d']/ Mpc)
plt.xlabel(r"$\Delta \theta$ / mas")
plt.ylabel(r"$D$ / Mpc")
plt.xscale('log')
plt.savefig(f"{PLOT_DIR}/u.png", bbox_inches='tight')

labs = {"df1": "KN",
        "df2": "KN+AG(e)",
        "df3": "KN+AG(e)+VLBI",
        "df6": "AG(p)",
        "df7": "AG(e)",
        "df8": "VLBI",
        "df9": "KNAG",
        "df10": "GRB",
        "df11": "KN+AG(e)+VLBI0",}

# for i, a in enumerate([df1, df2, df3]):
#     print(labs[f"df{i+1}"] + f": {len(a):7d}, {100 * float(len(a)) / len(df):.2f} (%)")
#     plt.figure()
#     plt.hist2d(a['d'] / Mpc, a['ci'], cmap="jet", range=[[0, 530], [0, 1]], bins=BINS)
#     plt.xlabel(r"$D$ [Mpc]")
#     plt.ylabel(r"$\cos i$")
#     plt.colorbar()
#     plt.title(labs[f"df{i+1}"])
#     plt.savefig(f"{PLOT_DIR}/{dLmax}_{energy}_" + labs[f"df{i+1}"] + ".pdf", bbox_inches='tight')

# edges for VLBI pdet interpolation
dc = 1./BINS
dd = 530./BINS

dedges = np.array([0] + [dd/2 + dd * n for n in range(BINS)] + [530])
cedges = np.array([0] + [dc/2 + dc * n for n in range(BINS)] + [1])

dgrid = np.array([dd * n for n in range(BINS+1)])
cgrid = np.array([dc * n for n in range(BINS+1)])
Htot, _, _ = np.histogram2d(df['d'] / Mpc, df['ci'], bins=(dedges, cedges))
H1, _, _ = np.histogram2d(df1['d'] / Mpc, df1['ci'], bins=(dedges, cedges))
H2, _, _ = np.histogram2d(df2['d'] / Mpc, df2['ci'], bins=(dedges, cedges))
H3, _, _ = np.histogram2d(df3['d'] / Mpc, df3['ci'], bins=(dedges, cedges))

plt.figure()
X, Y = np.meshgrid(dedges, cedges)
plt.pcolormesh(X, Y, Htot.T, cmap="jet")
plt.xlabel(r"$D$ [Mpc]")
plt.ylabel(r"$\cos i$")
plt.colorbar()
plt.title("Number of events per bin (total)")
plt.savefig(f"{PLOT_DIR}/{dLmax}_{energy}_Htot.pdf", bbox_inches='tight')

plt.figure()
plt.pcolormesh(X, Y, H1.T, cmap="jet")
plt.xlabel(r"$D$ [Mpc]")
plt.ylabel(r"$\cos i$")
plt.colorbar()
plt.title(f"Number of events per bin (KN {energy})")
plt.savefig(f"{PLOT_DIR}/{dLmax}_{energy}_H1.pdf", bbox_inches='tight')

plt.figure()
plt.pcolormesh(X, Y, H2.T, cmap="jet")
plt.xlabel(r"$D$ [Mpc]")
plt.ylabel(r"$\cos i$")
plt.colorbar()
plt.title(f"Number of events per bin (AG {energy})")
plt.savefig(f"{PLOT_DIR}/{dLmax}_{energy}_H2.pdf", bbox_inches='tight')

plt.figure()
plt.pcolormesh(X, Y, H3.T, cmap="jet")
plt.xlabel(r"$D$ [Mpc]")
plt.ylabel(r"$\cos i$")
plt.colorbar()
plt.title(f"Number of events per bin (VLBI {energy})")
plt.savefig(f"{PLOT_DIR}/{dLmax}_{energy}_H3.pdf", bbox_inches='tight')

plt.figure()
plt.pcolormesh(X, Y, H1.T/Htot.T, cmap="jet")
plt.xlabel(r"$D$ [Mpc]")
plt.ylabel(r"$\cos i$")
plt.colorbar()
plt.title(f"Pdet (KN {energy})")
plt.savefig(f"{PLOT_DIR}/{dLmax}_{energy}_pdetKN.pdf", bbox_inches='tight')

plt.figure()
plt.pcolormesh(X, Y, H2.T/Htot.T, cmap="jet")
plt.xlabel(r"$D$ [Mpc]")
plt.ylabel(r"$\cos i$")
plt.colorbar()
plt.title(f"Pdet (AG {energy})")
plt.savefig(f"{PLOT_DIR}/{dLmax}_{energy}_pdetAG.pdf", bbox_inches='tight')

plt.figure()
plt.pcolormesh(X, Y, H3.T/Htot.T, cmap="jet")
plt.xlabel(r"$D$ [Mpc]")
plt.ylabel(r"$\cos i$")
plt.colorbar()
plt.title(f"Pdet (VLBI {energy})")
plt.savefig(f"{PLOT_DIR}/{dLmax}_{energy}_pdetVLBI.pdf", bbox_inches='tight')

np.savetxt(f"{output}/{dLmax}_{energy}_pdet_VLBI.txt", H3.T/Htot.T, '%8.6f')
np.savetxt(f"{output}/{dLmax}_{energy}_pdet_AG.txt", H2.T/Htot.T, '%8.6f')
np.savetxt(f"{output}/{dLmax}_{energy}_pdet_KN.txt", H1.T/Htot.T, '%8.6f')
np.savetxt(f"{output}/530_distance_grid.txt", dgrid, '%8.6f')
np.savetxt(f"{output}/cosi_grid.txt", cgrid, '%8.6f')
# with open(f"{dLmax}_{energy}_kn.pickle", 'wb') as pickle_file:
#     pickle.dump({"dL": np.array(df1['d']) / Mpc,
#                  "cosi": np.array(cos(df1['tv']))}, pickle_file)
#
# with open(f"{dLmax}_{energy}_age.pickle", 'wb') as pickle_file:
#     pickle.dump({"dL": np.array(df2['d']) / Mpc,
#                  "cosi": np.array(cos(df2['tv']))}, pickle_file)
#
# with open(f"{dLmax}_{energy}_vlbi.pickle", 'wb') as pickle_file:
#     pickle.dump({"dL": np.array(df3['d']) / Mpc,
#                  "cosi": np.array(cos(df3['tv']))}, pickle_file)
