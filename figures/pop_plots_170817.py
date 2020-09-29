'''Plot the histogram of distances.'''

from sys import exit
import matplotlib.pyplot as plt
import pandas as pd
from numpy import log10 as log
from numpy import sin, cos
import numpy as np

DATA_FILE = "170817_like_1000.data"
PLOT_DIR = "plots/170817_like"
BINS = 20

GWH = 385.
vlas = 15.
skas = 3.
ngvlas = 1.
def gwd(d, tv):
    if GWH * np.sqrt(0.125 * (1. + 6. * cos(tv) ** 2 + cos(tv) ** 4)) > d:
        return 1
    else:
        return 0

# Read event data and create pandas data frame.
print(f"Collecting data from {DATA_FILE}")
all_events = pd.read_csv(DATA_FILE, sep=" ",
        names=['ti', 'd', 'n', 'e0', 'eb', 'tv', 'tj', 'pt', 'pf', 'ba'])
all_events = all_events.sort_values(by=['d'], ascending=[True])

# Add other columns to frame
all_events = all_events.assign(lpt=log(all_events['pt']))
all_events = all_events.assign(ln=log(all_events['n']))
all_events = all_events.assign(lpf=log(all_events['pf']))
all_events = all_events.assign(be=all_events['ba'] / (sin(all_events['tv']) + cos(all_events['tv']) * all_events['ba']))
all_events = all_events.assign(ga=1./np.sqrt(1. - all_events['be'] ** 2))
all_events = all_events.assign(ld=log(all_events['d']))
all_events = all_events.assign(ig=gwd(all_events['d'], all_events['tv']))
all_events = all_events.assign(iv = 1 if all_events['pf'] > vlas else 0)
all_events = all_events.assign(igv = all_events['ig'] * all_events['iv'])

# Collect event data according to detections in GW, AG, etc.
print("Creating derived dataframes.")
gw_vla = all_events.loc[all_events['igv'] == 1]
vla = all_events.loc[all_events['iv'] == 1]
gw = all_events.loc[all_events['ig'] == 1]
#ska = all_events.loc[all_events['is'] == 1]
ar = gw_vla.loc[np.abs(gw_vla['tv'] - 0.1) < 0.05] 

# Collect fractions
print("Calculating fractions")
distances = list(all_events['d'])
distances.sort()
print(len(distances))
c = np.array([float(i) + 1. for i in range(len(distances))])
gw_vlan = np.cumsum(all_events['igv'])
vlan = np.cumsum(all_events['iv'])
gwn = np.cumsum(all_events['ig'])

print("Fractions.")
plt.figure()
plt.plot(distances, gw_vlan / c, label='GW+VLA')
plt.plot(distances, gwn / c, label='GW')
plt.plot(distances, vlan / c, label='VLA')
plt.xlabel('D/Mpc')
plt.ylabel('fraction')
plt.legend()
plt.savefig(f"{PLOT_DIR}/fractions.png")

print(len(gw_vla), len(vla), len(gw), len(all_events))
for lab in ['d', 'ln', 'tv', 'lpt', 'lpf', 'ba', 'ga']:
    print(f"{lab}")
    plt.figure()
    plt.hist([all_events[lab], gw[lab], vla[lab], ska[lab], gw_vla[lab]],
            label=['All', 'GW', 'VLA', 'SKA', 'GW+VLA'],
            histtype='step', bins=BINS, density=True)
    plt.xlabel(f"{lab}")
    plt.ylabel(f"dN/d{lab}")
    plt.legend()
    plt.savefig(f"{PLOT_DIR}/d{lab}.png", bbox_inches='tight')

print("ld")
plt.figure()
plt.hist([all_events['ld'], gw['ld'], vla['ld'], gw_vla['ld']],
        label=['All', 'GW', 'VLA', 'GW+VLA'],
        histtype='step', bins=BINS, cumulative=True,
        log=True)
plt.xlabel("log(D/Mpc)")
plt.ylabel("N(<D)")
plt.legend()
plt.savefig(f"{PLOT_DIR}/ld.png", bbox_inches='tight')

print("dtvdd")
plt.figure()
plt.axes = gw_vla.plot.hexbin('d', 'tv', gridsize=BINS)
plt.savefig(f"{PLOT_DIR}/dtvddcp.png", bbox_inches='tight')
plt.figure()
plt.axes = gw.plot.hexbin('d', 'tv', gridsize=BINS)
plt.savefig(f"{PLOT_DIR}/dtvddgw.png", bbox_inches='tight')
plt.figure()
plt.axex = vla.plot.hexbin('d', 'tv', gridsize=BINS)
plt.savefig(f"{PLOT_DIR}/dtvddvla.png", bbox_inches='tight')

print("dptdpf")
plt.figure()
plt.axes = all_events.plot.hexbin('lpt', 'lpf', gridsize=BINS)
plt.savefig(f"{PLOT_DIR}/dtpdfp.png", bbox_inches='tight')
plt.figure()
plt.axes = gw.plot.hexbin('lpt', 'lpf', gridsize=BINS)
plt.savefig(f"{PLOT_DIR}/dtpdfpgw.png", bbox_inches='tight')
plt.figure()
plt.axex = vla.plot.hexbin('lpt', 'lpf', gridsize=BINS)
plt.savefig(f"{PLOT_DIR}/dtpdfpvla.png", bbox_inches='tight')

print("dndpf")
plt.figure()
plt.axes = all_events.plot.hexbin('ln', 'lpf', gridsize=BINS)
plt.savefig(f"{PLOT_DIR}/dndpf.png", bbox_inches='tight')
plt.figure()
plt.axes = gw.plot.hexbin('ln', 'lpf', gridsize=BINS)
plt.savefig(f"{PLOT_DIR}/dndpfgw.png", bbox_inches='tight')
plt.figure()
plt.axex = vla.plot.hexbin('ln', 'lpf', gridsize=BINS)
plt.savefig(f"{PLOT_DIR}/dndpfvla.png", bbox_inches='tight')

exit(0)

print("Peak times (cumulative).")
plt.figure()
plt.hist([log(vla['pt']), log(gw['pt']), log(all_events['pt'])],
    label=['Detection Time (AG counterpart to GW)', 'Peak Time (all GW)',
    'Peak Time (all events)'], cumulative=True, histtype='step', bins=BINS)
plt.xlabel("log(time/days)")
plt.ylabel("#")
plt.legend()
plt.savefig(f"{PLOT_DIR}/times.png", bbox_inches='tight')

print("Peak times (cumulative), counterpart only.")
plt.figure()
plt.hist(log(vla['pt']), label='Detection Time (AG counterpart to GW)',
    cumulative=True, histtype='step', bins=BINS)
plt.xlabel("log(time/days)")
plt.ylabel("#")
plt.legend()
plt.savefig(f"{PLOT_DIR}/timescp.png", bbox_inches='tight')

print("Peak times.")
plt.figure()
plt.hist([log(vla['pt']), log(gw['pt']), log(all_events['pt'])],
        label=['Detection Time (AG counterpart to GW)', 'Peak Time (all GW)',
        'Peak Time (all events)'], histtype='step', bins=BINS)
plt.xlabel("log(time/days)")
plt.ylabel("#")
plt.legend()
plt.savefig(f"{PLOT_DIR}/dtimes.png", bbox_inches='tight')

print("Peak times, counterpart only.")
plt.figure()
plt.hist(log(vla['pt']), label='Detection Time (AG counterpart to GW)',
    histtype='step', bins=BINS)
plt.xlabel("log(time/days)")
plt.ylabel("#")
plt.legend()
plt.savefig(f"{PLOT_DIR}/dtimescp.png", bbox_inches='tight')

print("Peak fluxes (cumulative).")
plt.figure()
plt.hist([log(gw_vla['pf']), log(gw['pf']), log(all_events['pt'])],
        label=['(AG counterpart to GW)', '(all GW)', '(all mergers)'],
        cumulative=True, histtype='step', bins=BINS)
plt.xlabel("log(Peak Flux/muJy)")
plt.ylabel("#")
plt.legend()
plt.savefig(f"{PLOT_DIR}/fluxes.png", bbox_inches='tight')

print("Peak betas.")
plt.figure()
plt.hist([vla['ba'], gw['ba'], all_events['ba']],
        label=['Peak apparent beta(AG counterpart to GW)', 'Peak apparent beta (all GW)',
        'Peak apparent beta (all events)'], histtype='step', bins=BINS)
plt.xlabel("beta_app")
plt.ylabel("#")
plt.legend()
plt.savefig(f"{PLOT_DIR}/dba.png", bbox_inches='tight')

print("Peak apparent betas (counterpart only).")
plt.figure()
plt.hist([gw_vla['ba']],
        label=['Peak apparent beta (AG counterpart to GW)'],
        histtype='step', bins=BINS)
plt.xlabel("beta_app")
plt.ylabel("#")
plt.legend()
plt.savefig(f"{PLOT_DIR}/dbacp.png", bbox_inches='tight')

print("Peak gammas.")
plt.figure()
plt.hist([vla['ga'], gw['ga'], all_events['ga']],
        label=['Peak gamma (vla)', 'Peak gamma (GW)', 'Peak gamma (all merger)'],
        histtype='step', bins=BINS, density=True, )
plt.xlabel("gamma")
plt.ylabel("#")
plt.legend()
plt.savefig(f"{PLOT_DIR}/dga.png", bbox_inches='tight')

print("Peak gammas (counterpart only).")
plt.figure()
plt.hist([gw_vla['ga']],
        label=['Peak gamma (AG counterpart to GW)'],
        histtype='step', bins=BINS)
plt.xlabel("gamma")
plt.ylabel("#")
plt.legend()
plt.savefig(f"{PLOT_DIR}/dgacp.png", bbox_inches='tight')

print("Peak fluxes (cumulative), GW events only.")
plt.figure()
plt.hist(log(gw['pf']), label='(all GW)',
        cumulative=True, histtype='step', bins=BINS)
plt.xlabel("log(Peak Flux/muJy)")
plt.ylabel("#")
plt.legend()
plt.savefig(f"{PLOT_DIR}/fluxesgw.png", bbox_inches='tight')

print("Peak fluxes.")
plt.figure()
plt.hist([log(gw_vla['pf']), log(gw['pf']), log(all_events['pt'])],
        label=['(AG counterpart to GW)', '(all GW)', '(all mergers)'],
        histtype='step', density=True, bins=BINS)
plt.xlabel("log(Peak Flux/muJy)")
plt.ylabel("#")
plt.legend()
plt.savefig(f"{PLOT_DIR}/dfluxes.png", bbox_inches='tight')

print("Distance (cumulative).")
plt.figure()
plt.hist([gw_vla['d'], gw['d'], all_events['d']],
        label=['(AG counterpart to GW)', '(all GW)', '(all mergers)'],
        cumulative=True, histtype='step', bins=BINS)
plt.xlabel("Distance/Mpc")
plt.ylabel("#")
plt.legend()
plt.savefig(f"{PLOT_DIR}/distances.png", bbox_inches='tight')

print("Distance.")
plt.figure()
plt.hist([gw_vla['d'], gw['d'], all_events['d']],
        label=['(AG counterpart to GW)', '(all GW)', '(all mergers)'],
        histtype='step', bins=BINS)
plt.xlabel("Distance/Mpc")
plt.ylabel("#")
plt.legend()
plt.savefig(f"{PLOT_DIR}/ddistances.png", bbox_inches='tight')

print("Distance, counterpart only.")
plt.figure()
plt.hist(gw_vla['d'], label='(AG counterpart to GW)', histtype='step', bins=BINS)
plt.xlabel("Distance/Mpc")
plt.ylabel("#")
plt.legend()
plt.savefig(f"{PLOT_DIR}/ddistancescp.png", bbox_inches='tight')

print("Density.")
plt.figure()
plt.hist([log(gw_vla['n']), log(gw['n']), log(all_events['n'])],
        label=['(AG counterpart to GW)', '(all GW)', '(all mergers)'],
        histtype='step', bins=BINS)
plt.xlabel("log(n/cm-3)")
plt.ylabel("#")
plt.legend()
plt.savefig(f"{PLOT_DIR}/dn.png", bbox_inches='tight')

print("Density, counterpart only.")
plt.figure()
plt.hist(log(gw_vla['n']), label='(AG counterpart to GW)',
        histtype='step', bins=BINS)
plt.xlabel("log(n/cm-3)")
plt.ylabel("#")
plt.legend()
plt.savefig(f"{PLOT_DIR}/dncp.png", bbox_inches='tight')

print("Viewing angle.")
plt.figure()
plt.hist([gw_vla['tv'], gw['tv'], all_events['tv']],
        label=['(AG counterpart to GW)', '(all GW)', '(all mergers)'],
        histtype='step', bins=BINS)
plt.xlabel("tv/rad")
plt.ylabel("#")
plt.legend()
plt.savefig(f"{PLOT_DIR}/dtv.png", bbox_inches='tight')

print("Viewing angle, counterpart only.")
plt.figure()
plt.hist(gw_vla['tv'], label='(AG counterpart to GW)', histtype='step', bins=BINS)
plt.xlabel("tv/rad")
plt.ylabel("#")
plt.legend()
plt.savefig(f"{PLOT_DIR}/dtvcp.png", bbox_inches='tight')

print("Around 0.1 radians")
plt.figure()
plt.hist(ar['ba'], label='Around 0.1rd, cp', histtype = 'step', bins=BINS)
plt.xlabel("ba")
plt.ylabel("#")
plt.legend()
plt.savefig(f"{PLOT_DIR}/arba.png")

print("Done.")
