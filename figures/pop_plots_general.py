'''Plot the histogram of distances.'''

import matplotlib.pyplot as plt
import pandas as pd
from numpy import log10 as log
DATA_FILE = "general_sample.data"
PLOT_DIR = "plots/general_sample"
BINS = 25
# Read event data and create pandas data frame.
all_events = pd.read_csv(DATA_FILE, sep=" ",
        names=['d', 'n', 'e0', 'eb', 'tv', 'tj', 'pt', 'pf', 'ba',
        'ig', 'iv', 'is', 'in'])

# Collect event data according to detections in GW, AG, etc.
gw_vla = all_events.loc[all_events['ig'] * all_events['iv'] == 1]
vla = all_events.loc[all_events['iv'] == 1]
gw = all_events.loc[all_events['ig'] == 1]

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

print("Peak betas.")
plt.figure()
plt.hist([vla['ba'], gw['ba'], all_events['ba']],
        label=['Peak apparent beta(AG counterpart to GW)', 'Peak apparent beta (all GW)',
        'Peak apparent beta (all events)'], histtype='step', bins=BINS)
plt.xlabel("beta_app")
plt.ylabel("#")
plt.legend()
plt.savefig(f"{PLOT_DIR}/dba.png", bbox_inches='tight')

print("Peak betas (AG only).")
plt.figure()
plt.hist([gw_vla['ba']],
        label=['Peak apparent beta(AG counterpart to GW)'],
        histtype='step', bins=BINS)
plt.xlabel("beta_app")
plt.ylabel("#")
plt.legend()
plt.savefig(f"{PLOT_DIR}/dbacp.png", bbox_inches='tight')

print("Peak fluxes (cumulative).")
plt.figure()
plt.hist([log(gw_vla['pf']), log(gw['pf']), log(all_events['pt'])],
        label=['(AG counterpart to GW)', '(all GW)', '(all mergers)'],
        cumulative=True, histtype='step', bins=BINS)
plt.xlabel("log(Peak Flux/muJy)")
plt.ylabel("#")
plt.legend()
plt.savefig(f"{PLOT_DIR}/fluxes.png", bbox_inches='tight')

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
        histtype='step', bins=BINS)
plt.xlabel("log(Peak Flux/muJy)")
plt.ylabel("#")
plt.legend()
plt.savefig(f"{PLOT_DIR}/dfluxes.png", bbox_inches='tight')

print("Energy.")
plt.figure()
plt.hist([log(gw_vla['e0']), log(gw['e0']), log(all_events['e0'])],
        label=['(AG counterpart to GW)', '(all GW)', '(all mergers)'],
        histtype='step', bins=BINS)
plt.xlabel("E_0/erg")
plt.ylabel("#")
plt.legend()
plt.savefig(f"{PLOT_DIR}/de0.png", bbox_inches='tight')

print("Energy, counterpart only.")
plt.figure()
plt.hist(log(gw_vla['e0']), label='Energy (AG counterpart to GW)',
    histtype='step', bins=BINS)
plt.xlabel("log(E_0/erg)")
plt.ylabel("#")
plt.legend()
plt.savefig(f"{PLOT_DIR}/de0cp.png", bbox_inches='tight')

print("Distance (Cumulative).")
plt.figure()
plt.hist([vla['d'], gw_vla['d'], gw['d'], all_events['d']],
        label=['vla', '(AG counterpart to GW)', '(all GW)', '(all mergers)'],
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

print("Done.")
