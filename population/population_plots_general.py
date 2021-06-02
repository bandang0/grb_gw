'''Plot the histogram of distances.'''

import matplotlib.pyplot as plt
import pandas as pd
from numpy import log10 as log
from sys import argv
from helpers import *
from scipy.optimize import fsolve

case = argv[1]
DATA_FILE = f"data/full_ag_{case}.data"
PLOT_DIR = f"plots/full_ag_{case}"
BINS = 300
p=2.2
ee = 0.1

# Read event data and create pandas data frame.
all_events = pd.read_csv(DATA_FILE, sep=" ",
        names=['ig', 'iv', 'd', 'n', 'e', 'eb',
               'tv', 'pt', 'pf', 'pa', 'pi', 'pc', 'pg', 'pb', 'pd', 'pr', 'pn'])




# Collect event data according to detections in GW, AG, etc.
gw_vla = all_events.loc[all_events['ig'] * all_events['iv'] == 1]
vla = all_events.loc[all_events['iv'] == 1]
gw = all_events.loc[all_events['ig'] == 1]

print(f"Maximum D in GW: {gw['d'].max()}, mean: {gw['d'].mean()}")
print(f"Mean viewing angle: {gw['tv'].mean() / Deg}")
print(f"Maximum nu_i: {gw_vla['pi'].max() / GHz} GHz")

all_events['yy'] = ee * (all_events['pi'] / all_events['pc']) ** (0.1) / (0.8 * all_events['eb'])

all_events['y'] = 0.5 * (sqrt(4 * all_events['yy']+1)-1)
all_events['yp'] = 0.
for i in range(len(all_events.y)):
    first_guess = all_events['y'].iloc[i]
    x = all_events['yy'].iloc[i]
    func = lambda z: z * (1 + z)**0.8 - x
    all_events['yp'].iloc[i] = fsolve(func, first_guess)

print(f"Median log-Compton (basic): {log(gw_vla['y']).median()}, 95%: {np.percentile(log(gw_vla['y']), 95)}")
print(f"Median log-Compton  (full): {log(gw_vla['yp']).median()}, 95%: {np.percentile(log(gw_vla['yp']), 95)}")
print(f"Minimum peak index: {gw_vla['pn'].min()}")

print("Peak times (cumulative).")
plt.figure()
plt.hist([log(vla['pt']), log(gw['pt']), log(all_events['pt'])],
    label=['Peak Time (AG counterpart to GW)', 'Peak Time (all GW)',
    'Peak Time (all events)'], cumulative=True, histtype='step', bins=BINS)
plt.xlabel(r"$\log( T_p/{\rm days})$")
plt.ylabel('N')
plt.legend()
plt.title(f"{case}")
plt.savefig(f"{PLOT_DIR}/pt.pdf", bbox_inches='tight')

print("Peak times (cumulative), counterpart only.")
plt.figure()
plt.hist(log(vla['pt']), label='Peak Time (AG counterpart to GW)',
    cumulative=True, histtype='step', bins=BINS)
plt.xlabel(r"$\log( T_p/{\rm days})$")
plt.ylabel('N')
plt.legend()
plt.title(f"{case}")
plt.savefig(f"{PLOT_DIR}/ptcp.pdf", bbox_inches='tight')

print("Peak times.")
plt.figure()
plt.hist([log(vla['pt']), log(gw['pt']), log(all_events['pt'])],
        label=['Peak Time (AG counterpart to GW)', 'Peak Time (all GW)',
        'Peak Time (all events)'], histtype='step', bins=BINS)
plt.xlabel(r"$\log( T_p/{\rm days})$")
plt.ylabel('N')
plt.legend()
plt.title(f"{case}")
plt.savefig(f"{PLOT_DIR}/dpt.pdf", bbox_inches='tight')

print("Peak times, counterpart only.")
plt.figure()
plt.hist(log(vla['pt']), label='Peak Time (AG counterpart to GW)',
    histtype='step', bins=BINS)
plt.xlabel(r"$\log( T_p/{\rm days})$")
plt.ylabel('N')
plt.legend()
plt.title(f"{case}")
plt.savefig(f"{PLOT_DIR}/dptcp.pdf", bbox_inches='tight')

print("Peak gammas.")
plt.figure()
plt.hist([vla['pg'], gw['pg'], all_events['pg']],
        label=['Peak gamma (AG counterpart to GW)', 'Peak gamma (all GW)',
        'Peak gamma (all events)'], histtype='step', bins=BINS)
plt.xlabel("$\Gamma$")
plt.ylabel('N')
plt.legend()
plt.title(f"{case}")
plt.savefig(f"{PLOT_DIR}/dpg.pdf", bbox_inches='tight')

print("Peak gammas (AG only).")
plt.figure()
plt.hist([gw_vla['pg']],
        label=['Peak gamma (AG counterpart to GW)'],
        histtype='step', bins=BINS)
plt.xlabel("$\Gamma$")
plt.ylabel('N')
plt.legend()
plt.title(f"{case}")
plt.savefig(f"{PLOT_DIR}/dpgcp.pdf", bbox_inches='tight')

print("Peak fluxes (cumulative).")
plt.figure()
plt.hist([log(gw_vla['pf']), log(gw['pf']), log(all_events['pf'])],
        label=['(AG counterpart to GW)', '(all GW)', '(all mergers)'],
        cumulative=True, histtype='step', bins=BINS)
plt.xlabel(r"$\log(F_p/\mu {\rm Jy})$")
plt.ylabel('N')
plt.legend()
plt.title(f"{case}")
plt.savefig(f"{PLOT_DIR}/f.pdf", bbox_inches='tight')

print("Peak fluxes (cumulative), GW events only.")
plt.figure()
plt.hist(log(gw['pf']), label='(all GW)',
        cumulative=True, histtype='step', bins=BINS)
plt.xlabel(r"$\log(F_p/\mu {\rm Jy})$")
plt.ylabel('N')
plt.legend()
plt.title(f"{case}")
plt.savefig(f"{PLOT_DIR}/fgw.pdf", bbox_inches='tight')

print("Peak fluxes.")
plt.figure()
plt.hist([log(gw_vla['pf']), log(gw['pf']), log(all_events['pf'])],
        label=['(AG counterpart to GW)', '(all GW)', '(all mergers)'],
        histtype='step', bins=BINS)
plt.xlabel(r"$\log(F_p/\mu {\rm Jy})$")
plt.ylabel('N')
plt.legend()
plt.title(f"{case}")
plt.savefig(f"{PLOT_DIR}/dpf.pdf", bbox_inches='tight')

print("Peak frequencies.")
plt.figure()
plt.hist([log(gw_vla['pa'] / GHz), log(gw_vla['pi'] / GHz), log(hPlanck * gw_vla['pc'] / keV)],
        label=['absorption', 'injection', 'cooling'],
        histtype='step', bins=BINS)
plt.xlabel(r"$\log (\nu_{a, i}/{\rm GHz}), \log (\nu_c / {\rm keV})$")
plt.ylabel('N')
plt.legend()
plt.title(f"{case}")
plt.savefig(f"{PLOT_DIR}/df.pdf", bbox_inches='tight')

print("Energy.")
plt.figure()
plt.hist([log(gw_vla['e']), log(gw['e']), log(all_events['e'])],
        label=['(AG counterpart to GW)', '(all GW)', '(all mergers)'],
        histtype='step', bins=BINS)
plt.xlabel(r"$\log E_{\rm iso}/{\rm erg}$")
plt.ylabel('N')
plt.legend()
plt.title(f"{case}")
plt.savefig(f"{PLOT_DIR}/de.pdf", bbox_inches='tight')

print("Energy, counterpart only.")
plt.figure()
plt.hist(log(gw_vla['e']),
    histtype='step', bins=BINS)
plt.xlabel(r"$\log E_{\rm iso}/{\rm erg}$")
plt.ylabel('N')
plt.title(f"{case}")
plt.savefig(f"{PLOT_DIR}/decp.pdf", bbox_inches='tight')

print("Distance (Cumulative).")
plt.figure()
plt.hist([vla['d'], gw_vla['d'], gw['d'], all_events['d']],
        label=['vla', '(AG counterpart to GW)', '(all GW)', '(all mergers)'],
        cumulative=True, histtype='step', bins=BINS)
plt.xlabel(r"$D /{\rm Mpc}$")
plt.ylabel('N')
plt.legend()
plt.title(f"{case}")
plt.savefig(f"{PLOT_DIR}/d.pdf", bbox_inches='tight')

print("Distance.")
plt.figure()
plt.hist([gw_vla['d'], gw['d'], all_events['d']],
        label=['(AG counterpart to GW)', '(all GW)', '(all mergers)'],
        histtype='step', bins=BINS)
plt.xlabel(r"$D/{\rm Mpc}$")
plt.ylabel('N')
plt.legend()
plt.title(f"{case}")
plt.savefig(f"{PLOT_DIR}/dd.pdf", bbox_inches='tight')

print("Distance, counterpart only.")
plt.figure()
plt.hist(gw_vla['d'], histtype='step', bins=BINS)
plt.xlabel(r"$D/{\rm Mpc}$")
plt.ylabel('N')
plt.legend()
plt.title(f"{case}")
plt.savefig(f"{PLOT_DIR}/ddcp.pdf", bbox_inches='tight')

print("Density.")
plt.figure()
plt.hist([log(gw_vla['n']), log(gw['n']), log(all_events['n'])],
        label=['(AG counterpart to GW)', '(all GW)', '(all mergers)'],
        histtype='step', bins=BINS)
plt.xlabel(r"$\log( n/{\rm cm}^{-3})$")
plt.ylabel('N')
plt.legend()
plt.title(f"{case}")
plt.savefig(f"{PLOT_DIR}/dn.pdf", bbox_inches='tight')

print("Density, counterpart only.")
plt.figure()
plt.hist(log(gw_vla['n']),
        histtype='step', bins=BINS)
plt.xlabel(r"$\log( n/{\rm cm}^{-3})$")
plt.ylabel('N')
plt.title(f"{case}")
plt.savefig(f"{PLOT_DIR}/dncp.pdf", bbox_inches='tight')

print("Viewing angle.")
plt.figure()
plt.hist([gw_vla['tv'] / Deg, gw['tv'] / Deg, all_events['tv'] / Deg],
        label=['(AG counterpart to GW)', '(all GW)', '(all mergers)'],
        histtype='step', bins=BINS)
plt.xlabel(r"$\theta_v/{\rm deg}$")
plt.ylabel('N')
plt.legend()
plt.title(f"{case}")
plt.savefig(f"{PLOT_DIR}/dtv.pdf", bbox_inches='tight')

print("Viewing angle, counterpart only.")
plt.figure()
plt.hist(gw_vla['tv'] / Deg, label='(AG counterpart to GW)', histtype='step', bins=BINS)
plt.xlabel(r"$\theta_v/{\rm deg}$")
plt.ylabel('N')
plt.legend()
plt.title(f"{case}")
plt.savefig(f"{PLOT_DIR}/dtvcp.pdf", bbox_inches='tight')

print("Viewing angle times gamma")
plt.figure()
plt.hist(gw_vla['tv'] * gw_vla['pg'], label=r'$\theta_v \times \Gamma$', histtype='step', bins=BINS)
plt.xlabel(r"$\theta_v \times \Gamma$")
plt.ylabel('N')
plt.legend()
plt.title(f"{case}")
plt.savefig(f"{PLOT_DIR}/dtvpg.pdf", bbox_inches='tight')
print("Done.")

print("Radius reached by jet")
plt.figure()
plt.hist(log(gw_vla['pr'] / pc), histtype='step', bins=BINS)
plt.xlabel(r"$\log R_{\rm max}$/ pc")
plt.ylabel('N')
plt.legend()
plt.title(f"{case}")
plt.savefig(f"{PLOT_DIR}/dr.pdf", bbox_inches='tight')
print("Done.")

print("Index of peak flux")
plt.figure()
plt.hist(gw_vla['pn'], histtype='step', bins=BINS)
plt.xlabel("$n_p$")
plt.ylabel('N')
plt.legend()
plt.title(f"{case}")
plt.savefig(f"{PLOT_DIR}/dnp.pdf", bbox_inches='tight')
print("Done.")

print("Injection vs. viewing angle")
plt.figure()
plt.hist2d(log(gw_vla['pi'] / GHz), gw_vla['tv'] / Deg, bins=100)
plt.xlabel(r"$\log \nu_i / {\rm GHz}$")
plt.ylabel(r'$\theta_v$ / deg')
plt.legend()
plt.title(f"{case}")
plt.savefig(f"{PLOT_DIR}/dnuidtv.pdf", bbox_inches='tight')
print("Done.")

print("nu_a vs. density")
plt.figure()
plt.hist2d(log(gw_vla['pa'] / GHz), log(gw_vla['n']), bins=100)
plt.xlabel(r"$\log \nu_a / {\rm GHz}$")
plt.ylabel(r'$\log n$')
plt.legend()
plt.title(f"{case}")
plt.savefig(f"{PLOT_DIR}/dnuadn.pdf", bbox_inches='tight')
print("Done.")

print("Viewing angle vs. gamma")
plt.figure()
plt.hist2d(gw_vla['tv'] / Deg, gw_vla['pg'], bins=100)
plt.xlabel(r"$\theta_v$ / deg")
plt.ylabel(r'$\Gamma$')
plt.legend()
plt.title(f"{case}")
plt.savefig(f"{PLOT_DIR}/dtvdg.pdf", bbox_inches='tight')
print("Done.")

print('Compton Y')
plt.figure()
plt.hist(log(gw_vla['y']), histtype='step', bins=BINS, label="basic Y")
plt.hist(log(gw_vla['yp']), histtype='step', bins=BINS, label="full Y")
plt.xlabel("$\log Y_c$")
plt.ylabel('N')
plt.legend()
plt.title(f"{case}")
plt.savefig(f"{PLOT_DIR}/dy.pdf", bbox_inches='tight')
print("Done.")

print("Viewing angle vs. Peak time")
plt.figure()
plt.hist2d(gw_vla['tv'] / Deg, log(gw_vla['pt']), bins=50)
plt.xlabel(r"$\theta_v / {\rm deg}$")
plt.ylabel(r'$\log T_p$ / day')
plt.legend()
plt.title(f"{case}")
plt.savefig(f"{PLOT_DIR}/dtvdpt.pdf", bbox_inches='tight')
print("Done.")

print("Y vs. eb")
plt.figure()
plt.hist2d(log(gw_vla['yp']), log(gw_vla['eb']), bins=100)
plt.xlabel(r"$\log Y_c / {\rm GHz}$")
plt.ylabel(r'$\epsilon_B$')
plt.legend()
plt.title(f"{case}")
plt.savefig(f"{PLOT_DIR}/dydeb.pdf", bbox_inches='tight')
print("Done.")
