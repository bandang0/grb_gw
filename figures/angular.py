#!/bin/env python3
"""angular plot."""

from matplotlib import pyplot as plt

from helpers import *

print("Plotting angular plot.")
plt.figure(4)
#plt.errorbar(OBS_TIMES, OBS_FLUXES, yerr = OBS_ERRORS, color = "grey", fmt="+")
plt.xscale('log')
plt.yscale('log')
plt.xlim([0.3, 1000.])
plt.ylim([2.e-7, 3e1])
plt.xlabel('Observer Time (days)')
plt.ylabel(r'$S_{\nu} \times (D/40Mpc)^{-2}$ (mJy)')
#plt.title('Light Curve (at 3 GHz)')
plt.grid(True, which='major')
lab = { 'jet28':r'$\theta_{\rm{obs}} = 28$ deg',
        'jet21':r'$\theta_{\rm{obs}} = 21$ deg',
        'jet14':r'$\theta_{\rm{obs}} = 14$ deg',
        'jet7': r'$\theta_{\rm{obs}} = 7$ deg',
        'jet0': r'$\theta_{\rm{obs}} = 0$ deg',
        'cocoon': 'rad. structure'}

col = { 'jet28':r'red',
        'jet21':r'blue',
        'jet14':r'orange',
        'jet7': r'grey',
        'jet0': r'green',
        'cocoon': 'black',
        '28':r'red',
        '21':r'blue',
        '14':r'orange',
        '07': r'grey',
        '00': r'green',
        'cocoon': 'black'}

for s in ['jet0', 'jet7', 'jet14', 'jet21', 'jet28', 'cocoon']:
    t_obs, F_nu_obs = data_list("data/daigne{}.out".format(s))
    t_obs_days = t_obs / day
    F_nu_obs_mJy = F_nu_obs * 1.e26
    plt.plot(t_obs_days, F_nu_obs_mJy, "-", label=lab[s], color = col[s])

for s in ['00', '07', '14', '21', '28', 'cocoon']:
    t_obs, F_nu_obs = daigne_data_list("daigne_output/duque_jet_01_{}.spn.lc".format(s))
    t_obs_days = t_obs / day
    F_nu_obs_mJy = F_nu_obs * 1.e3
    plt.plot(t_obs_days, F_nu_obs_mJy, "--", color = col[s])
plt.legend()
plt.savefig("plots/angular.pdf", bbox_inches='tight')
