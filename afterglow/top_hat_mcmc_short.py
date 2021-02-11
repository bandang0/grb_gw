"""Fit a jet afterglow to the decreasing phase of the 3GHz photometry."""

from os import system
from sys import exit
import numpy as np
import matplotlib.pyplot as plt
import lmfit

from helpers import *

p = lmfit.Parameters()
p.add('n', value=1.e-3, min=1.e-6, max=1.)
p.add('eb', value=4.e-5, min=1.e-6, max=1.)
p.add('E0', value=10**52.3, min=1.e47)
p.add('tv', value=0.32, min=0.01, max=1.57)
p.add('tj', value=0.08, min=0.02, max=1.1)

def residual(p):
    conf_dic = {'L': 0,
                'n': p['n'].value,
                'eb': p['eb'].value,
                'E0': p['E0'].value,
                'tv': p['tv'].value,
                'tj': p['tj'].value,
                'first_day': 0,
                'last_day': 350}
    write_conf_file(TMP_CONF_FILE, conf_dic)

    system("./remnant_light_curve.exe tmp > /dev/null")
    # Parse fname into the listed data
    t_obs, F_nu_obs = data_list("data/tmp.out")

    # get good units
    t_obs_days = t_obs / day
    F_nu_obs_muJy = 1.e29 * F_nu_obs

    # calculate the chi2 terms array for the decreasing phase
    err = list()
    time_ind = 0
    for k in range(len(OBS_TIMES_2_d)):
        while t_obs_days[time_ind] < OBS_TIMES_2_d[k]:
            time_ind = time_ind + 1
        err.append((F_nu_obs_muJy[time_ind] - OBS_FLUXES_2_muJy[k])\
                / OBS_ERRORS_2_muJy[k])

    return np.array(err)

# create Minimizer
mini = lmfit.Minimizer(residual, p, nan_policy='omit')

# first solve with Nelder-Mead
out = mini.minimize()
print(lmfit.fit_report(out))

#Confidence intervals
ci = lmfit.conf_interval(mini, out, sigmas=[1, 2], verbose=False)
lmfit.printfuncs.report_ci(ci)
