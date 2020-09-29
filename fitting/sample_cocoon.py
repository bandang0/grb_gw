"""Determine uncertainties on jet parameters for decreasing phase."""

from os import system
from sys import exit
import pickle
import numpy as np
import matplotlib.pyplot as plt
import emcee

from helpers import *

def lnlike(theta):
    """Log-likelihood (i.e chi2) from paramters."""
    logn, logeb, logE0, um, uM = theta

    conf_dic = {'L': 1,
                'n': 10 ** logn,
                'eb': 10 ** logeb,
                'E0': 10 ** logE0,
                'um': um,
                'uM': uM,
                'first_day':10,
                'last_day': 300}
    write_conf_file(TMP_CONF_FILE, conf_dic)

    system("./remnant_light_curve tmp > /dev/null")
    # Parse fname into the listed data
    t_obs, F_nu_obs = data_list("data/tmp.out")

    # get good units
    t_obs_days = t_obs / day
    F_nu_obs_muJy = 1.e29 * F_nu_obs

    # calculate the chi2 terms array for the whole light curve
    err = 0
    time_ind = 0
    for k in range(len(OBS_TIMES_d)):
        while t_obs_days[time_ind] < OBS_TIMES_d[k]:
            time_ind = time_ind + 1
        err = err - 0.5 * (F_nu_obs_muJy[time_ind] - OBS_FLUXES_muJy[k]) ** 2 / OBS_ERRORS_muJy[k] ** 2
    print(err)
    return err

def lnprior(theta):
    """Uniform priors on the paramters."""
    logn, logeb, logE0, um, uM = theta

    if -6 < logn < -1 and -6 < logeb < -1 and 47. < logE0 < 54. and 0. < um < 2. and 0. < uM < 4. and um < uM - 0.1:
        return 0.
    return -np.inf

def lnprob(theta):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta)

ndim = 5
nwalkers = 50

optimal_set = [-3., -3.3, 51, 1.5, 3.]
first_guess = [optimal_set + 1.e-2 * np.random.randn(ndim) * optimal_set for i in range(nwalkers)]

sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args = ())
sampler.run_mcmc(first_guess, 5 * nwalkers)

samples = sampler.chain[:, nwalkers:, :].reshape((-1, ndim))

with open("cocoon.pickle", 'wb') as pickle_file:
    pickle.dump(samples, pickle_file)
