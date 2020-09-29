"""Determine uncertainties on jet parameters for decreasing phase."""

from os import system
from sys import exit
import pickle
import numpy as np
import matplotlib.pyplot as plt
import emcee

from helpers import *

logn_ = -3.
logeb_ = -3.3
logE0_ = 52.3
tv_ = 0.32
tj_ = 0.08

def lnlike(theta):
    """Log-likelihood (i.e chi2) from paramters."""
    logn, logeb, logE0, tv, tj = theta

    conf_dic = {'L': 0,
                'n': 10 ** logn,
                'eb': 10 ** logeb,
                'E0': 10 ** logE0,
                'tv': tv,
                'tj': tj,
                'first_day':100,
                'last_day': 300}
    write_conf_file(TMP_CONF_FILE, conf_dic)

    system("./remnant_light_curve.exe tmp > /dev/null")
    # Parse fname into the listed data
    t_obs, F_nu_obs = data_list("data/tmp.out")

    # get good units
    t_obs_days = t_obs / day
    F_nu_obs_muJy = 1.e29 * F_nu_obs

    # calculate the chi2 terms array for the decreasing phase
    err = 0
    time_ind = 0
    for k in range(len(OBS_TIMES_2_d)):
        while t_obs_days[time_ind] < OBS_TIMES_2_d[k]:
            time_ind = time_ind + 1
        err = err - 0.5 * (F_nu_obs_muJy[time_ind] \
                - OBS_FLUXES_2_muJy[k]) ** 2 / OBS_ERRORS_2_muJy[k] ** 2

    return err

def lnprior(theta):
    """Gaussian priors on the parameters."""
    logn, logeb, logE0, tv, tj = theta

    return -0.5 * ((logn - logn_) ** 2 / 1.5 ** 2 \
            + (logeb - logeb_) ** 2 /1.5 ** 2 \
            + (logE0 - logE0_) ** 2 /1.5 ** 2 \
            + (tv - tv_) ** 2 / 0.2 ** 2 \
            + (tj - tj_) ** 2 / 0.2 ** 2) 
    

def lnprob(theta):
    "Simply Bayes theorem"
    _, _, _, tv, tj = theta
    if (tv < 0. or tv > 1.5 or tj < 0.02 or tj > 1.1):
        return -np.inf
    else:
        return lnlike(theta) + lnprior(theta)

ndim = 5
nwalkers = 10
nsteps = 10000

optimal_set = [logn_, logeb_, logE0_, tv_,tj_]
first_guess = [optimal_set + 1.e-3 * np.random.randn(ndim) * optimal_set for i in range(nwalkers)]

sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args = ())
sampler.run_mcmc(first_guess, nsteps)

samples = sampler.chain[:, nwalkers:, :].reshape((-1, ndim))

with open("jet_decrease.pickle", 'wb') as pickle_file:
    pickle.dump(samples, pickle_file)

with open("jet_decrease.pickle", "rb") as pickle_file:
    samples = pickle.load(pickle_file)

fig = corner.corner(samples, labels = ["log(n/cm-3)", "log(eb)", "log(E0/erg)", "tv", "tj"], color = "black", quantiles=[0.32, 0.5, 0.68], show_titles=True, title_fmt=".2f", title_kwargs={"fontsize": 12, "color":"red"})

fig.savefig("plots/corner_jet_decrease.png", bbox_inches="tight")

