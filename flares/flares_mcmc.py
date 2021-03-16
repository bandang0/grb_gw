from os import system
from sys import exit
import pickle
import numpy as np
import matplotlib.pyplot as plt
import emcee
from sys import argv
import pandas as pd
from helpers import *
import corner
from multiprocessing import Pool
case = argv[1]

a = -0.5
b = -2.0
Ep = 1. * keV
nup = Ep / hPlanck
nuobs = 5.0 * keV / hPlanck
nug = 75. * keV / hPlanck

t0 = 0
t1 = np.inf

ldl = 27

with open(f"data/{case}.best", 'r') as conf_file:
    lines = conf_file.readlines()
    s1_ = float(lines[0])
    s2_ = float(lines[1])
    lts_ = float(lines[2])
    lls_ = float(lines[3])
    le_ = float(lines[4])
    chic_ = float(lines[5])
    n_ = float(lines[6])
    tej_ = float(lines[7])
    tau_ = float(lines[8])
    g_ = float(lines[9])

if case == '050607':
    t0 = 200
    t1 = 2.e4
if case == '050724':
    t0 = 4.e2
if case == '060111A':
    t0 = 170
if case == '110801A':
    t0 = 150
    t1 = 300
if case == '110820A':
    t0 = 80
    t1 = 1.e5
if case == '050915A':
    t0 = 0
if case == '050820A':
    t1 = 8.e3
xrt_df = pd.read_csv(f"data/{case}.xrt", sep=" ", names=['t', 't_', 't__', 'f', 'ef', 'ef_'])

xrt_df = xrt_df.sort_values('t')
xrt_df['liso'] = 4 * Pi * (10 ** ldl) ** 2 * xrt_df['f']
xrt_df['lisoe'] = 4 * Pi * (10 ** ldl) ** 2 * np.abs(xrt_df['ef'])

def Liso_pl(t, s1, s2, lts, lls):
     if t < 10 ** lts:
         return 10 ** lls * (t / 10 ** lts) ** s1
     else:
        return 10 ** lls * (t / 10 ** lts) ** s2

def Liso_flare(t, le, chic, n, tej, tau, g):
    beta = np.sqrt(1 - 1 / g ** 2)
    r = 0.5 * n / g
    re = 2 * cLight * tau * g ** 2
    te = tej + re / (beta * cLight)
    tm = te - re * cos(max(0., chic - r)) / cLight
    tM = te - re * cos(min(Pi, chic + r)) / cLight
    if t > tM or t < tm:
        return 0
    return simple_hle_patch_band(t, 0, chic, r, te, re, g, 10 ** le, nup, a, b, XRT_0, XRT_1)

tobs_l = np.logspace(log(np.min(xrt_df['t']))-0.08, log(np.max(xrt_df['t']))+0.1, 300)
L_pl = np.array([Liso_pl(t, s1_, s2_, lts_, lls_) for t in tobs_l])
L_flare = np.array([Liso_flare(t, le_, chic_, n_, tej_, tau_, g_) for t in tobs_l])

plt.figure()
plt.scatter(xrt_df['t'], xrt_df['liso'], marker = '.', color='red')
plt.errorbar(xrt_df['t'], xrt_df['liso'], yerr=xrt_df['lisoe'], fmt='none', color='red', label='XRT')
plt.plot(tobs_l, L_pl, color='black', label='PL')
plt.plot(tobs_l, L_flare, color='black', label='HLE flare')
plt.plot(tobs_l, L_pl + L_flare, color='blue', label='Total')
plt.xscale('log')
plt.yscale('log')
plt.text(np.min(xrt_df['t']), np.min(xrt_df['liso']), r"$\chi$ = " + f"{chic_}\n"\
                                                + r"$N$ = " + f"{n_}\n"\
                                                + r"$t_{ej}$ = " + f"{tej_} s\n"\
                                                + r"$\tau$ = " + f"{tau_} s\n"\
                                                + r"$\Gamma$ = " + f"{g_}")
plt.xlabel('Time since burst [s]')
plt.ylabel('XRT-integrated luminosity [erg/s]')
if t0 != 0:
    plt.axvspan(np.min(tobs_l), t0, color="grey", alpha=0.5)
if t1 != np.inf:
    plt.axvspan(np.max(tobs_l), t1, color="grey", alpha=0.5)
plt.legend()
plt.title(f"GRB{case}")
plt.savefig(f"plots/flares/{case}_lc.pdf", bbox_inches="tight")

if not 'mcmc' in argv:
    exit(0)

def lnlike(theta):
    if 'bestchi' not in lnlike.__dict__:
        lnlike.bestchi = -np.inf
        lnlike.bestparam = None
    s1, s2, lts, lls, le, chic, n, tej, tau, g = theta
    err = 0
    for k, tobs in enumerate(xrt_df['t']):
        if tobs < t0 or t1 < tobs:
            continue
        err = err - 0.5 * (Liso_pl(tobs, s1, s2, lts, lls) + Liso_flare(tobs, le, chic, n, tej, tau, g) - xrt_df['liso'][k]) ** 2 / xrt_df['lisoe'][k] ** 2

    if err > lnlike.bestchi:
        lnlike.bestchi  = err
        lnlike.bestparam = theta


    return err

def lnprob(theta):
    "Simply Bayes theorem"
    s1, s2, lts, lls, le, chic, n, tej, tau, g= theta

    if s1 < -11 or s1 > 0 or s2 < -3 or s2 > 0 or lls > 46 or lls < 41 or lts < 2 or lts > 6 or le < 45  or le > 58 or chic < 0.03 or tej < -30 or tej > 400 or n < 1 or g < 75 or g > 1000 or n > 50 or tau < 0.01 or tau > 40 or chic > 0.3:
        return -np.inf
    return lnlike(theta)

ndim = 10
nwalkers = 22
nsteps = 20000


optimal_set = [s1_, s2_, lts_, lls_, le_, chic_, n_, tej_, tau_, g_]
first_guess = [optimal_set + 0.5 * np.random.randn(ndim) * optimal_set for i in range(nwalkers)]
print(f"Optimal set chi^2: {lnlike(optimal_set)}")

sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args = ())
burnin = sampler.run_mcmc(first_guess, 2000, progress=True, skip_initial_state_check=True)
print(f"Best model for {case} ({lnlike.bestchi}):")
for param in lnlike.bestparam:
    print(f"{param:.4g}")
sampler.reset()
sampler.run_mcmc(burnin, nsteps, progress=True, skip_initial_state_check=True)
samples = sampler.chain[:, nwalkers:, :].reshape((-1, ndim))

print(f"Best model for {case} ({lnlike.bestchi}):")
for param in lnlike.bestparam:
    print(f"{param:.4g}")

with open(f"data/{case}.pickle", 'wb') as pickle_file:
    pickle.dump(samples, pickle_file)

with open(f"data/{case}.pickle", "rb") as pickle_file:
    samples = pickle.load(pickle_file)

fig = corner.corner(samples, labels = ["s1", "s2", "lts", "lls", "le", "chic", "n", "tej", "tau", "g"], color = "black", quantiles=[0.32, 0.5, 0.68], show_titles=True, title_fmt=".2f", title_kwargs={"fontsize": 12, "color":"red"})

fig.savefig(f"plots/flares/{case}_corner.pdf", bbox_inches="tight")
