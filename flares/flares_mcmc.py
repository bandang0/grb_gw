from os import system
from sys import exit
import pickle
import numpy as np
import matplotlib.pyplot as plt
import emcee
import sys
from sys import argv
import pandas as pd
from helpers import *
import corner
from multiprocessing import Pool
from scipy.stats import kstest, norm
from astropy.cosmology import FlatLambdaCDM

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

case = argv[1]

if 't90' in argv:
    t90_mode = True
else:
    t90_mode = False

print(f"T90 mode = {t90_mode}")

a = -0.1
b = -1.2
Ep = 1. * keV
nup = Ep / hPlanck
nuobs = 5.0 * keV / hPlanck
nug = 75. * keV / hPlanck

#defauts
t0 = 0
t1 = np.inf
z = 1
Eiso_p = 0

with open(f"data/{case}.best", 'r') as conf_file:
    lines = conf_file.readlines()
    s1_ = float(lines[0])
    s2_ = float(lines[1])
    lts_ = float(lines[2])
    lls_ = float(lines[3])
    le_ = float(lines[4])
    ti_ = float(lines[5])
    n_ = float(lines[6])
    tej_ = float(lines[7])
    tau_ = float(lines[8])
    g_ = float(lines[9])

with open(f"data/{case}.bestt90", 'r') as conf_file:
    lines = conf_file.readlines()
    s1_t90 = float(lines[0])
    s2_t90 = float(lines[1])
    lts_t90 = float(lines[2])
    lls_t90 = float(lines[3])
    le_t90 = float(lines[4])
    ti_t90 = float(lines[5])
    n_t90 = float(lines[6])
    tej_t90 = float(lines[7])
    tau_t90 = float(lines[8])
    g_t90 = float(lines[9])

if case == '050607':
    t0 = 200
    t1 = 2.e4
if case == '050724':
    t0 = 4.e2
if case == '050820A':
    t1 = 8.e3
if case == '050915A':
    t0 = 0
if case == '060111A':
    t0 = 170
    t1 = 8.e3
if case == '060204B':
    t0 = 200
    t1 = 8.e3
    t90 = 139.4
    dt90 = 31.5
    z = 2.3393
    Eiso_p = 3.7e52
if case == '060413':
    t0 = 200
    t1 = 2.e4
    t90 = 147.7
    dt90 = 26.8
if case == '070704':
    t0 = 0
    t90 = 384.2
    dt90 = 24.5
if case == '081008':
    t0 = 200
    t1 = 6.e3
    t90 = 185.5
    dt90 = 40
    z = 1.9685
    Eiso_p = 9.5e52
if case == '081102':
    t0 = 300
    t90 = 62.7
    dt90 = 16.3
if case == '100619A':
    t0 = 270
    t1 = 1.e5
    t90 = 97.5
    dt90 = 1.5
if case == '110801A':
    t0 = 250
    t1 = 2.e3
if case == '110407A':
    t0 = 160
    t1 = 1.e4
    t90 = 145.5
    dt90 = 13
if case == '110820A':
    t0 = 80
    t1 = 1.e5
if case == '120728A':
    t0 = 0
    t90 = 22
    dt90 = 2.6

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
ldl = log(cosmo.luminosity_distance(z).value * Mpc)

xrt_df = pd.read_csv(f"data/{case}.xrt", sep=" ", names=['t', 't_', 't__', 'f', 'ef', 'ef_'])

xrt_df = xrt_df.sort_values('t')
xrt_df['liso'] = 4 * Pi * (10 ** ldl) ** 2 * xrt_df['f']
xrt_df['lisoe'] = 4 * Pi * (10 ** ldl) ** 2 * np.abs(xrt_df['ef'])
xrt_df['llisoe'] = log(xrt_df['liso'] + xrt_df['lisoe']) - log(xrt_df['liso'])
xrt_df['lliso'] = log(xrt_df['liso'])

times = list(xrt_df['t'])
fluxes = list(xrt_df['lliso'])
errors = list(xrt_df['llisoe'])
filtered_indices = [k for k in range(len(times)) if (times[k] > t0 and times[k] < t1)]
filtered_times = [times[k] for k in filtered_indices]
filtered_fluxes = [fluxes[k] for k in filtered_indices]
filtered_errors = [errors[k] for k in filtered_indices]
print(f"Maxmimum flux: {10 ** max(filtered_fluxes):.2e} erg/s at {filtered_times[filtered_fluxes.index(max(filtered_fluxes))]:.1f} s." )

def Liso_pl(t, s1, s2, lts, lls):
     if t < 10 ** lts:
         return 10 ** lls * (t / 10 ** lts) ** s1
     else:
        return 10 ** lls * (t / 10 ** lts) ** s2

def Liso_flare(t, le, ti, n, tej, tau, g):
    beta = np.sqrt(1 - 1 / g ** 2)
    r = 0.5 * n / g
    re = 2 * cLight * tau * g ** 2
    te = tej + re / (beta * cLight)
    tm = te - re * cos(max(0., ti)) / cLight
    tM = te - re * cos(min(Pi, ti + 2 * r)) / cLight
    if t > tM or t < tm:
        return 0
    return simple_hle_patch_band(t, 0, ti + r, r, te, re, g, 10 ** le, nup, a, b, XRT_0, XRT_1)

def lnlike(theta):
    if 'bestchi' not in lnlike.__dict__:
        lnlike.bestchi = -np.inf
        lnlike.bestparam = None
    s1, s2, lts, lls, le, ti, n, tej, tau, g = theta
    #s1, s2, lts, lls, ti = theta
    err = 0
    for k, tobs in enumerate(filtered_times):
        err = err - 1. * (log(Liso_pl(tobs, s1, s2, lts, lls) + Liso_flare(tobs, le, ti, n, tej, tau, g)) - filtered_fluxes[k]) ** 2 / filtered_errors[k] ** 2

    if err > lnlike.bestchi:
        lnlike.bestchi  = err
        lnlike.bestparam = theta
    return err

def lnprob(theta):
    "Simply Bayes theorem"
    s1, s2, lts, lls, le, ti, n, tej, tau, g = theta
    #s1, s2, lts, lls, ti = theta
    if t90_mode and tej > (t90 + dt90):
        return -np.inf
    if s1 < -11 or s1 > 0 \
    or s2 < -3 or s2 > 1 \
    or lls < 41 or lls > 49 \
    or lts < 2 or lts > 6 \
    or le < 45  or le > 54 \
    or ti < 0.01 or ti > 0.3\
    or tej < -30 \
    or n < 1 or n > 50\
    or g < 50 or g > 1000 \
    or tau < 0.01 or tau > 40:
        return -np.inf
    return lnlike(theta)

tobs_l = np.logspace(log(np.min(xrt_df['t']))-0.08, log(np.max(xrt_df['t']))+0.1, 300)
L_pl = np.array([Liso_pl(t, s1_, s2_, lts_, lls_) for t in tobs_l])
L_flare = np.array([Liso_flare(t, le_, ti_, n_, tej_, tau_, g_) for t in tobs_l])
L_pl_t90 = np.array([Liso_pl(t, s1_t90, s2_t90, lts_t90, lls_t90) for t in tobs_l])
L_flare_t90 = np.array([Liso_flare(t, le_t90, ti_t90, n_t90, tej_t90, tau_t90, g_t90) for t in tobs_l])

normalized_residuals = [(log(Liso_pl(tobs, s1_, s2_, lts_, lls_) + Liso_flare(tobs, le_, ti_, n_, tej_, tau_, g_)) - filtered_fluxes[k]) / filtered_errors[k] for k, tobs in enumerate(filtered_times)]
D, p = kstest(normalized_residuals, 'norm', mode='exact')

normalized_residualst90 = [(log(Liso_pl(tobs, s1_t90, s2_t90, lts_t90, lls_t90) + Liso_flare(tobs, le_t90, ti_t90, n_t90, tej_t90, tau_t90, g_t90)) - filtered_fluxes[k]) / filtered_errors[k] for k, tobs in enumerate(filtered_times)]
D_t90, p_t90 = kstest(normalized_residualst90, 'norm', mode='exact')
# print("Residual cumulative distribution function")
# x_l = np.linspace(min(normalized_residuals), max(normalized_residuals), 200)
# plt.figure()
# plt.plot(np.sort(normalized_residuals), np.linspace(0, 1, len(normalized_residuals), endpoint=False))
# plt.plot(x_l, norm.cdf(x_l))
# plt.xlabel("Normalized residual")
# plt.title(f"GRB{case}")
# plt.ylim([0, 1])
# plt.savefig(f"plots/flares/{case}_nr.pdf", bbox_inches="tight")

print("Plotting light-curve with fit")
plt.figure()
plt.scatter(xrt_df['t'], xrt_df['liso'], marker = '.', color='red')
plt.errorbar(xrt_df['t'], xrt_df['liso'], yerr=xrt_df['lisoe'], fmt='none', color='red', label='XRT')
#plt.plot(tobs_l, L_pl, color='black', label='PL')
#plt.plot(tobs_l, L_flare, color='black', label='HLE flare')
plt.plot(tobs_l, L_pl + L_flare, color='blue', label=f'Total w/o T90 (p = {p:.3g})')
plt.plot(tobs_l, L_pl_t90 + L_flare_t90, color='green', label=f'Total w/ T90 (p = {p_t90:.3g})')
plt.xscale('log')
plt.yscale('log')
plt.text(np.min(xrt_df['t']), np.min(xrt_df['liso']),
         "w/o T90 constraint:\n"\
        + r"$E_{\rm iso}$ = " + f"{10 ** le_:.1e}\n"\
        + r"$\Gamma$ = " + f"{g_:.3g}\n"
        + r"$N$ = " + f"{n_:.3g}\n"\
        + r"$\theta_i$ = " + f"{ti_:.3g}\n"\
        + r"$t_{ej}$ = " + f"{tej_:.3g} s\n"\
        + r"$\tau$ = " + f"{tau_:.3g} s")

plt.text(12 * np.min(xrt_df['t']), np.min(xrt_df['liso']),
         "w/ T90 constraint:\n"\
        + r"$E_{\rm iso}$ = " + f"{10 ** le_t90:.1e}\n"\
        + r"$\Gamma$ = " + f"{g_t90:.3g}\n"
        + r"$N$ = " + f"{n_t90:.3g}\n"\
        + r"$\theta_i$ = " + f"{ti_t90:.3g}\n"\
        + r"$t_{ej}$ = " + f"{tej_t90:.3g} s\n"\
        + r"$\tau$ = " + f"{tau_t90:.3g} s")
plt.xlabel('Time after trigger [s]')
plt.ylabel('XRT luminosity (0.3-30 keV) [erg/s]')
if t0 != 0:
    plt.axvspan(np.min(tobs_l), t0, color="grey", alpha=0.5)
if t1 != np.inf:
    plt.axvspan(t1, np.max(tobs_l), color="grey", alpha=0.5)
plt.legend()
plt.ylim([0.9 * min(xrt_df['liso']), 2 * max(xrt_df['liso'])])
if Eiso_p == 0:
    plt.title(f"GRB{case} (" + r"$T_{90}$" + f" = {t90}" + r"$\pm$" + f"{dt90} s)")
else:
    plt.title(f"GRB{case} (" + r"$T_{90}$" + f" = {t90}" + r"$\pm$" + f"{dt90} s, {Eiso_p} erg)")
plt.savefig(f"plots/flares/{case}_lc.pdf", bbox_inches="tight")

if not 'mcmc' in argv:
    exit(0)

ndim = 10
nwalkers = 22
nsteps = 20000

if t90_mode:
    optimal_set = [s1_t90, s2_t90, lts_t90, lls_t90, le_t90, ti_t90,  n_t90, tej_t90, tau_t90, g_t90]
    corner_file = f"plots/flares/{case}_cornert90.pdf"
else:
    corner_file = f"plots/flares/{case}_corner.pdf"
    optimal_set = [s1_, s2_, lts_, lls_, le_, ti_,  n_, tej_, tau_, g_]

initial_range = [1, 1,   1,    1,    0.5,   0.01, 1.5,  20,   0.5,    15]
first_guess = []
while len(first_guess) < nwalkers:
    new = optimal_set + 1 * np.random.randn(ndim) * initial_range
    if lnprob(new) > -np.inf:
        first_guess.append(new)

print([int(lnprob(new)) for new in first_guess].sort())
print(f"Optimal set chi^2 (T90 = {t90_mode}): {lnlike(optimal_set):.1f}")

sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args = ())
burnin = sampler.run_mcmc(first_guess, 1000, progress=True, skip_initial_state_check=True)
print(f"Best model for {case} (T90 = {t90_mode}) (chi^2 = {lnlike.bestchi:.1f}):")
for param in lnlike.bestparam:
    print(f"{param:.4g}")
sampler.reset()
sampler.run_mcmc(burnin, nsteps, progress=True, skip_initial_state_check=True)
samples = sampler.chain[:, nwalkers:, :].reshape((-1, ndim))

print(f"Best model for {case} (T90 = {t90_mode}) (chi^2 = {lnlike.bestchi:.1f}):")
for param in lnlike.bestparam:
    print(f"{param:.4g}")

#with open(f"data/{case}.pickle", 'wb') as pickle_file:
#    pickle.dump(samples, pickle_file)

#with open(f"data/{case}.pickle", "rb") as pickle_file:
#    samples = pickle.load(pickle_file)

fig = corner.corner(samples[:, 4:], labels = [r"$\log E_{iso}$", r"$\theta_i$", r"$n$", r"$t_{ej}$", r"$\tau$", "$\Gamma$"], color = "black", quantiles=[0.32, 0.5, 0.68], show_titles=True, title_fmt=".2f", title_kwargs={"fontsize": 12, "color":"red"})

fig.savefig(corner_file, bbox_inches="tight")
