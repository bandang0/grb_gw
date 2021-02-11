from sys import exit
import matplotlib.pyplot as plt
import pandas as pd
from numpy import sin, cos
import numpy as np
from helpers import *
from sys import argv
import corner
import emcee

DATA_FILE = "data/g16_f.data"
PLOT_DIR = f"plots/aa18_plots"
plt.style.use('palmerio')
BINS = 70

df = pd.read_csv("data/dainotti_2011_MNRAS_418.data", sep=" ", comment="#")

L = df.iloc[:,14]
dL = df.iloc[:,15]
T = df.iloc[:,12]
dT = df.iloc[:,13]
E = df.iloc[:,16]
dE = df.iloc[:,17]
A = L - E
dA = dL + dE

plt.figure()
m = -1.06
s = 50.06
sig = 0.69
X = np.linspace(1, 6, 20)
plt.errorbar(T, L, xerr=dT, yerr=dL, fmt="o", linewidth=0.5, label="Dai.+11-MNRAS418")
plt.plot(X, -1.23 * X + 51.9, ":", color="black", label="Mar.+13-MNRAS428", linewidth=0.5, )
plt.plot(X, m * X + s, color="black", label="max. likeli. lin. fit.")
plt.fill_between(X, (m) * X + s - sig, (m) * X + s + sig,
color="red", alpha=0.5, label="int. scatter")
plt.xlabel(r"$\log T_{{\rm p}}$ [s]")
plt.ylabel(r"$\log L_{X, {\rm p}}$ [erg/s]")
plt.legend(loc="lower left")
plt.xlim([1.2,5.2])
plt.ylim([43., 50])
plt.savefig("plots/plateaus/TL.pdf", bbox_inches="tight")

plt.figure()
m = 0.67
s = 11.85
sig = 0.90
X = np.linspace(48, 55, 20)
plt.errorbar(E, L, xerr=dE, yerr=dL, fmt="o", linewidth=0.5, label="Dai.+11-MNRAS418")
plt.plot(X, 1.06 * X + -8.85, ":", color="black", label="Mar.+13-MNRAS428", linewidth=0.5, )
plt.plot(X, m * X + s, color="black", label="max. likeli. lin. fit.")
plt.fill_between(X, (m) * X + s - sig, (m ) * X + s + sig,
color="red", alpha=0.5, label="int. scatter")
plt.xlabel(r"$\log  E_{\gamma, {\rm iso}}$ [erg]")
plt.ylabel(r"$\log  L_{X, {\rm p}}$ [erg/s]")
plt.legend(loc="upper left")
plt.xlim([48.2,54.6])
plt.ylim([43., 50])
plt.savefig("plots/plateaus/EL.pdf", bbox_inches="tight")

plt.figure()
m = -0.74
s = -3.41
sig = 0.63
X = np.linspace(1, 6, 20)
plt.errorbar(T, A, xerr=dT, yerr=dA, fmt="o", linewidth=0.5, label="Dai.+11-MNRAS418")
plt.plot(X, m * X + s, color="black", label="max. likeli. lin. fit.")
plt.fill_between(X, (m) * X + s - sig, (m) * X + s + sig,
color="red", alpha=0.5, label="int. scatter")
plt.xlabel(r"$\log  T_{{\rm p}}$ [s]")
plt.ylabel(r"$\log(L_{X, {\rm p}} /  E_{\gamma, {\rm iso}})$ [1/s]")
plt.legend(loc="lower left")
plt.xlim([1.2,5.2])
plt.ylim([-9, -2])
plt.savefig("plots/plateaus/TA.pdf", bbox_inches="tight")

exit()

def lnlike(theta):
    """Log-likelihood from parameters."""
    sig, a, b = theta

    return -0.5 * (np.sum(np.log(sig ** 2 + dA ** 2 + b ** 2 * dT ** 2))
                + np.sum((A - a - b * T) ** 2 / (sig ** 2 + dA ** 2 + b ** 2 * dT ** 2)))

def lnprior(theta):
    """Uniform priors on the parameters."""
    sig, a, b = theta

    if -30. < a < 30. and -5 < b < 1. and 0 < sig < 3.:
        return 0.
    return -np.inf

def lnprob(theta):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta)

ndim = 3
nwalkers = 12

optimal_set = [0.90, 1.5, -2.]
first_guess = [optimal_set + 0.1 * np.random.randn(ndim) * optimal_set for i in range(nwalkers)]

sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args = ())
sampler.run_mcmc(first_guess, 20000)

samples = sampler.chain[:, nwalkers:, :].reshape((-1, ndim))

fig = corner.corner(samples,
labels = [r"$\sigma$", r"$a$", r"$b$", ],
color = "blue", quantiles = [0.32, 0.5, 0.68], show_titles=True,
range=[0.95, 0.95, 0.95])

fig.savefig("plots/plateaus/emcee.pdf", bbox_inches="tight")
