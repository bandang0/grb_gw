#!/bin/env python3
"""Plot the after-glow data in a perty manner."""

from sys import argv
from matplotlib import pyplot as plt

from helpers import *

ldic = {'cocoon_best': "Sph. structure radiale",
        'shell_best': "Sph. coquille",
        'jet_best': "Jet hors-axe",
        'some_jet' : "Jet hors-axe"}
cdic = {'cocoon_best': "red",
        'shell_best': "orange",
        'jet_best': "green",
        'some_jet': "blue"}

def plot(plt, l_string):
    plt.figure(sum([len(string) for string in l_string]))
    plt.errorbar(OBS_TIMES, OBS_FLUXES, yerr = OBS_ERRORS, fmt="+", color = "grey")
    for string in l_string:
        t_obs, F_nu_obs = data_list("data/{}.out".format(string))
        t_obs_days = t_obs / day
        F_nu_obs_mJy = F_nu_obs * 1.e26
        plt.plot(t_obs_days, F_nu_obs_mJy, label=ldic[string], color=cdic[string])
    plt.legend()
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim([10., 600.])
    plt.ylim([0.01, 0.3])
    plt.xlabel('Temps observateur (jours)')
    plt.ylabel(r'$S_{\nu} \times (D/40Mpc)^{-2}$ (mJy)')
    plt.title('Courbe de lumière @ 3 GHz')
    plt.grid(True, which='major')
    plt.savefig("misc/{}.pdf".format(str(sum([len(string) for string in l_string]))), bbox_inches='tight')

if __name__ == "__main__":
    print("Plotting light curve.")
    plt.figure(0)
    plt.errorbar(OBS_TIMES, OBS_FLUXES, yerr = OBS_ERRORS, fmt="+", color = "grey")
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim([10., 600.])
    plt.ylim([0.01, 0.3])
    plt.xlabel('Temps observateur (jours)')
    plt.ylabel(r'$S_{\nu} \times (D/40Mpc)^{-2}$ (mJy)')
    plt.title('Courbe de lumière @ 3 GHz')
    plt.grid(True, which='major')
    plt.savefig("misc/empty.pdf", bbox_inches='tight')

    plot(plt, ["shell_best"])
    plot(plt, ["cocoon_best"])
    plot(plt, ["jet_best"])
    plot(plt, ["some_jet", "cocoon_best"])
    plot(plt, ["shell_best", "cocoon_best", "jet_best"])
