#!/bin/env python3
"""Plot Daigne output data."""

import numpy as np
from matplotlib import pyplot as plt

from  helpers import *

# Globals.
DATA_FILE = "data/testRD_1.dat_1_iso.spn.lc.simple"
CASE = "daigne"
PLOTS_DIR = "plots/"


conv_vector = np.ones(11)/11.
still_vector = np.array([0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0.])
#still_vector = np.array([0., 1., 0.])
smooth = lambda x: np.convolve(x, conv_vector, mode='valid')
still = lambda x: np.convolve(x, still_vector, mode='valid')

if __name__ == "__main__":
    # Initialize data lists
    t_obs = list()
    Radio_flux_Jy = list()
    R_flux_Jy = list()
    X_flux_Jy = list()
    obj_list = [t_obs, R_flux_Jy, X_flux_Jy, Radio_flux_Jy]

    # Parse fname into the listed data, and read t_obs_break
    with open(DATA_FILE, 'r') as data_file:
        lines = data_file.readlines()
        for l in range(2210, 4630):
            line = lines[l]
            if '0.00000E+00' in line:
                continue
            atoms = line.split()
            for i, obj in enumerate(obj_list):
                obj_list[i].append(float(atoms[i]))

    # np.arrize
    t_obs = np.array(t_obs)
    Radio_flux_Jy = np.array(Radio_flux_Jy)
    R_flux_Jy = np.array(R_flux_Jy)
    X_flux_Jy = np.array(X_flux_Jy)
    N_DATA_POINTS = len(t_obs)

    # Better units
    Radio_flux_mJy = 1.e3 * Radio_flux_Jy
    R_flux_mJy = 1.e3 * R_flux_Jy
    X_flux_mJy = 1.e3 * X_flux_Jy
    t_obs_days = t_obs / day

    #smooth
    #Radio_flux_mJy_s = np.convolve(Radio_flux_mJy, conv_vector, mode = 'valid')
    #R_flux_mJy_s = np.convolve(R_flux_mJy, conv_vector, mode = 'valid')
    #X_flux_mJy_s = np.convolve(X_flux_mJy, conv_vector, mode = 'valid')
    #t_obs_days_s = np.convolve(t_obs, still_vector, mode = 'valid')

    # Light curve
    print("Plotting light curve.")
    plt.figure(1)
    #plt.plot(t_obs_days, Radio_flux_mJy, color='black', label=r'$F_{Radio}$')
    plt.plot(still(t_obs_days), smooth(Radio_flux_mJy), 'g.--', label=r'$F_{Radio}$')
    plt.plot(still(t_obs_days), smooth(R_flux_mJy), 'r--', label=r'$F_{R}$')
    plt.plot(still(t_obs_days), smooth(X_flux_mJy), 'b--', label=r'$F_{X}$')

    plt.legend()
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Observer Time (days)')
    plt.ylabel(r'$S_{\nu} \times (D/40Mpc)^{-2}$ (mJy)')
    plt.title('Light Curve (at 5 GHz)')
    plt.grid(True, which='major')
    plt.savefig(PLOTS_DIR + CASE + "_light_curve.pdf", bbox_inches='tight')

    # Indices
    print("Plotting indices.")
    plt.figure(6)
    #plt.plot(t_obs_days, d_log_log(Radio_flux_mJy, t_obs_days), color = 'black', label = r'$d(log S_{Radio}) / d(log t_{obs})$')
    plt.plot(still(still(t_obs_days)), smooth(d_log_log(smooth(Radio_flux_mJy), smooth(t_obs_days))), label = r'$d(log S_{Radio}) / d(log t_{obs})$')
    plt.plot(still(still(t_obs_days)), smooth(d_log_log(smooth(R_flux_mJy), smooth(t_obs_days))), label = r'$d(log S_{R}) / d(log t_{obs})$')
    plt.plot(still(still(t_obs_days)), smooth(d_log_log(smooth(X_flux_mJy), smooth(t_obs_days))), label = r'$d(log S_{X}) / d(log t_{obs})$')

    plt.legend(loc="upper right", fontsize=8)
    plt.xscale('log')
    plt.yscale('linear')
    plt.xlabel('Observer Time (days)')
    plt.title('Indices')
    plt.grid(True, which='major')
    plt.savefig(PLOTS_DIR + CASE + "_indices.pdf", bbox_inches='tight')
    plt.show()
