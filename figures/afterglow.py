#!/bin/env python3
"""Plot the after-glow data in plots directory."""

from sys import argv, exit
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from os import mkdir
from os.path import exists
from helpers import *
import configparser as cp

# Globals.
CONF_FILE = "data/%s.conf" % argv[1]
DATA_FILE = "data/%s.out" % argv[1]
CASE = argv[1]
PLOTS_DIR = f"plots/{CASE}/"


# Initialize conf dict
conf_list = dict()

#Parse configuration data into conf_list
with open(CONF_FILE, 'r') as conf_file:
    for line in conf_file:
        if line[0] == '#':
            continue
        bits = line.split()
        conf_list[bits[0]] = float(bits[-1])
print("Configuration: " + str(conf_list))

df = pd.read_csv(DATA_FILE, sep=" ", names=['r', 't', 't_obs', 'gam', 'rho_s', 'eps_s',
'Niso', 'epsilon_e', 'epsilon_b', 'p', 'beta', 'nu_obs', 'F_nu_obs',
'nu_a', 'nu_i', 'nu_c', 'gam_i', 'gam_c', 'dopp'])

# Plot away!
# Light curve
print("Plotting light curve.")
plt.figure(4)
plt.plot(df['t_obs']/day, df['F_nu_obs']/microJy, 'g--', label=r'$F_{\nu_{obs}}$')
plt.errorbar(OBS_TIMES_d, OBS_FLUXES_muJy, yerr = OBS_ERRORS_muJy, fmt="b+")
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Observer Time (days)')
plt.ylabel(r'$S_{\nu} (\mu Jy)$')
plt.title('Light Curve (at %s GHz)' % str(1.e-9 * conf_list['nu_obs']))
plt.grid(True, which='both')
plt.savefig(PLOTS_DIR + CASE + "_light_curve.pdf", bbox_inches='tight')

# Gamma plot
print("Plotting Gamma(t).")
plt.figure(1)
plt.plot(df['t_obs']/day, df['gam'], 'r-', label=r'$\Gamma$')
plt.legend()
plt.xscale('log')
plt.xlabel('time (days)')
plt.ylabel('$\Gamma$')
plt.title('Lorentz Factor')
plt.grid(True, which='both')
plt.savefig(PLOTS_DIR + CASE + '_Gamma_t.pdf', bbox_inches='tight')

# Gamma plot
print("Plotting Gamma(r).")
plt.figure(45)
plt.plot(df['r'], df['gam'], 'r-', label=r'$\Gamma$')
plt.legend()
plt.xscale('log')
plt.xlabel('radius (cm)')
plt.ylabel('$\Gamma$')
plt.title('Lorentz Factor')
plt.grid(True, which='both')
plt.savefig(PLOTS_DIR + CASE + '_Gamma_r.pdf', bbox_inches='tight')

# Gamma plot
print("Plotting u(r).")
plt.figure(87)
plt.plot(df['r'], df['gam'] * df['beta'], 'r-', label=r'$\Gamma\beta$')
plt.legend()
plt.xscale('log')
plt.xlabel('radius (cm)')
plt.ylabel('$\Gamma \beta$')
plt.title('Lorentz Factor')
plt.grid(True, which='both')
plt.savefig(PLOTS_DIR + CASE + '_u_r.pdf', bbox_inches='tight')

# a(t) plot
print("Plotting dopp(t).")
plt.figure(56)
plt.plot(df['t_obs']/day, df['dopp'], 'r-', label=r'$\delta$')
plt.legend()
plt.xscale('log')
plt.xlabel('time (days)')
plt.ylabel(r'$\delta$')
plt.title('Doppler factor')
plt.grid(True, which='both')
plt.savefig(PLOTS_DIR + CASE + '_dopp.pdf', bbox_inches='tight')

# Beta plot
print("Plotting 1 - beta.")
plt.figure(2)
plt.plot(df['t_obs']/day, 1 - df['beta'], 'b-', label=r'\beta')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('time (days)')
plt.ylabel(r'$1 - \beta$')
plt.title('1 - Beta factor')
plt.grid(True, which='both')
plt.savefig(PLOTS_DIR + CASE + '_beta.pdf', bbox_inches='tight')

# Characteristic frequencies
print("Plotting frequencies.")
plt.figure(6)
plt.plot(df['t_obs']/day, df['nu_i'] / GHz, 'b--', label=r"$\nu_i'$")
plt.plot(df['t_obs']/day, df['nu_c'] / GHz, 'g--', label=r"$\nu_c'$")
plt.plot(df['t_obs']/day, df['nu_a'] / GHz, 'r--', label=r"$\nu_a'$")
plt.plot(df['t_obs']/day, df['nu_obs'] / GHz, color = 'grey', label=r"$\nu_{obs}$")
plt.plot(df['t_obs']/day, df['nu_obs'] / (df['dopp'] * GHz), color = 'black',
        label=r"$\nu_{obs}'$")
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Observer Time (days)')
plt.ylabel('Frequency (GHz)')
plt.title('Characteristic Frequencies')
plt.grid(True, which='both')
plt.savefig(PLOTS_DIR + CASE + "_frequencies.pdf", bbox_inches='tight')

# Characteristic Lorentz factors
print("Plotting synchroton characteristic Lorentz factors.")
plt.figure(7)
plt.plot(df['t_obs']/day, df['gam_i'], 'b--', label='$\gamma_i$')
plt.plot(df['t_obs']/day, df['gam_c'], 'g--', label="$\gamma_c$")
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Observer Time (days)')
plt.ylabel('Gamma')
plt.title('Characteristic Lorentz Factors')
plt.grid(True, which='both')
plt.savefig(PLOTS_DIR + CASE + "_gammas.pdf", bbox_inches='tight')

# Indices
print("Plotting indices.")
plt.figure(8)
plt.plot(df['t_obs']/day, d_log_log(df['F_nu_obs'], df['t_obs']),
        label = r'$d(log S_{\nu}) / d(log t_{obs})$')
plt.plot(df['t_obs']/day, d_log_log(df['gam'], df['r']),
        label = r'$d(log \Gamma) / d(log r)$')
plt.plot(df['t_obs']/day, d_log_log(df['beta'], df['r']),
        label = r'$d(log \beta) / d(log r)$')
plt.plot(df['t_obs']/day, d_log_log(df['nu_i'], df['t_obs']),
        label = r'$d(log \nu_i) / d(log t_{obs})$')
plt.plot(df['t_obs']/day, d_log_log(df['nu_c'], df['t_obs']),
        label = r'$d(log \nu_c) / d(log t_{obs})$')
plt.plot(df['t_obs']/day, d_log_log(df['r'], df['t_obs']),
        label = r'$d(log r) / d(log t_{obs})$')
plt.plot(df['t_obs']/day, d_log_log(df['gam'], df['t_obs']),
        label = r'$d(log \Gamma) / d(log t_{obs})$')
plt.legend(loc="upper right", fontsize=8)
plt.xscale('log')
plt.yscale('linear')
plt.xlabel('Observer Time (days)')
plt.title('Indices')
plt.grid(True, which='both')
plt.savefig(PLOTS_DIR + CASE + "_indices.pdf", bbox_inches='tight')
