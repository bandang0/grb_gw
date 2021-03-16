import json
import numpy as np
from helpers import *
from numpy import cos
from sys import exit
from matplotlib import pyplot as plt
from astropy.time import Time

GW170817_FILE = "data/GW170817.json"
theta_170817 = 15 * Deg
dm_170817 = 33.15
PLOT_DIR = f"plots/kn_plots"

with open(GW170817_FILE) as json_file:
    data = json.load(json_file)

colors_d = {'g': 'green', 'r': 'red', 'i': 'black', 'z': 'blue', 'y': 'grey', 'J': 'orange', 'H': 'pink', 'K': 'yellow'}
bands_d = ['g', 'r', 'i', 'z', 'y', 'J', 'H', 'K']
top_d = {'g': -12.1, 'r': -12.0, 'i': -12.3, 'z': -12.4, 'y': -12.4, 'J': -12.5, 'H': -13.9, 'K': -14.6}
bottom_d = {'g': -13.5, 'r': -13.3, 'i': -13.5, 'z': -13.5, 'y': -13.4, 'J': -13.3, 'H': -13.9, 'K': -14.6}
side_d = {'g': -6.4, 'r': -8.8, 'i': -10., 'z': -10.3, 'y': -10.6, 'J': -11.4, 'H': -13.2, 'K': -14.}

polar_d = {f: 0.5 * (top_d[f] + bottom_d[f]) for f in bands_d}
delta_d = {f: - polar_d[f] + side_d[f] for f in bands_d}

optical_photometry = [dic for dic in data['GW170817']['photometry'] if 'band' in dic and 'model' not in dic and 'upperlimit' not in dic]

light_curve = dict()
for f in bands_d:
    light_curve[f] = np.array([[float(dic['time']), float(dic['magnitude']), float(dic.get('e_magnitude', 0.2))] for dic in optical_photometry if dic['band'] == f])

min_mags = {f: np.min(light_curve[f][:,1]) for f in bands_d}
m_d = {f: min_mags[f] - 2 * delta_d[f] * (1 - cos(theta_170817)) - dm_170817 for f in bands_d}

for f in bands_d:
    print(f"{f}: {min_mags[f]:.5g}, {m_d[f]:.5g}, {delta_d[f]:.5g}")


merger_t = Time('2017-08-17 12:41:04.4', scale='utc')

plt.figure()
for f in bands_d:
    lc_t = Time(light_curve[f][:,0], format='mjd')
    plt.scatter(lc_t.mjd - merger_t.mjd, light_curve[f][:,1], marker = '.', color=colors_d[f])
    plt.errorbar(np.array(lc_t.mjd) - merger_t.mjd, light_curve[f][:,1], yerr=light_curve[f][:,2], fmt='none', color=colors_d[f], label=f)
plt.xlabel(r'time since merger [d]')
plt.ylabel('AB mag')
plt.legend()
plt.ylim([22.5, 16.5])
plt.xlim([0, 15])
plt.savefig(f"{PLOT_DIR}/AT2017_gfo.pdf", bbox_inches='tight')

theta_max = 60 * Deg
theta_l = np.linspace(0, 90, 200) * Deg
plt.figure()
for f in bands_d:
    Y = np.array([m_d[f] + delta_d[f] * ((1 - cos(min(x, theta_max))) / 0.5) for x in theta_l])
    plt.plot(theta_l / Deg, Y, color=colors_d[f], label=f)
plt.xlabel(r'$\theta_v$ [deg]')
plt.ylabel('peak absolute magnitude')
plt.legend()
plt.ylim([-9.5, -17])
plt.xlim([0, 90])
plt.savefig(f"{PLOT_DIR}/kn_model.pdf", bbox_inches='tight')
