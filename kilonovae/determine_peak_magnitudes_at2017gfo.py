import json
import numpy as np
from helpers import *
from numpy import cos

GW170817_FILE = "data/GW170817.json"
theta_170817 = 15 * Deg
dm_170817 = 33.15

with open(GW170817_FILE) as json_file:
    data = json.load(json_file)

bands_d = ['g', 'r', 'i', 'z']
delta_m_d = {'g': 7, 'r': 4, 'i': 3.5, 'z': 2.5}

optical_photometry = [dic for dic in data['GW170817']['photometry'] if 'band' in dic.keys() and 'model' not in dic.keys() and 'upperlimit' not in dic.keys()]

mags = dict()
for f in bands_d:
    mags[f] = np.array([float(dic['magnitude']) for dic in optical_photometry if dic['band'] == f])

max_mags = {f: np.min(mags[f]) for f in bands_d}

for f in bands_d:
    print(f"{f}: {np.min(mags[f]):.5g}")

for f in bands_d:
    print(f"M for {f}: {max_mags[f] - 2 * delta_m_d[f] * (1 - cos(theta_170817)) - dm_170817:.5g}")
