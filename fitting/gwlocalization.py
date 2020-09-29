import numpy as np
import matplotlib.pyplot as plt
from numpy import sin, cos, arccos, pi, sqrt
from numpy.random import random as u
import pandas as pd
from sys import exit

N = 409600
rho_0 = 8.
H = 200.
thetam = 4.

def theta2(i, ksi, theta, phi):
    fp = 0.5 * cos(2. * ksi) * (1. + cos(theta) ** 2) * cos(2. * phi)\
    - sin(2. * ksi) * cos(theta) * sin(2. * phi)
    fc = 0.5 * sin(2. * ksi) * (1 + cos(theta) ** 2) * cos(2. * phi)\
    + cos(2. * ksi) * cos(theta) * sin(2. * phi)
    return 4. * (fp ** 2 * (1. + cos(i) ** 2) ** 2 + 4. * fc ** 2 * cos(i) ** 2)

def theta2avg(i):
    return 4. * (1. + 6. * cos(i) ** 2 + cos(i) ** 4) / 5.

Hl = 130. * 2.26
Hh = 100. * 2.26
Hv = 50. * 2.26

N = 8192
M = 40
distances = np.linspace(Hl / M, Hl, M)
Nl = np.zeros(M)
Nh = np.zeros(M)

for d in range(len(distances)):
    for j in range(N):
        i = arccos(1 - 2 * u())
        theta = arccos(1 - 2 * u())
        phi = 2 * pi * u()
        ksi = 2 * pi * u()
        rh = rho_0 * Hh * sqrt(theta2(i, ksi, theta, phi)) / (distances[d] * thetam)
        rl = rho_0 * Hl * sqrt(theta2(i, ksi, theta, phi)) / (distances[d] * thetam)
        if rh > rho_0:
            Nh[d] = Nh[d] + 1
        if rl > rho_0:
            Nl[d] = Nl[d] + 1

plt.figure()
plt.plot(distances, Nh / N, label="Hanford")
plt.plot(distances, Nl / N, label="Livingston")
plt.plot(distances, (Nl - Nh) / N, label="Diff")

plt.legend()
plt.show()

exit(0)
d = list()
rhol2 = list()
rhov2 = list()
rhoh2 = list()
#theta = list()
#phi = list()
#ksi = list()
i = list()
s = list()

for j in range(N):
    d.append((1. / 1.58) * Hh * u() ** (1. / 3.))
    ia = arccos(u())
    i.append(ia)
    t =  theta2avg(ia)
    rh2 = rho_0 ** 2 * Hh ** 2 * t / (d[j] ** 2 * thetam ** 2)
    rl2 = rho_0 ** 2 * Hl ** 2 * t / (d[j] ** 2 * thetam ** 2)
    rv2 = rho_0 ** 2 * Hv ** 2 * t / (d[j] ** 2 * thetam ** 2)

    rhol2.append(rl2)
    rhov2.append(rv2)
    rhoh2.append(rh2)

    s.append(30 * 32.4 ** 2 / (rl2  + rv2 + rh2))

all = pd.DataFrame({'d': d,
                    'i': i,
                    'rl2': rhol2,
                    'rh2': rhoh2,
                    'rv2': rhov2,
                    's': s})

print(all)
gw = all.loc[all.rh2 > 8. ** 2]

plt.figure()
plt.hist([gw['d'], all['d']], label=["gw", 'all'],
density=True, histtype='step', bins=100)
plt.xlabel("D / Mpc")
plt.legend()

print(gw)
plt.figure()
plt.hist([gw['s']], label=["gw"],
cumulative=True, histtype='step', bins=70)
plt.xlabel("S / deg2")
plt.legend()
plt.show()
