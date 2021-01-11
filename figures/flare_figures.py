from sys import exit
import matplotlib.pyplot as plt
from matplotlib import colors
import pandas as pd
from numpy import sin, cos, sqrt
from numpy import log10 as log
import numpy as np
from helpers import *
from sys import argv
from astropy.io import fits
PLOT_DIR = f"plots/flares"
plt.style.use('presentation')
cmap = colors.LinearSegmentedColormap.from_list('custom',rainbow(np.linspace(0, 1, 10)))

# Fixed parameters
Eiso_p = 1.e53
Eiso_los = Eiso_p / 4
tj = 0.1
tv = 0.13
a = 0.
b = -1.5
Ep = 0.1 * keV
nup = Ep / hPlanck
nuobs = 5.0 * keV / hPlanck
nug = 80. * keV / hPlanck

# plot position of patches
chis = []
radiis = []

# LOS material parameters
Gamma_los = 125
beta_los = sqrt(1 - 1 / Gamma_los ** 2)
r_los = 100 / Gamma_los
re_los = 1.e15
te_los = re_los / (beta_los * cLight)
tobs_los = te_los - re_los / cLight
tm_los = te_los - re_los / cLight
tM_los = te_los - re_los * cos(r_los) / cLight
print(f"LOS: r = {r_los}, tm = {tm_los}, tM = {tM_los}, R/2G^2c = {re_los / (2 * Gamma_los**2 * cLight):3.2f}")

# Plot LOS with plateau and patches
FILE_NAME = f"R_1"
with fits.open('data/plateau_03.fits') as hdul:
    data = hdul[1].data
    tobs = data['tobs']

L_los = XRT_c * np.array([simple_hle_patch(t, nuobs, 0, r_los, te_los, re_los, Gamma_los, Eiso_los, nup, a, b)
      for t in tobs])
L_plateau = data['dtheta=0.03']
plt.plot(tobs - tobs_los, L_los + L_plateau, color="black")
plt.plot(tobs - tobs_los, L_los, color="grey", linestyle="--", linewidth=0.8)
plt.plot(tobs - tobs_los, L_plateau, color="grey", linestyle="--",  linewidth=0.8)



# Plot the off-axis patches
# P1
Gamma = 250
beta = sqrt(1 - 1 / Gamma ** 2)
re = 4 * re_los
t_delay = 15
dia = 1 / Gamma
r = dia / 2
chi = 0.25 * tv
S = (1 - beta * cos(chi))/(1 - beta)
print(f"P1: G = {Gamma}, dia={dia * Gamma}/Gamma, chi={chi / tv:4.3f}*tv, re={re:3.2e}, td={t_delay:3.2f}, S={S:3.2f}")
print(f"P1: R/2G^2c = {re / (2 *Gamma**2 * cLight):3.2f}, ton = {t_delay + re / (2 *Gamma**2 * cLight):3.2f}")
radiis.append(r)
te = t_delay + re / (beta * cLight)
tm = te - re * cos(max(0., chi - r)) / cLight
tM = te - re * cos(min(Pi, chi + r)) / cLight
print(f"P1: Delta t / t = {(tM - tm) / (tm - tobs_los):3.2f}")
T = 10 ** np.linspace(log(tm), log(tM), 2000)
L_p = XRT_c * np.array([simple_hle_patch(t, nuobs, chi, r, te, re, Gamma, Eiso_p, nup, a, b)
      for t in T])
plt.plot(T - tobs_los, L_p, linewidth=0.8, linestyle = ls_l[0],
         color=cmap(((chi / tv) - 0.25) / (1.1-0.25)),
         label=r"P1")
t1mp = te - re / cLight - t_delay
t1Mp = te - re *cos(r)/ cLight - t_delay
T1p = np.linspace(t1mp, t1Mp, 2000)
L1_pp = BAT_c * np.array([simple_hle_patch(t, nug, 0, r, te - t_delay, re, Gamma, Eiso_p, nup, a, b)
      for t in T1p])
dia = 3 / Gamma
r = dia / 2
te = t_delay + re / (beta * cLight)
tm = te - re * cos(max(0., chi - r)) / cLight
tM = te - re * cos(min(Pi, chi + r)) / cLight
T = 10 ** np.linspace(log(tm), log(tM), 2000)
L_p = XRT_c * np.array([simple_hle_patch(t, nuobs, chi, r, te, re, Gamma, Eiso_p, nup, a, b)
      for t in T])
plt.plot(T - tobs_los, L_p, linewidth=0.8, linestyle = ls_l[1],
         color=cmap(((chi / tv) - 0.25) / (1.1-0.25)))

# P2
re = 2 * re_los
t_delay = 10
dia = 1 / Gamma
r = dia / 2
chi = 0.5 * tv
S = (1 - beta * cos(chi))/(1 - beta)
print(f"P2: G = {Gamma}, dia={dia * Gamma}/Gamma, chi={chi / tv:4.3f}*tv, re={re:3.2e}, td={t_delay:3.2f}, S={S:3.2f}")
print(f"P2: R/2G^2c = {re / (2 *Gamma**2 * cLight):3.2f}, ton = {t_delay + re / (2 *Gamma**2 * cLight):3.2f}")
chis.append(chi)
radiis.append(r)
te = t_delay + re / (beta * cLight)
tm = te - re * cos(max(0., chi - r)) / cLight
tM = te - re * cos(min(Pi, chi + r)) / cLight
print(f"P2: Delta t / t = {(tM - tm) / (tm - tobs_los):3.2f}")
T = 10 ** np.linspace(log(tm), log(tM), 2000)
L_p = XRT_c * np.array([simple_hle_patch(t, nuobs, chi, r, te, re, Gamma, Eiso_p, nup, a, b)
      for t in T])
plt.plot(T - tobs_los, L_p, linewidth=0.8, linestyle = ls_l[0],
         color=cmap(((chi / tv) - 0.25) / (1.1-0.25)),
         label=r"P2")
t2mp = te - re / cLight - t_delay
t2Mp = te - re *cos(r)/ cLight - t_delay
T2p = np.linspace(t2mp, t2Mp, 2000)
L2_pp = BAT_c * np.array([simple_hle_patch(t, nug, 0, r, te - t_delay, re, Gamma, Eiso_p, nup, a, b)
      for t in T2p])
dia = 3 / Gamma
r = dia / 2
te = t_delay + re / (beta * cLight)
tm = te - re * cos(max(0., chi - r))/ cLight
tM = te - re * cos(min(Pi, chi + r)) / cLight
T = 10 ** np.linspace(log(tm), log(tM), 2000)
L_p = XRT_c * np.array([simple_hle_patch(t, nuobs, chi, r, te, re, Gamma, Eiso_p, nup, a, b)
      for t in T])
plt.plot(T - tobs_los, L_p, linewidth=0.8, linestyle = ls_l[1],
         color=cmap(((chi / tv) - 0.25) / (1.1-0.25)))

# P3
re = 0.5 * re_los
t_delay = 150
dia = 1 / Gamma
r = dia / 2
chi = 0.6 * tv
S = (1 - beta * cos(chi))/(1 - beta)
print(f"P3: G = {Gamma}, dia={dia * Gamma}/Gamma, chi={chi / tv:4.3f}*tv, re={re:3.2e}, td={t_delay:3.2f}, S={S:3.2f}")
print(f"P3: R/2G^2c = {re / (2 *Gamma**2 * cLight):3.2f}, ton = {t_delay + re / (2 *Gamma**2 * cLight):3.2f}")
chis.append(chi)
radiis.append(r)
te = t_delay + re / (beta * cLight)
tm = te - re * cos(max(0., chi - r)) / cLight
tM = te - re * cos(min(Pi, chi + r)) / cLight
print(f"P3: Delta t / t = {(tM - tm) / (tm - tobs_los):3.2f}")
T = 10 ** np.linspace(log(tm), log(tM), 2000)
L_p = XRT_c * np.array([simple_hle_patch(t, nuobs, chi, r, te, re, Gamma, Eiso_p, nup, a, b)
      for t in T])
plt.plot(T - tobs_los, L_p, linewidth=0.8, linestyle = ls_l[0],
         color=cmap(((chi/tv) - 0.25) / (1.1-0.25)),
         label=r"P3")
t3mp = te - re / cLight - t_delay
t3Mp = te - re *cos(r)/ cLight - t_delay
T3p = np.linspace(t3mp, t3Mp, 2000)
L3_pp = BAT_c * np.array([simple_hle_patch(t, nug, 0, r, te - t_delay, re, Gamma, Eiso_p, nup, a, b)
      for t in T3p])
dia = 3 / Gamma
r = dia / 2
te = t_delay + re / (beta * cLight)
tm = te - re * cos(max(0., chi - r)) / cLight
tM = te - re * cos(min(Pi, chi + r)) / cLight
T = 10 ** np.linspace(log(tm), log(tM), 2000)
L_p = XRT_c * np.array([simple_hle_patch(t, nuobs, chi, r, te, re, Gamma, Eiso_p, nup, a, b)
      for t in T])
plt.plot(T - tobs_los, L_p, linewidth=0.8, linestyle = ls_l[1],
         color=cmap(((chi/tv) - 0.25) / (1.1-0.25)))

# P4
Gamma = 250
beta = sqrt(1 - 1 / Gamma**2)
re = 10 * re_los
t_delay = 80
dia = 1 / Gamma
r = dia / 2
chi = 0.35*tv
S = (1 - beta * cos(chi))/(1 - beta)
print(f"P4: G = {Gamma}, dia={dia * Gamma}/Gamma, chi={chi / tv:4.3f}*tv, re={re:3.2e}, td={t_delay:3.2f}, S={S:3.2f}")
print(f"P4: R/2G^2c = {re / (2 *Gamma**2 * cLight):3.2f}, ton = {t_delay + re / (2 *Gamma**2 * cLight):3.2f}")
chis.append(chi)
radiis.append(r)
te = t_delay + re / (beta * cLight)
tm = te - re * cos(max(0., chi - r)) / cLight
tM = te - re * cos(min(Pi, chi + r)) / cLight
print(f"P4: Delta t / t = {(tM - tm) / (tm - tobs_los):3.2f}")
T = 10 ** np.linspace(log(tm), log(tM), 2000)
L_p = XRT_c * np.array([simple_hle_patch(t, nuobs, chi, r, te, re, Gamma, Eiso_p, nup, a, b)
      for t in T])
plt.plot(T - tobs_los, L_p, linewidth=0.8, linestyle = ls_l[0],
         color=cmap(((chi/tv) - 0.25) / (1.1-0.25)),
         label=r"P4")
t4mp = te - re / cLight - t_delay
t4Mp = te - re *cos(r)/ cLight - t_delay
T4p = np.linspace(t4mp, t4Mp, 2000)
L4_pp = BAT_c * np.array([simple_hle_patch(t, nug, 0, r, te - t_delay, re, Gamma, Eiso_p, nup, a, b)
      for t in T4p])
dia = 3 / Gamma
r = dia / 2
te = t_delay + re / (beta * cLight)
tm = te - re * cos(max(0., chi - r)) / cLight
tM = te - re * cos(min(Pi, chi + r)) / cLight
T = 10 ** np.linspace(log(tm), log(tM), 2000)
L_p = XRT_c * np.array([simple_hle_patch(t, nuobs, chi, r, te, re, Gamma, Eiso_p, nup, a, b)
      for t in T])
plt.plot(T - tobs_los, L_p, linewidth=0.8, linestyle = ls_l[1],
         color=cmap(((chi/tv) - 0.25) / (1.1-0.25)))

# P5
Gamma = 125
beta = sqrt(1 - 1 / Gamma**2)
re = 15 * re_los
t_delay = 80
dia = 1 / Gamma
r = dia / 2
chi = 0.7*tv
S = (1 - beta * cos(chi))/(1 - beta)
print(f"P5: G = {Gamma}, dia={dia * Gamma}/Gamma, chi={chi / tv:4.3f}*tv, re={re:3.2e}, td={t_delay:3.2f}, S={S:3.2f}")
print(f"P5: R/2G^2c = {re / (2 *Gamma**2 * cLight):3.2f}, ton = {t_delay + re / (2 *Gamma**2 * cLight):3.2f}")
chis.append(chi)
radiis.append(r)
te = t_delay + re / (beta * cLight)
tm = te - re * cos(max(0., chi - r)) / cLight
tM = te - re * cos(min(Pi, chi + r)) / cLight
print(f"P5: Delta t / t = {(tM - tm) / (tm - tobs_los):3.2f}")
T = 10 ** np.linspace(log(tm), log(tM), 2000)
L_p = XRT_c * np.array([simple_hle_patch(t, nuobs, chi, r, te, re, Gamma, Eiso_p, nup, a, b)
      for t in T])
plt.plot(T - tobs_los, L_p, linewidth=0.8, linestyle = ls_l[0],
         color=cmap(((chi/tv) - 0.25) / (1.1-0.25)),
         label=r"P5")
t5mp = te - re / cLight - t_delay
t5Mp = te - re *cos(r)/ cLight - t_delay
T5p = np.linspace(t5mp, t5Mp, 2000)
L5_pp = BAT_c * np.array([simple_hle_patch(t, nug, 0, r, te - t_delay, re, Gamma, Eiso_p, nup, a, b)
      for t in T5p])
dia = 3 / Gamma
r = dia / 2
te = t_delay + re / (beta * cLight)
tm = te - re * cos(max(0., chi - r)) / cLight
tM = te - re * cos(min(Pi, chi + r)) / cLight
T = 10 ** np.linspace(log(tm), log(tM), 2000)
L_p = XRT_c * np.array([simple_hle_patch(t, nuobs, chi, r, te, re, Gamma, Eiso_p, nup, a, b)
      for t in T])
plt.plot(T - tobs_los, L_p, linewidth=0.8, linestyle = ls_l[1],
         color=cmap(((chi/tv) - 0.25) / (1.1-0.25)))

# P6
Gamma = 250
beta = sqrt(1 - 1 / Gamma**2)
re = 0.5 * re_los
t_delay = 100
dia = 1 / Gamma
r = dia / 2
chi = 0.55 * tv
S = (1 - beta * cos(chi))/(1 - beta)
print(f"P6: G = {Gamma}, dia={dia * Gamma}/Gamma, chi={chi / tv:4.3f}*tv, re={re:3.2e}, td={t_delay:3.2f}, S={S:3.2f}")
print(f"P6: R/2G^2c = {re / (2 *Gamma**2 * cLight):3.2f}, ton = {t_delay + re / (2 *Gamma**2 * cLight):3.2f}")
chis.append(chi)
radiis.append(r)
te = t_delay + re / (beta * cLight)
tm = te - re * cos(max(0., chi - r)) / cLight
tM = te - re * cos(min(Pi, chi + r)) / cLight
print(f"P6: Delta t / t = {(tM - tm) / (tm - tobs_los):3.2f}")
T = 10 ** np.linspace(log(tm), log(tM), 2000)
L_p = XRT_c * np.array([simple_hle_patch(t, nuobs, chi, r, te, re, Gamma, Eiso_p, nup, a, b)
      for t in T])
plt.plot(T - tobs_los, L_p, linewidth=0.8, linestyle = ls_l[0],
         color=cmap(((chi/tv) - 0.25) / (1.1-0.25)),
         label=r"P6")
t6mp = te - re / cLight - t_delay
t6Mp = te - re *cos(r)/ cLight - t_delay
T6p = np.linspace(t6mp, t6Mp, 2000)
L6_pp = BAT_c * np.array([simple_hle_patch(t, nug, 0, r, te - t_delay, re, Gamma, Eiso_p, nup, a, b)
      for t in T6p])
dia = 3 / Gamma
r = dia / 2
te = t_delay + re / (beta * cLight)
tm = te - re * cos(max(0., chi - r)) / cLight
tM = te - re * cos(min(Pi, chi + r)) / cLight
T = 10 ** np.linspace(log(tm), log(tM), 2000)
L_p = XRT_c * np.array([simple_hle_patch(t, nuobs, chi, r, te, re, Gamma, Eiso_p, nup, a, b)
      for t in T])
plt.plot(T - tobs_los, L_p, linewidth=0.8, linestyle = ls_l[1],
         color=cmap(((chi/tv) - 0.25) / (1.1-0.25)))

# P7
Gamma = 125
beta = sqrt(1 - 1 / Gamma**2)
re = 3 * re_los
t_delay = 275 #580100
dia = 1 / Gamma
r = dia / 2
chi = 1.1 * tv
S = (1 - beta * cos(chi))/(1 - beta)
print(f"P7: G = {Gamma}, dia={dia * Gamma}/Gamma, chi={chi / tv:4.3f}*tv, re={re:3.2e}, td={t_delay:3.2f}, S={S:3.2f}")
print(f"P7: R/2G^2c = {re / (2 *Gamma**2 * cLight):3.2f}, ton = {t_delay + re / (2 *Gamma**2 * cLight):3.2f}")
chis.append(chi)
radiis.append(r)
te = t_delay + re / (beta * cLight)
tm = te - re * cos(max(0., chi - r)) / cLight
tM = te - re * cos(min(Pi, chi + r)) / cLight
print(f"P7: Delta t / t = {(tM - tm) / (tm - tobs_los):3.2f}")
T = 10 ** np.linspace(log(tm), log(tM), 2000)
L_p = XRT_c * np.array([simple_hle_patch(t, nuobs, chi, r, te, re, Gamma, Eiso_p, nup, a, b)
      for t in T])
plt.plot(T - tobs_los, L_p, linewidth=0.8, linestyle = ls_l[0],
         color=cmap(((chi/tv) - 0.25) / (1.1-0.25)),
         label=r"P7")
t7mp = te - re  / cLight - t_delay
t7Mp = te - re  *cos(r)/ cLight - t_delay
T7p = np.linspace(t7mp, t7Mp, 2000)
L7_pp = BAT_c * np.array([simple_hle_patch(t, nug, 0, r, te - t_delay, re, Gamma, Eiso_p, nup, a, b)
      for t in T7p])
dia = 3 / Gamma
r = dia / 2
te = t_delay + re / (beta * cLight)
tm = te - re * cos(max(0., chi - r)) / cLight
tM = te - re * cos(min(Pi, chi + r)) / cLight
T = 10 ** np.linspace(log(tm), log(tM), 2000)
L_p = XRT_c * np.array([simple_hle_patch(t, nuobs, chi, r, te, re, Gamma, Eiso_p, nup, a, b)
      for t in T])
plt.plot(T - tobs_los, L_p, linewidth=0.8, linestyle = ls_l[1],
         color=cmap(((chi/tv) - 0.25) / (1.1-0.25)))

plt.ylabel(f"Luminosity [0.3-30] keV (erg/s)")
plt.xlabel("Time since last prompt pulse (s)")
plt.xscale('log')
plt.yscale('log')
plt.ylim([1.e45, 5.e51])
plt.xlim([50, 40000])
plt.legend(loc='upper right')
plt.savefig(f"{PLOT_DIR}/{FILE_NAME}.pdf", bbox_inches='tight')

plt.figure()
fig = plt.gcf()
ax = fig.gca()
plt.gca().set_aspect('equal', adjustable='box')
ax.add_artist(plt.Circle((0, 0), tj, edgecolor='black', fill=False))
ax.add_artist(plt.Circle((0, tj / 3), 0.5 / 250, edgecolor='grey', fill=False))
ax.add_artist(plt.Circle((0, 2 * tj / 3), 0.5 / 125, edgecolor='grey', fill=False))

#plot patches
for j, ch in enumerate(chis):
    ax.add_artist(plt.Circle((ch - tv, 0), radiis[j], color=cmap(((ch/tv) - 0.25) / (1.1-0.25))))

plt.scatter([-tv], [0], marker = 'x')
plt.xlim([-1.1 * tv, tv])
plt.ylim([-tv, tv])
plt.savefig(f"{PLOT_DIR}/positions.pdf")

plt.figure()
plt.plot(T1p, L1_pp, label = "P1", color = cmap((0.25-0.25)/(1.1-0.25)))
plt.plot(T2p, L2_pp, label = "P2", color = cmap((0.5-0.25)/(1.1-0.25)))
plt.plot(T3p, L3_pp, label = "P3", color = cmap((0.6-0.25)/(1.1-0.25)))
plt.plot(T4p, L4_pp, label = "P4", color = cmap((0.35-0.25)/(1.1-0.25)))
plt.plot(T5p, L5_pp, label = "P5", color = cmap((0.7-0.25)/(1.1-0.25)))
plt.plot(T6p, L6_pp, label = "P6", color = cmap((0.55-0.25)/(1.1-0.25)))
plt.plot(T7p, L7_pp, label = "P7", color = cmap((1.1-0.25)/(1.1-0.25)))
plt.legend(loc="upper right")
plt.yscale("log")
plt.xlabel("Time (s)")
plt.ylabel("Luminosity [15-150] keV (erg/s)")
plt.savefig(f"{PLOT_DIR}/prompt.pdf")
