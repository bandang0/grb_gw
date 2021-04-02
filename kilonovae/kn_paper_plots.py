'''Make plots for kilonova paper.'''

from sys import exit, argv
import sys

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib import gridspec
import pandas as pd
from numpy import sin, cos
import numpy as np

from helpers import *

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

# Setup
DATA_FILE_O3 = "data/kn_O3.data"
DATA_FILE_O4 = "data/kn_O4.data"
DATA_FILE_O5 = "data/kn_O5.data"
DATA_FILE_MAX = "data/kn_max.data"
DATA_FILE_190425 = "data/kn_190425.data"
PLOT_DIR = f"plots/kn_plots"

if 'test' in argv:
    print("Testing mode")
    DATA_FILE_O3 = "data/kn_O3_test.data"
    DATA_FILE_O4 = "data/kn_O4_test.data"
    DATA_FILE_O5 = "data/kn_O5_test.data"
    DATA_FILE_MAX = "data/kn_max_test.data"
    DATA_FILE_190425 = "data/kn_190425_test.data"

# Statistical resolution
BINS = 100
thin = 1.


# Bands and printing constants
col_names = ['d', 'tv', 'g0', 'dg', 'r0', 'dr', 'i0', 'di', 'z0', 'dz', 'r1', 'r10']
dtype = np.float64
bands_d = ['g', 'r', 'i', 'z']
mag_lim_d = {'g': 21, 'r': 21, 'i': 21.5}
cmap_d = {'g': 'Greens', 'r': 'Reds', 'i': 'Greys'}
c_d = {'g': 'green', 'r': 'red', 'i': 'black', "z": 'blue'}
ls_d = {'g': '-', 'r': '-.', 'i': ":", 'z': '--'}
tex_d = {'g': r'$g$', 'r': r'$r$', 'i': r'$i$', 'z': r'$z$'}

# Constants
g0 = 100.
GWH_O3 = 157.
GWH_O4 = 229.
GWH_O5 = 472.
GWH_190425 = 181.
GWH_MAX = 715.
vlas = 15.
mag_170817 = {'r': 17.14, 'g': 17.28, 'i':16.99, 'z': 16.85}
mag_onaxis = {'g': -16.347,'r': -16.283, 'i': -16.405,'z': -16.466}

print(f"Collecting data from '{DATA_FILE_O3}'")
df_O3 = pd.read_csv(DATA_FILE_O3, sep=" ", names=col_names, dtype=dtype)
for band in bands_d:
    df_O3[band] = df_O3[band + '0'] + df_O3['d' + band]

print(f"Collecting data from '{DATA_FILE_O4}'")
df_O4 = pd.read_csv(DATA_FILE_O4, sep=" ", names=col_names, dtype=dtype)
for band in bands_d:
    df_O4[band] = df_O4[band + '0'] + df_O4['d' + band]

print(f"Collecting data from '{DATA_FILE_O5}'")
df_O5 = pd.read_csv(DATA_FILE_O5, sep=" ", names=col_names, dtype=dtype)
for band in bands_d:
    df_O5[band] = df_O5[band + '0'] + df_O5['d' + band]

print(f"Collecting data from '{DATA_FILE_MAX}'")
df_max = pd.read_csv(DATA_FILE_MAX, sep=" ", names=col_names, dtype=dtype)
for band in bands_d:
    df_max[band] = df_max[band + '0'] + df_max['d' + band]

print(f"Collecting data from '{DATA_FILE_190425}'")
df_190425 = pd.read_csv(DATA_FILE_190425, sep=" ", names=col_names, dtype=dtype)
for band in bands_d:
    df_190425[band] = df_190425[band + '0'] + df_190425['d' + band]

print(f"Total elements:\n"
f"O3:     {len(df_O3):3.2e}\n"
f"O4:     {len(df_O4):3.2e}\n"
f"O5:     {len(df_O5):3.2e}\n"
f"max:    {len(df_max):3.2e}\n"
f"190425: {len(df_190425):3.2e}")

gw_O5 = df_O5[1. + 6. * cos(df_O5.tv) ** 2 + cos(df_O5.tv) ** 4\
            > 8 * df_O5.d ** 2 / GWH_O5 ** 2]
gw_O4 = df_O4[1. + 6. * cos(df_O4.tv) ** 2 + cos(df_O4.tv) ** 4\
            > 8 * df_O4.d ** 2 / GWH_O4 ** 2]
gw_O3 = df_O3[1. + 6. * cos(df_O3.tv) ** 2 + cos(df_O3.tv) ** 4\
            > 8 * df_O3.d ** 2 / GWH_O3 ** 2]
gw_190425 = df_190425[1. + 6. * cos(df_190425.tv) ** 2 + cos(df_190425.tv) ** 4\
            > 8 * df_190425.d ** 2 / GWH_190425 ** 2]
r_max_O3 = gw_O3.r.max()
r_max_O4 = gw_O4.r.max()
r_max_O5 = gw_O5.r.max()
N_O3 = len(gw_O3)
N_O4 = len(gw_O4)
N_O5 = len(gw_O5)
print(f"Largest distance to event with r < 22:  {df_max[df_max.r < 22].d.max()} Mpc")
print(f"Largest distance to event with r0 < 22: {df_max[df_max.r0 < 22].d.max()} Mpc")

x_0 = 1 / sqrt(8)
x_m = gw_O5.d.median() / GWH_O5
t_m = gw_O5.tv.median() / Deg

print(f"Fraction within 10 < tv < 20: {len(gw_O3[(gw_O3.tv > 10 * Deg) & (gw_O3.tv < 20 * Deg)]) / N_O3:5.3g}")
print(f"Fraction within 5 < tv < 25: {len(gw_O3[(gw_O3.tv > 5 * Deg) & (gw_O3.tv < 25 * Deg)]) / N_O3:5.3g}")
print("Report fractions in all bands")
print("  : {:>5} {:>5} {:>5} | {:>5} {:>5} {:>5} | {:>5} {:>5} {:>5}".format("O3", "O4", "O5", "O3", "O4", "O5", "O3", "O4", "O5"))
for band in bands_d:
    print(f"{band} & {100 * len(gw_O3[gw_O3[band] < 18]) / N_O3:5.2g} & {100 * len(gw_O4[gw_O4[band] < 18]) / N_O4:5.2g} & {100 * len(gw_O5[gw_O5[band] < 18]) / N_O5:5.2g} & {100 * len(gw_O3[(gw_O3[band] > 18) & (gw_O3[band] < 20)]) / N_O3:5.2g} & {100 * len(gw_O4[(gw_O4[band] > 18) & (gw_O4[band] < 20)]) / N_O4:5.2g} & {100 * len(gw_O5[(gw_O5[band] > 18) & (gw_O5[band] < 20)]) / N_O5:5.2g} & {100 * len(gw_O3[gw_O3[band] > 20]) / N_O3:5.2g} & {100 * len(gw_O4[gw_O4[band] > 20]) / N_O4:5.2g} & {100 * len(gw_O5[gw_O5[band] > 20]) / N_O5:5.2g} \\\\")

print("figure 5")
histgw, edgesgw = np.histogram(gw_O5.tv / Deg, density=True, bins=BINS)
tmp = gw_190425[(gw_190425.d < 159 + 69) & (gw_190425.d > 159 - 71)]
all_obs = tmp[(tmp.r > 21) & (tmp['g'] > 21) & (tmp['i'] > 21.5)]
print(f"Under 190425, theta = {np.percentile(all_obs.tv, 5) / Deg:.2f}, {all_obs.tv.median() / Deg:.2f}, {np.percentile(all_obs.tv, 95) / Deg:.2f}.")
for band in ['g', 'r', 'i']:
    plt.figure(figsize=(4.3, 6.45))
    obs = tmp[tmp[band] > mag_lim_d[band]]
    gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1])

    # the fisrt subplot
    ax0 = plt.subplot(gs[0])
    ax0.hist2d(tmp.tv / Deg, tmp[band], range = [[0, 90], [10.95, 26]], bins=BINS, cmin=1,
                cmap = cmap_d[band])
    ax0.errorbar([15], mag_170817[band], xerr=[[2.46], [1.64]], color=c_d[band])
    ax0.text(20, mag_170817[band], "AT2017gfo")
    ax0.set_xlim(0, 90)
    ax0.set_ylim(26, 16)
    ax0.hlines(mag_lim_d[band], xmin=0, xmax=90, linestyle="--", linewidth=thin, color="grey")
    ax0.set_ylabel(tex_d[band])

    #the second subplot
    ax1 = plt.subplot(gs[1], sharex = ax0)
    ax1.plot(edgesgw[:-1], histgw / histgw.max(), color="grey", linewidth=thin)
    hist, edges = np.histogram(obs.tv / Deg, density=True, range=[0, 90], bins=BINS)
    ax1.plot(edges[:-1], hist / hist.max(), color='black')
    a = np.percentile(obs.tv, 5) / Deg
    m = np.percentile(obs.tv, 50) / Deg
    b = np.percentile(obs.tv, 95) / Deg
    ax1.vlines(a, ymin = 0, ymax=lsub(a, edges, hist)/hist.max(), color="black", linestyle=":", linewidth=thin)
    ax1.vlines(m, ymin = 0, ymax=lsub(m, edges, hist)/hist.max(), color="black", linestyle="--", linewidth=thin)
    ax1.vlines(b, ymin = 0, ymax=lsub(b, edges, hist)/ hist.max(), color="black", linestyle=":", linewidth=thin)
    ax1.set_xlabel(r"$\theta_v$ [deg]")
    if band == 'g':
        ax1.set_ylabel(r"${\rm d} N / {\rm d} \theta_v$")
    plt.setp(ax0.get_xticklabels(), visible=False)
    ax1.set_yticklabels([])
    ax1.set_ylim(0)

    # remove vertical gap between subplots
    plt.subplots_adjust(hspace=.0)
    plt.savefig(f"{PLOT_DIR}/fig5_{band}.pdf", bbox_inches='tight')


print("figure 4")
plt.close("All")
hmax = 1
for mag_lim in [22, 21, 20, 19]:
    print(f"r_lim = {mag_lim}")
    plt.figure()
    blue = df_max[(1. + 6. * cos(df_max.tv) ** 2 + cos(df_max.tv) ** 4\
                > 8 * df_max.d ** 2 / GWH_O4 ** 2) & (df_max.r0 > mag_lim)]
    green = df_max[(1. + 6. * cos(df_max.tv) ** 2 + cos(df_max.tv) ** 4\
                > 8 * df_max.d ** 2 / GWH_O4 ** 2) & (df_max.r0 < mag_lim)]
    red = df_max[(1. + 6. * cos(df_max.tv) ** 2 + cos(df_max.tv) ** 4\
                < 8 * df_max.d ** 2 / GWH_O4 ** 2) & (df_max.r0 < mag_lim)]

    tot1 = len(blue)
    tot2 = len(green)
    tot3 = len(red)

    r1 = 10 * tot1 / (tot1 + tot2)
    r2 = 10 * tot2 / (tot1 + tot2)
    r3 = 10 * tot3 / (tot1 + tot2)

    gw_kn_grb = len(df_max[(1. + 6. * cos(df_max.tv) ** 2 + cos(df_max.tv) ** 4\
                > 8 * df_max.d ** 2 / GWH_O4 ** 2) & (df_max.r < mag_lim) & (df_max.tv < 0.1)])
    kn_grb = len(df_max[(df_max.tv < 0.1) & (df_max.r < mag_lim)])
    rgw_kn_grb = 10 * gw_kn_grb / (tot1 + tot2)
    rkn_grb = 10 * kn_grb / (tot1 + tot2)
    print(f"  GW+KN+GRB: {rgw_kn_grb:5.3g} /yr")
    print(f"  KN+GRB:    {rkn_grb:5.3g} /yr")

    h1 = np.histogram2d(blue.d, blue.tv / Deg, range=[[0, 500], [0, 90]],bins=BINS)
    h2 = np.histogram2d(green.d, green.tv / Deg, range=[[0, 500], [0, 90]],bins=BINS)
    h3 = np.histogram2d(red.d, red.tv / Deg, range=[[0, 500], [0, 90]],bins=BINS)
    h1max = np.nanmax(h1[0])
    h2max = np.nanmax(h2[0])
    h3max = np.nanmax(h3[0])
    hmax = max(hmax, h1max, h2max, h3max)
    print(f"  hs: {h1max}, {h2max}, {h3max}, {hmax}")
    h1 = plt.hist2d(blue.d, blue.tv / Deg, range=[[0, 500], [0, 90]], cmin = 1,
        bins=BINS, cmap="Blues", norm=Normalize(0, hmax))
    h2 = plt.hist2d(green.d, green.tv / Deg, range=[[0, 500], [0, 90]], cmin = 1,
        bins=BINS, cmap="Greens", norm=Normalize(0, hmax))
    h3 = plt.hist2d(red.d, red.tv / Deg, range=[[0, 500], [0, 90]], cmin = 1,
        bins=BINS, cmap="Reds", norm=Normalize(0, hmax))

    plt.text(400, 80, r"$r_{\rm lim} = $" + f" {mag_lim}", fontsize="larger")
    plt.text(250, 60, f"GW w/o KN:", color="blue")
    plt.text(340, 60, f"{r1:4.2g} /yr")
    plt.text(250, 55, f"GW w/  KN:", color="green")
    plt.text(340, 55, f"{r2:4.2g} /yr")
    plt.text(250, 50, f"KN w/o GW:", color="red")
    plt.text(340, 50, f"{r3:4.2g} /yr")

    d_l = np.linspace(GWH_O4/sqrt(8), GWH_O4, 200)
    plt.plot(d_l, [np.arccos(sqrt(-3 + sqrt(8) * sqrt(1 + d**2 / GWH_O4**2))) / Deg for d in d_l],
            color = "black")
    d_m = 10 * pc * 10 ** ((mag_lim - mag_onaxis['r'] - 4) / 5) / Mpc
    d_l = np.linspace(d_m, GWH_MAX, 200)
    plt.plot(d_l, [np.arccos(0.125 * (8 - mag_lim + mag_onaxis['r'] + 5 * log(d * Mpc / (10 * pc)))) / Deg for d in d_l],
            color = "black")
    plt.vlines(d_m, ymin = 60, ymax = 90, color="black")

    plt.axhspan(0, 0.1 / Deg, color="grey", alpha=0.5)
    plt.text(30, 2, "on-axis", fontsize="larger")

    if mag_lim >= 21:
        plt.xlabel(r"$D$ [Mpc]")
    if mag_lim in [19, 21]:
        plt.ylabel(r"$\theta_v$ [deg]")
    plt.xlim(0, 500)
    plt.ylim(0, 90)
    plt.savefig(f"{PLOT_DIR}/fig4_{mag_lim}.pdf", bbox_inches='tight')


print("figure 1a")
plt.close("All")
plt.figure()
hist, edges = np.histogram(gw_O5.d / GWH_O5, density=True, bins=BINS)
plt.plot(edges[:-1], hist)
plt.vlines(x = x_0, ymin = 0, ymax=lsub(x_0, edges, hist), linestyle="--", linewidth=thin, color="grey")
plt.vlines(x = x_m, ymin = 0, ymax=lsub(x_m, edges, hist), linestyle="--", linewidth=thin, color="grey")
plt.text(x_0, 0.4, r"$D_0 = D_H/\sqrt{8}$", rotation=90)
plt.text(x_m, 0.4, r"median", rotation=90)
plt.xlim(0, 1)
plt.ylim(0, 2)
plt.tick_params(labelleft=False)
plt.xlabel(r"$D/D_H$")
plt.ylabel(r"${\rm d}N/{\rm d}(D/D_H)$")
plt.savefig(f"{PLOT_DIR}/fig1a.pdf", bbox_inches='tight')

print("figure 1b")
plt.figure()
hist = np.cumsum(hist) / BINS
plt.plot(edges[:-1], hist)
plt.vlines(x = x_0, ymin = 0, ymax=lsub(x_0, edges, hist), linestyle="--", linewidth=thin, color="grey")
plt.vlines(x = x_m, ymin = 0, ymax=lsub(x_m, edges, hist), linestyle="--", linewidth=thin, color="grey")
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.xlabel(r"$D/D_H$")
plt.savefig(f"{PLOT_DIR}/fig1b.pdf", bbox_inches='tight')

print("figure 1c")
plt.figure()
hist, edges = np.histogram(gw_O5.tv / Deg, density=True, bins=BINS)
plt.plot(edges[:-1], hist)
plt.vlines(x = t_m, ymin = 0, ymax=lsub(t_m, edges, hist), linestyle="--", linewidth=thin, color="grey")
plt.text(t_m, 0.0025, r"median", rotation=90)
plt.xlim(0, 90)
plt.ylim(0, 0.025)
plt.tick_params(labelleft=False)
plt.xlabel(r"$\theta_v$ [deg]")
plt.ylabel(r"${\rm d} N /{\rm d}\theta_v$")
plt.savefig(f"{PLOT_DIR}/fig1c.pdf", bbox_inches='tight')

print("figure 1d")
plt.figure()
hist = np.cumsum(hist) * 90 / BINS
plt.plot(edges[:-1], hist)
plt.vlines(x = t_m, ymin=0, ymax=lsub(t_m, edges, hist), linestyle="--", linewidth=thin, color="grey")
plt.xlim(0, 90)
plt.ylim(0, 1)
plt.xlabel(r"$\theta_v$ [deg]")
plt.savefig(f"{PLOT_DIR}/fig1d.pdf", bbox_inches='tight')

print("figure 1e")
plt.figure()
hist, edges = np.histogram(cos(gw_O5.tv), density=True, bins=BINS)
plt.plot(edges[:-1], hist)
plt.xlim(0, 1)
plt.tick_params(labelleft=False)
plt.xlabel(r"$\cos \theta_v$")
plt.ylabel(r"${\rm d} N /{\rm d}\cos\theta_v$")
plt.savefig(f"{PLOT_DIR}/fig1e.pdf", bbox_inches='tight')

print("figure 2a")
plt.close("All")
plt.figure()
print("Maximum and minimum magnitudes for O4 GW events")
for f in bands_d:
    print(f"{f}: {gw_O4[f].max():7.5g}, {gw_O4[f].min():7.5g}")
    hist, edges = np.histogram(gw_O4[f], density=True, bins=BINS, range=[12, 27.5])
    plt.plot(edges[:-1], hist, label=tex_d[f], color=c_d[f], linestyle=ls_d[f])
plt.xlim(28, 14)
plt.ylim(0)
plt.xlabel("AB magnitude")
plt.ylabel(r"${\rm d} N / {\rm d}m$")
plt.legend()
plt.tick_params(labelleft=False)
plt.savefig(f"{PLOT_DIR}/fig2a.pdf", bbox_inches='tight')

print("figure 2b")
plt.figure()
hist3, edges3 = np.histogram(gw_O3.r, bins=BINS)
hist4, edges4 = np.histogram(gw_O4.r, bins=BINS)
hist5, edges5 = np.histogram(gw_O5.r, bins=BINS)
hist3 = np.cumsum(hist3)
hist4 = np.cumsum(hist4)
hist5 = np.cumsum(hist5)
alpha = hist4[-1] / 10.
rate3 = hist3 * (GWH_O3 / GWH_O4) ** 3 / alpha
rate4 = hist4 / alpha
rate5 = hist5 * (GWH_O5 / GWH_O4) ** 3 / alpha
plt.plot(np.append(edges5[:-1], 27), np.append(rate5, rate5[-1]), "-", label="O5", color="black")
plt.plot(np.append(edges4[:-1], 27), np.append(rate4, rate4[-1]), "--", label="O4", color="black")
plt.plot(np.append(edges3[:-1], 27), np.append(rate3, rate3[-1]), "-.", label="O3", color="black")

plt.vlines(r_max_O3, ymin = 0.1, ymax=rate3[-1], linestyle=":", color="grey", linewidth=thin)
plt.vlines(r_max_O4, ymin=0.1, ymax=rate4[-1], linestyle=":", color="grey", linewidth=thin)
plt.vlines(r_max_O5, ymin=0.1, ymax=rate5[-1], linestyle=":", color="grey", linewidth=thin)
plt.text(r_max_O3, 0.2, r"$r_{{\rm max}}^{{\rm O3}}$", rotation=90)
plt.text(r_max_O4, 0.2, r"$r_{{\rm max}}^{{\rm O4}}$", rotation=90)
plt.text(r_max_O5, 0.2, r"$r_{{\rm max}}^{{\rm O5}}$", rotation=90)
plt.xlabel(r"$r$")
plt.ylabel(r"$\tau_{\rm KN}~[{\rm yr}^{-1}]$")
plt.xlim(27, 17)
plt.ylim(0.1, 300)
plt.yscale("log")
plt.legend()
plt.savefig(f"{PLOT_DIR}/fig2b.pdf", bbox_inches='tight')

print("figure 3")
plt.close("All")
plt.figure()
hist, edges = np.histogram(gw_O4.tv / Deg, bins = BINS)
plt.vlines(x = t_m, ymin = 0, ymax=lsub(t_m, edges, hist), linestyle="--", linewidth=thin, color="black")
plt.plot(edges[:-1], hist, linestyle=":", label="GW events", color="black")
plt.text(2, 6 * hist.max() / 7, "on-axis", rotation=90, fontsize="larger")
for mag_lim in [21, 20, 19, 18]:
    c = cmap((21 - mag_lim)/(21 - 18))
    if mag_lim == 19:
        c = "green"
    tmp = gw_O4[gw_O4.r < mag_lim]
    hist, edges = np.histogram(tmp.tv / Deg, bins=BINS)
    a = np.percentile(tmp.tv, 5) / Deg
    m = np.percentile(tmp.tv, 50) / Deg
    b = np.percentile(tmp.tv, 95) / Deg
    plt.plot(edges[:-1], hist, label=r"$r <$" + f" {mag_lim}", color=c, linestyle=ls_l[(2+mag_lim) % 4])
    plt.vlines(m, ymin = 0, ymax=lsub(m, edges, hist), color=c, linestyle="--", linewidth=thin)

plt.axvspan(0, 0.1 / Deg, color="grey", alpha=0.5)

plt.xlabel(r"$\theta_v$ [deg]")
plt.ylabel(r"${\rm d} N / {\rm d} \theta_v$")
plt.tick_params(labelleft=False)
plt.xlim(0, 90)
plt.ylim(0)
plt.legend()
plt.savefig(f"{PLOT_DIR}/fig3.pdf", bbox_inches='tight')

print("figure 6a")
plt.close("All")
plt.figure()
hist, edges = np.histogram(gw_O4.tv / Deg, bins = BINS)
plt.text(2, 6 * hist.max() / 7, "on-axis", rotation=90, fontsize="larger")
plt.plot(edges[:-1], hist, linestyle=":", label="GW events", color="grey")
print("phi = 1")
for mag_lim in [21, 20, 19, 18]:
    c = cmap((21 - mag_lim)/(21 - 18))
    if mag_lim == 19:
        c = "green"
    tmp = gw_O4[gw_O4.r < mag_lim]
    tmp2 = tmp[tmp.r1 > 3 * vlas]
    tmp3 = tmp[tmp.r1 > 1 * vlas]
    hist, edges = np.histogram(tmp.tv / Deg, bins=BINS)
    hist2, edges2 = np.histogram(tmp2.tv / Deg, bins=BINS)
    m = np.percentile(tmp2.tv, 50) / Deg
    plt.plot(edges[:-1], hist, color=c, linestyle=":")
    plt.plot(edges2[:-1], hist2, label=r"$r <$" + f" {mag_lim}", color=c, linestyle=ls_l[(2 + mag_lim) % 4])
    plt.vlines(m, ymin = 0, ymax=lsub(m, edges2, hist2), color=c, linestyle="--", linewidth=thin)
    print(f"  r_lim = {mag_lim}: ( > vlas) = {100 * len(tmp3) / len(tmp):4.2g}%, (> 3 x vlas) = {100 * len(tmp2) / len(tmp):4.2g}%")
    print(f"  r_lim = {mag_lim}: ( > vlas) = {10 * len(tmp3) / N_O4:4.2g},  (> 3 x vlas) = {10 * len(tmp2) / N_O4:4.2g}")

plt.axvspan(0, 0.1 / Deg, color="grey", alpha=0.5)
plt.xlabel(r"$\theta_v$ [deg]")
plt.ylabel(r"${\rm d} N / {\rm d} \theta_v$")
plt.tick_params(labelleft=False)
plt.xlim(0, 90)
plt.ylim(0)
plt.legend()
plt.savefig(f"{PLOT_DIR}/fig6a.pdf", bbox_inches='tight')

print("figure 6b")
plt.figure()
hist, edges = np.histogram(gw_O4.tv / Deg, bins = BINS)
plt.plot(edges[:-1], hist, linestyle=":", label="GW events", color="grey")
plt.text(2, 6 * hist.max() / 7, "on-axis", rotation=90, fontsize="larger")
print("phi = 10")
for mag_lim in [21, 20, 19, 18]:
    c = cmap((21 - mag_lim)/(21 - 18))
    if mag_lim == 19:
        c = "green"
    tmp = gw_O4[gw_O4.r < mag_lim]
    tmp2 = tmp[tmp.r10 > 3 * vlas]
    tmp3 = tmp[tmp.r10 > 1 * vlas]
    hist, edges = np.histogram(tmp.tv / Deg, bins=BINS)
    hist2, edges2 = np.histogram(tmp2.tv / Deg, bins=BINS)
    m = np.percentile(tmp2.tv, 50) / Deg
    plt.plot(edges[:-1], hist, color=c, linestyle="--")
    plt.plot(edges2[:-1], hist2, label=r"$r <$" + f" {mag_lim}", color=c, linestyle=ls_l[(2 + mag_lim) % 4])
    plt.vlines(m, ymin = 0, ymax=lsub(m, edges2, hist2), color=c, linestyle="--", linewidth=thin)
    print(f"  r_lim = {mag_lim}: ( > vlas) = {100 * len(tmp3) / len(tmp):4.2g}%, (> 3 x vlas) = {100 * len(tmp2) / len(tmp):4.2g}%")
    print(f"  r_lim = {mag_lim}: ( > vlas) = {10 * len(tmp3) / N_O4:4.2g},  (> 3 x vlas) = {10 * len(tmp2) / N_O4:4.2g}")

plt.axvspan(0, 0.1 / Deg, color="grey", alpha=0.5)
plt.xlabel(r"$\theta_v$ [deg]")
plt.tick_params(labelleft=False)
plt.xlim(0, 90)
plt.ylim(0)
plt.legend()
plt.savefig(f"{PLOT_DIR}/fig6b.pdf", bbox_inches='tight')
