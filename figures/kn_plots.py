'''Make plots for populations paper.'''

from sys import exit
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib import gridspec
import pandas as pd
from numpy import sin, cos
import numpy as np
from helpers import *
from sys import argv

# Setup
DATA_FILE_O3 = "data/kn_O3.data"
DATA_FILE_O4 = "data/kn_O4.data"
DATA_FILE_O5 = "data/kn_O5.data"
DATA_FILE_MAX = "data/kn_max.data"
DATA_FILE_190425 = "data/kn_190425.data"
PLOT_DIR = f"plots/kn_plots"
plt.style.use('presentation')
BINS = 80
thin = 1.

colors = {"g": "green",
            "r": "red",
            "i": "black",
            "z": "blue"}
# Constants
g0 = 100.
GWH_O3 = 157.
GWH_O4 = 229.
GWH_O5 = 472.
GWH_190425 = 181.
GWH_MAX = 650.
vlas = 15.
mag_170817={'r': 17.30, 'g': 17.18, 'i':17}

print(f"Collecting data from data file")
df_O3 = pd.read_csv(DATA_FILE_O3, sep=" ",
        names=['d', 'tv', 'g', 'r', 'i', 'z', 'r1', 'r10'])

df_O4 = pd.read_csv(DATA_FILE_O4, sep=" ",
        names=['d', 'tv', 'g', 'r', 'i', 'z', 'r1', 'r10'])

df_O5 = pd.read_csv(DATA_FILE_O5, sep=" ",
        names=['d', 'tv', 'g', 'r', 'i', 'z', 'r1', 'r10'])

df_max = pd.read_csv(DATA_FILE_MAX, sep=" ",
        names=['d', 'tv', 'g', 'r', 'i', 'z', 'r1', 'r10'])

df_190425 = pd.read_csv(DATA_FILE_190425, sep=" ",
        names=['d', 'tv', 'g', 'r', 'i', 'z', 'r1', 'r10'])

print(f"Total elements: {len(df_O3['d'])}.")
gw_O5 = df_O5.loc[1. + 6. * cos(df_O5['tv']) ** 2 + cos(df_O5['tv']) ** 4\
            > 8 * df_O5['d'] ** 2 / GWH_O5 ** 2]
gw_O4 = df_O4.loc[1. + 6. * cos(df_O4['tv']) ** 2 + cos(df_O4['tv']) ** 4\
            > 8 * df_O4['d'] ** 2 / GWH_O4 ** 2]
gw_O3 = df_O3.loc[1. + 6. * cos(df_O3['tv']) ** 2 + cos(df_O3['tv']) ** 4\
            > 8 * df_O3['d'] ** 2 / GWH_O3 ** 2]
gw_190425 = df_190425.loc[1. + 6. * cos(df_190425['tv']) ** 2 + cos(df_190425['tv']) ** 4\
            > 8 * df_190425['d'] ** 2 / GWH_190425 ** 2]
r_max_O3 = gw_O3['r'].max()
r_max_O4 = gw_O4['r'].max()
r_max_O5 = gw_O5['r'].max()

x_0 = 1 / sqrt(8)
x_m = gw_O5['d'].median() / GWH_O5
t_m = gw_O5['tv'].median() / Deg

print("figure 5")
histgw, edgesgw = np.histogram(gw_O5['tv'] / Deg, density=True, bins=BINS)
mag_lim_d = {'g': 21, 'r': 21, 'i': 21.5}
cmap_d = {'g': 'Greens', 'r': 'Reds', 'i': 'Greys'}
c_d = {'g': 'green', 'r': 'red', 'i': 'grey'}
tex_d = {'g': r'$g$', 'r': r'$r$', 'i': r'$i$'}
tmp = gw_190425.loc[(gw_190425['d'] < 159 + 69) & (gw_190425['d'] > 159 - 71)]
all_obs = tmp.loc[(tmp['r'] > 21) & (tmp['g'] > 21) & (tmp['i'] > 21.5)]
print(f"Under 190425, theta = {np.percentile(all_obs['tv'], 5) / Deg:.2f}, {all_obs['tv'].median() / Deg:.2f}, {np.percentile(all_obs['tv'], 95) / Deg:.2f}.")
for band in ['g', 'r', 'i']:
    plt.figure()
    obs = tmp.loc[tmp[band] > mag_lim_d[band]]
    gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1])

    # the fisrt subplot
    ax0 = plt.subplot(gs[0])
    ax0.hist2d(tmp['tv'] / Deg, tmp[band], range = [[0, 90], [10.95, 26]], bins=BINS, cmin=1,
                cmap = cmap_d[band])
    ax0.errorbar([15], mag_170817[band], xerr=[[2.46], [1.64]], color=c_d[band])
    ax0.text(20, mag_170817[band], "AT2017gfo")
    ax0.set_xlim(0, 90)
    ax0.set_ylim(26, 16)
    ax0.hlines(mag_lim_d[band], xmin=0, xmax=90, linestyle="--", linewidth=thin, color="grey")
    ax0.set_ylabel(tex_d[band])

    #the second subplot
    ax1 = plt.subplot(gs[1], sharex = ax0)
    ax1.plot(edgesgw[:-1], histgw / max(histgw), color="grey", linewidth=thin)
    hist, edges = np.histogram(obs['tv'] / Deg, density=True, range=[0, 90], bins=BINS)
    ax1.plot(edges[:-1], hist / max(hist), color='black')
    a = np.percentile(obs['tv'], 5) / Deg
    m = np.percentile(obs['tv'], 50) / Deg
    b = np.percentile(obs['tv'], 95) / Deg
    ax1.vlines(a, ymin = 0, ymax=lsub(a, edges, hist)/max(hist), color="black", linestyle=":", linewidth=thin)
    ax1.vlines(m, ymin = 0, ymax=lsub(m, edges, hist)/max(hist), color="black", linestyle="--", linewidth=thin)
    ax1.vlines(b, ymin = 0, ymax=lsub(b, edges, hist)/ max(hist), color="black", linestyle=":", linewidth=thin)
    ax1.set_xlabel(r"$\theta_v$ [deg]")
    plt.setp(ax0.get_xticklabels(), visible=False)
    ax1.set_yticklabels([])
    ax1.set_ylim(0)

    # remove vertical gap between subplots
    plt.subplots_adjust(hspace=.0)
    plt.savefig(f"{PLOT_DIR}/fig5_{band}.pdf", bbox_inches='tight')


print("figure 4")
plt.close("All")
hmax = 1
for mag_lim in [22, 19, 20, 21]:
    plt.figure()
    blue = df_max.loc[(1. + 6. * cos(df_max['tv']) ** 2 + cos(df_max['tv']) ** 4\
                > 8 * df_max['d'] ** 2 / GWH_O4 ** 2) & (df_max['r'] > mag_lim)]
    green = df_max.loc[(1. + 6. * cos(df_max['tv']) ** 2 + cos(df_max['tv']) ** 4\
                > 8 * df_max['d'] ** 2 / GWH_O4 ** 2) & (df_max['r'] < mag_lim)]
    red = df_max.loc[(1. + 6. * cos(df_max['tv']) ** 2 + cos(df_max['tv']) ** 4\
                < 8 * df_max['d'] ** 2 / GWH_O4 ** 2) & (df_max['r'] < mag_lim)]

    tot1 = len(blue['d'])
    tot2 = len(green['d'])
    tot3 = len(red['d'])

    r1 = 10 * tot1 / (tot1 + tot2)
    r2 = 10 * tot2 / (tot1 + tot2)
    r3 = 10 * tot3 / (tot1 + tot2)

    h1 = np.histogram2d(blue['d'], blue['tv'] / Deg, range=[[0, 650], [0, 90]],bins=BINS)
    h2 = np.histogram2d(green['d'], green['tv'] / Deg, range=[[0, 650], [0, 90]],bins=BINS)
    h3 = np.histogram2d(red['d'], red['tv'] / Deg, range=[[0, 650], [0, 90]],bins=BINS)
    h1max = np.nanmax(h1[0])
    h2max = np.nanmax(h2[0])
    h3max = np.nanmax(h3[0])
    hmax = max(hmax, h1max, h2max, h3max)
    print(f"Mag_lim = {mag_lim}: {h1max}, {h2max}, {h3max}, {hmax}")
    h1 = plt.hist2d(blue['d'], blue['tv'] / Deg, range=[[0, 650], [0, 90]], cmin = 1,
        bins=BINS, cmap="Blues", norm=Normalize(0, hmax))
    h2 = plt.hist2d(green['d'], green['tv'] / Deg, range=[[0, 650], [0, 90]], cmin = 1,
        bins=BINS, cmap="Greens", norm=Normalize(0, hmax))
    h3 = plt.hist2d(red['d'], red['tv'] / Deg, range=[[0, 650], [0, 90]], cmin = 1,
        bins=BINS, cmap="Reds", norm=Normalize(0, hmax))

    plt.text(350, 60, f"GW w/o KN: {r1:.2f} /yr")
    plt.text(350, 55, f"GW w/ KN: {r2:.2f} /yr")
    plt.text(350, 50, f"KN w/o GW: {r3:.2f} /yr")


    d_l = np.linspace(GWH_O4/sqrt(8), GWH_O4, 200)
    plt.plot(d_l, [np.arccos(sqrt(-3 + sqrt(8) * sqrt(1 + d**2 / GWH_O4**2))) / Deg for d in d_l],
            color = "black")
    d_m = 10 * pc * 10 ** ((mag_lim + 1 + 16.9 - 4) / 5) / Mpc
    d_l = np.linspace(d_m, GWH_MAX, 200)
    plt.plot(d_l, [np.arccos(0.125 * (8-mag_lim - 1 - 16.9 + 5 * log(d * Mpc / (10 * pc)))) / Deg for d in d_l],
            color = "black")
    plt.vlines(d_m, ymin = 60, ymax = 90, color="black")

    plt.axhspan(0, 0.1 / Deg, color="grey", alpha=0.5)
    plt.text(30, 2, "on-axis")


    plt.xlabel(r"$D$ [Mpc]")
    plt.ylabel(r"$\theta_v$ [deg]")
    plt.xlim(0, 650)
    plt.ylim(0, 90)
    plt.savefig(f"{PLOT_DIR}/fig4_{mag_lim}.pdf", bbox_inches='tight')

print("figure 1a")
plt.close("All")
plt.figure()
hist, edges = np.histogram(gw_O5['d'] / GWH_O5, density=True, bins=BINS)
plt.plot(edges[:-1], hist)
plt.vlines(x = x_0, ymin = 0, ymax=lsub(x_0, edges, hist), linestyle="--", linewidth=thin, color="grey")
plt.vlines(x = x_m, ymin = 0, ymax=lsub(x_m, edges, hist), linestyle="--", linewidth=thin, color="grey")
plt.text(x_0, 0.4, r"$D / D_H = 1/\sqrt{8}$", rotation=90)
plt.text(x_m, 0.4, r"median", rotation=90)
plt.xlim(0, 1)
plt.ylim(0, 2)
plt.tick_params(labelleft=False)
plt.xlabel(r"$D/D_H$")
plt.savefig(f"{PLOT_DIR}/fig1a.pdf", bbox_inches='tight')

print("figure 1b")
plt.figure()
hist = np.cumsum(hist) / BINS
plt.plot(edges[:-1], hist)
plt.vlines(x = x_0, ymin = 0, ymax=lsub(x_0, edges, hist), linestyle="--", linewidth=thin, color="grey")
plt.vlines(x = x_m, ymin = 0, ymax=lsub(x_m, edges, hist), linestyle="--", linewidth=thin, color="grey")
plt.text(x_0, 0.3, r"$D / D_H = 1/\sqrt{8}$", rotation=90)
plt.text(x_m, 0.3, r"median", rotation=90)
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.xlabel(r"$D/D_H$")
plt.savefig(f"{PLOT_DIR}/fig1b.pdf", bbox_inches='tight')

print("figure 1c")
plt.figure()
hist, edges = np.histogram(gw_O5['tv'] / Deg, density=True, bins=BINS)
plt.plot(edges[:-1], hist)
plt.vlines(x = t_m, ymin = 0, ymax=lsub(t_m, edges, hist), linestyle="--", linewidth=thin, color="grey")
plt.text(t_m, 0.0025, r"median", rotation=90)
plt.xlim(0, 90)
plt.ylim(0, 0.025)
plt.tick_params(labelleft=False)
plt.xlabel(r"$\theta_v$ [deg]")
plt.savefig(f"{PLOT_DIR}/fig1c.pdf", bbox_inches='tight')

print("figure 1d")
plt.figure()
hist = np.cumsum(hist) * 90 / BINS
plt.plot(edges[:-1], hist)
plt.vlines(x = t_m, ymin=0, ymax=lsub(t_m, edges, hist), linestyle="--", linewidth=thin, color="grey")
plt.text(t_m, 1/10, r"median", rotation=90)
plt.xlim(0, 90)
plt.ylim(0, 1)
plt.xlabel(r"$\theta_v$ [deg]")
plt.savefig(f"{PLOT_DIR}/fig1d.pdf", bbox_inches='tight')

print("figure 2a")
plt.close("All")
plt.figure()
for f in ['g', 'r', 'i', 'z']:
    print(f"{f}: {np.max(gw_O4[f])}, {np.min(gw_O4[f])}")
    hist, edges = np.histogram(gw_O4[f], density=True, bins=BINS, range=[10.95, 26])
    plt.plot(edges[:-1], hist, label=f, color=colors[f])
plt.xlim(26, 16)
plt.ylim(0)
plt.xlabel("AB magnitude")
plt.legend()
plt.tick_params(labelleft=False)
plt.savefig(f"{PLOT_DIR}/fig2a.pdf", bbox_inches='tight')

print("figure 2b")
plt.figure()
hist3, edges3 = np.histogram(gw_O3['r'], bins=BINS)
hist4, edges4 = np.histogram(gw_O4['r'], bins=BINS)
hist5, edges5 = np.histogram(gw_O5['r'], bins=BINS)
hist3 = np.cumsum(hist3)
alpha = hist3[-1] / 5
hist4 = np.cumsum(hist4)
hist5 = np.cumsum(hist5)
rate3 = hist3 / alpha
rate4 = hist4 * (GWH_O4 / GWH_O3)**3 / alpha
rate5 = hist5 * (GWH_O5 / GWH_O3)**3 / alpha
plt.plot(np.append(edges5[:-1], 26), np.append(rate5, rate5[-1]), "-", label="O5", color="black")
plt.plot(np.append(edges4[:-1], 26), np.append(rate4, rate4[-1]), "--", label="O4", color="black")
plt.plot(np.append(edges3[:-1], 26), np.append(rate3, rate3[-1]), "-.", label="O3", color="black")

plt.vlines(r_max_O3, ymin = 0.1, ymax=rate3[-1], linestyle=":", color="grey", linewidth=thin)
plt.vlines(r_max_O4, ymin=0.1, ymax=rate4[-1], linestyle=":", color="grey", linewidth=thin)
plt.vlines(r_max_O5, ymin=0.1, ymax=rate5[-1], linestyle=":", color="grey", linewidth=thin)
plt.text(r_max_O3, 0.2, r"$r_{{\rm max}}^{{\rm O3}}$", rotation=90)
plt.text(r_max_O4, 0.2, r"$r_{{\rm max}}^{{\rm O4}}$", rotation=90)
plt.text(r_max_O5, 0.2, r"$r_{{\rm max}}^{{\rm O5}}$", rotation=90)
plt.xlabel(r"$r$")
plt.ylabel(r"$\tau_{\rm KN}~[{\rm yr}^{-1}]$")
plt.xlim(26, 16)
plt.ylim(0.1, 300)
plt.yscale("log")
plt.legend()
plt.savefig(f"{PLOT_DIR}/fig2b.pdf", bbox_inches='tight')

print("figure 3")
plt.close("All")
plt.figure()
hist, edges = np.histogram(gw_O4['tv'] / Deg, bins = BINS)
plt.plot(edges[:-1], hist, linestyle=":", label="GW events", color="black")
plt.text(2, 6 * max(hist) / 7, "on-axis", rotation=90)
for mag_lim in [21, 20, 19, 18]:
    c = cmap((21 - mag_lim)/(21 - 18))
    tmp = gw_O4.loc[gw_O4['r'] < mag_lim]
    hist, edges = np.histogram(tmp['tv'] / Deg, bins=BINS)
    a = np.percentile(tmp['tv'], 5) / Deg
    m = np.percentile(tmp['tv'], 50) / Deg
    b = np.percentile(tmp['tv'], 95) / Deg
    plt.plot(edges[:-1], hist, label=r"$r <$" + f" {mag_lim}", color=c)
    plt.vlines(a, ymin = 0, ymax=lsub(a, edges, hist), color=c, linestyle=":", linewidth=thin)
    plt.vlines(m, ymin = 0, ymax=lsub(m, edges, hist), color=c, linestyle="--", linewidth=thin)
    plt.vlines(b, ymin = 0, ymax=lsub(b, edges, hist), color=c, linestyle=":", linewidth=thin)

plt.axvspan(0, 0.1 / Deg, color="grey", alpha=0.5)

plt.xlabel(r"$\theta_v$ [deg]")
plt.tick_params(labelleft=False)
plt.xlim(0, 90)
plt.ylim(0)
plt.legend()
plt.savefig(f"{PLOT_DIR}/fig3.pdf", bbox_inches='tight')

print("figure 6a")
plt.close("All")
plt.figure()
hist, edges = np.histogram(gw_O4['tv'] / Deg, bins = BINS)
plt.text(2, 6 * max(hist) / 7, "on-axis", rotation=90)
plt.plot(edges[:-1], hist, linestyle=":", label="GW events", color="grey")
for mag_lim in [21, 20, 19, 18]:
    c = cmap((21 - mag_lim)/(21 - 18))
    tmp = gw_O4.loc[gw_O4['r'] < mag_lim]
    tmp2= tmp.loc[tmp['r1'] > 45]
    hist, edges = np.histogram(tmp['tv'] / Deg, bins=BINS)
    hist2, edges2 = np.histogram(tmp2['tv'] / Deg, bins=BINS)
    plt.plot(edges[:-1], hist, color=c, linestyle="--")
    plt.plot(edges2[:-1], hist2, label=r"$r <$" + f" {mag_lim}", color=c, linestyle="-")


plt.axvspan(0, 0.1 / Deg, color="grey", alpha=0.5)
plt.xlabel(r"$\theta_v$ [deg]")
plt.tick_params(labelleft=False)
plt.xlim(0, 90)
plt.ylim(0)
plt.legend()
plt.savefig(f"{PLOT_DIR}/fig6a.pdf", bbox_inches='tight')

print("figure 6b")
plt.figure()
hist, edges = np.histogram(gw_O4['tv'] / Deg, bins = BINS)
plt.plot(edges[:-1], hist, linestyle=":", label="GW events", color="grey")
plt.text(2, 6 * max(hist) / 7, "on-axis", rotation=90)
for mag_lim in [21, 20, 19, 18]:
    c = cmap((21 - mag_lim)/(21 - 18))
    tmp = gw_O4.loc[gw_O4['r'] < mag_lim]
    tmp2= tmp.loc[tmp['r10'] > 45]
    hist, edges = np.histogram(tmp['tv'] / Deg, bins=BINS)
    hist2, edges2 = np.histogram(tmp2['tv'] / Deg, bins=BINS)
    plt.plot(edges[:-1], hist, color=c, linestyle="--")
    plt.plot(edges2[:-1], hist2, label=r"$r <$" + f" {mag_lim}", color=c, linestyle="-")


plt.axvspan(0, 0.1 / Deg, color="grey", alpha=0.5)
plt.xlabel(r"$\theta_v$ [deg]")
plt.tick_params(labelleft=False)
plt.xlim(0, 90)
plt.ylim(0)
plt.legend()
plt.savefig(f"{PLOT_DIR}/fig6b.pdf", bbox_inches='tight')
