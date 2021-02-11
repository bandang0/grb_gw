'''Make plots for populations paper.'''

from sys import exit
import matplotlib.pyplot as plt
import pandas as pd
import corner
from numpy import sin, cos, sqrt
import numpy as np
from helpers import *
from sys import argv
from matplotlib.ticker import ScalarFormatter

DATA_FILE = "data/g16_f_pl.data"
PLOT_DIR = f"plots/aa18_plots_referee"
plt.style.use('palmerio')
BINS = 70

# Constants
g0 = 100.
GWH_O2 = 135.6
GWH_O3 = 226.
GWH_DESIGN = 429.4
GWH_MAX = 500.0
AH_O2 = GWH_O2 / 1.58
AH_O3 = GWH_O3 / 1.58
AH_DESIGN = GWH_DESIGN / 1.58
AH_MAX = GWH_MAX / 1.58
GWH_MIN = 50.
LS_MIN = -1.0
LS_MAX = 1.5
vlas = 15.
ska1s = 3.
ska2s = 0.3
fnx_fnr = 1.80e-5
fx_fnr = 4.23e-19
thin=0.4
thick=0.7

# Read event data and create pandas data frame.
print(f"Collecting data from {DATA_FILE}")
df = pd.read_csv(DATA_FILE, sep=" ",
        names=['d', 'n', 'e0', 'eb',
        'tv', 'tj', 'pte', 'ptne', 'pf', 'pm'],
        dtype={'ti': float, 'd': float, 'n': float, 'e0': float,
        'eb': float, 'tv': float, 'tj': float, 'pte': float, 'ptne': float,
        'pf': float, 'pm': float})
df = df.loc[df['eb'] < 0.8]
# Add other columns to frame
df['lpte'] = log(df['pte'])
df['lptne'] = log(df['ptne'])
df['ln'] = log(df['n'])
df['lpf'] = log(df['pf'])
df['le0'] = log(df['e0'])
df['ld'] = log(df['d'])

# Collect event data according to detections in GW, AG, etc.
print("Creating derived dataframes.")
gw_O2 = df.loc[1. + 6. * cos(df['tv']) ** 2 + cos(df['tv']) ** 4\
            > 8 * df['d'] ** 2 / AH_O2 ** 2]
gw_O3 = df.loc[1. + 6. * cos(df['tv']) ** 2 + cos(df['tv']) ** 4\
            > 8 * df['d'] ** 2 / AH_O3 ** 2]
gw_design = df.loc[1. + 6. * cos(df['tv']) ** 2 + cos(df['tv']) ** 4\
            > 8 * df['d'] ** 2 / AH_DESIGN ** 2]
gw_max = df.loc[1. + 6. * cos(df['tv']) ** 2 + cos(df['tv']) ** 4\
            > 8 * df['d'] ** 2 / AH_MAX ** 2]

gw_vla_O2 = gw_O2.loc[gw_O2['pf'] > vlas]
gw_vla_O3 = gw_O3.loc[gw_O3['pf'] > vlas]
gw_vla_design = gw_design.loc[gw_design['pf'] > vlas]
gw_ska1_design = gw_design.loc[gw_design['pf'] > ska1s]
gw_ska2_design = gw_design.loc[gw_design['pf'] > ska2s]
gw_ska1_max = gw_max.loc[gw_max['pf'] > ska1s]

print('tv')
print(np.mean(gw_vla_O2['tv']) * 180 / Pi)
print(np.mean(gw_vla_O3['tv']) * 180 / Pi)
print(np.mean(gw_vla_design['tv']) * 180 / Pi)
print(np.mean(gw_ska1_design['tv']) * 180 / Pi)
print(np.mean(gw_ska2_design['tv']) * 180 / Pi)

print('frac')
f_vla_O2 = float(len(gw_vla_O2))/len(gw_O2)
f_vla_O3 = float(len(gw_vla_O3))/len(gw_O3)
f_vla_design = float(len(gw_vla_design))/len(gw_design)
f_ska1_design = float(len(gw_ska1_design))/len(gw_design)
f_ska2_design = float(len(gw_ska2_design))/len(gw_design)
print(f"O2+VLA {f_vla_O2}")
print(f"O3+VLA {f_vla_O3}")
print(f"D+VLA {f_vla_design}")
print(f"D+SKA1 {f_ska1_design}")
print(f"D+SKA2 {f_ska2_design}")

# Fractions as function of horizon
gwh_l = np.linspace(GWH_MIN, GWH_MAX, 10)
sens_l = np.linspace(LS_MIN, LS_MAX, 10)
X, Y = np.meshgrid(gwh_l, sens_l)
mean_tv = np.zeros_like(X)
std_tv = np.zeros_like(X)
frac = np.zeros_like(X)
frac_vla = np.zeros_like(gwh_l)
frac_ska1 = np.zeros_like(gwh_l)
frac_ska2 = np.zeros_like(gwh_l)
num_vla = np.zeros_like(gwh_l)
num_ska1 = np.zeros_like(gwh_l)
num_ska2 = np.zeros_like(gwh_l)
mean_lpf = np.zeros_like(X)
n_joint = np.zeros_like(X)
frac_onaxis = np.zeros_like(X)
n_jointO3 = len(gw_vla_O3)
n_jointO2 = len(gw_vla_O2)
e0l = 10 ** np.array([50., 50.5, 51., 51.5, 52.])
e0ls = ['500', '505', '510', '515', '520']
frac ={"005":np.zeros_like(e0l),
        "010":np.zeros_like(e0l),
        "015": np.zeros_like(e0l)}
ang ={"005":np.zeros_like(e0l),
        "010":np.zeros_like(e0l),
        "015": np.zeros_like(e0l)}
f150 ={"005":np.zeros_like(e0l),
        "010":np.zeros_like(e0l),
        "015": np.zeros_like(e0l)}
g16f150 = dict()
g16ang = dict()
g16frac = dict()

print("Functions of energy.")
for i, e0 in enumerate(e0ls):
    for tj in ["005", "010", "015"]:
        df_e0 = pd.read_csv(f"data/wp_{e0}_{tj}.data", sep=" ",
            names=['d', 'n', 'e0', 'eb',
            'tv', 'tj', 'pte', 'ptne', 'pf', 'pm'],
            dtype={'ti': float, 'd': float, 'n': float, 'e0': float,
            'eb': float, 'tv': float, 'tj': float, 'pte': float, 'ptne': float,
            'pf': float, 'pm': float})
        tmp_gw = df_e0.loc[1. + 6 * cos(df_e0['tv']) ** 2 \
                + cos(df_e0['tv']) ** 4\
                > 8 * df_e0['d'] ** 2 / AH_O3 ** 2]
        tmp_joint = tmp_gw.loc[tmp_gw['pf'] > vlas]
        frac[tj][i] = 100. * float(len(tmp_joint))/len(tmp_gw)
        ang[tj][i] = tmp_joint['tv'].mean() * 180. / Pi
        f150[tj][i] = 100. * float(len(tmp_joint.loc[tmp_joint['pte'] < 150.]))/len(tmp_joint)

for tj in ["005", "010", "015"]:
    df_e0 = pd.read_csv(f"data/g16_{tj}.data", sep=" ",
        names=['d', 'n', 'e0', 'eb',
        'tv', 'tj', 'pte', 'ptne', 'pf', 'pm'],
        dtype={'ti': float, 'd': float, 'n': float, 'e0': float,
        'eb': float, 'tv': float, 'tj': float, 'pte': float, 'ptne': float,
        'pf': float, 'pm': float})
    tmp_gw = df_e0.loc[1. + 6 * cos(df_e0['tv']) ** 2 \
            + cos(df_e0['tv']) ** 4\
            > 8 * df_e0['d'] ** 2 / AH_O3 ** 2]
    tmp_joint = tmp_gw.loc[tmp_gw['pf'] > vlas]
    g16frac[tj] = 100. * float(len(tmp_joint))/len(tmp_gw)
    g16ang[tj] = tmp_joint['tv'].mean() * 180. / Pi
    g16f150[tj] = 100. * float(len(tmp_joint.loc[tmp_joint['pte'] < 150.]))/len(tmp_joint)
print(g16frac, g16ang, g16f150)

print("Functions of horizon.")
for i, gwh in enumerate(gwh_l):
    ah = gwh / 1.58
    total_events = len(df.loc[df['d'] < ah])
    vla_df = df[(1. + 6 * cos(df['tv']) ** 2 + cos(df['tv']) ** 4\
        > 8 * df['d'] ** 2 / ah ** 2) & (df['pf'] > vlas)]
    ska1_df = df[(1. + 6 * cos(df['tv']) ** 2 + cos(df['tv']) ** 4\
        > 8 * df['d'] ** 2 / ah ** 2) & (df['pf'] > ska1s)]
    ska2_df = df[(1. + 6 * cos(df['tv']) ** 2 + cos(df['tv']) ** 4\
        > 8 * df['d'] ** 2 / ah ** 2) & (df['pf'] > ska2s)]
    frac_vla[i] = float(len(vla_df)) / total_events
    frac_ska1[i] = float(len(ska1_df)) / total_events
    frac_ska2[i] = float(len(ska2_df)) / total_events
    num_vla[i] = float(len(vla_df))
    num_ska1[i] = float(len(ska1_df))
    num_ska2[i] = float(len(ska2_df))
    for j, s in enumerate(sens_l):
        tmp_df = df[(1. + 6 * cos(df['tv']) ** 2 + cos(df['tv']) ** 4\
            > 8 * df['d'] ** 2 / ah ** 2) & (df['lpf'] > s)]
        on_axis_df = tmp_df.loc[tmp_df['tv'] < 0.1]
        mean_tv[j, i] = tmp_df['tv'].mean()
        std_tv[j, i] = np.std(tmp_df['tv'])
        frac[j, i] = float(len(tmp_df)) / total_events
        n_joint[j, i] = len(tmp_df)
        frac_onaxis[j, i] = float(len(on_axis_df))/len(tmp_df)



print("fig2a")
plt.figure()
plt.plot(gwh_l, 100 * frac_vla / 0.29, label="VLA", linestyle="--")
plt.plot(gwh_l, 100 * frac_ska1 / 0.29, label="SKA1", linestyle="-")
plt.plot(gwh_l, 100 * frac_ska2 / 0.29, label="SKA2", linestyle="-.")
plt.xlabel(r"$H/{\rm Mpc}$")
plt.ylabel(r"$N_{\rm joint} / N_{\rm GW}$ (\%)")
plt.xlim([50, 500])
plt.ylim([0, 100])
plt.vlines([GWH_O2, GWH_O3, GWH_DESIGN],
        0, 100, linestyles='dashed', linewidth=thin)
plt.text(GWH_O2-15, 5, "O2", rotation=90)
plt.text(GWH_O3-15, 5, "O3", rotation=90)
plt.text(GWH_DESIGN-15, 5, "Design", rotation=90)

plt.legend()
plt.savefig(f"{PLOT_DIR}/fig2a.pdf", bbox_inches='tight')
plt.close("All")

print("fig2b")
plt.figure()
f = 1.
plt.plot(gwh_l, f * num_vla / len(gw_vla_O3), label="VLA", linestyle="--")
plt.plot(gwh_l, f * num_ska1 / len(gw_vla_O3), label="SKA1", linestyle="-")
plt.plot(gwh_l, f * num_ska2 / len(gw_vla_O3), label="SKA2", linestyle="-.")
plt.xlabel(r"$H/{\rm Mpc}$")
plt.ylabel(r"$N_{\rm joint}$ (norm. to O3+VLA)")
plt.xlim([50, 500])
plt.xscale('log')
plt.yscale('log')
plt.ylim([4.e-2, 30])
plt.vlines([GWH_O2, GWH_O3, GWH_DESIGN],
        0, 100, linestyles='dashed', linewidth=thin)
plt.text(GWH_O2-15, 6.e-2, "O2", rotation=90)
plt.text(GWH_O3-15, 6.e-2, "O3", rotation=90)
plt.text(GWH_DESIGN-15, 6.e-2, "Design", rotation=90)
plt.legend()
plt.savefig(f"{PLOT_DIR}/fig2b.pdf", bbox_inches='tight')
plt.close("All")

print("fig3a")
plt.figure()
xbins= np.linspace(0., 1., BINS)
counts, bin_edges = np.histogram(gw_O3['d'] / AH_O3, bins=xbins, density=True)
plt.plot(xbins[:-1], counts , label="GW-only",
linewidth=thick, color="black", linestyle="-")
counts2, bin_edges2 = np.histogram(gw_vla_O2['d'] / AH_O2, bins=xbins, density=True)
plt.plot(xbins[:-1], counts2 * f_vla_O2, label="O2+VLA",
linewidth=thick, color="blue", linestyle="--")
counts3, bin_edges3 = np.histogram(gw_vla_O3['d'] / AH_O3, bins=xbins, density=True)
plt.plot(xbins[:-1], counts3 * f_vla_O3, label="O3+VLA",
linewidth=thick, color="red", linestyle="-.")
plt.xlabel(r"$D/\bar{H}$")
plt.ylabel(r"$1/N_{\rm GW}~{\rm d}N_{\rm joint}/{\rm d}D/\bar{H}$")
plt.legend(loc="upper left")
plt.xlim([0,1])
plt.ylim([0,1.78])
plt.savefig(f"{PLOT_DIR}/fig3a.pdf", bbox_inches='tight')
plt.close("All")

print("fig3b")
plt.figure()
xbins= np.linspace(0., 1., BINS)
counts, bin_edges = np.histogram(gw_O3['d'] / AH_O3, bins=xbins, density=True)
plt.plot(xbins[:-1], np.cumsum(counts) / BINS , label="GW-only",
linewidth=thick, color="black", linestyle="-")
counts2, bin_edges2 = np.histogram(gw_vla_O2['d'] / AH_O2, bins=xbins, density=True)
plt.plot(xbins[:-1], np.cumsum(counts2) * f_vla_O2 / BINS, label="O2+VLA",
linewidth=thick, color="blue", linestyle="--")
counts3, bin_edges3 = np.histogram(gw_vla_O3['d'] / AH_O3, bins=xbins, density=True)
plt.plot(xbins[:-1], np.cumsum(counts3) * f_vla_O3 / BINS, label="O3+VLA",
linewidth=thick, color="red", linestyle="-.")
plt.xlabel(r"$D/\bar{H}$")
plt.ylabel(r"$N_{\rm joint}(< D/\bar{H})/N_{\rm GW}$")
plt.legend(loc="upper left")
plt.xlim([0,1])
plt.ylim([0,1])
plt.savefig(f"{PLOT_DIR}/fig3b.pdf", bbox_inches='tight')
plt.close("All")

print("fig4a")
plt.figure()
xbins= np.linspace(0., 90., BINS)
counts, bin_edges = np.histogram(gw_O3['tv'] * 180 / Pi, bins=xbins, density=True)
plt.plot(xbins[:-1], counts , label="GW-only",
linewidth=thick, color="black", linestyle="-")
counts2, bin_edges2 = np.histogram(gw_vla_O2['tv'] * 180 / Pi, bins=xbins, density=True)
plt.plot(xbins[:-1], counts2, label="O2+VLA",
linewidth=thick, color="blue", linestyle="--")
counts3, bin_edges3 = np.histogram(gw_vla_O3['tv'] * 180 / Pi, bins=xbins, density=True)
plt.plot(xbins[:-1], counts3, label="O3+VLA",
linewidth=thick, color="red", linestyle="-.")
plt.xlabel(r"$\theta_v / {\rm deg}$")
plt.ylabel(r"$1/N~{\rm d}N/{\rm d}\theta_v$")
plt.legend(loc="upper right")
plt.xlim([0,90])
plt.ylim([0,0.035])
plt.savefig(f"{PLOT_DIR}/fig4a.pdf", bbox_inches='tight')
plt.close("All")

print("fig4b")
plt.figure()
xbins= np.linspace(0., 90., BINS)
counts, bin_edges = np.histogram(gw_O3['tv'] * 180 / Pi, bins=xbins, density=True)
plt.plot(xbins[:-1], np.cumsum(counts) * 90. / BINS, label="GW-only",
linewidth=thick, color="black", linestyle="-")
counts2, bin_edges2 = np.histogram(gw_vla_O2['tv'] * 180 / Pi, bins=xbins, density=True)
plt.plot(xbins[:-1], np.cumsum(counts2)  * 90. / BINS, label="O2+VLA",
linewidth=thick, color="blue", linestyle="--")
counts3, bin_edges3 = np.histogram(gw_vla_O3['tv'] * 180 / Pi, bins=xbins, density=True)
plt.plot(xbins[:-1], np.cumsum(counts3) * 90. / BINS, label="O3+VLA",
linewidth=thick, color="red", linestyle="-.")
plt.xlabel(r"$\theta_v / {\rm deg}$")
plt.ylabel(r"$N_{\rm joint}(< \theta_v)/N_{\rm joint}$")
plt.legend(loc="lower right")
plt.xlim([0,90])
plt.ylim([0,1])
plt.savefig(f"{PLOT_DIR}/fig4b.pdf", bbox_inches='tight')
plt.close("All")

print("fig5a")
plt.figure()
plt.xlabel("$H$ / Mpc")
plt.ylabel(r"$\log(s~/~\mu {\rm Jy})$")
plt.vlines(GWH_O2, LS_MIN, LS_MAX, linestyles='dashed', linewidth=thin)
plt.vlines(GWH_O3, LS_MIN, LS_MAX, linestyles='dashed', linewidth=thin)
plt.vlines(GWH_DESIGN, LS_MIN, LS_MAX, linestyles='dashed', linewidth=thin)
plt.hlines(log(vlas), GWH_MIN, GWH_MAX, linestyles='dashed', linewidth=thin)
plt.hlines(log(ska1s), GWH_MIN, GWH_MAX, linestyles='dashed', linewidth=thin)
plt.hlines(log(ska2s), GWH_MIN, GWH_MAX, linestyles='dashed', linewidth=thin)
plt.contourf(X, Y, mean_tv * 180 / Pi, cmap='jet', levels=11)
cbar = plt.colorbar()
cbar.set_label(r'mean observed $\theta_v$ (deg)')
plt.text(GWH_O2-15, -0.8, "O2", rotation=90)
plt.text(GWH_O3-15, -0.8, "O3", rotation=90)
plt.text(GWH_DESIGN-15, -0.8, "Design", rotation=90)
plt.text(80, log(vlas)+0.1, "VLA", rotation=0)
plt.text(80, log(ska1s)+0.1, "SKA1", rotation=0)
plt.text(80, log(ska2s)+0.1, "SKA2/ngVLA", rotation=0)
plt.savefig(f"{PLOT_DIR}/fig5a.pdf", bbox_inches = 'tight')

print("fig5b")
plt.figure()
plt.xlabel("$H$/Mpc")
plt.ylabel(r"$\log(s/\mu {\rm Jy})$")
plt.vlines(GWH_O2, LS_MIN, LS_MAX, linestyles='dashed')
plt.vlines(GWH_O3, LS_MIN, LS_MAX, linestyles='dashed')
plt.vlines(GWH_DESIGN, LS_MIN, LS_MAX, linestyles='dashed')
plt.hlines(log(vlas), GWH_MIN, GWH_MAX, linestyles='dashed', linewidth=thin)
plt.hlines(log(ska1s), GWH_MIN, GWH_MAX, linestyles='dashed', linewidth=thin)
plt.hlines(log(ska2s), GWH_MIN, GWH_MAX, linestyles='dashed', linewidth=thin)
plt.contourf(X, Y, 100 * frac_onaxis, cmap='jet', levels=11)
cbar = plt.colorbar()
cbar.set_label(r'fraction on axis ($\theta_v \leq \theta_j$) (\%)')
plt.text(GWH_O2-15, -0.8, "O2", rotation=90)
plt.text(GWH_O3-15, -0.8, "O3", rotation=90)
plt.text(GWH_DESIGN-15, -0.8, "Design", rotation=90)
plt.text(80, log(vlas)+0.1, "VLA", rotation=0)
plt.text(80, log(ska1s)+0.1, "SKA1", rotation=0)
plt.text(80, log(ska2s)+0.1, "SKA2/ngVLA", rotation=0)
plt.savefig(f"{PLOT_DIR}/fig5b.pdf", bbox_inches = 'tight')

print("fig6a")
plt.figure()
xbins= np.linspace(0., 600., BINS)
counts2, bin_edges2 = np.histogram(gw_vla_O3['pte'], bins=xbins, density=True)
plt.plot(xbins[:-1], np.cumsum(counts2) * 600. / BINS, label="O3+VLA (ex.)",
linewidth=thick, color="black", linestyle="-")
counts3, bin_edges3 = np.histogram(gw_vla_O3['ptne'], bins=xbins, density=True)
plt.plot(xbins[:-1], np.cumsum(counts3) * 600 / BINS, label="O3+VLA (no ex.)",
linewidth=thick, color="black", linestyle="--")
plt.xlabel(r"$t_p / {\rm day}$")
plt.ylabel(r"$N_{\rm joint}(< t_p)/N_{\rm joint}$")
plt.legend(loc="lower right")
plt.xlim([0., 500.])
plt.ylim([0,1])
plt.savefig(f"{PLOT_DIR}/fig6a.pdf", bbox_inches='tight')
plt.close("All")

print("fig6b")
plt.figure()
xbins= np.linspace(-5, 5, BINS)
counts, bin_edges = np.histogram(gw_O2['lpf'], bins=xbins, density=True)
plt.plot(xbins[:-1], counts , label="O2",
linewidth=thick, color="blue", linestyle="-")
counts2, bin_edges2 = np.histogram(gw_O3['lpf'], bins=xbins, density=True)
plt.plot(xbins[:-1], counts2, label="O3",
linewidth=thick, color="red", linestyle="--")
counts3, bin_edges3 = np.histogram(gw_design['lpf'], bins=xbins, density=True)
plt.plot(xbins[:-1], counts3, label="Design",
linewidth=thick, color="black", linestyle="-.")
plt.vlines(log(vlas), 0, 1, linestyles='dashed', label="VLA", linewidth=thin)
plt.vlines(log(ska1s), 0, 1, linestyles='dashed', label="SKA1", linewidth=thin)
plt.vlines(log(ska2s), 0, 1, linestyles='dashed', label="SKA2/ngVLA", linewidth=thin)
plt.xlabel(r"$\log F_p / {\rm \mu Jy}$")
plt.ylabel(r"$1/N_{\rm GW}~{\rm d}N_{\rm GW}/{\rm d}\log F_p$")
plt.legend(loc="upper left")
plt.xlim([-4., 4.])
plt.ylim([0, 0.4])
plt.text(log(vlas)-0.2, 0.025, "VLA", rotation=90)
plt.text(log(ska1s)-0.2, 0.025, "SKA1", rotation=90)
plt.text(log(ska2s)-0.2, 0.025, "SKA2", rotation=90)
plt.savefig(f"{PLOT_DIR}/fig6b.pdf", bbox_inches='tight')
plt.close("All")

print("fig6c")
plt.figure()
xbins= np.linspace(0., 1, BINS)
counts, bin_edges = np.histogram(gw_vla_O2['pm'], bins=xbins, density=True)
plt.plot(xbins[1:-1], counts[1:] , label=r"O2 ($s_{VLBI} = 15\mu{\rm Jy}$)",
linewidth=thick, color="blue", linestyle="-")
counts2, bin_edges2 = np.histogram(gw_vla_O3['pm'], bins=xbins, density=True)
plt.plot(xbins[1:-1], counts2[1:], label="O3",
linewidth=thick, color="red", linestyle="--")
counts3, bin_edges3 = np.histogram(gw_vla_design['pm'], bins=xbins, density=True)
plt.plot(xbins[1:-1], counts3[1:], label="Design",
linewidth=thick, color="black", linestyle="-.")
plt.xlabel(r"$\mu_{\rm max} / {\rm mas.month}^{-1}$")
plt.ylabel(r"$1/N_{\rm joint}~{\rm d}N_{\rm joint}/{\rm d}\mu_{\rm max}$")
plt.legend(loc="upper right")
plt.xlim([0., 0.5])
plt.ylim([0, 10])
plt.vlines(0.1, 0, 20, linestyles='dashed', linewidth=thin)
plt.savefig(f"{PLOT_DIR}/fig6c.pdf", bbox_inches='tight')
plt.close("All")

N = 200
M = int(f_vla_O3 * N)
r = gw_O3
s = gw_vla_O3

print("fig7a")
plt.figure()
f = lambda x: sqrt(1 - x ** 2) * sin(x)/(1 - sqrt(1 - x ** 2) * cos(x))
xl = np.linspace(0, 1, 200)
plt.plot(xl * 180 / Pi, 50 * f(xl), color="blue", linewidth=thin, linestyle="dashed")
plt.scatter(r['tv'][4:4 + N] * 180 / Pi, r['d'][4:4 + N], color="black", alpha=0.5,
    marker=".", s=80)
plt.scatter(s['tv'][-4 - M: -4] * 180 / Pi, s['d'][-4 - M: -4], color="red",
    marker="+")
plt.xlabel(r"$\theta_v / {\rm deg}$")
plt.ylabel(r"$D / {\rm Mpc}$")
plt.ylim([0, 143])
plt.savefig(f"{PLOT_DIR}/fig7a.pdf", bbox_inches='tight')

print("fig7b")
plt.figure()
plt.scatter(r['tv'][4:4 + N] * 180 / Pi, r['lpf'][4:4 + N], color="black", alpha=0.5,
    marker=".", s=80)
plt.scatter(s['tv'][-4 - M: -4] * 180 / Pi, s['lpf'][-4 - M: -4], color="red",
    marker="+")
plt.xlabel(r"$\theta_v / {\rm deg}$")
plt.ylabel(r"$\log(F_p / \mu{\rm Jy})$")
plt.savefig(f"{PLOT_DIR}/fig7b.pdf", bbox_inches='tight')

print("fig7c")
plt.figure()
plt.scatter(r['tv'][4:4 + N] * 180 / Pi, r['pte'][4:4 + N], color="black", alpha=0.5,
    marker=".", s=80)
plt.scatter(s['tv'][-4 - M: -4] * 180 / Pi, s['pte'][-4 - M: -4], color="red",
    marker="+")
plt.xlabel(r"$\theta_v / {\rm deg}$")
plt.ylabel(r"$t_p / {\rm day}$")
plt.savefig(f"{PLOT_DIR}/fig7c.pdf", bbox_inches='tight')

print("fig7d")
plt.figure()
plt.scatter(r['d'][4:4 + N], 1.e-6 + r['pm'][4:4 + N], color="black", alpha=0.5,
    marker=".", s=80)
plt.scatter(s['d'][-4 - M: -4], 1.e-6 + s['pm'][-4 - M: -4], color="red",
    marker="+")
plt.hlines(0.1, 0, 143, linestyle="dashed", color="blue", linewidth=thin)
plt.yscale('log')
plt.ylim([0.01, 10])
plt.xlim([0, 143])
plt.ylabel(r"$\mu_{\rm max} / {\rm mas/month}$")
plt.xlabel(r"$D / {\rm Mpc}$")
plt.savefig(f"{PLOT_DIR}/fig7d.pdf", bbox_inches='tight')

print(f"fig8a")
P = len(gw_vla_O3['ln'])
Q = len(gw_O3['ln'])
R = len(df['ln'])
plt.figure()
xbins=10 ** np.linspace(50, 53, BINS)
counts, bin_edges = np.histogram(gw_O3['e0'], bins=xbins)
plt.plot(xbins[:-1], counts / (Q * bin_edges[:-1]), label="GW-only",
linewidth=thick, color="black")
counts2, bin_edges2 = np.histogram(gw_vla_O3['e0'], bins=xbins)
plt.plot(xbins[:-1], counts2 / (Q * bin_edges2[:-1]), label="O3+VLA",
linewidth=thick, color="red")
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r"$E_{\rm c, iso}/{\rm erg}$")
plt.ylabel(r"$1/N_{\rm GW}~{\rm d}N/{\rm d}E_{\rm c,iso}$")
plt.legend(loc="lower left")
plt.xlim([bin_edges[0], bin_edges[-1]])
plt.savefig(f"{PLOT_DIR}/fig8a.pdf", bbox_inches='tight')
plt.close()

print("fig8b")
plt.figure()
xbins= np.linspace(-5.5, -0.5, BINS)
counts, bin_edges = np.histogram(gw_O3['ln'], bins=xbins, density=True)
plt.plot(xbins[:-1], counts , label="GW-only",
linewidth=thick, color="black", linestyle="-")
counts3, bin_edges3 = np.histogram(gw_vla_O3['ln'], bins=xbins, density=True)
plt.plot(xbins[:-1], counts3 * f_vla_O3, label="O3+VLA",
linewidth=thick, color="red", linestyle="-.")
plt.xlabel(r"$\log n/{\rm cm}^{-3}$")
plt.ylabel(r"$1/N_{\rm GW}~{\rm d}N/{\rm d}\log n$")
plt.legend(loc="upper left")
plt.xlim([-5,-1])
plt.ylim([0,0.6])
plt.savefig(f"{PLOT_DIR}/fig8b.pdf", bbox_inches='tight')
plt.close("All")

print("fig9a")
plt.figure()
plt.plot(e0l, frac["005"], linestyle="-", color="black", label=r"$\theta_j = 0.05$")
plt.plot(e0l, frac["010"], linestyle="-", color="blue", label=r"$\theta_j = 0.1$")
plt.plot(e0l, frac["015"], linestyle="-", color="red", label=r"$\theta_j = 0.15$")
plt.hlines(g16frac["005"], 10**50, 10**52, linestyle="--", color="black", linewidth=thin)
plt.hlines(g16frac["010"], 10**50, 10**52, linestyle="--", color="blue", linewidth=thin)
plt.hlines(g16frac["015"], 10**50, 10**52, linestyle="--", color="red", linewidth=thin)
plt.xscale('log')
plt.xlabel(r"$E_m / {\rm erg}$")
plt.ylabel(r"$N_{\rm joint} / N_{\rm GW}$")
plt.xlim([10*50, 10* 52])
plt.savefig(f"{PLOT_DIR}/fig9a.pdf", bbox_inches='tight')

print("fig9b")
plt.figure()
plt.plot(e0l, ang["005"], linestyle="-", color="black", label=r"$\theta_j = 0.05$")
plt.plot(e0l, ang["010"], linestyle="-", color="blue", label=r"$\theta_j = 0.1$")
plt.plot(e0l, ang["015"], linestyle="-", color="red", label=r"$\theta_j = 0.15$")
plt.hlines(g16ang["005"], 10**50, 10**52, linestyle="--", color="black", linewidth=thin)
plt.hlines(g16ang["010"], 10**50, 10**52, linestyle="--", color="blue", linewidth=thin)
plt.hlines(g16ang["015"], 10**50, 10**52, linestyle="--", color="red", linewidth=thin)
plt.xscale('log')
plt.xlabel(r"$E_m / {\rm erg}$")
plt.ylabel(r"$<\theta_v>_{\rm obs}$")
plt.xlim([10*50, 10* 52])
plt.savefig(f"{PLOT_DIR}/fig9b.pdf", bbox_inches='tight')

print("fig9c")
plt.figure()
plt.plot(e0l, f150["005"], linestyle="-", color="black", label=r"$\theta_j = 0.05$")
plt.plot(e0l, f150["010"], linestyle="-", color="blue", label=r"$\theta_j = 0.1$")
plt.plot(e0l, f150["015"], linestyle="-", color="red", label=r"$\theta_j = 0.15$")
plt.hlines(g16f150["005"], 10**50, 10**52, linestyle="--", color="black", linewidth=thin)
plt.hlines(g16f150["010"], 10**50, 10**52, linestyle="--", color="blue", linewidth=thin)
plt.hlines(g16f150["015"], 10**50, 10**52, linestyle="--", color="red", linewidth=thin)
plt.xscale('log')
plt.xlabel(r"$E_m / {\rm erg}$")
plt.ylabel(r"$N_{\rm joint}(t_p < 150~{\rm day}) / N_{\rm joint}$")
plt.xlim([10*50, 10* 52])
plt.legend(loc='lower left')
plt.savefig(f"{PLOT_DIR}/fig9c.pdf", bbox_inches='tight')

exit(0)
