'''Make plots for populations paper.'''

from sys import exit
import matplotlib.pyplot as plt
import pandas as pd
import corner
from numpy import sin, cos
import numpy as np
from helpers import *
from sys import argv

DATA_FILE = "data/g16_f.data"
PLOT_DIR = f"plots/aa18_plots"
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
GWH_MIN = 100.
LS_MIN = -0.5
LS_MAX = 1.5
vlas = 15.
skas = 3.
ska2s = 0.3
rad10s = 10.
ngvlas = 1.
fnx_fnr = 1.80e-5
fx_fnr = 4.23e-19

print("fig9a")
lw=3.5
wple=[50, 50.5, 51, 51.5, 52]
frac1 = [2.5, 5, 7.5, 11, 15]
frac2 = [3.5, 7, 9.5, 14, 20]
frac3 = [5.5, 10, 13.5, 17, 26]

plt.figure()
plt.plot(wple, frac1, "-", label=r"$\theta_j = 3~{\rm deg}$", color="blue")
plt.plot(wple, frac2, "-", label=r"$\theta_j = 6~{\rm deg}$", color="black")
plt.plot(wple, frac3, "-", label=r"$\theta_j = 12~{\rm deg}$", color="red")
plt.axhline(y=20, linestyle="-", color="blue", linewidth=lw)
plt.axhline(y=30, linestyle="-", color="black", linewidth=lw)
plt.axhline(y=40, linestyle="-", color="red", linewidth=lw)
plt.xlabel(r"${\rm Log}(E_m/{\rm erg})$")
plt.ylabel(r"$N_{\rm joint} / N_{\rm GW}$ (\%)")
plt.legend()
plt.savefig(f"{PLOT_DIR}/fig9a.pdf", bbox_inches = 'tight')

ta1 = [2.5, 5, 7.5, 11, 15]
ta2 = [3.5, 7, 9.5, 14, 20]
ta3 = [5.5, 10, 13.5, 17, 26]
plt.figure()
plt.plot(wple, ta1, "-", label=r"$\theta_j = 3~{\rm deg}$", color="blue")
plt.plot(wple, ta2, "-", label=r"$\theta_j = 6~{\rm deg}$", color="black")
plt.plot(wple, ta3, "-", label=r"$\theta_j = 12~{\rm deg}$", color="red")
plt.axhline(y=20, linestyle="-", color="blue", linewidth=lw)
plt.axhline(y=30, linestyle="-", color="black", linewidth=lw)
plt.axhline(y=40, linestyle="-", color="red", linewidth=lw)
plt.xlabel(r"${\rm Log}(E_m/{\rm erg})$")
plt.ylabel(r"$\langle \theta_v \rangle$ (deg)")
plt.legend()
plt.savefig(f"{PLOT_DIR}/fig9b.pdf", bbox_inches = 'tight')

plt.figure()
na1 = [2.5, 5, 7.5, 11, 15]
na2 = [3.5, 7, 9.5, 14, 20]
na3 = [5.5, 10, 13.5, 17, 26]
nao1 = [2, 5.8, 7.2, 11.6, 17]
nao2 = [3, 7.8, 9.2, 14.6, 22]
nao3 = [5, 10.8, 13.2, 17.6, 28]
plt.plot(wple, na1, "-", label=r"$\theta_j = 3~{\rm deg}$", color="blue")
plt.plot(wple, na2, "-", label=r"$\theta_j = 6~{\rm deg}$", color="black")
plt.plot(wple, na3, "-", label=r"$\theta_j = 12~{\rm deg}$", color="red")
plt.axhline(y=20, linestyle="-", color="blue", linewidth=lw)
plt.axhline(y=30, linestyle="-", color="black", linewidth=lw)
plt.axhline(y=40, linestyle="-", color="red", linewidth=lw)
plt.plot(wple, nao1, ":",color="blue")
plt.plot(wple, nao2, ":",color="black")
plt.plot(wple, nao3, ":", color="red")
plt.axhline(y=24, linestyle=":", color="blue", linewidth=lw)
plt.axhline(y=34, linestyle=":", color="black", linewidth=lw)
plt.axhline(y=44, linestyle=":", color="red", linewidth=lw)
plt.xlabel(r"${\rm Log}(E_m/{\rm erg})$")
plt.ylabel(r"${\rm f}(t_p \leq 100~{\rm days})$")
plt.legend()
plt.savefig(f"{PLOT_DIR}/fig9c.pdf", bbox_inches = 'tight')

# Read event data and create pandas data frame.
print(f"Collecting data from {DATA_FILE}")
df = pd.read_csv(DATA_FILE, sep=" ",
        names=['ti', 'd', 'n', 'e0', 'eb',
        'tv', 'tj', 'pt', 'pf'],
        dtype={'ti': int, 'd': float, 'n': float, 'e0': float,
        'eb': float, 'tv': float, 'tj': float, 'pt': float, 'pf': float})
df = df.loc[df['eb'] < 0.8]
# Add other columns to frame
df['lpt'] = log(df['pt'])
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
gw_ska_design = gw_design.loc[gw_design['pf'] > skas]
gw_ska2_design = gw_design.loc[gw_design['pf'] > ska2s]
gw_ska_max = gw_max.loc[gw_max['pf'] > skas]

print('tv')
print(np.mean(gw_vla_O2['tv']) * 180 / Pi)
print(np.mean(gw_vla_O3['tv']) * 180 / Pi)
print(np.mean(gw_vla_design['tv']) * 180 / Pi)
print(np.mean(gw_ska_design['tv']) * 180 / Pi)
print(np.mean(gw_ska2_design['tv']) * 180 / Pi)

print('frac')
print(float(len(gw_vla_O2))/len(gw_O2))
print(float(len(gw_vla_O3))/len(gw_O3))
print(float(len(gw_vla_design))/len(gw_design))
print(float(len(gw_ska_design))/len(gw_design))
print(float(len(gw_ska2_design))/len(gw_design))

print(f"fig5")
P = len(gw_vla_O3['ln'])
Q = len(gw_O3['ln'])
R = len(df['ln'])
plt.figure()
plt.hist([df['ln']],
        #label=['GW-only'],
        histtype='step', bins=BINS, weights=np.ones(R) / float(Q))
plt.hist([gw_vla_O3['ln']],
        #label=['O3+VLA'],
        histtype='step', bins=BINS,
        weights= float(R) * np.ones(P) / float(Q * Q))
plt.xlabel(r"${\rm Log}(n/{\rm cm}^{-3})$")
plt.ylabel(r"$1/N_{\rm GW}~{\rm d}N/{\rm d}{\rm Log}~n$")
plt.xlim([-6,0])
plt.savefig(f"{PLOT_DIR}/fig5.pdf", bbox_inches='tight')
plt.close()

print(f"fig6")
plt.figure()
xbins=10 ** np.linspace(50, 53, BINS)
counts, bin_edges = np.histogram(gw_O3['e0'], bins=xbins)
plt.plot(xbins[:-1], counts / (Q * bin_edges[:-1]), label="GW-only",
drawstyle="steps", linewidth=0.5)
counts2, bin_edges2 = np.histogram(gw_vla_O3['e0'], bins=xbins)
plt.plot(xbins[:-1], counts2 / (Q * bin_edges2[:-1]), label="O3+VLA",
drawstyle="steps", linewidth=0.5)
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r"$E_{\rm c, iso}/{\rm erg}$")
plt.ylabel(r"$1/N_{\rm GW}~{\rm d}N/{\rm d}E_{\rm c,iso}$")
plt.legend(loc="lower left")
#plt.ylim([1.e-2, 2.1])

plt.savefig(f"{PLOT_DIR}/fig6.pdf", bbox_inches='tight')
plt.close()

# Fractions as function of horizon
print("Functions of horizon.")
gwh_l = np.linspace(GWH_MIN, GWH_MAX, 35)
sens_l = np.linspace(LS_MIN, LS_MAX, 25)
X, Y = np.meshgrid(gwh_l, sens_l)
mean_tv = np.zeros_like(X)
std_tv = np.zeros_like(X)
frac = np.zeros_like(X)
mean_lpf = np.zeros_like(X)
n_joint = np.zeros_like(X)
frac_onaxis = np.zeros_like(X)
n_jointO3 = len(gw_vla_O3)
n_jointO2 = len(gw_vla_O2)
e0_l = ['500', '505', '510', '515', '520']
mean_fluxes = dict([[a, np.zeros_like(gwh_l)] for a in e0_l])
dev_fluxes = dict([[a, np.zeros_like(gwh_l)] for a in e0_l])
fraction_l = dict([[a, np.zeros_like(gwh_l)] for a in e0_l])
for i, gwh in enumerate(gwh_l):
    ah = gwh / 1.58
    total_events = len(df.loc[df['d'] < ah])
    for j, s in enumerate(sens_l):
        tmp_df = df[(1. + 6 * cos(df['tv']) ** 2 + cos(df['tv']) ** 4\
            > 8 * df['d'] ** 2 / ah ** 2) & (df['lpf'] > s)]
        on_axis_df = tmp_df.loc[tmp_df['tv'] < 0.1]
        mean_tv[j, i] = tmp_df['tv'].mean()
        std_tv[j, i] = np.std(tmp_df['tv'])
        frac[j, i] = float(len(tmp_df)) / total_events
        n_joint[j, i] = len(tmp_df)
        frac_onaxis[j, i] = float(len(on_axis_df))/len(tmp_df)

    for e0 in e0_l:
        df_e0 = pd.read_csv(f"data/wp15_{e0}.data", sep=" ",
            names=['ti', 'd', 'n', 'e0', 'eb',
            'tv', 'tj', 'pt', 'pf'],
            dtype={'ti': int, 'd': float, 'n': float, 'e0': float,
            'eb': float, 'tv': float, 'tj': float, 'pt': float, 'pf': float})
        df_e0['lpf'] = log(df_e0['pf'])
        total_events = len(df_e0.loc[df_e0['d'] < ah])
        tmp_gw = df_e0.loc[1. + 6 * cos(df_e0['tv']) ** 2 \
                + cos(df_e0['tv']) ** 4\
                > 8 * df_e0['d'] ** 2 / ah ** 2]
        joint_events = float(len(tmp_gw.loc[tmp_gw['pf'] > vlas]))
        mean_fluxes[e0][i] = tmp_gw['lpf'].mean()
        dev_fluxes[e0][i] = np.std(tmp_gw['lpf'])
        fraction_l[e0][i] = joint_events / total_events

print("fig1")
plt.figure()
plt.xlabel("$H$ / Mpc")
plt.ylabel(r"$\log(s/\mu {\rm Jy})$")
plt.vlines(GWH_O2, LS_MIN, LS_MAX, linestyles='dashed')
plt.vlines(GWH_O3, LS_MIN, LS_MAX, linestyles='dashed')
plt.vlines(GWH_DESIGN, LS_MIN, LS_MAX, linestyles='dashed')
plt.hlines(log(vlas), GWH_MIN, GWH_MAX, linestyles='dashed')
plt.hlines(log(skas), GWH_MIN, GWH_MAX, linestyles='dashed')
plt.hlines(log(ngvlas), GWH_MIN, GWH_MAX, linestyles='dashed')
plt.contourf(X, Y, mean_tv * 180 / Pi, cmap='jet')
cbar = plt.colorbar()
cbar.set_label(r'mean observed $\theta_v$ (deg)')
plt.savefig(f"{PLOT_DIR}/fig1.pdf", bbox_inches = 'tight')

print("fig45")
plt.figure()
plt.xlabel("$H$ / Mpc")
plt.ylabel(r"$\log(s/\mu {\rm Jy})$")
plt.vlines(GWH_O2, LS_MIN, LS_MAX, linestyles='dashed')
plt.vlines(GWH_O3, LS_MIN, LS_MAX, linestyles='dashed')
plt.vlines(GWH_DESIGN, LS_MIN, LS_MAX, linestyles='dashed')
plt.hlines(log(vlas), GWH_MIN, GWH_MAX, linestyles='dashed')
plt.hlines(log(skas), GWH_MIN, GWH_MAX, linestyles='dashed')
plt.hlines(log(ngvlas), GWH_MIN, GWH_MAX, linestyles='dashed')
plt.contourf(X, Y, 100 * frac_onaxis, cmap='jet')
cbar = plt.colorbar()
cbar.set_label(r'fraction on axis ($\theta_v \leq \theta_j$) (\%)')
plt.savefig(f"{PLOT_DIR}/fig45.pdf", bbox_inches = 'tight')

print("fig2")
plt.figure()
plt.xlabel("Horizon / Mpc")
plt.ylabel("log(Radio sensitivity/muJy)")
plt.vlines(GWH_O2, LS_MIN, LS_MAX, linestyles='dashed')
plt.vlines(GWH_O3, LS_MIN, LS_MAX, linestyles='dashed')
plt.vlines(GWH_DESIGN, LS_MIN, LS_MAX, linestyles='dashed')
plt.hlines(log(vlas), GWH_MIN, GWH_MAX, linestyles='dashed')
plt.hlines(log(skas), GWH_MIN, GWH_MAX, linestyles='dashed')
plt.hlines(log(ngvlas), GWH_MIN, GWH_MAX, linestyles='dashed')
plt.contourf(X, Y, std_tv * 180 / Pi, cmap='jet')
cbar = plt.colorbar()
cbar.set_label(r'observed $\theta_v$ width (deg)')
plt.savefig(f"{PLOT_DIR}/fig2.pdf", bbox_inches = 'tight')

print("fig3")
plt.figure()
plt.xlabel("Horizon / Mpc")
plt.ylabel("log(Radio sensitivity/muJy)")
plt.vlines(GWH_O2, LS_MIN, LS_MAX, linestyles='dashed')
plt.vlines(GWH_O3, LS_MIN, LS_MAX, linestyles='dashed')
plt.vlines(GWH_DESIGN, LS_MIN, LS_MAX, linestyles='dashed')
plt.hlines(log(vlas), GWH_MIN, GWH_MAX, linestyles='dashed')
plt.hlines(log(skas), GWH_MIN, GWH_MAX, linestyles='dashed')
plt.hlines(log(ngvlas), GWH_MIN, GWH_MAX, linestyles='dashed')
plt.contourf(X, Y, frac * 100, cmap='jet')
cbar = plt.colorbar()
cbar.set_label("fraction of joint detections in total events (%)")
plt.savefig(f"{PLOT_DIR}/fig3.pdf", bbox_inches='tight')

print("fig4")
plt.figure()
plt.xlabel("D / Mpc")
plt.ylabel("$\log(F_p/\mu Jy)$")
plt.hlines(log(vlas), 0, AH_O3, linestyles='dashed')
plt.hlines(log(skas), 0, AH_O3, linestyles='dashed')
plt.hlines(log(ngvlas), 0, AH_O3, linestyles='dashed')
plt.hist2d(gw_O3['d'], gw_O3['lpf'], cmap="jet", bins=BINS)
plt.colorbar()
plt.savefig(f"{PLOT_DIR}/fig4.pdf", bbox_inches='tight')


print("fig8")
plt.figure()
lab = {'500': '50.', "505": '50.5', '510': '51.', "515": '51.5', '520': "52.0"}
for e0 in e0_l:
    plt.plot(gwh_l, mean_fluxes[e0], label=r"$\log(E_m/{\rm erg})$ = " + lab[e0])
plt.xlabel("$H$ / Mpc")
plt.ylabel("Mean $F_p$ of GW events ($\mu$Jy)")
plt.hlines(log(vlas), GWH_MIN, GWH_MAX, linestyles='dashed')
plt.hlines(log(skas), GWH_MIN, GWH_MAX, linestyles='dashed')
plt.hlines(log(ngvlas), GWH_MIN, GWH_MAX, linestyles='dashed')
plt.legend()
plt.savefig(f"{PLOT_DIR}/fig8.pdf", bbox_inches='tight')
plt.close("All")

print(f"fig9")
plt.figure()
plt.hist([gw_max['d'] / AH_MAX, gw_vla_O3['d'] / AH_O3,
    gw_ska_design['d'] / AH_DESIGN],
        label=['GW only', 'Joint (O3+VLA)', 'Joint (Design+SKA)'],
        histtype='step', bins=BINS, density=True)
plt.xlabel(r"$D / \bar{H}$")
plt.ylabel(r"d$N$/d($D/\bar{H}$)")
plt.legend(loc="upper left")
plt.savefig(f"{PLOT_DIR}/fig9.pdf", bbox_inches='tight')
plt.close()

print("fig10")
plt.figure()
plt.xlabel("$H$ / Mpc")
plt.ylabel(r"$\log(s/\mu {\rm Jy})$")
plt.vlines(GWH_O2, LS_MIN, LS_MAX, linestyles='dashed')
plt.vlines(GWH_O3, LS_MIN, LS_MAX, linestyles='dashed')
plt.vlines(GWH_DESIGN, LS_MIN, LS_MAX, linestyles='dashed')
plt.hlines(log(vlas), GWH_MIN, GWH_MAX, linestyles='dashed')
plt.hlines(log(skas), GWH_MIN, GWH_MAX, linestyles='dashed')
plt.hlines(log(ngvlas), GWH_MIN, GWH_MAX, linestyles='dashed')
plt.contourf(X, Y, log(3. * n_joint / n_jointO3) , cmap='jet')
cbar = plt.colorbar()
cbar.set_label("$log(N_J)$")
plt.savefig(f"{PLOT_DIR}/fig10.pdf", bbox_inches='tight')

#print(f"fig11")
#plt.figure()
#frame = pd.read_csv("pf_tv.data", sep=" ",
#        names=['tv', 'ts', 'fs', 'ta', 'fa'])
#plt.plot(180 * frame['tv'] / Pi, frame['fs'], label="simplified")
#plt.plot(180 * frame['tv'] / Pi, frame['fa'], label="analytical")
#frame = pd.read_csv("data/structuredjet_daigne_better.data",
#        sep=" ", names=['tv', 'a', 'b', 'c', 'd', 'e', 'fc'])
#plt.plot(frame['tv'], 1000 * frame['fc'], label="complete")
#plt.xlabel(r"$\theta_v$ (deg)")
#plt.ylabel("$F_P$ $\mu Jy$")
#plt.yscale("log")
#plt.xscale("log")
#plt.legend()
#plt.savefig(f"{PLOT_DIR}/fig11.pdf", bbox_inches='tight')
#plt.close()

print(f"fig12")
plt.figure()
plt.hist([gw_O2['lpf'], gw_O3['lpf'], gw_design['lpf']],
        label=['O2', 'O3', 'Design'],
        histtype='step', bins=BINS, density=True)
plt.xlabel("$\log(F_P / \mu Jy)$")
plt.ylabel("$dN/d\log(F_P)$")
plt.legend()
plt.savefig(f"{PLOT_DIR}/fig12.pdf", bbox_inches='tight')
plt.close()

print(f"fig13")
plt.figure()
plt.hist([180 * df['tv'] / Pi,
        180 * gw_O2['tv'] / Pi,
        180 * gw_vla_O3['tv']/Pi],
        label=['All events', 'GW only', ' Joint (O3+VLA)'],
        histtype='step', bins=BINS, density=True, linewidth=1.2)
plt.xlabel(r"$\theta_v$ (deg)")
plt.ylabel(r"d$N$/d$\theta_v$")
plt.legend()
plt.savefig(f"{PLOT_DIR}/fig13.pdf", bbox_inches='tight')
plt.close()

print("fig14")
plt.figure()
plt.xlabel("D / Mpc")
plt.ylabel(r"$\theta_v$")
plt.hist2d(gw_vla_O3['d'], 180 * gw_vla_O3['tv'] / Pi, cmap="jet", bins=BINS,
        label="$d^2 N/dDd\theta_v$")
plt.colorbar()
plt.savefig(f"{PLOT_DIR}/fig14.pdf", bbox_inches='tight')

# print("fig15")
# plt.figure()
# plt.xlabel("D / Mpc")
# plt.ylabel("log(PM / (mas/d))")
# plt.hist2d(gw_vla_O3['d'], gw_vla_O3['laas'], cmap="jet", bins=BINS,
#         label = "$d^2 N / dDdPM$")
# plt.colorbar()
# plt.savefig(f"{PLOT_DIR}/fig15.pdf", bbox_inches='tight')

print("fig16")
plt.figure()
for e0 in e0_l:
    plt.plot(gwh_l, 100 * fraction_l[e0], label="$\log(E_m) = $" + f"{e0}")
plt.xlabel("H/Mpc")
plt.ylabel("\% fraction of joint detections with VLA")
plt.vlines([GWH_O2, GWH_O3, GWH_DESIGN],
        0, 15, linestyles='dashed')
plt.legend()
plt.savefig(f"{PLOT_DIR}/fig16.pdf", bbox_inches='tight')
plt.close("All")

# print("fig17")
# fig = plt.figure()
# ax1 = fig.add_subplot(111)
# ax2 = ax1.twiny()
# ax1.set_xlabel("$\log(F_X/\mu $Jy$)$")
# ax1.hist([gw_O2['lpfx'], gw_O3['lpfx'],
#         gw_design['lpfx']],
#         label=['O2', 'O3', 'Design'],
#         histtype='step', bins=BINS, density=True)
# new_tick_locations = ax1.get_xticks()

# def tick_function(X):
#     V = X - 10.63
#     return ["%.1f" % z for z in V]
#
# ax2.set_xlim(ax1.get_xlim())
# ax2.set_xticks(new_tick_locations)
# ax2.set_xticklabels(tick_function(new_tick_locations))
# ax2.set_xlabel(r"$\log(F_{[0.3-10]keV}~/~{\rm erg/s/cm}^2)$")
# ax1.legend()
# plt.savefig(f"{PLOT_DIR}/fig17.pdf", bbox_inches='tight')

#
# print("fig20")
# plt.figure()
#
# plt.plot([40.], [2.],
#         "s", label="O2+VLA")
# plt.plot(gw_vla_O3['d'][2243:2245], gw_vla_O3['lpf'][2243:2245],
#         "s", label="O3+VLA")
# plt.plot(gw_ska_max['d'][2243:2264], gw_ska_max['lpf'][2243:2264],
#         "s", label="700Mpc horizon + SKA")
# plt.xlabel("Distance (Mpc)")
# plt.ylabel("$\log(F_p / \mu$Jy$)$")
# plt.legend()
# plt.savefig(f"{PLOT_DIR}/fig20.pdf", bbox_inches='tight')
