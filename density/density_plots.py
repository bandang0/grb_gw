'''Make plots for density paper.'''

from sys import exit
import matplotlib.pyplot as plt
import pandas as pd
from numpy import sin, cos
import numpy as np
from scipy.interpolate import interp1d
from helpers import *
from sys import argv
import corner

PLOT_DIR = "plots/aa19_plots"
plt.style.use('presentation')
BINS = 70

N = 100000
p2_l = np.linspace(0., 1., N)

def ptilde(p2, frac1, frac2):
    return (1 - p2) * frac1 / ((1 - p2) * frac1 + p2 * frac2)

def posterior(frac1, frac2, n1, n2, Nt, NHD, lis):
    ret = np.array([ptilde(p2, frac1, frac2) ** (Nt - NHD) \
            * (1 - ptilde(p2, frac1, frac2)) ** (NHD) \
            for p2 in lis])
    return len(ret) * ret / np.sum(ret)

print("fig1")
plt.figure()
plt.xlabel(r"$i$ [deg]")
plt.ylabel(r"$n~[{\rm cm}^{-3}]$")
plt.yscale('log')

P = 200
x_l = np.linspace(0.1, 60., P)
x_l2 = np.linspace(24, 28, P)
x_l3 = np.linspace(10, 45, 10)
U = 1.e-3 * 10 ** (-4.02 + 5.78 * log(x_l / 5.73) + 1.5)
M = 1.e-3 * 10 ** (-4.02 + 5.78 * log(x_l / 5.73))
D = 1.e-3 * 10 ** (-4.02 + 5.78 * log(x_l / 5.73) - 1.5)
D2 = 1.e-3 * 10 ** (-4.02 + 5.78 * log(x_l2 / 5.73) - 1.5)
plt.plot(x_l, D, '--', color = "black")
plt.plot(x_l, M, '-', color = "black", label="JAG")
plt.plot(x_l, U, '--', color = "black")
plt.fill_between(x_l, D, U, alpha = 0.3, color="black")

plt.axvline(x=32-8, linestyle="dashed", color = "red")
plt.axvline(x=32, linestyle="solid", color = "red", label="GW")
plt.axvline(x=32+6, linestyle="dashed", color = "red")
plt.axvspan(32-8, 32+6, alpha = 0.3, color = "red")

plt.axvline(x=14, linestyle="dashed", color = "blue")
plt.axvline(x=21, linestyle="solid", color = "blue", label="VLBI")
plt.axvline(x=28, linestyle="dashed", color = "blue")
plt.axvspan(14, 28, alpha = 0.3, color = "blue")

#plt.plot(x_l, 1.e-4 * np.ones(P), '--', color = 'green')
plt.plot(x_l3, 10**(-3) * np.ones(10), 'v', color = 'green', label="KNAG")
#plt.plot(x_l, 1.e-2 * np.ones(P), '--', color = 'green')
plt.fill_between(x_l, 1.e-8 * np.ones(P),
    10**(-3) * np.ones(P), alpha = 0.3, color = "green")
plt.fill_between(x_l2, D2, [1.e-3 for a in x_l2], alpha = 1., color="purple")

plt.xlim([10., 45.])
plt.ylim([1.e-6, 1.e-1])
plt.legend()
plt.savefig(f"{PLOT_DIR}/fig1.pdf", bbox_inches = 'tight')

print(f"Collecting data")
g16 = pd.read_csv("data/g16_frac_n.data", sep=" ",
        names=['logn0', 'a', 'b', 'c', 'd', 'e', 'radio', 'r', 'x', 'm'])
wp15 = pd.read_csv("data/wp15_frac_n.data", sep=" ",
        names=['logn0', 'a', 'b', 'c', 'd', 'e', 'radio', 'r', 'x', 'm'])
fracg16 = interp1d(g16['logn0'], g16['m'], kind='cubic')
fracwp = interp1d(wp15['logn0'], wp15['m'], kind='cubic')

print(f"Collecting data (aligned)")
g16_aligned = pd.read_csv("data/g16_frac_n_aligned.data", sep=" ",
        names=['logn0', 'a', 'b', 'c', 'd', 'e', 'radio', 'r', 'x', 'm'])
fracg16_aligned = interp1d(g16_aligned['logn0'], g16_aligned['m'], kind='cubic')

print(f"Collecting data (aligned far)")
g16_far = pd.read_csv("data/g16_frac_n_aligned_far.data", sep=" ",
        names=['logn0', 'a', 'b', 'c', 'd', 'e', 'radio', 'r', 'x', 'm'])
fracg16_far = interp1d(g16_far['logn0'], g16_far['m'], kind='cubic')


print("fig4a")
plt.figure()
plt.xlabel(r"$\log n~[{\rm cm}^{-3}]$")
plt.ylabel(r"$r(n)$ [\%]")
plt.plot(g16['logn0'], 100 * g16['radio'],color="blue",
            label=r"Radio")
plt.plot(g16['logn0'], 100 * g16['r'], color="green",
            label=r"$r$")
plt.plot(g16['logn0'], 100 * g16['x'], color="orange",
            label=r"X")
plt.plot(g16['logn0'], 100 * g16['m'],color="black",
            label=r"Multi-$\lambda$")
plt.legend()
plt.title("Recovery of GW+AG (G16)")
plt.xlim([-5, 3.5])
plt.ylim([0, 100])

plt.savefig(f"{PLOT_DIR}/fig4a.pdf", bbox_inches = 'tight')

print("fig4a_aligned")
plt.figure()
plt.xlabel(r"$\log n~[{\rm cm}^{-3}]$")
plt.ylabel(r"$r(n)$ [\%]")
plt.plot(g16_aligned['logn0'], 100 * g16_aligned['radio'],color="blue",
            label=r"Radio")
plt.plot(g16_aligned['logn0'], 100 * g16_aligned['r'], color="green",
            label=r"$r$")
plt.plot(g16_aligned['logn0'], 100 * g16_aligned['x'], color="orange",
            label=r"X")
plt.plot(g16_aligned['logn0'], 100 * g16_aligned['m'],color="black",
            label=r"Multi-$\lambda$")
plt.legend()
plt.title(r"Recovery of GW+AG, with $\theta_v = 0$ (G16)")
plt.xlim([-5, 3.5])
plt.ylim([0, 100])
plt.savefig(f"{PLOT_DIR}/fig4a_aligned.pdf", bbox_inches = 'tight')

print("fig4a_far")
plt.figure()
plt.xlabel(r"$\log n~[{\rm cm}^{-3}]$")
plt.ylabel(r"$r(n)$ [\%]")
plt.plot(g16_far['logn0'], 100 * g16_far['radio'],color="blue",
            label=r"Radio")
plt.plot(g16_far['logn0'], 100 * g16_far['r'], color="green",
            label=r"$r$")
plt.plot(g16_far['logn0'], 100 * g16_far['x'], color="orange",
            label=r"X")
plt.plot(g16_far['logn0'], 100 * g16_far['m'],color="black",
            label=r"Multi-$\lambda$")
plt.legend()
plt.title(r"Recovery of AG only, with $\theta_v = 0$" + "\n" + r"and $D_L <$ 6 Gpc (G16)")
plt.xlim([-5, 3.5])
plt.ylim([0, 100])

plt.savefig(f"{PLOT_DIR}/fig4a_far.pdf", bbox_inches = 'tight')
print("fig4b")
plt.figure()
plt.xlabel(r"$\log n~[{\rm cm}^{-3}]$")
plt.ylabel(r"$r(n)$ [\%]")
plt.plot(wp15['logn0'], 100 * wp15['radio'], color="blue",
            label=r"Radio")
plt.plot(wp15['logn0'], 100 * wp15['r'], color="green",
            label=r"$r$")
plt.plot(wp15['logn0'], 100 * wp15['x'], color="orange",
            label=r"X")
plt.plot(wp15['logn0'], 100 * wp15['m'], color="black",
            label=r"Multi-$\lambda$")
plt.xlim([-5, 3.5])
plt.ylim([0, 100])
plt.legend()
plt.title("Recovery of GW+AG (WP15)")
plt.savefig(f"{PLOT_DIR}/fig4b.pdf", bbox_inches = 'tight')

AH_O3 = 143.
vlas = 15.

print("df0")
df0 = pd.read_csv("data/g16_0.data", sep=" ",
        names=['ig', 'iv', 'd', 'lpf', 'lpt'])
gw0_vla_O3 = df0.loc[df0['ig'] * df0['iv'] == 1]

print("df3")
df3 = pd.read_csv("data/g16_3.data", sep=" ",
        names=['ig', 'iv', 'd', 'lpf', 'lpt'])
gw3_vla_O3 = df3.loc[df3['ig'] * df3['iv'] == 1]

print("fig2")
col_to_plot = ['d','lpf', 'lpt']
fig = corner.corner(gw0_vla_O3[col_to_plot],
                    labels=['$D$ [Mpc]', r'$\log F_{\rm p,3GHz}$ [$\mu$Jy]',
                            r'$\log T_{\rm p}$ [d]'],
                    label_kwargs={'fontsize':'large'},
                    contour_kwargs={'linewidths':0.9,'levels':6},
                    range=[(0, 143), (1, 5), (-0.5, 4)],
                    bins=BINS,
                    color="blue",
                    plot_density=False,
                    plot_datapoints=False,
                    plot_contours=True, smooth=True,
                    hist_kwargs={"density": True, "linewidth":1.})

corner.corner(gw3_vla_O3[col_to_plot],
                    labels=['$D$ [Mpc]', r'$\log F_{\rm p,3GHz}$ [$\mu$Jy]',
                            r'$\log T_{\rm p}$ [d]'],
                    label_kwargs={'fontsize':'large'},
                    contour_kwargs={'linewidths':0.9, 'levels': 6},
                    range=[(0, 143), (log(15.), 5), (-0.5, 4)],
                    bins=BINS,
                    color="orange",
                    plot_density=False,
                    plot_datapoints=False,
                    plot_contours=True,fig=fig, smooth=True,
                    hist_kwargs={"density": True, "linewidth":1.})

plt.savefig(f"{PLOT_DIR}/fig2.pdf", bbox_inches = 'tight')

exit(0)
n1 = -3.
n2 = 0.
f1g16 = fracg16(n1)
f2g16 = fracg16(n2)
f1wp = fracwp(n1)
f2wp = fracwp(n2)

print("fig5a")
plt.figure()
plt.xlabel(r"$f_{\rm HD}$")
plt.ylabel(r"$p(f_{\rm HD} | f_{\rm HD}^{\rm obs})$")
#plt.xscale('log')
plt.plot(p2_l, posterior(0.5, 0.5, n1, n2, 10, 5, p2_l), ":",
color = "blue", linewidth=1.2)
plt.plot(p2_l, posterior(f1g16, f2g16, n1, n2, 10, 1, p2_l), "-",
color="green", linewidth=1.2, label=r"$f_{\rm HD}^{\rm obs}$ = 1/10")
plt.plot(p2_l, posterior(f1g16, f2g16, n1, n2, 10, 3, p2_l), "-",
color="red", linewidth=1.2, label="3/10")
plt.plot(p2_l, posterior(f1g16, f2g16, n1, n2, 10, 5, p2_l), "-",
color= "blue", linewidth=1.2, label="5/10")
plt.plot(p2_l, posterior(f1g16, f2g16, n1, n2, 1, 0, p2_l), "--",
color="grey", label="0/1", linewidth=1.2)
plt.xlim([0., .5])
plt.ylim([0., 12])
plt.legend(loc="upper left")

plt.savefig(f"{PLOT_DIR}/fig5a.pdf", bbox_inches = 'tight')

print("fig5b")
plt.figure()
plt.xlabel(r"$f_{\rm HD}$")
plt.ylabel(r"$p(f_{\rm HD} | f_{\rm HD}^{\rm obs})$")
#plt.xscale('log')
plt.plot(p2_l, posterior(0.5, 0.5, n1, n2, 10, 5, p2_l), ":",
color = "blue", linewidth=1.2)
plt.plot(p2_l, posterior(f1wp, f2wp, n1, n2, 10, 1, p2_l), "-",
color="green", label=r"$f_{HD}^{\rm obs} = 1/10$", linewidth=1.2)
plt.plot(p2_l, posterior(f1wp, f2wp, n1, n2, 10, 3, p2_l), "-",
color="red", label="3/10", linewidth=1.2)
plt.plot(p2_l, posterior(f1wp, f2wp, n1, n2, 10, 5, p2_l), "-",
color= "blue", label="5/10", linewidth=1.2)
plt.plot(p2_l, posterior(f1wp, f2wp, n1, n2, 1, 0, p2_l), "--",
color="grey", label="0/1", linewidth=1.2)
plt.xlim([0., .5])
plt.ylim([0., 12])

plt.savefig(f"{PLOT_DIR}/fig5b.pdf", bbox_inches = 'tight')



print("fig3a")
pp1 = 1 - ptilde(0.1, f1g16, f2g16)
pp3 = 1 - ptilde(0.3, f1g16, f2g16)
print(pp1, pp3)
plt.figure()
xbins=np.linspace(0,143, BINS)
counts0, bin_edges0 = np.histogram(gw0_vla_O3['d'], bins=xbins, density=True)
counts3, bin_edges3 = np.histogram(gw3_vla_O3['d'], bins=xbins, density=True)
plt.plot(xbins[:-1], (1 - pp1) * counts3 + pp1 * counts0, label=r"$f_{HD} = 10\%$",
color="blue", linestyle="-")
plt.plot(xbins[:-1], (1 - pp3) * counts3 + pp3 * counts0, label=r"$f_{HD} = 30\%$",
color="red", linestyle="-")
plt.xlabel(r"$D$/Mpc")
plt.ylabel(r"$dN/dD$")
plt.legend(loc="upper left")

plt.savefig(f"{PLOT_DIR}/fig3a.pdf", bbox_inches = 'tight')

print("fig3b")
plt.figure()
xbins=np.linspace(log(15.),7, BINS)
counts0, bin_edges0 = np.histogram(gw0_vla_O3['lpf'], bins=xbins, density=True)
counts3, bin_edges3 = np.histogram(gw3_vla_O3['lpf'], bins=xbins, density=True)
plt.plot(xbins[:-1], (1 - pp1) * counts3 + pp1 * counts0, label=r"$f_{HD} = 10\%$",
color="blue", linestyle="-")
plt.plot(xbins[:-1], (1 - pp3) * counts3 + pp3 * counts0, label=r"$f_{HD} = 30\%$",
color="red", linestyle="-")
plt.xlabel(r"$\log F_p$/${\rm \mu Jy}$")
plt.ylabel(r"$dN/d\log F_p$")
#plt.legend(loc="upper left")

plt.savefig(f"{PLOT_DIR}/fig3b.pdf", bbox_inches = 'tight')

print("fig3c")
plt.figure()
xbins=np.linspace(-1.5,3, BINS)
counts0, bin_edges0 = np.histogram(gw0_vla_O3['lpt'], bins=xbins, density=True)
counts3, bin_edges3 = np.histogram(gw3_vla_O3['lpt'], bins=xbins, density=True)
plt.plot(xbins[:-1], (1 - pp1) * counts3 + pp1 * counts0, label=r"$f_{HD} = 10\%$",
color="blue", linestyle="-")
plt.plot(xbins[:-1], (1 - pp3) * counts3 + pp3 * counts0, label=r"$f_{HD} = 30\%$",
color="red", linestyle="-")
plt.xlabel(r"$\log T_p$/${\rm d}$")
plt.ylabel(r"$dN/d\log T_p$")
#plt.legend(loc="upper left")

plt.savefig(f"{PLOT_DIR}/fig3c.pdf", bbox_inches = 'tight')
