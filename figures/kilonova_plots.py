'''Make plots for kilonova paper.'''

from sys import exit
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from numpy import cos, sin
from helpers import *
from sys import argv

plt.style.use('palmerio')
PLOT_DIR = f"plots/aa20_plots"
BINS = 150
AFILE = "data/pif.data"
TK = 0.35
AH = 210.
f = pd.read_csv(AFILE, sep=" ", header=None)

def ai(tv, tk):
    return f.values[min(99, int(tk * 200 / Pi))][min(99, int(tv * 200 / Pi))] / Pi

ai = np.vectorize(ai)
# Read event data and create pandas data frame.
df = pd.read_csv("data/uh.data", sep=" ",
        names=['t', 'ub', 'ur', 'hb', 'hr'])

print("fig1")
a = ai(0., TK)
b = ai(Pi / 8, TK)
c = ai(Pi / 4, TK)
d = ai(Pi / 2, TK)
plt.figure()
plt.xlabel("time / day")
plt.ylabel("magU")
plt.plot(df['t'] / day, magAB(df['ub']), "--", color="blue", label="U (Blue)")
plt.plot(df['t'] / day, magAB(df['ur']), "--", color="red", label="U (Red)")
plt.plot(df['t'] / day, magAB(df['ub'] * a + (1 - a) * df['ur']), label=r"U ($\theta_v = 0.$)" )
plt.plot(df['t'] / day, magAB(df['ub'] * b + (1 - b) * df['ur']), label=r"U ($\theta_v = 0.4$)" )
plt.plot(df['t'] / day, magAB(df['ub'] * c + (1 - c) * df['ur']), label=r"U ($\theta_v = 0.8$)" )
plt.plot(df['t'] / day, magAB(df['ub'] * d + (1 - d) * df['ur']), label=r"U ($\theta_v = 1.5$)" )
plt.xlim([0., 10])
plt.ylim([17, 40])
plt.legend()
plt.gca().invert_yaxis()
plt.savefig(f"{PLOT_DIR}/fig1.pdf", bbox_inches = 'tight')

print("fig2")
plt.figure()
plt.xlabel("time / day")
plt.ylabel("magH")
plt.plot(df['t'] / day, magAB(df['ub']), "--", color ="blue", label="U (Blue)")
plt.plot(df['t'] / day, magAB(df['ur']), "--", color="red", label="U (Red)")
plt.plot(df['t'] / day, magAB(df['ub'] * ai(0., TK) + (1 - ai(0., TK)) * df['ur']), label=r"U ($\theta_v = 0.$)" )
plt.plot(df['t'] / day, magAB(df['ub'] * ai(0.3, TK) + (1 - ai(0.3, TK)) * df['ur']), label=r"U ($\theta_v = 0.3$)" )
plt.plot(df['t'] / day, magAB(df['ub'] * ai(0.6, TK) + (1 - ai(0.6, TK)) * df['ur']), label=r"U ($\theta_v = 0.6$)" )
plt.plot(df['t'] / day, magAB(df['ub'] * ai(0.9, TK) + (1 - ai(0.9, TK)) * df['ur']), label=r"U ($\theta_v = 0.9$)" )
plt.plot(df['t'] / day, magAB(df['ub'] * ai(1.2, TK) + (1 - ai(1.2, TK)) * df['ur']), label=r"U ($\theta_v = 1.2$)" )
plt.plot(df['t'] / day, magAB(df['ub'] * ai(1.5, TK) + (1 - ai(1.5, TK)) * df['ur']), label=r"U ($\theta_v = 1.5$)" )

plt.plot(df['t'] / day, 9 + magAB(df['ur']), "--", color="red", label="U (Red)")
plt.plot(df['t'] / day, 9 + magAB(df['ub']), "--", color ="blue", label="U (Blue)")
plt.plot(df['t'] / day, 9 + magAB(df['ub'] * ai(0., 2 * TK) + (1 - ai(0.,2 * TK)) * df['ur']), label=r"U ($\theta_v = 0.$)" )
plt.plot(df['t'] / day, 9 + magAB(df['ub'] * ai(0.3, 2 * TK) + (1 - ai(0.3, 2 * TK)) * df['ur']), label=r"U ($\theta_v = 0.3$)" )
plt.plot(df['t'] / day, 9 + magAB(df['ub'] * ai(0.6, 2 * TK) + (1 - ai(0.6, 2 * TK)) * df['ur']), label=r"U ($\theta_v = 0.6$)" )
plt.plot(df['t'] / day, 9 + magAB(df['ub'] * ai(0.9, 2 * TK) + (1 - ai(0.9, 2 * TK)) * df['ur']), label=r"U ($\theta_v = 0.9$)" )
plt.plot(df['t'] / day, 9 + magAB(df['ub'] * ai(1.2, 2 * TK) + (1 - ai(1.2, 2 * TK)) * df['ur']), label=r"U ($\theta_v = 1.2$)" )
plt.plot(df['t'] / day, 9 + magAB(df['ub'] * ai(1.5, 2 * TK) + (1 - ai(1.5, 2 * TK)) * df['ur']), label=r"U ($\theta_v = 1.5$)" )


plt.xlim([0.0, 4])
plt.ylim([19, 40])
plt.legend()
plt.gca().invert_yaxis()
plt.savefig(f"{PLOT_DIR}/fig2.pdf", bbox_inches = 'tight')

print("fig2b")
t = np.linspace(0, Pi / 2, 1000)
plt.figure()
plt.xlabel("theta obs / deg")
plt.plot(t * 180 / Pi, magAB(ai(t, 0.1) * 3.63e-32 + (1 - ai(t, 0.1)) * 2.29e-29), "-", label="0.1 (U)")
plt.plot(t * 180 / Pi, magAB(ai(t, 0.3) * 3.63e-32 + (1 - ai(t, 0.3)) * 2.29e-29), "-", label="0.3 (U)")
plt.plot(t * 180 / Pi, magAB(ai(t, 0.5) * 3.63e-32 + (1 - ai(t, 0.5)) * 2.29e-29), "-", label="0.5 (U)")
plt.plot(t * 180 / Pi, magAB(ai(t, 0.7) * 3.63e-32 + (1 - ai(t, 0.7)) * 2.29e-29), "-", label="0.7 (U)")

#plt.plot(t * 180 / Pi, magAB(ai(t, 0.1) * 2.29e-28 + (1 - ai(t, 0.1)) * 3.63e-27), ".", label="0.1 (H)")
#plt.plot(t * 180 / Pi, magAB(ai(t, 0.3) * 2.29e-28 + (1 - ai(t, 0.3)) * 3.63e-27), ".", label="0.3 (H)")
#plt.plot(t * 180 / Pi, magAB(ai(t, 0.5) * 2.29e-28 + (1 - ai(t, 0.5)) * 3.63e-27), ".", label="0.5 (H)")
#plt.plot(t * 180 / Pi, magAB(ai(t, 0.7) * 2.29e-28 + (1 - ai(t, 0.7)) * 3.63e-27), ".", label="0.7 (H)")
plt.ylabel("mag")
plt.legend()
plt.gca().invert_yaxis()
plt.savefig(f"{PLOT_DIR}/fig2b.pdf", bbox_inches = 'tight')
exit(0)

pop_df = pd.read_csv("data/kilonova_pop.data", sep=" ",
        names=['i', 'd', 'tv', 'bb', 'br', 'hb', 'hr'])#, 'ub2', 'ur2', 'ub3', 'ur3',
                          #'bb0', 'br0', 'bb1', 'br1', 'bb2', 'br2', 'bb3', 'br3',
                          #'ib0', 'ir0', 'ib1', 'ir1', 'ib2', 'ir2', 'ib3', 'ir3',
                          #'hb0', 'hr0', 'hb1', 'hr1', 'hb2', 'hr2', 'hb3', 'hr3'])
pop_joint = pop_df.loc[1. + 6. * cos(pop_df['tv']) ** 2 + cos(pop_df['tv']) ** 4\
            > 8 * pop_df['d'] ** 2 / AH ** 2]
pop_joint['tk'] = TK

pop_joint['f'] = ai(pop_joint['tv'], pop_joint['tk'])

pop_joint['b'] = magAB(pop_joint['bb'] * pop_joint['f'] \
                        + pop_joint['br'] * (1 - pop_joint['f']))

pop_joint['h'] = magAB(pop_joint['hb'] * pop_joint['f'] \
                        + pop_joint['hr'] * (1 - pop_joint['f']))

print("fig3")
plt.figure()
plt.hist([pop_joint['b'], pop_joint['h'],
magAB(pop_joint['bb']), magAB(pop_joint['br']),
magAB(pop_joint['hb']), magAB(pop_joint['hr'])],
label=["B", "H", "B (Blue)", "B (Red)", "H (Blue)", "H (Red)"],
density=True, histtype='step', bins=BINS)
plt.xlabel("mag")
plt.ylabel("dN/dmag")
plt.gca().invert_xaxis()
plt.legend()
plt.savefig(f"{PLOT_DIR}/fig3.pdf", bbox_inches='tight')
plt.close("all")

print("fig4")
plt.figure()
plt.hist2d(pop_joint['d'], pop_joint['b'],
cmap='jet', bins=BINS)
plt.xlabel("D / Mpc")
plt.ylabel("B")
plt.gca().invert_yaxis()
plt.savefig(f"{PLOT_DIR}/fig4.png", bbox_inches='tight')
plt.close("all")

print("fig5")
plt.figure()
plt.hist2d(pop_joint['d'], pop_joint['h'],
cmap='jet', bins=BINS)
plt.xlabel("D / Mpc")
plt.ylabel("H")
plt.gca().invert_yaxis()
plt.savefig(f"{PLOT_DIR}/fig5.png", bbox_inches='tight')
plt.close("all")

print("fig6")
plt.figure()
plt.hist2d(180 * pop_joint['tv'] / Pi, pop_joint['b'],
cmap='jet', bins=BINS)
plt.xlabel("tv / deg")
plt.ylabel("B")
plt.gca().invert_yaxis()
plt.savefig(f"{PLOT_DIR}/fig6.png", bbox_inches='tight')
plt.close("all")

print("fig7")
plt.figure()
plt.hist2d(180 * pop_joint['tv'] / Pi, pop_joint['h'],
cmap='jet', bins=BINS)
plt.xlabel("tv / deg")
plt.ylabel("H")
plt.gca().invert_yaxis()
plt.savefig(f"{PLOT_DIR}/fig7.png", bbox_inches='tight')
plt.close("all")

print("fig8")
plt.figure()
plt.hist([180 * pop_joint['tv'] / Pi],
label=["tv"], density=True, histtype='step', bins=BINS)
plt.xlabel("tv")
plt.ylabel("dN/dtv")
plt.legend()
plt.savefig(f"{PLOT_DIR}/fig8.pdf", bbox_inches='tight')
plt.close("all")

print("fig9")
plt.figure()
plt.hist([pop_joint['d']],
label=["D / Mpc"], density=True, histtype='step', bins=BINS)
plt.xlabel("D")
plt.ylabel("dN/dD")
plt.legend()
plt.savefig(f"{PLOT_DIR}/fig9.pdf", bbox_inches='tight')
plt.close("all")
