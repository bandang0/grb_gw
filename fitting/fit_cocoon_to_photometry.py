"""Fit a cocoon to the photometry points."""

from os import system
import numpy as np
import matplotlib.pyplot as plt
import lmfit

from helpers import *

p = lmfit.Parameters()
p.add('n', value=1.e-3, min = 1.e-6, max = 1.)
p.add('eb', value=1.e-3, min = 1.e-6, max = 1.)
p.add('E0', value=5.e51, min = 1.e46, max = 1.e55)
p.add('um', value=1.48899627, max=3.19)
#p.add('uM', value=10**BEST_LOGGB_MAX, min=0.4)

def residual(p):
    conf_dic = {'L': 1,
                'n': p['n'].value,
                'eb': p['eb'].value,
                'E0': p['E0'].value,
                'um': p['um'].value,
                'uM': 3.19}
    print(conf_dic)
    write_conf_file(TMP_CONF_FILE, conf_dic)

    system("./remnant_light_curve tmp > /dev/null")

    return chi2_array()

# create Minimizer
mini = lmfit.Minimizer(residual, p, nan_policy='omit')

#optimal parameters
out = mini.minimize(method='nelder')
print(lmfit.fit_report(out))

#confidence intervals
ci, trace = lmfit.conf_interval(mini, out, sigmas=[1, 2],
                                trace=True, verbose=False)
lmfit.printfuncs.report_ci(ci)

#plots
plot_type = 1

if plot_type == 0:
    plt.plot(x, y)
    plt.plot(x, residual(out2.params) + y)

elif plot_type == 1:
    cx, cy, grid = lmfit.conf_interval2d(mini, out2, 'n', 'eb', 30, 30)
    plt.contourf(cx, cy, grid, np.linspace(0, 1, 11))
    plt.xlabel('n')
    plt.colorbar()
    plt.ylabel('eb')

elif plot_type == 2:
    cx, cy, grid = lmfit.conf_interval2d(mini, out2, 'a1', 't2', 30, 30)
    plt.contourf(cx, cy, grid, np.linspace(0, 1, 11))
    plt.xlabel('a1')
    plt.colorbar()
    plt.ylabel('t2')

elif plot_type == 3:
    cx1, cy1, prob = trace['a1']['a1'], trace['a1']['t2'], trace['a1']['prob']
    cx2, cy2, prob2 = trace['t2']['t2'], trace['t2']['a1'], trace['t2']['prob']
    plt.scatter(cx1, cy1, c=prob, s=30)
    plt.scatter(cx2, cy2, c=prob2, s=30)
    plt.gca().set_xlim((2.5, 3.5))
    plt.gca().set_ylim((11, 13))
    plt.xlabel('a1')
    plt.ylabel('t2')

if plot_type > 0:
    plt.show()
