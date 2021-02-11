''''Kilonova dynamics and light curves.'''
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

from helpers import *

def 
# Our case
D0 = 120.e24
MvR = 0.04 * Msun
MvB = 0.012 * Msun
kR = 10
kB = 0.5
vR = 0.15 * cLight
vB = 0.3 * cLight
nuobs = cLight / 0.70e-6
nuir = cLight / 2.2e-6
time = day * 10 ** np.linspace(-8, 1, 100000, dtype=np.float64)
_, RR, TR, _, _ = shell_dyn(time, MvR, vR, kR)
_, RB, TB, _, _ = shell_dyn(time, MvB, vB, kB)
FR = fBB(nuobs, TR, RR, D0)
FB = fBB(nuobs, TB, RB, D0)
f = ac(0, 0.4) / Pi
F = f * FB + (1 - f) * FR
plt.xscale("log")
plt.yscale("log")
plt.plot(time / day, FR, color="red", label="Red (B)")
plt.plot(time / day, FB, color="blue", label="Blue (B)")
plt.plot(time / day, F, color="black", label="Total")
plt.xlabel("time (day)")
plt.ylabel("mag")
plt.legend()
plt.show()
