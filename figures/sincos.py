'''GRB/GW detection probability.'''

from scipy import integrate
import numpy as np
import matplotlib.pyplot as plt

def func(y, x, theta):
    return np.abs(0.5 * (1 + np.cos(theta) ** 2) * np.cos(2 * x) * np.cos(2 * y) - np.cos(theta) * np.sin(2 * x) * np.sin(2 * y))
N = 30

theta_g = np.linspace(0., np.pi/2, N)

def p(theta):
    print(theta)
    i = integrate.dblquad(func, 0, 2*np.pi, lambda x: 0, lambda x: 2*np.pi, args=[theta], epsrel=1e-3)[0]
    print(i)
    return i
y = np.array([p(t) for t in theta_g])
y = y / (sum(y) * np.pi * 0.5 / N)
plt.plot(theta_g, y, label = r"$p(\theta)$")
plt.show()
