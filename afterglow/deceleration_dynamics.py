from helpers import *
import numpy as np
from sys import argv
from matplotlib import pyplot as plt

case = argv[1]
PLOT_DIR = "plots/deceleration"
def dgamma(eps, gamma, m, mdot, mej):
    return - (mdot/mej) * (gamma ** 2 * (1 + eps) - eps * gamma) / (1 + 2 * (m/mej) * gamma)
n = 1
A = 5.e11
Eiso = 1.e53
gamma_0 = 400
mej = Eiso / (gamma_0 * cLight ** 2)
if case == "ism":
    r_dec = (3 * Eiso / (4 * np.pi * n * mProton * gamma_0 ** 2 * cLight ** 2)) ** (1 / 3)
else:
    r_dec = mej / (4 * np.pi * gamma_0 * A)

print(f"Eiso = {mej * gamma_0 * cLight ** 2:.3g} erg")
print(f"R_dec = {r_dec / pc:.3g} pc")

N = 200000

r_l = np.linspace(0, 500 * r_dec, N)
dr_l = [(r_l[i+1] - r_l[i]) for i in range(len(r_l)-1)]
if case == 'ism':
    m_l = [n * mProton * r ** 3 for r in r_l]
    eta_l = [(r / r_dec) ** 3 for r in r_l]
else:
    m_l = [4 * np.pi * A * r for r in r_l]
    eta_l = [(r / r_dec) for r in r_l]
mdot_l = [(m_l[i+1] - m_l[i])/dr_l[i] for i in range(len(r_l)-1)]
eta_l = [(r / r_dec) ** 3 for r in r_l]

fig = plt.figure()
gs = fig.add_gridspec(2, hspace=0, height_ratios=[1, 3])
axs = gs.subplots(sharex=True)


print("adiabatic")
gamma_a = [gamma_0 * (-1 + np.sqrt(1 + 4 * eta + 4 * eta ** 2 / gamma_0 ** 2)) / (2 * eta) for eta in eta_l]
dldl_a = [(log(gamma_a[i+1]) - log(gamma_a[i]))/(log(r_l[i+1]) - log(r_l[i])) for i in range(len(r_l)-1)]
dldl_a.append(dldl_a[-1])
axs[1].plot(r_l / r_dec, gamma_a, label ="abiabatic", linestyle="-", color="brown")
axs[0].plot(r_l / r_dec, dldl_a, color="brown", linestyle="-")

print("epsilon")
for k, eps in enumerate([0.001, 0.01, 0.1]):
    gamma_eq = [gamma_0]
    integral = 0.
    for i in range(1, len(r_l)-1):
        integral = integral + dr_l[i-1] * mdot_l[i-1] * gamma_eq[i-1] * (gamma_eq[i-1] - 1) / mej
        gamma_eq.append(gamma_0 * (-1 + np.sqrt(1 + 4 * eta_l[i] + 4 * eta_l[i] ** 2 / gamma_0 ** 2 - 4 * eta_l[i] * eps * integral)) / (2 * eta_l[i]))
    gamma_eq.append(gamma_eq[-1])
    dldl_eq = [(log(gamma_eq[i+1]) - log(gamma_eq[i]))/(log(r_l[i+1]) - log(r_l[i])) for i in range(len(r_l)-1)]
    dldl_eq.append(dldl_eq[-1])
    axs[1].plot(r_l / r_dec, gamma_eq, label=r"$\epsilon$" + f" = {eps}", color=colors_l[k])
    axs[0].plot(r_l / r_dec, dldl_eq, color=colors_l[k])

print("radiative")
eps = 1
gamma_rad = [gamma_0]
integral = 0.
for i in range(1, len(r_l)-1):
    integral = integral + dr_l[i-1] * mdot_l[i-1] * gamma_rad[i-1] * (gamma_rad[i-1] - 1)
    gamma_rad.append((gamma_0 * mej + m_l[i] - integral)/ (mej + m_l[i]))
gamma_rad.append(gamma_rad[-1])
dldl_rad = [(log(gamma_rad[i+1]) - log(gamma_rad[i]))/(log(r_l[i+1]) - log(r_l[i])) for i in range(len(r_l)-1)]
dldl_rad.append(dldl_rad[-1])
axs[1].plot(r_l / r_dec, gamma_rad, label=f"radiative", linestyle="--", color="red")
axs[0].plot(r_l / r_dec, dldl_rad, linestyle="--", color="red")

print("balistic")
gamma_balistic = [(gamma_0 * mej + m) / (mej + m) for m in m_l]
dldl_balistic = [(log(gamma_balistic[i+1]) - log(gamma_balistic[i]))/(log(r_l[i+1]) - log(r_l[i])) for i in range(len(r_l)-1)]
dldl_balistic.append(dldl_balistic[-1])
axs[1].plot(r_l / r_dec, gamma_balistic, label ="balistic", linestyle="-.", color="blue")
axs[0].plot(r_l / r_dec, dldl_balistic, color="blue", linestyle="-.")

axs[1].set_xlabel(r"$R/R_{\rm dec}$")
axs[1].set_ylabel(r"$\Gamma$")
axs[0].set_ylabel(r"d$\log\Gamma$/d$\log R$")
#axs[0].set_yscale('linear')
axs[1].set_xscale('log')
axs[1].set_yscale('log')
if case == 'ism':
    axs[1].set_xlim([0.1, 80])
else:
    axs[1].set_xlim([0.01, 500])
axs[0].set_ylim([-3,0])
#axs[1].set_ylim([0, 400])
plt.legend()
plt.suptitle("ISM" if case == 'ism' else 'Wind')
plt.savefig(f"{PLOT_DIR}/gamma_{case}.pdf", bbox_inches='tight')
