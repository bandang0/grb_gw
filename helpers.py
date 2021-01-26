"""Helper stuff for all the fitting stuff."""

import numpy as np
from os import system
from numpy import cos, sin, tan, arccos, sqrt
from numpy import log10 as log
from matplotlib import colors
from matplotlib.cm import rainbow


Pi = 3.141592653589793238
Deg         = (1.74532925e-2)    # rad
mas         = (4.8481368e-9)     # rad
muas        = (4.8481368e-12)     # rad
cLight      = (2.99792458e10)    # cm/s
hPlanck     = (6.6260755e-27)    # erg.s
hbar        = (1.05457266e-27)   # erg.s=h/2pi)
eElectron   = (4.8032068e-10)    # esu
mElectron   = (9.1093897e-28)    # g= (i.e. 510.99906 keV/c2)
mProton     = (1.6726231e-24)    # g= (i.e. 938.27231 MeV/c2)
uamu        = (1.6605402e-24)    # g= (unite de mass atomique)
SigmaThomson= (6.6524616e-25)    # cm2
GNewton     = (6.67259e-8)       # cgs
NAvogadro   = (6.0221367e23)     # mol
kBoltzmann  = (1.380658e-16)     # erg/K
SigmaStefan = (5.67051e-5)       # erg/s/cm2/K4
aStefan     = (7.56591e-15)      # erg/cm3/K4= (recalcule = 4Sigma/c)
au          = (1.4959787066e13)  # cm
km          = (1.0000000000e5)  # cm
pc          = (3.0856775807e18)  # cm
kpc         = (3.0856775807e21)  # cm
Mpc         = (3.0856775807e24)  # cm
Gpc         = (3.0856775807e27)  # cm
Rsun        = (6.96e10)          # cm
Msun        = (1.98892e33)       # g
Lsun        = (3.846e33)         # erg/s
day         = (86400.0e0)                 # s      = (24 h)
month       = (2.63e6)            # s
yr          = (3.156e7)          # s
ly          = (0.9461e18)        # cm
eV          = (1.60217733e-12)   # erg
keV         = (1.60217733e-9)    # erg
MeV         = (1.60217733e-6)    # erg
GeV         = (1.60217733e-3)    # erg
TeV         = (1.60217733e0)     # erg
erg         = (6.24150636e5)     # MeV
MeVperc2    = (1.78266270e-27)   # g
Jy          = (1.0e-23)           # erg/cm2/s/Hz
mJy         = (1.0e-26)           # erg/cm2/s/Hz
microJy     = (1.0e-29)           # erg/cm2/s/Hz
micron      = (1.0e-4)            # cm
XRT_c = (30. - 0.3) * keV / hPlanck
BAT_c = (150. - 15.) * keV / hPlanck
cmap = colors.LinearSegmentedColormap.from_list('custom',rainbow(np.linspace(0, 1, 10)))
ls_l = ['-', ':', '--', '-.']

def binsum(a, b):
    return np.cumsum(a * np.append(np.diff(b), 0.))

def magAB(F):
    return - 2.5 * log(F) - 48.60

def data_list(infile):
    """Parse t_obs and F_nu_obs into the listed data"""
    a = list()
    b = list()
    with open(infile, 'r') as data_file:
        lines = data_file.readlines()
        for l in range(2, len(lines) - 1):
            atoms = lines[l].split()
            a.append(float(atoms[2]))
            b.append(float(atoms[13]))
    return np.array(a), np.array(b)

def d_log_log(f, x):
    """Returns d(log f)/d(log x) as a list with the same length as x and f."""
    lf = log(f)
    lx = log(x)
    ret = [(lf[i + 1] - lf[i]) / (lx[i + 1] - lx[i]) for i in range(len(f) - 1)]
    ret.append(ret[-1])
    return ret

def write_conf_file(fname, dic):
    """write configuration file fname with configuration read in dic"""
    fmt = "{:.8e}\n"
    k = dic.keys()
    with open(fname, 'w') as tmp_conf:
        tmp_conf.write(fmt.format(default_L \
                if 'L' not in k else dic['L']))
        tmp_conf.write(fmt.format(default_n \
                if 'n' not in k else dic['n']))
        tmp_conf.write(fmt.format(default_E_0 \
                if 'E0' not in k else dic['E0']))
        tmp_conf.write(fmt.format(default_gamma_0 \
                if 'gamma_0' not in k else dic['gamma_0']))
        tmp_conf.write(fmt.format(default_epsilon_e \
                if 'ee' not in k else dic['ee']))
        tmp_conf.write(fmt.format(default_epsilon_B \
                if 'eb' not in k else dic['eb']))
        tmp_conf.write(fmt.format(default_zeta))
        tmp_conf.write(fmt.format(default_p if 'p' not in k else dic['p']))
        tmp_conf.write(fmt.format(default_nu_obs \
                if 'nu_obs' not in k else dic['nu_obs']))
        tmp_conf.write(fmt.format(default_D))
        tmp_conf.write("{}\n".format(default_N_POINTS))
        tmp_conf.write(fmt.format(default_alpha))
        tmp_conf.write(fmt.format(default_gb_min \
                if 'um' not in k else dic['um']))
        tmp_conf.write(fmt.format(default_gb_max \
                if 'uM' not in k else dic['uM']))
        tmp_conf.write(fmt.format(default_theta_obs \
                if 'tv' not in k else dic['tv']))
        tmp_conf.write(fmt.format(default_theta_jet \
                if 'tj' not in k else dic['tj']))
        tmp_conf.write(fmt.format(default_first_day \
                if 'first_day' not in k else dic['first_day']))
        tmp_conf.write(fmt.format(default_last_day \
                if 'last_day' not in k else dic['last_day']))


def STR(l):
    ret = ""
    for a in l:
        ret = ret + ("%s " % str(a))
    return ret

def lsub(a, x, y):
    """Return the value of y closest to the result of a inter polation of y over x"""
    return y[min(range(len(x)), key=lambda t:np.abs(a - x[t]))]

def fit_err(dic):
    """Calculate the fitting chi2 with configuration read in dic"""
    write_conf_file(TMP_CONF_FILE, dic)
    system("./remnant_light_curve tmp > /dev/null")

    #Initialize data lists
    t_obs = list()
    F_nu_obs = list()

    # Parse fname into the listed data
    t_obs, F_nu_obs = data_list("data/tmp.out")

    # get good units
    t_obs_days = t_obs / day
    F_nu_obs_muJy = 1.e29 * F_nu_obs

    # calculate the error
    err = 0.
    time_ind = 0
    for k in range(len(OBS_TIMES)):
        while t_obs_days[time_ind] < OBS_TIMES_d[k]:
            time_ind = time_ind + 1
        err = err + np.abs(F_nu_obs_muJy[time_ind] - OBS_FLUXES_muJy[k]) ** 2 \
                / OBS_ERRORS_muJy[k] ** 2

    return err

def calc_err():
    """Run tmp conf file and return the chi2."""
    system("./remnant_light_curve tmp > /dev/null")

    # Parse fname into the listed data
    t_obs, F_nu_obs = data_list("data/tmp.out")

    # get good units
    t_obs_days =  t_obs / day
    F_nu_obs_muJy = 1.e29 * F_nu_obs

    # calculate the summed chi2
    err = 0.
    time_ind = 0
    for k in range(len(OBS_TIMES)):
        while t_obs_days[time_ind] < OBS_TIMES_d[k]:
            time_ind = time_ind + 1
        err = err + np.abs(F_nu_obs_muJy[time_ind] - OBS_FLUXES_muJy[k]) ** 2 \
                / OBS_ERRORS_muJy[k] ** 2

    return err

def chi2_array():
    """Return the residual array of chi2 terms from tmp data file."""
    # Parse fname into the listed data
    t_obs, F_nu_obs = data_list("data/tmp.out")

    # get good units
    t_obs_days = t_obs / day
    F_nu_obs_muJy = 1.e29 * F_nu_obs

    # calculate the chi2 terms array array
    err = list()
    time_ind = 0
    for k in range(len(OBS_TIMES_d)):
        while t_obs_days[time_ind] < OBS_TIMES_d[k]:
            time_ind = time_ind + 1
        err.append((F_nu_obs_muJy[time_ind] - OBS_FLUXES_muJy[k]) \
                / OBS_ERRORS_muJy[k])

    return np.array(err)

def hide_jet_err(dic):
    """Error for not hiding and not not fitting the radio points."""
    write_conf_file(TMP_CONF_FILE, dic)
    #print("cal")
    system("./remnant_light_curve tmp > /dev/null")

    # Parse fname into the listed data
    t_obs, F_nu_obs = data_list("data/tmp.out")

    # get good units
    t_obs_days = t_obs / day
    F_nu_obs_muJy = F_nu_obs * 1.e29

    # calculate the error
    err = 0.
    time_ind = 0
    for k in range(len(OBS_TIMES)):
        while t_obs_days[time_ind] < OBS_TIMES_d[k]:
            time_ind = time_ind + 1
        err = err + errf((OBS_FLUXES_muJy[k] - F_nu_obs_muJy[time_ind] \
                - OBS_ERRORS_muJy[k])/OBS_ERRORS_muJy[k])

    return err

def t_max_F_max(dic):
    """Return the maximum flux and time calculated from dic parameters."""
    write_conf_file(TMP_CONF_FILE, dic)
    system("./remnant_light_curve tmp > /dev/null")

    # Parse fname into the listed data
    t_obs, F_nu_obs = data_list("data/tmp.out")

    # get good units
    t_obs_days = t_obs / day
    F_nu_obs_muJy = F_nu_obs * 1.e29

    # determine the maximum flux and maximum flux time
    max_ind = F_nu_obs_muJy.argmax()

    return t_obs_days[max_ind], F_nu_obs_muJy[max_ind]

def polyfit(x, y, sigma):
    '''Fit the (x, y) points in a straight line with errors sigma.'''

    T = sum([1 / s ** 2 for s in sigma])
    bar = lambda l: sum([l[i] / sigma[i] ** 2 for i in range(len(x))]) / T

    m = (bar(x * y) - bar(x) * bar(y))/(bar(x * x) - bar(x) ** 2)
    c = bar(y) - m * bar(x)
    delta_m = 1 / (T * (bar(x * x) - bar(x) ** 2))
    delta_c = bar(x * x) * delta_m

    return m, c, sqrt(delta_m), sqrt(delta_c)

def add_plot(plt, file_name, **kwargs):
    a, b = data_list("data/{}.out".format(file_name))
    plt.plot(a / day, b * 1.e29, **kwargs)

def deltaphi(chic, r, chi):
    if chic < r: #on axis case
        if chi <= r - chic:
            return 2 * Pi
        if chi >= chic + r:
            return 0.
        else:
            return 2 * arccos((cos(r) - cos(chic) * cos(chi))/(sin(chic) * sin(chi)))
    else: # off-axis case
        if (chi <= chic - r
            or chi >= chic + r
            or np.abs((cos(r) - cos(chic) * cos(chi))/(sin(chic) * sin(chi))) >= 1.):
            return 0.
        else:
            return 2 * arccos((cos(r) - cos(chic) * cos(chi))/(sin(chic) * sin(chi)))
def bpl(nu, nup, a, b):
    if nu < nup:
        return (nu/nup) ** a / (nup * (1/(a+1) - 1/(b+1)))
    else:
        return (nu/nup) ** b / (nup * (1/(a+1) - 1/(b+1)))

def simple_hle_patch(t, nuobs, chic, r, te, re, g, Eiso, nup, a, b):
    beta = sqrt(1 - 1 / g ** 2)
    return ((Eiso / 4 * Pi) * (cLight / re) * deltaphi(chic, r, np.arccos((te - t) * cLight / re))
            * bpl(nuobs * g * (1 - beta * cLight * (te - t) / re), nup, a, b)
            / (g * (1 - beta * cLight * (te - t) / re)) ** 2)

def peak_patch_luminosity(nuobs, chic, r, te, re, g, Eiso, nup, a, b):
    beta = sqrt(1 - 1 / g ** 2)
    if chic - r <0.:
        return ((Eiso * cLight) / (2 * re)
                 * bpl(nuobs * g * (1 - beta), nup, a, b)
                 / (g * (1 - beta))**2)
    else:
        theta_peak = 1 * (chic - 0.5 * r)
        return ((Eiso * cLight) / (4 * Pi * re) * deltaphi(chic, r, theta_peak)
                 * bpl(nuobs * g * (1 - beta * cos(theta_peak)), nup, a, b)
                 / (g * (1 - beta * cos(theta_peak)))**2)

def toy_esd(t):
    # in erg/s in the XRT band
    return 1.e48 * (t/1.e2) ** (-3)

def toy_plateau(t):
    maxi = 1.e46 / ((30 - 0.3) * keV / hPlanck)
    if t < 1.e3:
        return maxi * (t/1.e3)**(-3)
    if t > 1.e5:
        return maxi * (t/1.e5)**(-1)
    else:
        return maxi
