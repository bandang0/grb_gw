"""Helper stuff for all the fitting stuff."""

import numpy as np
from os import system

TMP_CONF_FILE = "data/tmp.conf"
TMP_DATA_FILE = "data/tmp.out"
MAX_TRIES = 256
EXPLO_EXTENT = 4

log = np.log10
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

# Observations
OBS_TIMES_d = np.array([16.42,    17.39,  18.33,  22.36,  24.26,  31.22,
    46.26,  54.27,  57.22,  93.13,  115.05, 162.89, 196.79, 216.91,
    220., 256.76, 267., 294.])
OBS_FLUXES_muJy = np.array([18.7, 15.1,   14.5,   22.5,   25.6,   34.0,
    44.0,   48.,    61.0,   70.0,   89.05,  98.0,   78.9,   68.,
    64.7, 55.,    40.3, 31.2])
OBS_ERRORS_muJy = np.array([6.3,  3.9,    3.7,    3.4,    2.9,    3.6,
    4.,     6.,     9.,     5.7,    20.,    22.5,   9.,     21.,
    2.7,  12.,    2.7,  3.6])

OBS_FLUXES_mJy = 1.e-3 * OBS_FLUXES_muJy
OBS_ERRORS_mJy = 1.e-3 * OBS_ERRORS_muJy

#increase
OBS_TIMES_1_d = OBS_TIMES_d[:12]
OBS_FLUXES_1_muJy = OBS_FLUXES_muJy[:12]
OBS_ERRORS_1_muJy = OBS_ERRORS_muJy[:12]

#decrease
OBS_TIMES_2_d = OBS_TIMES_d[10:]
OBS_FLUXES_2_muJy = OBS_FLUXES_muJy[10:]
OBS_ERRORS_2_muJy = OBS_ERRORS_muJy[10:]

default_L = 0.
default_gamma_0 = 100.
default_E_0 = 5.26e50
default_epsilon_e = 0.1
default_epsilon_B = 1.11e-3
default_n = 8.87e-4
default_zeta = 1.
default_p = 2.2
default_gb_min = 1.51
default_gb_max = 3.10
default_nu_obs = 3.e9
default_N_POINTS = 1024
default_D = 1.29e26
default_alpha = 5.0
default_theta_obs = 0.
default_theta_jet = 2.0
default_first_day = 10
default_last_day = 1000

BEST_LOGN, BEST_LOGE, BEST_LOGEB, BEST_LOGGB_MIN, BEST_LOGGB_MAX \
        = -3.09, 50.75, -2.93, 0.18, 0.50

# Dictonaries for input parameters
jet_dic = lambda x: {'E_0': x[0],
                        'n': x[1],
                        'epsilon_B': x[2],
                        'theta_obs': x[3],
                        'theta_jet': 2.,
                        'L': 0.0}

cocoon_dic = lambda x: {'n': x[0],
                        'E_0': x[1],
                        'epsilon_B': x[2],
                        'gb_min': x[3],
                        'gb_max': x[4],
                        'L': 1.0}

shell_dic = lambda x: {'n': x[0],
                        'E_0': x[1],
                        'gamma_0': x[2],
                        'epsilon_B': x[3],
                        'L': 0.0}

errf = lambda x: np.exp(- x * np.sqrt(2)) - 1 if x < 0. else 0

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

    return m, c, np.sqrt(delta_m), np.sqrt(delta_c)

def add_plot(plt, file_name, **kwargs):
    a, b = data_list("data/{}.out".format(file_name))
    plt.plot(a / day, b * 1.e29, **kwargs)
