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
GHz         = (1.e9)              # s^{-1}  
XRT_0 = 0.3 * keV / hPlanck
XRT_1 = 30. * keV / hPlanck
BAT_0 = 15. * keV / hPlanck
BAT_1 = 150. * keV / hPlanck
XRT_c = (30. - 0.3) * keV / hPlanck
BAT_c = (150. - 15.) * keV / hPlanck


cmap = colors.LinearSegmentedColormap.from_list('custom',rainbow(np.linspace(0, 1, 10)))
ls_l = ['-', ':', '--', '-.']
colors_l = [u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', u'#8c564b', u'#e377c2', u'#7f7f7f', u'#bcbd22', u'#17becf']

def magAB(F):
    '''AB magnitude from spectral flux. F in erg/s/cm2/Hz'''
    return - 2.5 * log(F) - 48.60


def lsub(a, x, y):
    """Return the value of y closest to the result of an interpolation of y over x"""
    return y[min(range(len(x)), key=lambda t:np.abs(a - x[t]))]

def polyfit(x, y, sigma):
    '''Fit the (x, y) points in a straight line with errors sigma.'''

    T = sum([1 / s ** 2 for s in sigma])
    bar = lambda l: sum([l[i] / sigma[i] ** 2 for i in range(len(x))]) / T

    m = (bar(x * y) - bar(x) * bar(y))/(bar(x * x) - bar(x) ** 2)
    c = bar(y) - m * bar(x)
    delta_m = 1 / (T * (bar(x * x) - bar(x) ** 2))
    delta_c = bar(x * x) * delta_m

    return m, c, sqrt(delta_m), sqrt(delta_c)

def deltaphi(chic, r, chi):
    '''Length of phi interval of intersection of patch with parallel at chi.'''
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
    '''Normalized broken-power law function.'''
    k = 1 / (nup * (1/(a+1) - 1/(b+1)))
    if nu < nup:
        return k * (nu/nup) ** a
    else:
        return k * (nu/nup) ** b

def ipl(nu, nup, a, b):
    '''A primitive of bpl, such that \int_x^y bpl = ipl(y) - ipl(x)'''
    k = 1 / (1 / (a + 1) - 1 / (b + 1))
    if nu < nup:
        return k * (nu / nup) ** (a + 1) / (a + 1)
    else:
        return k * (1 / (a + 1) + ((nu / nup) ** (b + 1) - 1) / (b + 1))

def simple_hle_patch_spectral(t, nuobs, chic, r, te, re, g, Eiso, nup, a, b):
    '''HLE from circular patch centered on chic, radius r'''
    beta = sqrt(1 - 1 / g ** 2)
    alpha = np.arccos((te - t) * cLight / re)
    doppler = 1 / (g * (1 - beta * cos(alpha)))
    return ((Eiso * cLight / (4 * Pi * g * re)) * deltaphi(chic, r, alpha)
            * bpl(nuobs / doppler, nup, a, b)
            * doppler ** 2)

def simple_hle_patch_band(t, nuobs, chic, r, te, re, g, Eiso, nup, a, b, nu0, nu1):
    '''HLE from circular patch centered on chic, radius r, integrated on [nu0, nu1]'''
    beta = sqrt(1 - 1 / g ** 2)
    alpha = np.arccos((te - t) * cLight / re)
    doppler = 1 / (g * (1 - beta * cos(alpha)))
    return ((Eiso * cLight / (4 * Pi * g * re)) * deltaphi(chic, r, alpha)
            * (ipl(nu1 / doppler, nup, a, b) - ipl(nu0 / doppler, nup, a, b))
            * doppler ** 3)

def simple_hle_patch_band_zero(t, nuobs, chic, r, te, re, g, Eiso, nup, a, b, nu0, nu1):
    tm = te - re * cos(max(0., chic - r)) / cLight
    tM = te - re * cos(min(Pi, chic + r)) / cLight
    if t < tm or t > tM:
        return 0
    else:
        return simple_hle_patch_band(t, nuobs, chic, r, te, re, g, Eiso, nup, a, b, nu0, nu1)

def peak_patch_luminosity_spectral(nuobs, chic, r, te, re, g, Eiso, nup, a, b):
    '''Approximate peak spectral luminosity of emission from patch'''
    beta = sqrt(1 - 1 / g ** 2)
    if chic - r <0.:
        return ((Eiso * cLight / (2 * g * re)) * bpl(nuobs * g * (1 - beta), nup, a, b)
                 / (g * (1 - beta))**2)
    else:
        tm = te - re * cos(max(0., chic - r)) / cLight
        tM = te - re * cos(min(Pi, chic + r)) / cLight
        T = 10 ** np.linspace(log(tm), log(tM), 7)
        return np.max([simple_hle_patch_spectral(t, nuobs, chic, r, te, re, g, Eiso, nup, a, b) for t in T])

def peak_patch_luminosity_band(nuobs, chic, r, te, re, g, Eiso, nup, a, b, nu0, nu1):
    '''Approximate peak band luminosity of emission from patch'''
    beta = sqrt(1 - 1 / g ** 2)
    tm = te - re * cos(max(0., chic - r)) / cLight
    tM = te - re * cos(min(Pi, chic + r)) / cLight
    T = 10 ** np.linspace(log(tm), log(tM), 7)
    return np.max([simple_hle_patch_band(t, nuobs, chic, r, te, re, g, Eiso, nup, a, b, nu0, nu1) for t in T])

def toy_afterglow_wind(t):
    '''Toy XRT band afterglow model with ESD, plateau and regular decay'''
    maxi = 1.e45 / ((30 - 0.3) * keV / hPlanck)
    if t < 1.e3:
        return maxi * (t/1.e3)**(-3)
    if t > 1.e6:
        return maxi * (t/1.e6)**(-1)
    else:
        return maxi

def toy_esd(t):
    return 1.e48 * (t/100)**(-3.) / XRT_c
