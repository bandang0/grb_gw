'''Determine the scalings of t_obs_max and F_nu_max with parameters.'''

from sys import exit
import numpy as np
import matplotlib.pyplot as plt
import subprocess
from helpers import *

num_str = lambda x: "{" + "".join(["{},".format(k) for k in range(x)]) + "{}".format(x) + "}"
around = lambda x: np.linspace(x - 1., x + 1., 60)

logE_d = 52.
logn_d = -3.
logeb_d = -3.
logtv_d = log(0.38)
logtj_d = -1

logE_ = around(logE_d)
logn_ = around(logn_d)
logeb_ = around(logeb_d)
logtv_ = np.linspace(-5, log(Pi / 2), 60)
logtj_ = np.linspace(-5, 0, 60)


print("tv:")
fluxes = list()
times = list()
for logeb in logeb_:
    write_conf_file(TMP_CONF_FILE, {'L': 0,
                'n': 10. ** logn_d,
                'eb': 10. ** logeb,
                'E0': 10. ** logE_d,
                'tv': 10. ** logtv_d,
                'tj': 10. ** logtj_d})
    _, t1, f1 = tuple(map(float, subprocess.check_output(["./remnant_detection.exe", "tmp"]).split()))
    fluxes.append(f1)
    times.append(t1)

times = np.array(times)
fluxes = np.array(fluxes)
plt.figure()
plt.plot(logeb_, log(times), label='time')
plt.plot(logeb_, log(fluxes), label='flux')
plt.legend()
plt.show()
