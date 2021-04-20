from sys import argv
from numpy import log10 as log
from helpers import Gpc

case = argv[1]

ldl_before = 27
ldl_after = log(6.6 * Gpc)
delta = ldl_after - ldl_before

with open(f"data/{case}.best", 'r') as conf_file:
    lines = conf_file.readlines()
    s1_ = float(lines[0])
    s2_ = float(lines[1])
    lts_ = float(lines[2])
    lls_ = float(lines[3])
    le_ = float(lines[4])
    ti_ = float(lines[5])
    n_ = float(lines[6])
    tej_ = float(lines[7])
    tau_ = float(lines[8])
    g_ = float(lines[9])

with open(f"data/{case}.best", 'w') as conf_file:
    conf_file.write(f"{s1_:4g}\n")
    conf_file.write(f"{s2_:4g}\n")
    conf_file.write(f"{lts_:4g}\n")
    conf_file.write(f"{lls_ + 2 * delta:4g}\n")
    conf_file.write(f"{le_ + 2 * delta:4g}\n")
    conf_file.write(f"{ti_:4g}\n")
    conf_file.write(f"{n_:4g}\n")
    conf_file.write(f"{tej_:4g}\n")
    conf_file.write(f"{tau_:4g}\n")
    conf_file.write(f"{g_:4g}\n")
