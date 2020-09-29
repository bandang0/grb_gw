'''Hide the jet in the GRB70817 data.'''

from os import system
from itertools import product

from helpers import *

num_str = lambda x: "{" + "".join(["{},".format(k) for k in range(x)]) + "{}".format(x) + "}"

E_ = [50., 50.5, 51., 51.5, 52., 53]
n_ = [-2.5]
eb_ = [-4,] #-3.75, -3.5, -3.25, -3, -2.75, -2.5]
theta_obs_ = [ -.46]
theta_jet_ = [ -1]
N = 0
for k, param in enumerate(product(E_, n_, eb_, theta_obs_, theta_jet_)):
    err = hide_jet_err(jet_dic(param))
    system('cp data/tmp.out data/tmp{}.out'.format(k))
    print("{}, err = {}".format(jet_dic(param), err))
    N = k

system("./over_plot.py tmp{}".format(num_str(N)))
