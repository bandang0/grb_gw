'''Fit the GRB70817 data with a jet using the amoeba optimizing algorithm (NR 10.4).'''

import numpy as np
from matplotlib import pyplot as plt

from helpers import *


initial_simplex = [ np.array([-3.8, 51.25, 0.1, -3, -1.6, -1.6]),
                    np.array([-3.6, 51.5, 0.5, -2.8, -1.4, -1.4]),
                    np.array([-3.4, 51.75, 0.8, -2.6, -1.2, -1.2]),
                    np.array([-3.2, 52., 1, -2.2, -1, -1]),
                    np.array([-3, 52.25, 1.3, -2., -0.8, -0.8]),
                    np.array([-2.8, 52.5, 1.6, -1.8, -0.6, -0.6]),
                    np.array([-2.6, 52.75, 2., -1.6, -0.4, -0.4])]

current_err = 1.e34

a_alpha = 1.
a_gamma = 1.
a_rho = .5
a_sigma = .5

simplex = initial_simplex

best_params = simplex[0]
best_err = current_err
n_try = 0

while(best_err > 7.e-3 and n_try < MAX_TRIES):
    n_try = n_try + 1
    maps = [(param, fit_err(jet_dic(param))) for param in simplex]
    sorted_simplex = sorted(maps, key = lambda x: x[1]) # order according to error
    current_err = sorted_simplex[-1][1]
    best = sorted_simplex[0]
    if best[1] < best_err:
        best_err = best[1]
        best_params = best[0]
    worst = sorted_simplex[-1]
    second_worst = sorted_simplex[-2]
    print("{}, err = {:.8e}".format(jet_dic(worst[0]), worst[1]))

    centroid = sum([a[0] for a in sorted_simplex[:-1]])/(len(sorted_simplex) - 1)
    reflected = centroid + a_alpha * (centroid - worst[0])
    fr = fit_err(jet_dic(reflected))
    if best[1] <= fr and fr < second_worst[1]:
        simplex = [a[0] for a in sorted_simplex[:-1]]
        simplex.append(reflected)
        continue

    if fr < best[1]:

        expanded = centroid + a_gamma * (reflected - centroid)
        fe = fit_err(jet_dic(expanded))
        if fe < fr:
            simplex = [a[0] for a in sorted_simplex[:-1]]
            simplex.append(expanded)
        else:
            simplex = [a[0] for a in sorted_simplex[:-1]]
            simplex.append(reflected)
        continue

    contracted = centroid + a_rho * (worst[0] - centroid)
    fc = fit_err(jet_dic(contracted))
    if fc < worst[1]:
        simplex = [a[0] for a in sorted_simplex[:-1]]
        simplex.append(contracted)
        continue
    simplex = [x[0] + a_sigma * (x[0] - best[0]) for x in sorted_simplex[1:]]
    simplex.append(best[0])

print("best : {}, err = {:.8e}".format(best_params, best_err))

write_conf_file("data/jet_best.conf", jet_dic(best_params))
