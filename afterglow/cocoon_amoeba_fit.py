'''Fit the GRB70817 data with a cocoon using the amoeba optimizing algorithm (NR 10.4).'''

import numpy as np

from helpers import *

initial_simplex = [ np.array([-3.3, 49., -3.2, 0.05, 0.1]),
                    np.array([-3.2, 49.5, -3.1, 0.1, 0.5]),
                    np.array([-3.1, 50., -3., 0.15, 0.55]),
                    np.array([-3., 50.5, -2.9, 0.2, 0.6]),
                    np.array([-2.9, 51., -2.8, 0.25, 0.65]),
                    np.array([-2.8, 51, -2.7, 0.3, 0.8])]

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
    maps = [(param, fit_err(cocoon_dic(param))) for param in simplex]
    sorted_simplex = sorted(maps, key = lambda x: x[1]) # order according to error
    current_err = sorted_simplex[-1][1]
    best = sorted_simplex[0]
    if best[1] < best_err:
        best_err = best[1]
        best_params = best[0]
    worst = sorted_simplex[-1]
    second_worst = sorted_simplex[-2]
    print("{}, err = {:.8e}".format(cocoon_dic(worst[0]), worst[1]))

    centroid = sum([a[0] for a in sorted_simplex[:-1]])/(len(sorted_simplex) - 1)
    reflected = centroid + a_alpha * (centroid - worst[0])
    fr = fit_err(cocoon_dic(reflected))
    if best[1] <= fr and fr < second_worst[1]:
        simplex = [a[0] for a in sorted_simplex[:-1]]
        simplex.append(reflected)
        continue

    if fr < best[1]:
        expanded = centroid + a_gamma * (reflected - centroid)
        fe = fit_err(cocoon_dic(expanded))
        if fe < fr:
            simplex = [a[0] for a in sorted_simplex[:-1]]
            simplex.append(expanded)
        else:
            simplex = [a[0] for a in sorted_simplex[:-1]]
            simplex.append(reflected)
        continue

    contracted = centroid + a_rho * (worst[0] - centroid)
    fc = fit_err(cocoon_dic(contracted))
    if fc < worst[1]:
        simplex = [a[0] for a in sorted_simplex[:-1]]
        simplex.append(contracted)
        continue

    simplex = [x[0] + a_sigma * (x[0] - best[0]) for x in sorted_simplex[1:]]
    simplex.append(best[0])

print("best : {}, err = {:.8e}".format(best_params, best_err))

write_conf_file("data/cocoon_best.conf", cocoon_dic(best_params))
