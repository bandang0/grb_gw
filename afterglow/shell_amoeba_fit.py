'''Fit the GRB70817 data with a shell using the amoeba optimizing algorithm (NR 10.4).'''

import numpy as np
from helpers import *

initial_simplex = [ np.array([-4.84, 50.95, 0.45, -2.05]),
                    np.array([-4.86, 51., 0.5, -2.0]),
                    np.array([-4.88, 51.05, 0.55, -1.95]),
                    np.array([-4.92, 51.1, 0.6, -1.9]),
                    np.array([-4.94, 51.15, 0.65, -1.85])]

current_err = 1.e34

a_alpha = 1.
a_gamma = 1.
a_rho = .5
a_sigma = .5

simplex = initial_simplex
#print(simplex)
best_params = simplex[0]
best_err = current_err
n_try = 0

while(best_err > 7.e-3 and n_try < MAX_TRIES):
    #print(simplex)
    #print("err = {:.2e}".format(current_err))
    n_try = n_try + 1
    maps = [(param, fit_err(shell_dic(param))) for param in simplex]
    sorted_simplex = sorted(maps, key = lambda x: x[1]) # order according to error
    current_err = sorted_simplex[-1][1]
    best = sorted_simplex[0]
    if best[1] < best_err:
        best_err = best[1]
        best_params = best[0]
    worst = sorted_simplex[-1]
    second_worst = sorted_simplex[-2]
    print("{}, err = {:.8e}".format(shell_dic(worst[0]), worst[1]))

    centroid = sum([a[0] for a in sorted_simplex[:-1]])/(len(sorted_simplex) - 1)
    reflected = centroid + a_alpha * (centroid - worst[0])
    fr = fit_err(shell_dic(reflected))
    if best[1] <= fr and fr < second_worst[1]:
        #print("case 1")
        simplex = [a[0] for a in sorted_simplex[:-1]]
        simplex.append(reflected)
        continue

    if fr < best[1]:
        #print("case 2")
        expanded = centroid + a_gamma * (reflected - centroid)
        fe = fit_err(shell_dic(expanded))
        if fe < fr:
            #print("case 21")
            simplex = [a[0] for a in sorted_simplex[:-1]]
            simplex.append(expanded)
        else:
            #print("case 22")
            simplex = [a[0] for a in sorted_simplex[:-1]]
            simplex.append(reflected)
        continue

    contracted = centroid + a_rho * (worst[0] - centroid)
    fc = fit_err(shell_dic(contracted))
    if fc < worst[1]:
        #print("case 3")
        #print(simplex)
        #print(simplex[:-1])
        simplex = [a[0] for a in sorted_simplex[:-1]]
        simplex.append(contracted)
        #print(simplex)
        continue

    #print("case 4")
    simplex = [x[0] + a_sigma * (x[0] - best[0]) for x in sorted_simplex[1:]]
    simplex.append(best[0])

print("best : {}, err = {:.8e}".format(best_params, best_err))

write_conf_file("data/shell_best.conf", shell_dic(best_params))
