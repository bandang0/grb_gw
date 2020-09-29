#!/bin/env python3
"""Plot the after-glow data with two slope model."""

from helpers import *

F1 = lambda x: 0.0776 * 1.e-26 * (x / (24 * 3600 * 100)) ** 0.79
F2 = lambda x: 0.0426 * 1.e-26 * (x / (24 * 3600 * 300)) ** -1.3

if __name__ == "__main__":

    print("writing lines")
    with open("data/tmp.out", 'r') as conf_file:
        with open("data/slope1.out", 'w') as s1:
            with open("data/slope2.out", 'w') as s2:
                for line in conf_file.readlines()[2:-3]:
                    atoms = line.split()
                    #print(atoms)
                    t = float(atoms[2])
                    if t < 160 * 3600 * 24.:
                        s1.write("0 0 " + str(t) + " 0 0 0 0 0 0 0 0 0 0 " + str(F1(t)) + "\n")
                    if t > 140 * 24 * 3600.:
                        s2.write("0 0 " + str(t) + " 0 0 0 0 0 0 0 0 0 0 " + str(F2(t)) + "\n")
