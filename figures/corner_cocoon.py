"""Generate the corner plot for the cocoon fitting."""

import pickle
import corner

from helpers import *

with open("cocoon.pickle", "rb") as pickle_file:
    samples = pickle.load(pickle_file)

fig = corner.corner(samples, labels = [r"$\log(n/cm^{{-3}})$", r"$\log(\epsilon_b)$", r"$\log(E_0/$erg)", r"$u_m$", r"$u_M$"], color = "blue", quantiles = [0.32, 0.5, 0.68], show_titles=True)

fig.savefig("plots/corner_cocoon.pdf", bbox_inches="tight")
