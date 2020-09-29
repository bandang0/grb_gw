"""Generate the corner plot for the jet decrease fitting."""

import pickle
import corner

from helpers import *

with open("jet_decrease.pickle", "rb") as pickle_file:
    samples = pickle.load(pickle_file)

fig = corner.corner(samples, labels = ["log(n/cm-3)", "log(eb)", "log(E0/erg)", "tv", "tj"], color = "black", quantiles=[0.32, 0.5, 0.68], show_titles=True, title_fmt=".2f", title_kwargs={"fontsize": 12, "color":"red"})

fig.savefig("plots/corner_jet_decrease.png", bbox_inches="tight")
