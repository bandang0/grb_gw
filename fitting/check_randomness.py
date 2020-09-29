'''Check the performance of the PRNG.'''

from sys import exit
import pandas as pd
from scipy.stats import moment
from numpy import mean
import matplotlib.pyplot as plt

# before you should do this:
# sed 's/^ *//g' | awk '/^[0-9]+ [0-9]\.[0-9]{5}E.[0-9]{2}$/{print $0}' > rn.data
#

DATA_FILE = 'random_numbers.data'

df = pd.read_csv(DATA_FILE, sep=" ", names=['ti', 'u'], header=None, dtype={'ti':int, 'u': float})

for i in range(16):
    par = df.loc[df['ti'] == i]['u']
    print(f"Proc: {i}, mean: {mean(par)}, var: {moment(par, 2)}, skew: {moment(par, 3)}, kurt: {moment(par, 4)}")

par = df['u']
print(f"Total: mean: {mean(par)}, var: {moment(par, 2)}, skew: {moment(par, 3)}, kurt: {moment(par, 4)}")
