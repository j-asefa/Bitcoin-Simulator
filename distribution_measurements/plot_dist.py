import sys
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
from scipy import stats, integrate

from math import floor

import seaborn as sns
sns.set(color_codes=True)

interval = 5 # seconds

plt.ylim(0, 0.02)
plt.xlim(0, 1400)

plt.xlabel("Number of invs per " + str(interval) + ' seconds')
plt.ylabel('Density')


filename = sys.argv[1]
raw_data = []

plt.title(filename)

with open(filename, 'r') as f:
    for line in f:
	raw_data = [int(x) for x in line.split(';')[:-1]]
	break

min_value = min(raw_data)
max_value = max(raw_data)
nbins = (max_value - min_value) / interval + 1
data = nbins  * [0]

for val in raw_data:
    point = val - min_value
    data[int(floor(point / interval))]+=1


sns_plot = sns.distplot(data,kde=False,norm_hist=True,bins=nbins/32).get_figure()


sns_plot.savefig("output_" + filename)
