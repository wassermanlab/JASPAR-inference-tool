#!/usr/bin/env python


python regression.py -p ./0.05/pairwise.pickle -o 0.05/ -v
python regression.py -p ./0.01/pairwise.pickle -o 0.01/ -v
python regression.py -p ./0.001/pairwise.pickle -o 0.001/ -v



# Defaults for plotting
import matplotlib.pyplot as plt 
import pandas as pd
import seaborn as sns
plt.rc("font", size=12)
sns.set(style="white")
sns.set(style="whitegrid", color_codes=True)