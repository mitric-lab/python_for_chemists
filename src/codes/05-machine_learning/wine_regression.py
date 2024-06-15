#!/usr/bin/env python

### ANCHOR: imports
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
### ANCHOR_END: imports

### ANCHOR: load_data
url = 'https://archive.ics.uci.edu/ml/machine-learning-databases/wine-quality/winequality-red.csv'
df = pd.read_csv(url, sep=';')
print(df.shape) # (1599, 12)
print(df.head())
### ANCHOR_END: load_data

#print(df.head().to_markdown())

### ANCHOR: split_data
#df.