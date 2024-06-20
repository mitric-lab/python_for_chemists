#!/usr/bin/env python

### ANCHOR: imports
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
### ANCHOR_END: imports

### ANCHOR: load_data
url = 'https://archive.ics.uci.edu/ml/machine-learning-databases/wine-quality/winequality-red.csv'
df = pd.read_csv(url, sep=';')

print(type(df)) # <class 'pandas.core.frame.DataFrame'>
print(df.shape) # (1599, 12)
print(df.head())
### ANCHOR_END: load_data

#print(df.head().to_markdown())

### ANCHOR: correlation
# Plot correlation matrix for the first four features
features = df.columns[:4]
pd.plotting.scatter_matrix(df[features], figsize=(10, 10))
### ANCHOR_END: correlation

#plt.tight_layout()
#plt.savefig('../../assets/figures/05-machine_learning/wine_correlation.svg')

### ANCHOR: preprocess_data
# Define data matrix and labels
X = df[["alcohol", "volatile acidity"]].to_numpy()
y = df["quality"].to_numpy()

# Add a column of ones to the data matrix
X = np.hstack([np.ones((X.shape[0], 1)), X])
### ANCHOR_END: preprocess_data

### ANCHOR: linear_regression
# Perform linear regression
theta = np.linalg.inv(X.T @ X) @ X.T @ y
print(theta)
### ANCHOR_END: linear_regression

### ANCHOR: plot_results
# Make a 3D plot
fig, ax = plt.subplots(subplot_kw={'projection': '3d'})

# Plot alcohol and volatile acidity against quality
ax.scatter(df['alcohol'], df['volatile acidity'], y, c='r', marker='o')

# Plot the linear regression
x = np.linspace(8, 15, 100) # alcohol
y = np.linspace(0, 1.6, 100) # volatile acidity
X, Y = np.meshgrid(x, y)
Z = theta[0] + theta[1] * X + theta[2] * Y
ax.plot_surface(X, Y, Z, alpha=0.5)

# Set labels
ax.set_xlabel('Alcohol')
ax.set_ylabel('Volatile acidity')
ax.set_zlabel('Quality')

fig.tight_layout()

plt.show()
### ANCHOR_END: plot_results

#fig.savefig('../../assets/figures/05-machine_learning/wine_regression.svg')