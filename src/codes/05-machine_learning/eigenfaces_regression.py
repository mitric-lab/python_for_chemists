#!/usr/bin/env python

### ANCHOR: eigenfaces_regression
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Load the data
path = './eigenfaces_pca.csv'
df = pd.read_csv(path, sep=';')

print(df.shape) # (1599, 12)
print(df.head())

# Define data matrix and labels
X = df[["pca1", "pca2"]].to_numpy()
y = df["label"].to_numpy()

# Add a column of ones to the data matrix
X = np.hstack([np.ones((X.shape[0], 1)), X])

# Perform linear regression
theta = np.linalg.inv(X.T @ X) @ X.T @ y
print(theta)

# Make a 3D plot
fig, ax = plt.subplots(subplot_kw={'projection': '3d'})

# Plot the data points, color-coded by the labels
ax.scatter(X[:,1][y == 1], X[:,2][y == 1], y[y == 1], color='blue', label='Person 5')
ax.scatter(X[:,1][y == -1], X[:,2][y == -1], y[y == -1], color='red', label='Person 7')

# Plot the linear regression
x = np.linspace(-5000, 5000, 1000)
y = np.linspace(-5000, 5000, 1000)
X, Y = np.meshgrid(x, y)
Z = theta[0] + theta[1] * X + theta[2] * Y
ax.plot_surface(X, Y, Z, alpha=0.5)

# Plot decision boundary, where Z = 0
ax.contour(X, Y, Z, levels=[0], colors='black', linestyles='dashed')

# Set labels
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_xlabel('PCA 5')
ax.set_ylabel('PCA 6')
ax.set_zlabel('Label')
ax.view_init(elev=30, azim=-45)

fig.tight_layout()

plt.show()
### ANCHOR_END: eigenfaces_regression

#fig.savefig('../../assets/figures/05-machine_learning/eigenfaces_regression.svg')