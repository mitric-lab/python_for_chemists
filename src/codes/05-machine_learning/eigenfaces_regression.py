#!/usr/bin/env python

### ANCHOR: eigenfaces_regression
import numpy as np
import matplotlib.pyplot as plt

# Load the data
pca1, pca2, labels = np.loadtxt('./eigenfaces_pca.csv', unpack=True, delimiter=',')

# Define data matrix and labels
X = np.vstack([np.ones_like(pca1), pca1, pca2]).T
y = labels

# Perform linear regression
theta = np.linalg.inv(X.T @ X) @ X.T @ y
print(theta)

# Make a 3D plot
fig, ax = plt.subplots(subplot_kw={'projection': '3d'})

# Plot the data points, color-coded by the labels
ax.scatter(pca1[labels == 1], pca2[labels == 1], labels[labels == 1], color='blue', label='Person 5')
ax.scatter(pca1[labels == -1], pca2[labels == -1], labels[labels == -1], color='red', label='Person 7')

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

fig.tight_layout()

plt.show()
### ANCHOR_END: eigenfaces_regression

#fig.savefig('../../assets/figures/05-machine_learning/eigenfaces_regression.svg')