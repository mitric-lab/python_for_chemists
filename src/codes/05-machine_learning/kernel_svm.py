#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from svm import SupportVectorMachine


def phi(x_1, x_2):
    return np.array([
        x_1, 
        x_2, 
        x_1**2 + x_2**2,
    ])

path = './circles.csv'
df = pd.read_csv(path, sep=';')
# Define data matrix and labels
X = df[['x_1', 'x_2']].to_numpy()
y = df['label'].to_numpy()

X = phi(X[:,0], X[:,1]).T


# Define hyperparameters
tau = 0.005
lam = 10.0
epochs = 200
dim = X.shape[1]

#np.random.seed(42)

# Initialize the perceptron
classifier = SupportVectorMachine(dim=dim, tau=tau, lam=lam, epochs=epochs)

# Fit the perceptron
classifier.fit(X, y)
### ANCHOR_END: fit 

### ANCHOR: plot_results
# Make plot
fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(8, 4))

# Plot the data points, color-coded by the labels
ax1.scatter(X[:,0][y == 1], X[:,1][y == 1], color='red', label='+1')
ax1.scatter(X[:,0][y == -1], X[:,1][y == -1], color='blue', label='-1')
ax1.legend()

# Plot the final decision boundary
x = np.linspace(-8, 8, 1000)
y = np.linspace(-8, 8, 1000)
X, Y = np.meshgrid(x, y)
X, Y, Z = phi(X, Y)
labels = classifier.w[0] * X + classifier.w[1] * Y \
    + classifier.w[2] * Z + classifier.b
ax1.contour(X, Y, labels, levels=[0], colors='black', linestyles='dashed')

# Set labels
ax1.set_xticks([])
ax1.set_yticks([])
ax1.set_xlabel('x_1')
ax1.set_ylabel('x_2')

ax2.set_yscale('log')
ax2.plot(range(1, epochs + 1), classifier.losses, marker='o')
ax2.set_xlabel('Epoch')
ax2.set_ylabel('Loss')

fig.tight_layout()

plt.show()
### ANCHOR_END: plot_results

