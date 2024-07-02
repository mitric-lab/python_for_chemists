#!/usr/bin/env python

### ANCHOR: eigenfaces_logistic_regression
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Define model class
class Logistic_Regression:
    def __init__(self, dim=2, tau=0.1, epochs=100):
        self.tau = tau
        self.epochs = epochs
        self.weights = np.random.rand(dim)
        self.bias = 0

    def sigmoid(self, x):
        return 1 / (1 + np.exp(-x)) 

    def fit(self, X, y):
        N = X.shape[0]
        for e in range(self.epochs):
            print(f"Epoch {e + 1}/{self.epochs}")
            
            for xi, yi in zip(X, y):
                self.weights += self.tau / N * (yi - self.sigmoid(self.net_output(xi))) * xi
                self.bias += self.tau / N * (yi - self.sigmoid(self.net_output(xi)))

    def net_output(self, x):
        return np.dot(x, self.weights) + self.bias

    def predict(self, x):
        return np.where(self.net_output >= 0, 1, 0)

# Load the data
path = './eigenfaces_pca.csv'
df = pd.read_csv(path, sep=';')

# Define data matrix and labels
X = df[["pca2"]].to_numpy()
y = df["label"].to_numpy()

# Normalize data
X = (X - X.mean()) / X.std()
y[y == -1] = 0

# Define hyperparameters
dim = X.shape[1]
tau = 0.1
epochs = 100

# Instantiate the model
LR = Logistic_Regression(dim=dim, tau=tau, epochs=epochs)

# Fit the model
LR.fit(X, y)

# Make predictions
x_grid = np.linspace(-3, 3, 100)
decision_line = np.array([LR.net_output([x]) for x in x_grid])

# Plot the data and decision line
fig, ax = plt.subplots(figsize=(6, 5))

ax.plot(X[y == 0], y[y == 0], 'o', color='b', label='Data')
ax.plot(X[y == 1], y[y == 1], 'o', color='r', label='Data')
ax.plot(x_grid, decision_line, 'g', label=r'$w x + b$')
ax.plot(x_grid, LR.sigmoid(decision_line), 'm', label=f'$\sigma(w x + b)$')

ax.set_ylim(-0.1, 1.1)
ax.set_xticklabels([])
ax.set_xlabel('PCA2')
plt.legend()

fig.tight_layout()

plt.show()
### ANCHOR_END: eigenfaces_logistic_regression

#fig.savefig('./eigenfaces_logistic_regression.png')