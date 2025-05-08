#!/usr/bin/env python

### ANCHOR: logistic_regression
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
### ANCHOR_END: logistic_regression

### ANCHOR: eigenfaces_logistic_regression
# Load the data
path = './eigenfaces_pca.csv'
df = pd.read_csv(path, sep=';')

# Define data matrix and labels
X = df[["pca2"]].to_numpy() # take only one feature for simplicity
y = df["label"].to_numpy()

# Normalize data
X = (X - X.mean()) / X.std()
y[y == -1] = 0 # change -1 to 0

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

# Perform simple linear regression for comparison
beta_1, beta_0 = np.polyfit(X.flatten(), y, 1)

# Plot the data and decision line
fig, ax = plt.subplots(figsize=(6, 5))

ax.plot(X[y == 0], y[y == 0], 'o', color='b', label='Data')
ax.plot(X[y == 1], y[y == 1], 'o', color='r', label='Data')
ax.plot(x_grid, beta_1 * x_grid + beta_0, 'grey', linestyle='--', label='Linear regression')
ax.plot(x_grid, decision_line, 'g', label=r'$w x + b$')
ax.plot(x_grid, LR.sigmoid(decision_line), 'm', label=f'$\sigma(w x + b)$')

ax.set_ylim(-0.1, 1.1)
ax.set_xticklabels([])
ax.set_xlabel('PCA2')
plt.legend()

fig.tight_layout()

plt.show()
### ANCHOR_END: eigenfaces_logistic_regression

#fig.savefig('../../assets/figures/06-neural_networks/eigenfaces_logistic_regression.svg')

# Define data matrix and labels
X = df[["pca1", "pca2"]].to_numpy()
y = df["label"].to_numpy()

# Normalize data
X = (X - X.mean()) / X.std()
y[y == -1] = 0 # change -1 to 0

# Define hyperparameters
dim = X.shape[1]
tau = 0.1
epochs = 1000

# Instantiate the model
LR = Logistic_Regression(dim=dim, tau=tau, epochs=epochs)

# Fit the model
LR.fit(X, y)

# Plot class probabilities in 2D
x1_grid = np.linspace(-3, 3, 100)
x2_grid = np.linspace(-3, 3, 100)
X1, X2 = np.meshgrid(x1_grid, x2_grid)
Z = np.zeros_like(X1)
for i in range(X1.shape[0]):
    for j in range(X1.shape[1]):
        Z[i, j] = LR.sigmoid(LR.net_output([X1[i, j], X2[i, j]]))


fig, ax = plt.subplots(figsize=(6, 5), subplot_kw={'projection': '3d'})

ax.plot_surface(X1, X2, Z, cmap='coolwarm', alpha=0.5)
ax.scatter(X[y == 0, 0], X[y == 0, 1], y[y == 0], color='b', label='Data')
ax.scatter(X[y == 1, 0], X[y == 1, 1], y[y == 1], color='r', label='Data')
ax.contour(X1, X2, Z, levels=[0.5], colors='black', linestyles='dashed')
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_zticklabels([])
ax.set_xlabel('PCA1')
ax.set_ylabel('PCA2')
ax.set_zlabel('Label')
ax.view_init(elev=30, azim=-45)

fig.tight_layout()

plt.show()

#fig.savefig('../../assets/figures/06-neural_networks/eigenfaces_logistic_regression_2d.svg')