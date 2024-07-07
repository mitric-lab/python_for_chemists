#!/usr/bin/env python

### ANCHOR: sigmoid
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Define activation function class
class Sigmoid:
    def __call__(self, x):
        return 1 / (1 + np.exp(-x))

    def gradient(self, x):
        return self(x) * (1 - self(x))
### ANCHOR_END: sigmoid

# Plot sigmoid function
x_grid = np.linspace(-5, 5, 100)
sigmoid = Sigmoid()
fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(x_grid, sigmoid(x_grid), label='Sigmoid')
ax.set_xlabel('x')
ax.set_ylabel(r'$\sigma$(x)')

fig.tight_layout()

plt.show()

#fig.savefig('../../assets/figures/06-neural_networks/sigmoid.svg')

### ANCHOR: slp_feedforward
# Define model class
class SLP:
    def __init__(self, dim=2, hidden_size=2, activation='Sigmoid', epochs=100, tau=0.1, batch_size=5):
        self.weights = np.random.randn(dim, hidden_size)
        self.bias = np.random.randn(hidden_size)
        self.linear_weights = np.random.randn(hidden_size)
        if activation == "Sigmoid":
            self.activation = Sigmoid()
        else:
            raise NotImplementedError(f"Activation function not implemented.")
        self.epochs = epochs
        self.tau = tau
        self.batch_size = batch_size
        self.losses = []
    
    def feedforward(self, x): 
        z = np.dot(self.weights.T, x) + self.bias
        return np.dot(self.linear_weights.T, self.activation(z))
### ANCHOR_END: slp_feedforward
    
### ANCHOR: slp_train
    def train(self, X, y):        
        N = X.shape[0]
        for e in range(self.epochs):
            print(f"Epoch {e + 1}/{self.epochs}")
            
            # Shuffle data
            indices = np.arange(N)
            np.random.shuffle(indices)

            # Iterate over batches
            loss = 0
            for i in range(0, N, self.batch_size):

                # Define batch
                batch_indices = indices[i:i + self.batch_size]
                X_batch = X[batch_indices]
                y_batch = y[batch_indices]

                # Initialize gradients
                gradient_w = np.zeros_like(self.weights)
                gradient_b = np.zeros_like(self.bias)
                gradient_lw = np.zeros_like(self.linear_weights)

                # Accumulate gradients over the batch
                for xi, yi in zip(X_batch, y_batch):
            
                    zi = np.dot(self.weights.T, xi) + self.bias
                    d_inner = self.linear_weights * self.activation.gradient(zi)
                    residue = self.feedforward(xi) - yi
                    loss += residue ** 2

                    # Compute gradients
                    gradient_w += residue * np.outer(d_inner, xi).T
                    gradient_b += residue * d_inner
                    gradient_lw += residue * self.activation(zi)

                # Update parameters after each batch
                self.weights -= self.tau / self.batch_size * gradient_w
                self.bias -= self.tau / self.batch_size * gradient_b
                self.linear_weights -= self.tau / self.batch_size * gradient_lw
            
            self.losses.append(loss / N)
### ANCHOR_END: slp_train

### ANCHOR: slp_example
# Load the data
path = './circles.csv'
df = pd.read_csv(path, sep=';')

# Define data matrix and labels
X = df[['x_1', 'x_2']].to_numpy()
y = df['label'].to_numpy()

# Set hyperparameters
hidden_size = 50
tau = 0.01
dim = X.shape[1]
epochs = 100
batch_size = 12

# Instantiate the model
f_hat = SLP(dim=dim, hidden_size=hidden_size, tau=tau, epochs=epochs, batch_size=batch_size)

# Train the model
f_hat.train(X, y)
### ANCHOR_END: slp_example

### ANCHOR: slp_plot
# Make plot
fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(12, 6))

# Plot the data points, color-coded by the labels
ax1.scatter(X[y == 1, 0], X[y == 1, 1], color='blue', label='Class 0')
ax1.scatter(X[y == -1, 0], X[y == -1, 1], color='red', label='Class 1')

# Plot the decision boundary
x_grid = np.linspace(-8, 8, 100)
y_grid = np.linspace(-8, 8, 100)
X_grid, Y_grid = np.meshgrid(x_grid, y_grid)
Z = np.zeros_like(X_grid)
for i, x in enumerate(x_grid):
    for j, y in enumerate(y_grid):
        Z[i, j] = f_hat.feedforward([x, y])

ax1.contour(X_grid, Y_grid, Z, levels=[0.0], colors='black', linestyles='dashed')
ax1.contourf(X_grid, Y_grid, Z, levels=[-10.0, 0.0, 10.0], colors=['red', 'blue'], alpha=0.2)

# Plot the loss over epochs
ax2.plot(f_hat.losses)
ax2.set_xlabel('Epoch')
ax2.set_ylabel('Loss')

fig.tight_layout()

plt.show()
### ANCHOR_END: slp_plot

#fig.savefig('../../assets/figures/06-neural_networks/slp_circles.svg')