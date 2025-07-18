
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

# Plot sigmoid function
x_grid = np.linspace(-5, 5, 100)
sigmoid = Sigmoid()
fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(x_grid, sigmoid(x_grid), label='Sigmoid')
ax.set_xlabel('x')
ax.set_ylabel(r'$\sigma$(x)')

fig.tight_layout()

plt.show()
### ANCHOR_END: sigmoid

#fig.savefig('../../assets/figures/06-neural_networks/sigmoid.svg')

np.random.seed(42)

### ANCHOR: slp_init_predict
# Define model class
class SLP:
    def __init__(self, dim=3, hidden_size=5, activation='Sigmoid', epochs=100, tau=0.1):
        self.weights = np.random.randn(hidden_size, dim)  # dim includes bias column
        self.linear_weights = np.random.randn(hidden_size)
        if activation == "Sigmoid":
            self.activation = Sigmoid()
        else:
            raise NotImplementedError(f"Activation function not implemented.")
        self.epochs = epochs
        self.tau = tau
        self.losses = []
    
    def predict(self, x): 
        z = np.dot(self.weights, x)
        return np.dot(self.linear_weights, self.activation(z))
### ANCHOR_END: slp_init_predict

### ANCHOR: slp_fit
    def fit(self, X, y):        
        
        # Shuffle data
        N = X.shape[0]
        indices = np.random.permutation(N)
        X = X[indices]
        y = y[indices]

        # Training loop
        for e in range(self.epochs):
            print(f"Epoch {e + 1}/{self.epochs}")
            loss = 0

            # Iterate over all data points and update after each one (SGD)
            for xi, yi in zip(X, y):
        
                zi = np.dot(self.weights, xi)
                d_inner = self.linear_weights * self.activation.gradient(zi)
                residue = self.predict(xi) - yi
                loss += residue ** 2

                # Compute gradients for this single data point
                gradient_w = residue * np.outer(d_inner, xi)
                gradient_lw = residue * self.activation(zi)

                # Update parameters after each data point
                self.weights -= self.tau * gradient_w
                self.linear_weights -= self.tau * gradient_lw
            
            # Append loss after each epoch
            self.losses.append(loss / N)
### ANCHOR_END: slp_fit

### ANCHOR: load_data
# Load the data
path = 'aptamer_xor_data.csv'
df = pd.read_csv(path)
print(df.head())

# Define data matrix and labels
X = df.drop('labels', axis=1).values
X = np.hstack([X, np.ones((X.shape[0], 1))])
y = df['labels'].values
### ANCHOR_END: load_data

### ANCHOR: train_model
# Set hyperparameters
hidden_size = 5
tau = 0.01
dim = X.shape[1]
epochs = 100

# Instantiate the model
model = SLP(dim=dim, hidden_size=hidden_size, tau=tau, epochs=epochs)

# Train the model
model.fit(X, y)
### ANCHOR_END: train_model

### ANCHOR: plot_results
# Make plot
fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(12, 6))

# Plot the data points, color-coded by the labels
ax1.scatter(X[y == 1, 0], X[y == 1, 1], color='blue', label='Class 0')
ax1.scatter(X[y == -1, 0], X[y == -1, 1], color='red', label='Class 1')

# Plot the decision boundary
x1_grid = np.linspace(-8, 8, 100)
x2_grid = np.linspace(-8, 8, 100)
X1_grid, X2_grid = np.meshgrid(x1_grid, x2_grid, indexing='ij')
Y = np.zeros_like(X1_grid)
for i, x1 in enumerate(x1_grid):
    for j, x2 in enumerate(x2_grid):
        Y[i, j] = model.predict([x1, x2, 1])  # Add bias term

# Assign colors to the grid points based on the predicted class
ax1.contour(X1_grid, X2_grid, Y, levels=[0.0], colors='black', linestyles='dashed')
ax1.contourf(X1_grid, X2_grid, Y, levels=[-10.0, 0.0, 10.0], colors=['red', 'blue'], alpha=0.2)

# Plot the loss over epochs
ax2.plot(model.losses)
ax2.set_xlabel('Epoch')
ax2.set_ylabel('Loss')

fig.tight_layout()

plt.show()
### ANCHOR_END: plot_results

# fig.savefig('../../assets/figures/06-neural_networks/slp_xor_classification.svg')