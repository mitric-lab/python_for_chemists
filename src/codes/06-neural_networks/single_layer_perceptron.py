#!/usr/bin/env python

### ANCHOR: single_layer_perceptron
import numpy as np
import matplotlib.pyplot as plt

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

#fig.savefig('../../assets/figures/06-neural_networks/sigmoid.svg')

# Define model class
class SLP:
    def __init__(self, dim=2, hidden_size=2, activation='Sigmoid', epochs=100, tau=0.1, batch_size=5):
        self.weights = np.random.randn(dim, hidden_size) # (d, n)
        self.bias = np.random.randn(hidden_size) # (n,)
        self.linear_weights = np.random.randn(hidden_size) # (n,)
        if activation == "Sigmoid":
            self.activation = Sigmoid()    
        self.epochs = epochs
        self.tau = tau
        self.batch_size = batch_size
        self.losses = []
    
    def feedforward(self, x): 
        # x: (d,)
        z = np.dot(self.weights.T, x) + self.bias # (n,)
        return np.dot(self.linear_weights.T, self.activation(z)) # (1,)
    
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
            
                    zi = np.dot(self.weights.T, xi) + self.bias # (n,)
                    d_inner = self.linear_weights * self.activation.gradient(zi) # (n,)
                    residue = self.feedforward(xi) - yi # (1,)
                    loss += residue ** 2

                    gradient_w += residue * np.outer(d_inner, xi).T
                    gradient_b += residue * d_inner
                    gradient_lw += residue * self.activation(zi)

                # Update parameters after each batch
                self.linear_weights -= self.tau / N * gradient_lw
                self.weights -= self.tau / N * gradient_w
                self.bias -= self.tau / N * gradient_b
                
            self.losses.append(loss / N)


hidden_size = 50
tau = 0.1
dim = 1
N = 20
epochs = 10000
batch_size = 5

# initialize network
f_hat = SLP(dim=dim, hidden_size=hidden_size, tau=tau, epochs=epochs, batch_size=batch_size)
f = lambda x: np.sin(4*x) + np.cos(6*x)

# plot data points
X = np.linspace(-1.5, 1.5, N)
X = X.reshape(-1, 1) # (N, n)
y = f(X).flatten() # (N,)

# train network
f_hat.train(X, y)

# plot loss
plt.plot(f_hat.losses)
plt.show()

# plot data points and network output
plt.scatter(X, y, color='blue')
x_grid = np.linspace(-1.5, 1.5, 100)
predictions = np.array([f_hat.feedforward([x]) for x in x_grid])
plt.plot(x_grid, predictions, color='red')
plt.plot(x_grid, f(x_grid), color='green')
plt.show()