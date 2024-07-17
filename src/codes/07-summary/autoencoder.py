#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

class Sigmoid():
    def __call__(self, x):
        return 1 / (1 + np.exp(-x))
    
    def gradient(self, x):
        return self(x) * (1 - self(x))

class ReLU():
    def __call__(self, x):
        return np.maximum(x, 0)
    
    def gradient(self, x):
        return np.heaviside(x, 0) 

class Identity():
    def __call__(self, x):
        return x
    
    def gradient(self, x):
        return np.ones_like(x)

class Softmax():
    def __call__(self, x):
        exps = np.exp(x - np.max(x))
        return exps / np.sum(exps)

class MLP():

    def __init__(self, sizes, tau=0.1, batch_size=5):
        # Initialize the network's parameters
        self.sizes = sizes['sizes']
        self.activations = sizes['activations']
        self.num_layers = len(self.sizes)
        self.weights = [np.random.randn(x, y) for x, y in zip(self.sizes[:-1], self.sizes[1:])]
        self.biases = [np.random.randn(y) for y in self.sizes[1:]]
        self.tau = tau
        self.batch_size = batch_size
            
    def feedforward(self, xi):
        # Compute the output of the network for input xi
        self.h_list = [xi]
        h_l = xi
        for b_l, w_l, activation_l in zip(self.biases, self.weights, self.activations[1:]):
            h_l = activation_l(np.dot(w_l.T, h_l) + b_l)
            self.h_list.append(h_l)
        return h_l

    def train_step(self, X, y):
        # Train the network using mini-batch gradient descent

        # Shuffle the training data
        N = X.shape[0]
        indices = np.arange(N)
        np.random.shuffle(indices)

        # Iterate over mini-batches
        for i in range(0, N, self.batch_size):
            batch_indices = indices[i:i+self.batch_size]
            self.update_mini_batch(X[batch_indices], y[batch_indices])
                
    def update_mini_batch(self, X_batch, y_batch):
        # Update the network's parameters using gradient descent

        # Initialize gradients for this minibatch
        gradient_weights = [np.zeros_like(w) for w in self.weights]
        gradient_bias = [np.zeros_like(b) for b in self.biases]

        # Iterate over all training pairs in this minibatch
        for xi, yi in zip(X_batch, y_batch):
            
            # Compute gradients for this training pair
            gradient_weights_i, gradient_bias_i = self.backprop(xi, yi)

            # Update the gradients for this training pair
            gradient_weights = [dw + dw_i for dw, dw_i in zip(gradient_weights, gradient_weights_i)]
            gradient_bias = [db + db_i for db, db_i in zip(gradient_bias, gradient_bias_i)]

        # Update the network's parameters using the gradients
        self.weights = [w - (self.tau / self.batch_size) * dw for w, dw in zip(self.weights, gradient_weights)]
        self.biases = [b - (self.tau / self.batch_size) * db for b, db in zip(self.biases, gradient_bias)]

    def backprop(self, xi, yi):
        # Compute the gradients of the network's parameters for input xi and target yi

        # Initialize gradients for this input
        gradient_weights_i = [np.zeros_like(w) for w in self.weights]
        gradient_bias_i = [np.zeros_like(b) for b in self.biases]
        
        ### Feedforward pass ###
        
        # First layer is just the input
        h_l = xi
        h_vectors = [h_l] # store hidden layer outputs
        a_vectors = [] # store activations
        
        # Perform forward pass through the network and store activations and outputs
        for w_l, b_l, activation_l in zip(self.weights, self.biases, self.activations[1:]):
            a_l = np.dot(w_l.T, h_l) + b_l
            a_vectors.append(a_l)
            h_l = activation_l(a_l)
            h_vectors.append(h_l)
        
        ### Backward pass ###

        # Compute delta for output layer
        delta_l = (h_vectors[-1] - yi) * self.activations[-1].gradient(a_vectors[-1])
        
        # Compute derivatives of parameters in output layer
        gradient_weights_i[-1] = np.outer(h_vectors[-2], delta_l)
        gradient_bias_i[-1] = delta_l
                
        # Iterate over hidden layers
        for l in range(2, self.num_layers):

            # Compute delta for this layer
            delta_l = self.activations[-l].gradient(a_vectors[-l]) * np.dot(delta_l, self.weights[-(l-1)].T)

            # Compute derivatives of parameters in this layer
            gradient_weights_i[-l] = np.outer(h_vectors[-(l+1)], delta_l)
            gradient_bias_i[-l] = delta_l
        
        return  gradient_weights_i, gradient_bias_i

    def train(self, X, epochs=10, validate=False):
        # Train or validate the network for a number of epochs

        for e in range(epochs):
            
            # Train the network
            if not validate:
                self.train_step(X, X)
            
            # Compute reconstruction loss
            N = X.shape[0]
            total_loss = 0
            for xi in X:
                total_loss += np.linalg.norm(xi - self.feedforward(xi))**2 / N
            
            print(f"Epoch: {e+1}, Loss: {total_loss:.3f}")

import mnist_loader as loader

# Load the MNIST dataset
train_data, valid_data, test_data = loader.load_MNIST('./mnist.pkl.gz')

tau = 3.0
batch_size = 100
epochs = 10
            
sizes = {'sizes': [784, 256, 128, 2, 128, 256, 784],
        'activations': [None, Sigmoid(), Sigmoid(), Identity(), Sigmoid(),Sigmoid(), Sigmoid()]}

autoencoder = MLP(sizes, tau=tau, batch_size=batch_size)

N = len(train_data)
X_train = np.array([train_data[i][0] for i in range(N)])

autoencoder.train(X_train, epochs=epochs, validate=False)

N_valid = len(valid_data)
X_valid = np.array([valid_data[i][0] for i in range(N_valid)])

autoencoder.train(X_valid, epochs=1, validate=True)

# Recover latent space representation
latent_vectors = []
reconstructed_images = []
for xi in X_train:
    x_hat = autoencoder.feedforward(xi)
    z = autoencoder.h_list[3]
    latent_vectors.append(z)
    reconstructed_images.append((x_hat, xi))

# Plot original and reconstructed images
fig, ax = plt.subplots(2, 5, figsize=(10, 4))

for i in range(5):
    ax[0, i].imshow(reconstructed_images[i][0].reshape(28, 28), cmap='gray')
    ax[0, i].axis('off')
    ax[1, i].imshow(reconstructed_images[i][1].reshape(28, 28), cmap='gray')
    ax[1, i].axis('off')

fig.tight_layout()

plt.show()

# Plot latent space
fig, ax = plt.subplots(1, 1, figsize=(6, 6))

latent_vectors = np.array(latent_vectors)
ax.scatter(latent_vectors[:, 0], latent_vectors[:, 1], c='b', marker='o')

fig.tight_layout()

plt.show()

