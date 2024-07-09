#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

class ReLU():
    def __call__(self, x):
        return np.maximum(x, 0)
    
    def gradient(self, x):
        return np.heaviside(x, 0) 
    
class Sigmoid():
    def __call__(self, x):
        return 1 / (1 + np.exp(-x))
    
    def gradient(self, x):
        return self(x) * (1 - self(x))

class Identity():
    def __call__(self, x):
        return x
    
    def gradient(self, x):
        return np.ones_like(x)

# define a mutli layer Network class
class MLP():

    # class constructor
    def __init__(self, sizes, tau=0.1, batch_size=5):        

        self.sizes = list(sizes.keys())
        self.activations = list(sizes.values())
        self.num_layers = len(self.sizes)
        
        self.biases = [np.random.randn(y) for y in self.sizes[1:]]
        self.weights = [np.random.randn(x, y) for x, y in zip(self.sizes[:-1], self.sizes[1:])]
        self.tau = tau
        self.batch_size = batch_size
            
            
    def feedforward(self, x):

        h_l = x
        for b_l, w_l, activation_l in zip(self.biases, self.weights, self.activations[1:]):
            h_l = activation_l(np.dot(w_l.T, h_l) + b_l)
        return h_l
        
    def train_step(self, X, y):

        N = X.shape[0]
        indices = np.arange(N)
        np.random.shuffle(indices)

        for i in range(0, N, self.batch_size):
            batch_indices = indices[i:i + self.batch_size]
            self.update_mini_batch(X[batch_indices], y[batch_indices])
                
    def update_mini_batch(self, X_batch, y_batch):

        print(X_batch.shape, y_batch.shape)
        gradient_bias = [np.zeros_like(b) for b in self.biases]
        gradient_weights = [np.zeros_like(w) for w in self.weights]

        # iterate over all training pairs in mini batch
        for xi, yi in zip(X_batch, y_batch):
            print(xi.shape, yi.shape)
            
            # compute gradient for this training pair
            gradient_bias_i, gradient_weights_i = self.backprop(xi, yi)

            # add gradient to already computed gradients for this minibatch
            gradient_bias = [db + db_i for db, db_i in zip(gradient_bias, gradient_bias_i)]
            gradient_weights = [dw + dw_i for dw, dw_i in zip(gradient_weights, gradient_weights_i)]

        # update free parameters using gradient descent step with average gradient of this minibatch
        self.weights = [w - (self.tau / self.batch_size) * dw for w, dw in zip(self.weights, gradient_weights)]
        self.biases = [b - (self.tau / self.batch_size) * db for b, db in zip(self.biases, gradient_bias)]

        
    def backprop(self, xi, yi):

        # initialize partial derivatives w.r.t. free parameters as zero
        gradient_bias_i = [np.zeros_like(b) for b in self.biases]
        gradient_weights_i = [np.zeros_like(w) for w in self.weights]
        
        ### Feedforward step ###
        
        # save input data as first output of input layer
        h_l = xi
        h_vectors = [h_l] # list to store all the outputs, layer by layer
        a_vectors = [] # list to store all the activations a, layer by layer
        
        # perform a feedforward step for second to last layer
        # with input x and save all intermediate results
        for b_l, w_l, activation_l in zip(self.biases, self.weights, self.activations[1:]):
            a_l = np.dot(w_l.T, h_l) + b_l
            print('a_l', a_l.shape)
            a_vectors.append(a_l)
            h_l = activation_l(a_l)
            print('h_l', h_l.shape)
            h_vectors.append(h_l)
        
        ### Backward pass ###
        #print(len(h_vectors), len(a_vectors))
        
        # compute derivative of the networks output layer
        delta_l = (h_vectors[-1] - y) * self.activations[-1].gradient(a_vectors[-1])
        print('delta_L', delta_l.shape)
        
        # compute partial derivatives of free parameters in output layer
        gradient_bias_i[-1] = delta_l
        gradient_weights_i[-1] = np.outer(h_vectors[-2], delta_l)
                
        # iterate from second last layer to input layer (using negative indices)
        for l in range(2, self.num_layers):

            delta_l = self.activations[-l].gradient(a_vectors[-l]) * np.dot(delta_l, self.weights[-(l-1)].T)
            print('delta_l', delta_l.shape)

            # compute derivatives of parameters in this layer using delta and z_vectors
            gradient_bias_i[-l] = delta_l
            gradient_weights_i[-l] = np.outer(h_vectors[-(l+1)], delta_l)
        
        return gradient_bias_i, gradient_weights_i

import mnist_loader as loader

tau = 0.1
batch_size = 10
            
sizes = {784: None,
        30: Sigmoid(),
        10: Sigmoid()}

f_hat = MLP(sizes, tau=tau, batch_size=batch_size)

# load MNIST handwritten digits database
train_data, valid_data, test_data = loader.load_MNIST('./mnist.pkl.gz')
N = len(train_data)
X = np.array([train_data[i][0] for i in range(N)])
y = np.array([train_data[i][1] for i in range(N)])
print(X.shape, y.shape)

epochs = 20
for e in range(epochs):
    f_hat.train_step(X, y)   
    if np.mod(e,1) == 0:
        accuracy = 0
        for i in range(N):
            x = train_data[i][0]
            y = train_data[i][1]
            if np.argmax(f_hat.feedforward(x)) == np.argmax(y):
                accuracy += 1
        accuracy /= N
        print("Epoch: {}, Train Accuracy: {}".format(e+1, accuracy))
        
N_valid = len(valid_data)

valid_accuracy = 0
for i in range(N_valid):
    x = valid_data[i][0]
    y = valid_data[i][1]
    if np.argmax(f_hat.feedforward(x)) == np.argmax(y):
        valid_accuracy += 1
valid_accuracy /= N_valid

print("Validation Accuracy: {}".format(valid_accuracy))