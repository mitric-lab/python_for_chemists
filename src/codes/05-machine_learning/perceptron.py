#!/usr/bin/env python

### ANCHOR: rosenblatt_perceptron
import numpy as np
import matplotlib.pyplot as plt

class Perceptron:
    """Perceptron classifier.
    
    Parameters
    ------------
    tau : float
      Learning rate (between 0.0 and 1.0)
    epochs : int
      Passes over the training dataset.
    
    Attributes
    -----------
    w : 1d-array
      Weights after fitting.
    b : Scalar
      Bias unit after fitting.
    w_list : list
      Weights in every epoch.
    b_list : list
      Bias units in every epoch.
    errors : list
      Number of misclassifications (updates) in each epoch.

    """
    def __init__(self, tau=.1, epochs=100):
        self.tau = tau
        self.epochs = epochs
        self.w = np.random.randn(X.shape[-1])
        self.b = 0.0
        self.w_list = [self.w.copy()]
        self.b_list = [self.b]
        self.errors = []
        
    def fit(self, X, y):
        '''Fit training data'''
        
        N = X.shape[0]
        for _ in range(epochs):
            
            errors = 0
            for xi, yi in zip(X, y):
                
                # yi if wrongly classified, zero if correct
                update = yi if yi * self.net_input(xi) < 0.0 else 0.0
                
                # update parameters
                self.w += self.tau / N * update * xi
                self.b += self.tau / N * update
                
                # count wrong classifications
                errors += 1 if yi * self.net_input(xi) < 0.0 else 0
            
            # save parameters and errors after epoch
            self.w_list.append(self.w.copy())
            self.b_list.append(self.b)
            self.errors.append(errors)
            
        return self
            
    def net_input(self, x):
        """Calculate net input"""
        return np.dot(x, self.w) + self.b
    
    def predict(self, x):
        """Return class label after unit step"""
        return np.where(self.net_input(x) >= 0.0, 1, -1)
### ANCHOR_END: rosenblatt_perceptron