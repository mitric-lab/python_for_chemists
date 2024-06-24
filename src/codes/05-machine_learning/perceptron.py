#!/usr/bin/env python

### ANCHOR: perceptron_init
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

class Perceptron:
    """Perceptron classifier.
    
    Parameters
    ------------
    dim : int
      Dimension of the input data.
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
    def __init__(self, dim=2, tau=.1, epochs=100):
        self.tau = tau
        self.epochs = epochs
        self.w = np.random.randn(dim)
        self.b = 0.0
        self.w_list = [self.w.copy()] # need to copy to avoid reference
        self.b_list = [self.b] # no need to copy, scalar
        self.errors = []
### ANCHOR_END: perceptron_init

### ANCHOR: perceptron_fit        
    def fit(self, X, y):
        '''Fit training data'''
        
        N = X.shape[0]
        for e in range(epochs):
            print(f"Epoch {e+1}/{epochs}")
            
            errors = 0
            for xi, yi in zip(X, y):
                
                # yi if wrongly classified, zero if correct
                update = yi if yi * self.net_input(xi) < 0.0 else 0.0

                # count misclassifications
                errors += 1 if update != 0.0 else 0
                
                # update parameters
                self.w += self.tau / N * update * xi
                self.b += self.tau / N * update
            
            # save parameters and errors after epoch
            self.w_list.append(self.w.copy())
            self.b_list.append(self.b)
            self.errors.append(errors)
            
        return self
            
    def net_input(self, x):
        """Calculate net input"""
        return np.dot(x, self.w) + self.b
### ANCHOR_END: perceptron_fit
    
### ANCHOR: perceptron_predict
    def predict(self, x):
        """Return class label after unit step"""
        return np.where(self.net_input(x) >= 0.0, 1, -1)
### ANCHOR_END: perceptron_predict

### ANCHOR: fit
# Load the data
path = './eigenfaces_pca.csv'
df = pd.read_csv(path, sep=';')

# Define data matrix and labels
X = df[["pca1", "pca2"]].to_numpy()
y = df["label"].to_numpy()

# Define hyperparameters
tau = 0.001
epochs = 10
dim = X.shape[1]

# Initialize the perceptron
classifier = Perceptron(dim=dim, tau=tau, epochs=epochs)

# Fit the perceptron
classifier.fit(X, y)
### ANCHOR_END: fit 

### ANCHOR: plot_results
# Make plot
fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(12, 6))

# Plot the data points, color-coded by the labels
ax1.scatter(X[:,0][y == 1], X[:,1][y == 1], color='blue', label='Person 5')
ax1.scatter(X[:,0][y == -1], X[:,1][y == -1], color='red', label='Person 7')

# Plot the final decision boundary as a line
x = np.linspace(-5000, 5000, 1000)
y = np.linspace(-5000, 5000, 1000)
X, Y = np.meshgrid(x, y)
Z = classifier.w[0] * X + classifier.w[1] * Y + classifier.b
ax1.contour(X, Y, Z, levels=[0], colors='black', linestyles='dashed')

# Set labels
ax1.set_xticks([])
ax1.set_yticks([])
ax1.set_xlabel('PCA 5')
ax1.set_ylabel('PCA 6')

# Plot the error rate
ax2.plot(range(1, epochs+1), classifier.errors, marker='o')

# Set labels
ax2.set_xlabel('Epoch')
ax2.set_ylabel('Number of misclassifications')

fig.tight_layout()

plt.show()
### ANCHOR_END: plot_results

#fig.savefig('../../assets/figures/05-machine_learning/perceptron.svg')