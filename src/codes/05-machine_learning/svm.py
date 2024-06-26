#!/usr/bin/env python

### ANCHOR: svm_init
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

class SupportVectorMachine(object):
    """Support Vector Machine classifier. 
    
    Parameters
    ------------
    dim : int
      Dimension of the input data.
    tau : float
      Learning rate (between 0.0 and 1.0)
    lam : float
      Weight for maximising the margin.
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
    margins : list
      Width of the margin in each epoch.
    """

    def __init__(self, dim=2, tau=0.1, lam=1.0, epochs=100):
        self.tau = tau
        self.lam = lam
        self.epochs = epochs
        self.w = np.random.randn(dim)
        self.b = 0.0
        self.w_list = [self.w.copy()] # need to copy to avoid reference
        self.b_list = [self.b] # no need to copy, scalar
        self.errors = []
        self.losses = []
        self.margins = []
### ANCHOR_END: svm_init

    def fit(self, X, y):
        """Fit training data"""
        
        N = X.shape[0]
        for e in range(self.epochs):
            print(f'Epoch {e + 1}/{self.epochs}')
            
            error = 0.0
            loss = 0.0
            for xi, yi in zip(X, y):
                
                # yi if wrongly classified, zero if correct
                update = yi if yi * self.net_input(xi) < 0.0 else 0.0

                # count misclassifications
                error += 1 if update != 0.0 else 0

                # calculate loss
                loss += 0.5 * self.lam * np.dot(self.w, self.w) \
                    - update * (np.dot(self.w, xi) + self.b)
                
                # update parameters
                self.w += self.tau / N * (update * xi - self.lam * self.w)
                self.b += self.tau / N * update
            
            # save parameters and errors after epoch
            self.w_list.append(self.w.copy())
            self.b_list.append(self.b)
            self.errors.append(error)
            self.losses.append(loss)
            self.margins.append(2.0 / np.linalg.norm(self.w))
            
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

if __name__ == '__main__':

    ### ANCHOR: fit
    # Load the data
    path = './eigenfaces_pca.csv'
    df = pd.read_csv(path, sep=';')
    
    # Define data matrix and labels
    X = df[["pca1", "pca2"]].to_numpy()
    y = df["label"].to_numpy()
    
    # Define hyperparameters
    tau = 0.01
    lam = 10.0
    epochs = 50
    dim = X.shape[1]
    
    np.random.seed(42)
    
    # Initialize the perceptron
    classifier = SupportVectorMachine(dim=dim, tau=tau, lam=lam, epochs=epochs)
    
    # Fit the perceptron
    classifier.fit(X, y)
    ### ANCHOR_END: fit 
    
    ### ANCHOR: plot_results
    # Make plot
    fig, [[ax1, ax2], [ax3, ax4]] = plt.subplots(2, 2, figsize=(8.4, 7.2))
    
    # Plot the data points, color-coded by the labels
    ax1.scatter(X[:,0][y == 1], X[:,1][y == 1], color='blue', label='Person 5')
    ax1.scatter(X[:,0][y == -1], X[:,1][y == -1], color='red', label='Person 7')
    ax1.legend()
    
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
    ax2.plot(range(1, epochs + 1), classifier.errors, marker='o')
    
    # Set labels
    ax2.set_xlabel('Epoch')
    ax2.set_ylabel('Number of misclassifications')
    
    ax3.set_yscale('log')
    ax3.plot(range(1, epochs + 1), classifier.losses, marker='o')
    ax3.set_xlabel('Epoch')
    ax3.set_ylabel('Loss')
    
    ax4.set_yscale('log')
    ax4.plot(range(1, epochs + 1), classifier.margins, marker='o')
    ax4.set_xlabel('Epoch')
    ax4.set_ylabel('Margin width')
    
    fig.tight_layout()
    
    plt.show()
    ### ANCHOR_END: plot_results
    
    #fig.savefig('../../assets/figures/05-machine_learning/perceptron.svg')

