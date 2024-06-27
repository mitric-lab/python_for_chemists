#!/usr/bin/env python

### ANCHOR: pca_init
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

class PCA:
    def __init__(self, n_components=2):
        self.n_components = n_components
        self.components = None
        self.explained_variance = None
### ANCHOR_END: pca_init

### ANCHOR: pca_fit    
    def fit(self, X):
        # Center the data
        X_centered = X - np.mean(X, axis=0)
        
        # Compute the covariance matrix
        cov_matrix = np.cov(X_centered, rowvar=False)
        
        # Compute the eigenvalues and eigenvectors
        eigenvalues, eigenvectors = np.linalg.eigh(cov_matrix)
        
        # Sort the eigenvalues and corresponding eigenvectors in descending order
        sorted_indices = np.argsort(eigenvalues)[::-1]
        eigenvalues = eigenvalues[sorted_indices]
        eigenvectors = eigenvectors[:, sorted_indices]
        
        # Store the top n_components eigenvectors (principal components)
        self.components = eigenvectors[:, :self.n_components]
        
        # Calculate the explained variances
        self.explained_variance = eigenvalues[:self.n_components] / np.sum(eigenvalues)
### ANCHOR_END: pca_fit
    
### ANCHOR: pca_transform
    def transform(self, X):
        # Project the data onto the principal components
        X_centered = X - np.mean(X, axis=0)
        return np.dot(X_centered, self.components)
    
    def fit_transform(self, X):
        # Fit the model and return the transformed data
        self.fit(X)
        return self.transform(X)
### ANCHOR_END: pca_transform

### ANCHOR: pca_example
# Import Iris dataset
csv_url = 'https://archive.ics.uci.edu/ml/machine-learning-databases/iris/iris.data'
col_names = ['Sepal_Length', 'Sepal_Width', 'Petal_Length', 'Petal_Width', 'Class']
df =  pd.read_csv(csv_url, names=col_names)

# Show the first few rows
print(df.head())

# Convert class labels to integers
df['Class'] = df['Class'].astype('category').cat.codes

# Define data matrix and labels
X = df.drop('Class', axis=1).to_numpy()
y = df['Class'].to_numpy()

# Perform PCA
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X)

# Plot projected data and color-code by class
fig, ax = plt.subplots(figsize=(7, 5))

ax.scatter(X_pca[y == 0, 0], X_pca[y == 0, 1], color='blue', label='Iris-setosa')
ax.scatter(X_pca[y == 1, 0], X_pca[y == 1, 1], color='red', label='Iris-versicolor')
ax.scatter(X_pca[y == 2, 0], X_pca[y == 2, 1], color='green', label='Iris-virginica')

ax.set_xlabel('Principal Component 1')
ax.set_ylabel('Principal Component 2')
ax.legend()

fig.tight_layout()

plt.show()
### ANCHOR_END: pca_example

#fig.savefig('../../assets/figures/05-machine_learning/pca_iris.svg')