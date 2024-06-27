#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pca import PCA

### ANCHOR: kmeans_init
class KMeans:
    def __init__(self, n_clusters=3, num_iter=300):
        self.n_clusters = n_clusters
        self.num_iter = num_iter
        self.centroids = None # array of shape (n_clusters, n_features)
        self.labels = None # array of shape (n_points)
### ANCHOR_END: kmeans_init
    
### ANCHOR: kmeans_fit
    def fit(self, X):
        # Randomly initialize centroids
        random_indices = np.random.choice(X.shape[0], self.n_clusters, replace=False)
        self.centroids = X[random_indices]
        
        for i in range(self.num_iter):
            # Assign labels based on closest centroid
            self.labels = self.assign_labels(X)
            # Calculate new centroids from the means of the points
            self.centroids = self.compute_centroids(X)
### ANCHOR_END: kmeans_fit

### ANCHOR: kmeans_assign_labels
    def assign_labels(self, X):
        # Calculate the distance between each point and each centroid
        distances = np.linalg.norm(X[:, None, :] - self.centroids, axis=2)
        # Assign the nearest centroid to each point
        return np.argmin(distances, axis=1)
### ANCHOR_END: kmeans_assign_labels
    
### ANCHOR: kmeans_compute_centroids
    def compute_centroids(self, X):
        # Calculate new centroids as the mean of all points assigned to each centroid
        return np.array([np.mean(X[self.labels == i], axis=0) for i in range (self.n_clusters)])
    
    def predict(self, X):
        # Assign labels to new data points based on the current centroids
        return self.assign_labels(X)
### ANCHOR_END: kmeans_compute_centroids

# import Iris dataset
csv_url = 'https://archive.ics.uci.edu/ml/machine-learning-databases/iris/iris.data'
col_names = ['Sepal_Length', 'Sepal_Width', 'Petal_Length', 'Petal_Width', 'Class']
df =  pd.read_csv(csv_url, names=col_names)

# Convert class labels to integers
df['Class'] = df['Class'].astype('category').cat.codes

# Define data matrix and labels
X = df.drop('Class', axis=1).to_numpy()
y = df['Class'].to_numpy()

# Perform PCA
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X)

### ANCHOR: kmeans_example
# Define hyperparameters
n_clusters = 3
num_iter = 100

# Perform K-means clustering
kmeans = KMeans(n_clusters=n_clusters, num_iter=num_iter)
kmeans.fit(X_pca)

# Extract the predicted labels
y_pred = kmeans.predict(X_pca)

# Plot the data points, color-coded by the predicted labels
fig, ax = plt.subplots(figsize=(7, 5))

ax.scatter(X_pca[y_pred == 0, 0], X_pca[y_pred == 0, 1], color='blue')
ax.scatter(X_pca[y_pred == 1, 0], X_pca[y_pred == 1, 1], color='red')
ax.scatter(X_pca[y_pred == 2, 0], X_pca[y_pred == 2, 1], color='green')

ax.set_xlabel('Principal Component 1')
ax.set_ylabel('Principal Component 2')

fig.tight_layout()

plt.show()
### ANCHOR_END: kmeans_example

#fig.savefig('../../assets/figures/05-machine_learning/k_means_iris.svg')
