#!/usr/bin/env python

SEED = 0

### ANCHOR: kmeans_init
class KMeans:
    def __init__(self, n_clusters=3, num_iter=50):
        self.n_clusters = n_clusters
        self.num_iter = num_iter
        self.centroids = None # array of shape (n_clusters, n_features)
        self.labels = None # array of shape (n_points)
### ANCHOR_END: kmeans_init
    
### ANCHOR: kmeans_fit
    def fit(self, X):
        # Randomly initialize centroids
        # Randomly select indices of data points to serve as the initial centroids.
        # np.random.choice selects 'self.n_clusters' unique indices from the range of available data points (X.shape[0]).
        # 'replace=False' ensures that the same data point is not selected more than once.
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

### ANCHOR: load_data_from_csv
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

path_to_csv = "aptamer_classification_data.csv"
df = pd.read_csv(path_to_csv)
X = df[["PC1", "PC2"]].values
print(df.head())
### ANCHOR_END: load_data_from_csv

np.random.seed(SEED)

### ANCHOR: kmeans_fit_and_predict
kmeans = KMeans(n_clusters=4)
kmeans.fit(X)
cluster_labels = kmeans.predict(X)
### ANCHOR_END: kmeans_fit_and_predict

### ANCHOR: plot_kmeans_result
fig, ax = plt.subplots(figsize=(7, 6))

# Plot the data points colored by the cluster labels
ax.scatter(X[:, 0], X[:, 1], c=cluster_labels)

# Plot the centroids as red crosses
ax.scatter(kmeans.centroids[:, 0], kmeans.centroids[:, 1], c='red', s=100, marker='x')

# Set the labels
ax.set_xlabel("PC1")
ax.set_ylabel("PC2")

plt.show()
### ANCHOR_END: plot_kmeans_result

# fig.savefig('../../assets/figures/05-machine_learning/k_means_aptamers.svg')

### ANCHOR: kmeans_with_history
class KMeansWithHistory:
    def __init__(self, n_clusters=3, num_iter=50):
        self.n_clusters = n_clusters
        self.num_iter = num_iter
        self.centroids = None
        self.labels = None
        # Store history of centroids and labels during optimization
        self.centroids_history = []
        self.labels_history = []
        self.step_descriptions = []  # Describe what each step represents
    
    def fit(self, X):
        # Randomly initialize centroids
        random_indices = np.random.choice(X.shape[0], self.n_clusters, replace=False)
        self.centroids = X[random_indices]
        
        # Store initial state
        self.centroids_history.append(self.centroids.copy())
        self.labels = self.assign_labels(X)
        self.labels_history.append(self.labels.copy())
        self.step_descriptions.append("Initial state")
        
        for i in range(self.num_iter):
            # Step 1: Assign labels based on current centroids
            self.labels = self.assign_labels(X)
            self.labels_history.append(self.labels.copy())
            self.centroids_history.append(self.centroids.copy())  # Keep same centroids
            self.step_descriptions.append(f"Iteration {i+1}: After label assignment")
            
            # Step 2: Calculate new centroids from the means of the points
            self.centroids = self.compute_centroids(X)
            self.centroids_history.append(self.centroids.copy())
            self.labels_history.append(self.labels.copy())  # Keep same labels
            self.step_descriptions.append(f"Iteration {i+1}: After centroid update")
    
    def assign_labels(self, X):
        # Calculate the distance between each point and each centroid
        distances = np.linalg.norm(X[:, None, :] - self.centroids, axis=2)
        # Assign the nearest centroid to each point
        return np.argmin(distances, axis=1)
    
    def compute_centroids(self, X):
        # Calculate new centroids as the mean of all points assigned to each centroid
        return np.array([np.mean(X[self.labels == i], axis=0) for i in range(self.n_clusters)])
    
    def predict(self, X):
        # Assign labels to new data points based on the current centroids
        return self.assign_labels(X)
### ANCHOR_END: kmeans_with_history

import os

# Create directory for iteration plots
os.makedirs('../../assets/figures/05-machine_learning/kmeans_iterations', exist_ok=True)

np.random.seed(SEED)

# Fit k-means with history tracking
kmeans_hist = KMeansWithHistory(n_clusters=4, num_iter=10)  # Reduced iterations for clarity
kmeans_hist.fit(X)

# Create plots for each step (both label assignment and centroid update)
for i in range(len(kmeans_hist.centroids_history)):
    fig, ax = plt.subplots(figsize=(7, 6))
    
    # Get current labels and centroids
    current_labels = kmeans_hist.labels_history[i]
    current_centroids = kmeans_hist.centroids_history[i]
    
    # Plot the data points colored by current cluster labels
    ax.scatter(X[:, 0], X[:, 1], c=current_labels)
    
    # Plot the current centroids as red crosses
    ax.scatter(current_centroids[:, 0], current_centroids[:, 1], c='red', s=100, marker='x')
    
    # Set the labels and simple title with iteration number
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    
    # Calculate iteration number from step index
    if i == 0:
        ax.set_title("K-Means Clustering - Initial State")
    else:
        iteration_num = (i + 1) // 2
        ax.set_title(f"K-Means Clustering - Iteration {iteration_num}")
    
    # Save the plot
    # fig.savefig(f'../../assets/figures/05-machine_learning/kmeans_iterations/step_{i:02d}.svg', dpi=300, bbox_inches='tight')
    plt.close(fig)  # Close to save memory

print(f"Saved {len(kmeans_hist.centroids_history)} step plots")

# Also save the final result
fig, ax = plt.subplots(figsize=(7, 6))
final_labels = kmeans_hist.labels_history[-1]
final_centroids = kmeans_hist.centroids_history[-1]

ax.scatter(X[:, 0], X[:, 1], c=final_labels)
ax.scatter(final_centroids[:, 0], final_centroids[:, 1], c='red', s=100, marker='x')

ax.set_xlabel("PC1")
ax.set_ylabel("PC2")
ax.set_title("K-Means Clustering - Final Result")

# fig.savefig('../../assets/figures/05-machine_learning/k_means_final_with_history.svg', dpi=300, bbox_inches='tight')
plt.close(fig)

print("Final result saved as k_means_final_with_history.svg")
### ANCHOR_END: kmeans_optimization_visualization
