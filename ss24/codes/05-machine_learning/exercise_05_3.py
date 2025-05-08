#!/usr/bin/env python

### ANCHOR: svm_init
import numpy as np


class KMeans:
    
    def __init__(self, n_clusters=3, num_iter=300):
        self.n_clusters = n_clusters
        self.num_iter = num_iter
        self.centroids = None # array of shape (n_clusters, n_features)
        self.labels = None # array of shape (n_points)
        self.energies = []

    def fit(self, X):
        # Randomly initialize centroids
        random_indices = np.random.choice(X.shape[0], self.n_clusters, replace=False)
        self.centroids = X[random_indices]
        self.labels = self.assign_labels(X, self.centroids)

        for i in range(self.num_iter):
            print(f'Iteration {i + 1}/{self.num_iter}: {self.n_clusters} clusters')
            rand = np.random.rand()
            centroids = np.copy(self.centroids)
            current_energy = self.objective(X, centroids, self.labels)
            self.energies.append(current_energy)
            if rand < 0.5 and self.n_clusters > 1:
                # combine two clusters
                cluster1, cluster2 = np.random.choice(
                    self.n_clusters, 2, replace=False,
                )
                n1 = len(self.labels[self.labels == cluster1])
                n2 = len(self.labels[self.labels == cluster2])
                centroids[cluster1] = (1.0 / (n1 + n2)) \
                    * (n1 * centroids[cluster1] + n2 * centroids[cluster2])
                centroids = np.delete(centroids, cluster2, axis=0)
                labels = self.assign_labels(X, centroids)
                new_energy = self.objective(X, centroids, labels)
                if new_energy < current_energy:
                    self.centroids = centroids
                    self.labels = labels
                    self.n_clusters -= 1
            else:
                # split a cluster
                cluster = np.random.choice(self.n_clusters)
                n = len(self.labels[self.labels == cluster])
                if n < 2:
                    continue
                centroids[cluster] = np.mean(X[self.labels == cluster][0:n//2], axis=0)
                centroids = np.concatenate(
                    [
                        centroids, 
                        np.mean(X[self.labels == cluster][n//2:], axis=0)[None, :],
                    ],
                    axis=0,
                )
                labels = self.assign_labels(X, centroids)
                new_energy = self.objective(X, centroids, labels)
                if new_energy < current_energy:
                    self.centroids = centroids
                    self.labels = labels
                    self.n_clusters += 1

    def objective(self, X, centroids, labels):
        distances = np.linalg.norm(X[:, None, :] - centroids, axis=2)
        energy = 0.0 + 10.0 * len(centroids)**2
        for k in range(0, len(centroids)):
            energy += np.sum(distances[labels == k, k]**2)
        return energy

    def assign_labels(self, X, centroids):
        # Calculate the distance between each point and each centroid
        distances = np.linalg.norm(X[:, None, :] - centroids, axis=2)
        # Assign the nearest centroid to each point
        return np.argmin(distances, axis=1)

    def compute_centroids(self, X):
        # Calculate new centroids as the mean of all points assigned to each centroid
        return np.array([
            np.mean(X[self.labels == i], axis=0) 
            for i in range (self.n_clusters)
        ])
    
    def predict(self, X):
        # Assign labels to new data points based on the current centroids
        return self.assign_labels(X, self.centroids)



if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from matplotlib import colormaps
    import pandas as pd
    from pca import PCA
    
    # np.random.seed(0)
    cm = colormaps['jet']

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
    X_pca = np.unique(X_pca, axis=0)
    print(X_pca.shape)
    
    n_clusters = 3
    num_iter = 1000

    # Perform K-means clustering
    kmeans = KMeans(n_clusters=n_clusters, num_iter=num_iter)
    kmeans.fit(X_pca)
    
    # Extract the predicted labels
    y_pred = kmeans.predict(X_pca)
    
    colors = cm(np.linspace(0, 1, kmeans.n_clusters))

    # Plot the data points, color-coded by the predicted labels
    fig, ax = plt.subplots(figsize=(7, 5))
    
    for k, c in zip(range(kmeans.n_clusters), colors):
        ax.scatter(
            X_pca[y_pred == k, 0], X_pca[y_pred == k, 1], 
            color=c, label=f'Cluster {k}'
        )
    
    ax.set_xlabel('Principal Component 1')
    ax.set_ylabel('Principal Component 2')
    ax.legend()
    
    fig.tight_layout()
    

    fig2, ax2 = plt.subplots(figsize=(7, 5))
    ax2.plot(kmeans.energies)

    plt.show()
    
