#!/usr/bin/env python

## ANCHOR: load_data_from_csv
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

path_to_csv = "aptamer_classification_data.csv"
df = pd.read_csv(path_to_csv)
print(df.head())
### ANCHOR_END: load_data_from_csv

### ANCHOR: process_data
target_column = "lambda_abs_class"

X = df[["PC1", "PC2"]].values # Select the PC1 and PC2 columns as features
X = np.hstack([X, np.ones((X.shape[0], 1))]) # Add a column of ones to the data matrix
y = df[target_column].values # Target column

print(X.shape)
print(y.shape)
### ANCHOR_END: process_data

### ANCHOR: plot_data
fig, ax = plt.subplots(figsize=(7, 6))
ax.scatter(X[:, 0], X[:, 1], c=y, cmap='coolwarm', alpha=0.7)
ax.set_xlabel('PC1')
ax.set_ylabel('PC2')
plt.show()
### ANCHOR_END: plot_data

fig.savefig('../../assets/figures/05-machine_learning/classification_data.svg')

### ANCHOR: classification_class
class RosenblattPerceptron:
    def __init__(self, learning_rate=0.1, n_iterations=100):
        self.learning_rate = learning_rate
        self.n_iterations = n_iterations
        self.weights = None
    
    def fit(self, X, y):
        self.weights = np.zeros(X.shape[1])

        for _ in range(self.n_iterations):
            for x_i, y_i in zip(X, y):
                loss = - y_i * np.dot(self.weights, x_i)
                if loss >= 0:
                    self.weights += self.learning_rate * y_i * x_i

    def predict(self, X):
        return np.sign(np.dot(X, self.weights))
### ANCHOR_END: classification_class

### ANCHOR: fit_model
model = RosenblattPerceptron()
model.fit(X, y)
y_pred = model.predict(X)
### ANCHOR_END: fit_model

### ANCHOR: calculate_accuracy
accuracy = np.mean(y_pred == y)
print(f"Accuracy: {accuracy}")
### ANCHOR_END: calculate_accuracy

### ANCHOR: plot_decision_boundary
fig, ax = plt.subplots(figsize=(7, 6))
ax.scatter(X[:, 0], X[:, 1], c=y, cmap='coolwarm', alpha=0.7)

# Plot decision boundary (where w·x = 0)
# For weights [w0, w1, w2] and features [PC1, PC2, 1]:
# w0*PC1 + w1*PC2 + w2 = 0
# Solving for PC2: PC2 = -(w0*PC1 + w2) / w1

if model.weights[1] != 0:  # Avoid division by zero
    pc1_min, pc1_max = X[:, 0].min(), X[:, 0].max()
    pc1_range = np.linspace(pc1_min, pc1_max, 100)
    pc2_boundary = -(model.weights[0] * pc1_range + model.weights[2]) / model.weights[1]
    ax.plot(pc1_range, pc2_boundary, 'k-', linewidth=2, label='Decision Boundary')
    ax.legend()

ax.set_xlabel('PC1')
ax.set_ylabel('PC2')
plt.show()
### ANCHOR_END: plot_decision_boundary

fig.savefig('../../assets/figures/05-machine_learning/classification_decision_boundary.svg')
