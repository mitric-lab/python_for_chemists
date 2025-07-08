#!/usr/bin/env python

SEED = 1234

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

# fig.savefig('../../assets/figures/05-machine_learning/classification_data.svg')

### ANCHOR: classification_class_init
class RosenblattPerceptron:
    def __init__(self, learning_rate=0.01, n_iterations=50):
        self.learning_rate = learning_rate
        self.n_iterations = n_iterations
        self.weights = None
### ANCHOR_END: classification_class_init
    
### ANCHOR: classification_class_fit
    def fit(self, X, y):
        self.weights = np.random.randn(X.shape[1])

        for _ in range(self.n_iterations):
            for x_i, y_i in zip(X, y):
                loss = - y_i * np.dot(self.weights, x_i)
                if loss > 0:
                    self.weights += self.learning_rate * y_i * x_i
### ANCHOR_END: classification_class_fit

### ANCHOR: classification_class_predict
    def predict(self, X):
        return np.sign(np.dot(X, self.weights))
### ANCHOR_END: classification_class_predict

np.random.seed(SEED)
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

# Define the decision boundary as a function of the first feature
x1_range = np.linspace(X[:, 0].min(), X[:, 0].max(), 100)
x2_boundary = -(model.weights[0] * x1_range + model.weights[2]) / model.weights[1]

# Plot the decision boundary
ax.plot(x1_range, x2_boundary, 'k--', linewidth=2, label='Decision Boundary')
ax.legend(loc='upper right')

ax.set_xlabel('PC1')
ax.set_ylabel('PC2')
ax.set_xlim(X[:, 0].min()-0.1, X[:, 0].max()+0.1)
ax.set_ylim(X[:, 1].min()-0.1, X[:, 1].max()+0.1)
plt.show()
### ANCHOR_END: plot_decision_boundary

# fig.savefig('../../assets/figures/05-machine_learning/classification_decision_boundary.svg')

class RosenblattPerceptronWithHistory:
    def __init__(self, learning_rate=0.01, n_iterations=50, save_interval=10):
        self.learning_rate = learning_rate
        self.n_iterations = n_iterations
        self.save_interval = save_interval
        self.weights = None
        self.weight_history = []
    
    def fit(self, X, y):
        self.weights = np.random.randn(X.shape[1])
        self.weight_history = [self.weights.copy()]
        
        for iteration in range(self.n_iterations):
            for x_i, y_i in zip(X, y):
                loss = - y_i * np.dot(self.weights, x_i)
                if loss > 0:
                    self.weights += self.learning_rate * y_i * x_i
            
            self.weight_history.append(self.weights.copy())
            
            if iteration % self.save_interval == 0 or iteration == self.n_iterations - 1:
                self.plot_and_save_boundary(X, y, iteration)
    
    def predict(self, X):
        return np.sign(np.dot(X, self.weights))
    
    def plot_and_save_boundary(self, X, y, iteration):
        fig, ax = plt.subplots(figsize=(7, 6))
        ax.scatter(X[:, 0], X[:, 1], c=y, cmap='coolwarm', alpha=0.7)
        
        x1_range = np.linspace(X[:, 0].min(), X[:, 0].max(), 100)
        
        if self.weights[1] != 0:
            x2_boundary = -(self.weights[0] * x1_range + self.weights[2]) / self.weights[1]
            ax.plot(x1_range, x2_boundary, 'k--', linewidth=2, label=f'Decision Boundary')
        
        ax.set_title(f'Iteration {iteration+1}')
        ax.legend(loc='upper right')
        ax.set_xlabel('PC1')
        ax.set_ylabel('PC2')
        ax.set_xlim(X[:, 0].min()-0.1, X[:, 0].max()+0.1)
        ax.set_ylim(X[:, 1].min()-0.1, X[:, 1].max()+0.1)
        
        plt.savefig(f'../../assets/figures/05-machine_learning/classification_boundary/iter_{iteration:03d}.png', dpi=150, bbox_inches='tight')
        plt.close()

np.random.seed(SEED)
model_with_history = RosenblattPerceptronWithHistory(learning_rate=0.01, n_iterations=50, save_interval=1)
model_with_history.fit(X, y)

