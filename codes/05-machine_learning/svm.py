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

# fig.savefig('../../assets/figures/05-machine_learning/classification_data.svg')

### ANCHOR: svm_init
class SupportVectorMachine:
    def __init__(self, learning_rate=0.01, n_iterations=50, lam=10.0):
        self.learning_rate = learning_rate
        self.n_iterations = n_iterations
        self.lam = lam  # regularization parameter lambda
        self.weights = None
        self.losses = []  # store loss values for each epoch
        self.margins = []  # store margin values (2 / ||w||) for each epoch
### ANCHOR_END: svm_init
    
### ANCHOR: svm_fit
    def fit(self, X, y):
        n_samples, n_features = X.shape
        self.weights = np.random.randn(n_features)
        
        for epoch in range(self.n_iterations):
            epoch_loss = 0
            
            for i, (x_i, y_i) in enumerate(zip(X, y)):
                # Calculate prediction
                prediction = np.dot(x_i, self.weights)
                
                # Calculate hinge loss for this sample
                hinge_loss = max(0, 1 - y_i * prediction)
                
                # Update weights based on whether point is misclassified
                if y_i * prediction < 1:  # misclassified or within margin
                    # Gradient of hinge loss + regularization
                    self.weights = (1 - self.learning_rate * self.lam) * self.weights + self.learning_rate * y_i * x_i
                else:  # correctly classified
                    # Only regularization term
                    self.weights = (1 - self.learning_rate * self.lam) * self.weights
                
                # Accumulate loss for this epoch
                epoch_loss += hinge_loss
            
            # Calculate total loss for this epoch (hinge loss + regularization)
            regularization_loss = 0.5 * self.lam * np.dot(self.weights, self.weights)
            total_loss = epoch_loss / n_samples + regularization_loss
            self.losses.append(total_loss)
            
            # Calculate margin (2 / ||w||)
            weight_norm = np.linalg.norm(self.weights)
            margin = 2 / weight_norm if weight_norm > 0 else 0
            self.margins.append(margin)
            
            if epoch % 10 == 0:
                print(f"Epoch {epoch}, Loss: {total_loss:.4f}, Margin: {margin:.4f}")
### ANCHOR_END: svm_fit

### ANCHOR: svm_predict
    def predict(self, X):
        return np.sign(np.dot(X, self.weights))
### ANCHOR_END: svm_predict

np.random.seed(SEED)
### ANCHOR: fit_svm_model
svm_model = SupportVectorMachine(learning_rate=0.01, n_iterations=100, lam=0.1)
svm_model.fit(X, y)
y_pred_svm = svm_model.predict(X)
### ANCHOR_END: fit_svm_model

### ANCHOR: calculate_svm_accuracy
accuracy_svm = np.mean(y_pred_svm == y)
print(f"SVM Accuracy: {accuracy_svm}")
### ANCHOR_END: calculate_svm_accuracy

### ANCHOR: plot_svm_decision_boundary
fig, ax = plt.subplots(figsize=(7, 6))

# Plot decision boundary
ax.scatter(X[:, 0], X[:, 1], c=y, cmap='coolwarm', alpha=0.7)

# Define the decision boundary as a function of the first feature
x1_range = np.linspace(X[:, 0].min(), X[:, 0].max(), 100)
if svm_model.weights[1] != 0:
    x2_boundary = -(svm_model.weights[0] * x1_range + svm_model.weights[2]) / svm_model.weights[1]
    ax.plot(x1_range, x2_boundary, 'k--', linewidth=2, label='SVM Decision Boundary')

ax.legend(loc='upper right')
ax.set_xlabel('PC1')
ax.set_ylabel('PC2')
ax.set_title('SVM Decision Boundary')
ax.set_xlim(X[:, 0].min()-0.1, X[:, 0].max()+0.1)
ax.set_ylim(X[:, 1].min()-0.1, X[:, 1].max()+0.1)

plt.show()
### ANCHOR_END: plot_svm_decision_boundary


