#! /usr/bin/env python

class LinearRegression:
    def __init__(self):
        self.weights = None
    
    def fit(self, X, y):
        self.weights = np.linalg.pinv(X) @ y
    
    def predict(self, X):
        return X @ self.weights

## ANCHOR: exercise_01
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def split_data(X, y, test_size=0.2):
    """Split the data into training and test sets."""
    n_samples = X.shape[0]
    n_test = int(n_samples * test_size)
    indices = np.random.permutation(n_samples)

    X_train = X[indices[:-n_test]]
    y_train = y[indices[:-n_test]]
    X_test = X[indices[-n_test:]]
    y_test = y[-n_test:]

    return X_train, y_train, X_test, y_test
### ANCHOR_END: exercise_01

np.random.seed(42)

### ANCHOR: exercise_02_a
path_to_csv = "aptamer_fingerprints_regression_data.csv"
df = pd.read_csv(path_to_csv)
print(df.head())

target_column = "fl_int"

X = df.drop(target_column, axis=1).values # Ignore the target column
X = np.hstack([np.ones((X.shape[0], 1)), X]) # Add a column of ones to the data matrix
y = df[target_column].values # Target column

X_train, y_train, X_test, y_test = split_data(X, y)

try:
    model = LinearRegression() # Create a model
    model.fit(X_train, y_train) # Fit the model to the training data
except Exception as e:
    print(f"Error: {e}")
### ANCHOR_END: exercise_02_a

### ANCHOR: exercise_02_b
class LinearRegression:
    def __init__(self):
        self.weights = None
    
    def fit(self, X, y):
        self.weights = np.linalg.pinv(X) @ y
    
    def predict(self, X):
        return X @ self.weights
    
model = LinearRegression() # Create a model
model.fit(X_train, y_train) # Fit the model to the training data
y_train_pred = model.predict(X_train) # Predict the target values for the training data
y_test_pred = model.predict(X_test) # Predict the target values for the test data

mae_train = np.mean(np.abs(y_train - y_train_pred))
mae_test = np.mean(np.abs(y_test - y_test_pred))
print(f"MAE train: {mae_train}")
print(f"MAE test: {mae_test}")

fig, ax = plt.subplots(figsize=(6, 6))
ax.scatter(y_train, y_train_pred, color='red', label='Train Data')
ax.scatter(y_test, y_test_pred, color='blue', label='Test Data')

ax.set_xlabel(f"True Target")
ax.set_ylabel(f"Predicted Target")
ax.set_title(f"True vs Predicted Targets")

# Get the min and max across both axes to set equal limits
min_val = min(min(y_test), min(y_test_pred)) - 0.1
max_val = max(max(y_test), max(y_test_pred)) + 0.1
ax.set_xlim(min_val, max_val)
ax.set_ylim(min_val, max_val)

# Add diagonal line representing perfect predictions
ax.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.5, label='Perfect Prediction')

ax.legend()
plt.show()
### ANCHOR_END: exercise_02_b

# fig.savefig('../../assets/figures/05-machine_learning/aptamer_fingerprints_linear_regression.svg')

### ANCHOR: exercise_02_c
class RidgeRegression:
    def __init__(self, lambda_=0.1):
        self.weights = None
        self.lambda_ = lambda_
    
    def fit(self, X, y):
        self.weights = np.linalg.pinv(X.T @ X + self.lambda_ * np.eye(X.shape[1])) @ X.T @ y # Ridge regression
    
    def predict(self, X):
        return X @ self.weights

model = RidgeRegression(lambda_=1.0) # Create a model
model.fit(X_train, y_train) # Fit the model to the training data
y_train_pred = model.predict(X_train) # Predict the target values for the training data
y_test_pred = model.predict(X_test) # Predict the target values for the test data

mae_train = np.mean(np.abs(y_train - y_train_pred))
mae_test = np.mean(np.abs(y_test - y_test_pred))
print(f"MAE train: {mae_train}")
print(f"MAE test: {mae_test}")

fig, ax = plt.subplots(figsize=(6, 6))
ax.scatter(y_train, y_train_pred, color='red', label='Train Data')
ax.scatter(y_test, y_test_pred, color='blue', label='Test Data')
ax.set_xlabel(f"True Target")
ax.set_ylabel(f"Predicted Target")
ax.set_title(f"True vs Predicted Targets")

# Get the min and max across both axes to set equal limits
min_val = min(min(y_test), min(y_test_pred)) - 0.1
max_val = max(max(y_test), max(y_test_pred)) + 0.1
ax.set_xlim(min_val, max_val)
ax.set_ylim(min_val, max_val)

# Add diagonal line representing perfect predictions
ax.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.5, label='Perfect Prediction')

ax.legend()
plt.show()
### ANCHOR_END: exercise_02_c

# fig.savefig('../../assets/figures/05-machine_learning/aptamer_fingerprints_ridge_regression.svg')