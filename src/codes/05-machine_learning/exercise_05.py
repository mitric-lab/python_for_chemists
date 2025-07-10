#! /usr/bin/env python

import numpy as np

class LinearRegression:
    def __init__(self):
        self.weights = None
    
    def fit(self, X, y):
        # self.weights = np.linalg.pinv(X) @ y
        self.weights = np.linalg.inv(X.T @ X) @ X.T @ y
    
    def predict(self, X):
        return X @ self.weights
    
### ANCHOR: exercise_01_a
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

path_to_csv = "aptamer_fingerprints_regression_data.csv"
df = pd.read_csv(path_to_csv)
df = df.loc[:, df.std() != 0]
print(df.head())

target_column = "fl_int"

X = df.drop(target_column, axis=1).values # Ignore the target column
X = np.hstack([np.ones((X.shape[0], 1)), X]) # Add a column of ones to the data matrix
y = df[target_column].values # Target column

print(np.linalg.matrix_rank(X))
print(X.shape)
print(y.shape)

model = LinearRegression() # Create a model
model.fit(X, y) # Fit the model to the data
y_pred = model.predict(X) # Predict the target values

mae = np.mean(np.abs(y - y_pred))
print(f"MAE: {mae}")

fig, ax = plt.subplots(figsize=(6, 6))
ax.scatter(y, y_pred, color='blue', label='Data')
ax.set_xlabel(f"True Target")
ax.set_ylabel(f"Predicted Target")
ax.set_title(f"True vs Predicted Targets")

# Get the min and max across both axes to set equal limits
min_val = min(min(y), min(y_pred)) - 0.1
max_val = max(max(y), max(y_pred)) + 0.1
ax.set_xlim(min_val, max_val)
ax.set_ylim(min_val, max_val)

# Add diagonal line representing perfect predictions
ax.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.5, label='Perfect Prediction')

ax.legend()
plt.show()
### ANCHOR_END: exercise_01_a

### ANCHOR: exercise_01_b
class RidgeRegression:
    def __init__(self, lambda_=0.1):
        self.weights = None
        self.lambda_ = lambda_
    
    def fit(self, X, y):
        self.weights = np.linalg.pinv(X.T @ X + self.lambda_ * np.eye(X.shape[1])) @ X.T @ y # Ridge regression
    
    def predict(self, X):
        return X @ self.weights

model = RidgeRegression() # Create a model
model.fit(X, y) # Fit the model to the data
y_pred = model.predict(X) # Predict the target values

mae = np.mean(np.abs(y - y_pred))
print(f"MAE: {mae}")

fig, ax = plt.subplots(figsize=(6, 6))
ax.scatter(y, y_pred, color='blue', label='Data')
ax.set_xlabel(f"True Target")
ax.set_ylabel(f"Predicted Target")
ax.set_title(f"True vs Predicted Targets")

# Get the min and max across both axes to set equal limits
min_val = min(min(y), min(y_pred)) - 0.1
max_val = max(max(y), max(y_pred)) + 0.1
ax.set_xlim(min_val, max_val)
ax.set_ylim(min_val, max_val)

# Add diagonal line representing perfect predictions
ax.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.5, label='Perfect Prediction')

ax.legend()
plt.show()