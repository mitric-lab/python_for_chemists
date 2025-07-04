#!/usr/bin/env python

### ANCHOR: load_data_from_csv
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

path_to_csv = "aptamer_regression_data.csv"
df = pd.read_csv(path_to_csv)
print(df.head())
### ANCHOR_END: load_data_from_csv

### ANCHOR: process_data
target_column = "fl_int"

X = df.drop(target_column, axis=1).values # Ignore the target column
X = np.hstack([np.ones((X.shape[0], 1)), X]) # Add a column of ones to the data matrix
y = df[target_column].values # Target column

print(X.shape)
print(y.shape)
### ANCHOR_END: process_data

### ANCHOR: regression_class
class LinearRegression:
    def __init__(self):
        self.weights = None
    
    def fit(self, X, y):
        self.weights = np.linalg.inv(X.T @ X) @ X.T @ y
    
    def predict(self, X):
        return X @ self.weights
### ANCHOR_END: regression_class

### ANCHOR: fit_model
model = LinearRegression() # Create a model
model.fit(X, y) # Fit the model to the data
y_pred = model.predict(X) # Predict the target values
### ANCHOR_END: fit_model

### ANCHOR: calculate_mae
mae = np.mean(np.abs(y - y_pred))
print(f"MAE: {mae}")
### ANCHOR_END: calculate_mae

### ANCHOR: plot_predictions
# Plot the true and predicted values
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
### ANCHOR_END: plot_predictions

# fig.savefig('../../assets/figures/05-machine_learning/regression_predictions.svg')

### ANCHOR: plot_model_weights
feature_names = df.drop(target_column, axis=1).columns
feature_weights = model.weights[1:] # Exclude the bias term (which is the first weight)

fig, ax = plt.subplots(figsize=(7, 6))
ax.bar(feature_names, feature_weights)
ax.set_ylabel("Weight")
ax.set_title("Learned Feature Weights")
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()
### ANCHOR_END: plot_model_weights

# fig.savefig('../../assets/figures/05-machine_learning/regression_model_weights.svg')





