#!/usr/bin/env python

import numpy as np
import pandas as pd


def multilinear_regression(X, y):
    X = np.hstack([np.ones((X.shape[0], 1)), X])
    theta = np.linalg.inv(X.T @ X) @ X.T @ y
    return theta


def mean_squared_error(y_true, y_pred):
    return np.mean((y_true - y_pred) ** 2)


def coefficient_of_determination(y_true, y_pred):
    y_mean = np.mean(y_true)
    numerator = np.sum(
        (y_true - y_mean) * (y_pred - y_mean)
    )**2
    denominator = np.sum(
        (y_true - y_mean)**2
    ) * np.sum(
        (y_pred - y_mean)**2
    )
    return numerator / denominator


url = 'https://archive.ics.uci.edu/ml/machine-learning-databases/wine-quality/winequality-red.csv'
df = pd.read_csv(url, sep=';')

# alcohol and volatile acidity
X = df[['alcohol', 'volatile acidity']].to_numpy()
y = df['quality'].to_numpy()
theta = multilinear_regression(X, y)
y_pred = np.hstack([np.ones((X.shape[0], 1)), X]) @ theta
mse = mean_squared_error(y, y_pred)
r2 = coefficient_of_determination(y, y_pred)
print(f'Mean squared error (2 features): {mse}')
print(f'Coefficient of determination (2 features): {r2}')

# all features
X = df.drop('quality', axis=1).to_numpy()
y = df['quality'].to_numpy()
theta = multilinear_regression(X, y)
y_pred = np.hstack([np.ones((X.shape[0], 1)), X]) @ theta
mse = mean_squared_error(y, y_pred)
r2 = coefficient_of_determination(y, y_pred)
print(f'Mean squared error (all features): {mse}')
print(f'Coefficient of determination (all features): {r2}')

