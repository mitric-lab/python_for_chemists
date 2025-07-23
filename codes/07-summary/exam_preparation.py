#!/usr/bin/env python

import numpy as np
from typing import Callable

### ANCHOR: euler_errors
def dydx(x: float, y: float) -> float:
    k = 0.0039022970
    return -k * y

def euler_step(
    x_n: float,
    y_n: float,
    h: float,
    dydx: Callable[[float, float], float],
) -> float:
    return x_n + h * dydx(x_n, y_n)

def euler_method(
    x0: float, 
    y0: float, 
    h: float, 
    dydx: Callable[[float, float], float], 
    nsteps: int,
) -> np.ndarray:
    x = x0 * np.arange(0, nsteps) * h
    y = np.zeros(nsteps + 1)
    y[1] = y0

    for i in range(0, nsteps):
        y[i] = euler_step(x[i], y[i], h, dydx)

    return x, y
### ANCHOR_END: euler_errors

### ANCHOR: gradient_descent_errors
def double_well_gradient(x):  # x: float
    grad = 6 * x**5 - 2 * x - 0.25
    return grad

def gradient_descent(func_grad, x0, tau=0.01, maxgrad=1e-6, maxiter=500)
    x = 0
    converged = False
    
    for _ in range(0, maxiter):
        grad = func_grad(x)
        x = x - tau * grad
        if grad < maxgrad:
            converged = True
            break
    
        return x, converged

x_opt, converged = gradient_descent(double_well_gradient, -1.2)
print("Ein lokales Minimum liegt bei x = " , x_opt)  # x_opt = -0.8872
### ANCHOR_END: gradient_descent_errors

### ANCHOR: gradient_descent_correct
def double_well_gradient(x):  # x: float
    grad = 6 * x**5 - 4 * x - 0.25  # Error 1: -2*x should be -4*x (derivative of -2x^2 is -4x)
    return grad

def gradient_descent(func_grad, x0, tau=0.01, maxgrad=1e-6, maxiter=500):
    x = x0  # Error 2: Should initialize x with x0, not 0
    converged = False

    for _ in range(0, maxiter):
        grad = func_grad(x)
        x = x - tau * grad
        if abs(grad) < maxgrad:  # Error 3: Should use abs(grad) for convergence check
            converged = True
            break

    return x, converged  # Error 4: return must be outside the for loop

x_opt, converged = gradient_descent(double_well_gradient, -1.2)
print("Ein lokales Minimum liegt bei x = " , x_opt)  # x_opt = -0.8872
### ANCHOR_END: gradient_descent_correct

### ANCHOR: knn_incomplete
import numpy as np
import matplotlib.pyplot as plt

class kNNClassifier:
    def __init__(self, k):
        self.k = k

    def predict(self, X, y, xi):
        # TODO: Implement the kNN label prediction for xi

        






        y_pred = ...
        return y_pred
### ANCHOR_END: knn_incomplete

### ANCHOR: knn_weighted_incomplete
class kNNWeightedClassifier
    def __init__(self, k, sigma=1.0):
        self.k = k
        self.sigma = sigma

    def predict(self, X, y, xi):
        # TODO: Implement the weighted kNN label prediction for xi








        y_pred = ...
        return y_pred
### ANCHOR_END: knn_weighted_incomplete

### ANCHOR: knn_complete
import numpy as np
import matplotlib.pyplot as plt

class kNN_Classifier:
    def __init__(self, k):
        self.k = k

    def predict(self, X, y, xi):
        distances = np.linalg.norm(X - xi, axis=1)
        nearest = np.argsort(distances)[:self.k]
        y_nearest = y[nearest]
        unique_labels, label_counts = np.unique(y_nearest, return_counts=True)
        y_pred = unique_labels[np.argmax(label_counts)]
        return y_pred
### ANCHOR_END: knn_complete

### ANCHOR: knn_weighted_complete
class kNNWeightedClassifier:
    def __init__(self, k, sigma=1.0):
        self.k = k
        self.sigma = sigma

    def predict(self, X, y, xi):
        distances = np.linalg.norm(X - xi, axis=1)
        weights = np.exp(-distances**2 / (2 * self.sigma**2))
        nearest = np.argsort(distances)[:self.k]
        y_nearest = y[nearest]
        unique_labels, label_counts = np.unique(y_nearest, return_counts=True)
        y_pred = unique_labels[np.argmax(label_counts)]
        return y_pred
### ANCHOR_END: knn_weighted_complete

### ANCHOR: knn_example
N = 20
X = np.random.randn(N, 2)
y = np.hstack((np.zeros(N//2), np.ones(N//2)))

knn = kNN_Classifier(k=3)

xi = np.array([0, 0])
y_pred = knn.predict(X, y, xi)

print(y_pred)  # Output: 0 or 1
### ANCHOR_END: knn_example

fig, ax = plt.subplots(figsize=(3, 3))

ax.scatter(X[y == 0][:, 0], X[y == 0][:, 1], c='red')
ax.scatter(X[y == 1][:, 0], X[y == 1][:, 1], c='blue')
ax.scatter(xi[0], xi[1], c='green')

circle = plt.Circle(xi, 0.5, color='black', fill=False)
ax.add_artist(circle)
circle = plt.Circle(xi, 1.0, color='black', fill=False, linestyle='dashed')
ax.add_artist(circle)

ax.set_xlim(-2, 2)
ax.set_ylim(-2, 2)

ax.set_xticks([])
ax.set_yticks([])

fig.tight_layout()

plt.show()

#fig.savefig('../../assets/figures/07-summary/kNN.svg')
