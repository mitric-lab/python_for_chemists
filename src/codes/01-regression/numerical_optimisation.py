#!/usr/bin/env python

### ANCHOR: import_numpy
import numpy as np
### ANCHOR_END: import_numpy


### ANCHOR: finite_difference
def finite_difference(func, beta, h=1e-5, args=()):
    n = len(beta)
    grad = np.zeros(n)

    for i in range(0, n):
        e = np.zeros(n)
        e[i] = 1
        grad[i] = (
            func(beta + h * e, *args) - func(beta - h * e, *args)
        ) / (2 * h)

    return grad
### ANCHOR_END: finite_difference


### ANCHOR: least_squares_loss
def least_squares_loss(beta, x, y):
    return np.sum((y - (beta[0] + beta[1] * x))**2)
### ANCHOR_END: least_squares_loss


### ANCHOR: gradient_descent
def gradient_descent(func, beta0, alpha=0.001, 
                     max_norm=1e-6, max_iter=10000, args=()):
    beta = np.copy(beta0)
    for niter in range(0, max_iter):
        grad = finite_difference(func, beta, args=args)
        beta = beta - alpha * grad
        if np.linalg.norm(grad) < max_norm:
            break
    if niter == max_iter - 1:
        print('Warning: Maximum iterations reached. '
              'Result may not be reliable.')
    
    return beta, niter
### ANCHOR_END: gradient_descent


### ANCHOR: data_array
concentrations = [
    2.125, 4.250, 6.375, 8.500, 10.63, 12.75, 14.88, 17.00, 19.13, 21.25,
    23.38, 25.50, 27.63, 29.75, 31.88, 34.00, 36.13, 38.25, 40.38, 42.50,
]
absorbances = [
    0.0572, 0.1391, 0.2049, 0.2754, 0.3420, 
    0.4139, 0.4956, 0.5815, 0.6806, 0.7481,
    0.8242, 0.9130, 1.0043, 1.0809, 1.1511,
    1.2483, 1.3373, 1.4027, 1.4927, 1.5853,
]
concentrations = np.array(concentrations)
absorbances = np.array(absorbances)
### ANCHOR_END: data_array

### ANCHOR: gradient_descent_call
beta_guess = np.array([1.0, 1.0])
beta_opt, niter = gradient_descent(
    least_squares_loss, 
    beta_guess,
    alpha=0.00005,
    max_iter=100000,
    args=(concentrations, absorbances),
)

assert np.isclose(beta_opt[0], -0.04907034)
assert np.isclose(beta_opt[1], 0.03800109)
### ANCHOR_END: gradient_descent_call

print(beta_opt)
print(niter)

### ANCHOR: scipy_minimize
from scipy.optimize import minimize

res = minimize(
    least_squares_loss,
    beta_guess,
    args=(concentrations, absorbances),
    method='CG',
    jac=lambda beta, x, y: finite_difference(least_squares_loss, beta, args=(x, y)),
    options={'maxiter': 10000, 'gtol': 1e-6},
)

beta_opt = res.x
niter = res.nit
assert np.isclose(beta_opt[0], -0.04907034)
assert np.isclose(beta_opt[1], 0.03800109)
### ANCHOR_END: scipy_minimize

print(beta_opt)
print(niter)

