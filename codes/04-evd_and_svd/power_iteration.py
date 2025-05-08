#!/usr/bin/env python

### ANCHOR: import
import numpy as np
### ANCHOR_END: import


### ANCHOR: power_iteration_function
def power_iteration(
    a: np.ndarray, eps: float, maxiter: int = 1000,
) -> (float, np.ndarray):
    n = a.shape[0]
    
    x = np.random.rand(n)
    x /= np.linalg.norm(x)
    x_prev = np.zeros(n)
    
    for i in range(0, maxiter):
        x_prev = x
        x = np.dot(a, x)
        x /= np.linalg.norm(x)
        if np.linalg.norm(x - x_prev) < eps:
            break
    if i == maxiter - 1:
        print('Power iteration did not converge.')
        print('Residual norm:', np.linalg.norm(x - x_prev))

    lambda_ = np.dot(np.dot(x, a), x)

    return lambda_, x
### ANCHOR_END: power_iteration_function


### ANCHOR: test_power_iteration
NDIM = 10

a_mat = np.random.rand(NDIM, NDIM)
a_mat = a_mat + a_mat.T

lambda_, x = power_iteration(a_mat, 1e-8)
### ANCHOR_END: test_power_iteration

### ANCHOR: verify_results
eigvals, eigvecs = np.linalg.eigh(a_mat)
lambda_ref = eigvals[-1]
x_ref = eigvecs[:, -1]

print(lambda_, lambda_ref)
print(x, x_ref)
assert np.isclose(lambda_, lambda_ref)
assert np.allclose(np.abs(np.dot(x, x_ref)), 1.0)
### ANCHOR_END: verify_results

