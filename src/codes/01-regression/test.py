#!/usr/bin/env python

### ANCHOR: import_numpy
import numpy as np
from typing import Callable, Any
### ANCHOR_END: import_numpy


### ANCHOR: finite_difference
def finite_difference(
    func: Callable[[np.ndarray, Any], float],
    x0: np.ndarray,
    h: float = 1e-5,
    args: tuple = (),
) -> np.ndarray:
    n: int = len(x0)
    grad: np.ndarray = np.zeros(n)

    for i in range(0, n):
        e: np.ndarray = np.zeros(n)
        e[i] = 1
        grad[i] = (
            func(x0 + h * e, *args) - func(x0 - h * e, *args)
        ) / (2 * h)

    return grad
### ANCHOR_END: finite_difference


### ANCHOR: objective_function
def objective_function(
    beta: np.ndarray,
    *args: np.ndarray,
) -> float:
    concentrations: np.ndarray = args[0]
    absorbances: np.ndarray = args[1]
    return np.sum((absorbances - (beta[0] + beta[1] * concentrations))**2)
### ANCHOR_END: objective_function


### ANCHOR: objective_function_gradient
def objective_function_gradient(
    beta: np.ndarray,
    *args: np.ndarray,
) -> np.ndarray:
    concentrations: np.ndarray = args[0]
    absorbances: np.ndarray = args[1]
    grad: np.ndarray = finite_difference(
        objective_function, 
        beta, 
        args=(concentrations, absorbances),
    )
    return grad
### ANCHOR_END: objective_function_gradient


### ANCHOR: gradient_descent
def gradient_descent(
    func_grad: Callable[[np.ndarray, Any], np.ndarray],
    x0: np.ndarray,
    alpha: float = 0.001,
    max_norm: float = 1e-6,
    max_iter: int = 10000,
    args: tuple = (),
) -> tuple[np.ndarray, int]:
    x: np.ndarray = x0
    for niter in range(0, max_iter):
        grad: np.ndarray = func_grad(x, *args)
        x = x - alpha * grad
        if np.linalg.norm(grad) < max_norm:
            break
    if niter == max_iter - 1:
        print('Warning: Maximum iterations reached. '
              'Result may not be reliable.')
    
    return x, niter
### ANCHOR_END: gradient_descent


### ANCHOR: data_array
concentrations = [
    0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
    1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
]
absorbances = [
    0.0522, 0.1010, 0.1518, 0.2017, 0.2453, 
    0.2920, 0.3404, 0.3831, 0.4240, 0.4685,
    0.5080, 0.5486, 0.5899, 0.6245, 0.6630,
    0.7020, 0.7382, 0.7766, 0.8091, 0.8424,
]
concentrations: np.ndarray = np.array(concentrations)
absorbances: np.ndarray = np.array(absorbances)
### ANCHOR_END: data_array

### ANCHOR: gradient_descent_call
beta_guess = np.array([1.0, 1.0])
beta_opt, niter = gradient_descent(
    objective_function_gradient, 
    beta_guess,
    alpha=0.001,
    args=(concentrations, absorbances),
)

beta0, beta1 = beta_opt
assert np.isclose(beta0, 0.03735158)
assert np.isclose(beta1, 0.41501278)

print(beta_opt)
print(niter)
### ANCHOR_END: gradient_descent_call


### ANCHOR: scipy_minimize
from scipy.optimize import minimize

res = minimize(
    objective_function,
    beta_guess,
    args=(concentrations, absorbances),
    method='CG',
    options={'maxiter': 10000, 'gtol': 1e-6},
)

beta0, beta1 = res.x
niter = res.nit
assert np.isclose(beta0, 0.03735158)
assert np.isclose(beta1, 0.41501278)
### ANCHOR_END: scipy_minimize

print(res.x)
print(niter)
