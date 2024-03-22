#!/usr/bin/env python

### ANCHOR: import
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
### ANCHOR_END: import

### ANCHOR: read_data
time, absorbance = np.loadtxt('mn_decay.txt', unpack=True)
### ANCHOR_END: read_data


### ANCHOR: exp_model
def exp_decay(t: np.ndarray, a0: float, k: float) -> np.ndarray:
    return a0 * np.exp(-k * t)
### ANCHOR_END: exp_model


### ANCHOR: objective_function
def objective_function(beta, *args):
    a0, k = beta
    time, absorbance = args
    return np.sum((absorbance - exp_decay(time, a0, k))**2)
### ANCHOR_END: objective_function


### ANCHOR: optimise
beta_guess = (1.0, 0.01)
res = minimize(
    objective_function, beta_guess, method='Nelder-Mead', 
    args=(time, absorbance),
)
a0, k = res.x

print(f'k = {k} s^-1')
print(f'A0 = {a0}')
### ANCHOR_END: optimise

### ANCHOR: verification
assert np.isclose(k, 0.0039022970)
assert np.isclose(a0, 1.6475263)
### ANCHOR_END: verification

### ANCHOR: plot
fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(time, absorbance, 'o', label='data')
ax.plot(time, exp_decay(time, a0, k), label='exponential fit')

ax.set_xlabel('time / s')
ax.set_ylabel('absorbance')
ax.legend()

fig.tight_layout()
plt.show()
### ANCHOR_END: plot

fig.savefig('../../assets/figures/01-regression/nonlinreg_mn.svg')

beta1, beta0 = np.polyfit(time, np.log(absorbance), 1)
a0_lin = np.exp(beta0)
k_lin = -beta1
print(f'k_lin = {k_lin} s^-1')
print(f'A0_lin = {a0_lin}')

### ANCHOR: verification_lin
assert np.isclose(k_lin, 0.0041675912)
assert np.isclose(a0_lin, 1.7633374)
### ANCHOR_END: verification_lin

ax.plot(time, a0_lin * np.exp(-k_lin * time), label='linearised fit')
ax.legend()
plt.show()

fig.savefig('../../assets/figures/01-regression/nonlinreg_mn_wlin.svg')

