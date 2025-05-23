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
def exp_decay(t, a0, k):
    return a0 * np.exp(-k * t)
### ANCHOR_END: exp_model


### ANCHOR: objective_function
def objective_function(beta, x, y):
    # x: time; y: absorbance
    a0, k = beta
    return np.sum((y - exp_decay(x, a0, k))**2)
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
time_interp = np.linspace(time.min(), time.max(), 1000)

ax.plot(time, absorbance, 'o', label='data')
ax.plot(time_interp, exp_decay(time_interp, a0, k), label='exponential fit')

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

ax.plot(time_interp, a0_lin * np.exp(-k_lin * time_interp), label='linearised fit')
ax.legend()
plt.show()

fig.savefig('../../assets/figures/01-regression/nonlinreg_mn_wlin.svg')

