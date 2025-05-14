#!/usr/bin/env python

### ANCHOR: import
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
### ANCHOR_END: import


### ANCHOR: titration_model
def titration_sasb_model(v_b, c0_b, c0_a, v0):
    k_w = 1e-14
    c_a = c0_a * v0 / (v0 + v_b)
    c_b = c0_b * v_b / (v0 + v_b)
    delta = c_a - c_b
    c_h = 0.5 * (delta + np.sqrt(delta**2 + 4 * k_w))
    ph = -np.log10(c_h)
    return ph
### ANCHOR_END: titration_model


### ANCHOR: objective_function
def objective_function(beta, x, y, c0_b):
    # x: v_b; y: ph
    c0_a, v0 = beta
    ph_fit = titration_sasb_model(x, c0_b, c0_a, v0)
    return np.sum((y - ph_fit)**2)
### ANCHOR_END: objective_function

### ANCHOR: read_data
C0_B = 0.1  # mol/l
v_b, ph = np.loadtxt('titration_sasb.txt', unpack=True)
### ANCHOR_END: read_data


### ANCHOR: optimise
beta_guess = (0.01, 100)
res = minimize(
    objective_function, 
    beta_guess, 
    args=(v_b, ph, C0_B),
    method='Nelder-Mead',
    options={'maxiter': 1000},
)
print(res.message)

c0_a, v0 = res.x
print(f'c0_a = {c0_a} mol/l')
print(f'v0 = {v0} ml')
print(f'n0_a = {c0_a * v0} mmol')
### ANCHOR_END: optimise


### ANCHOR: plot
fig, ax = plt.subplots(figsize=(8, 6))
v_b_interp = np.linspace(v_b.min(), v_b.max(), 1000)

ax.plot(v_b, ph, 'o', label='data')
ax.plot(
    v_b_interp, 
    titration_sasb_model(v_b_interp, C0_B, c0_a, v0), 
    label='pH fit',
)

ax.set_xlabel('volume of base / ml')
ax.set_ylabel('pH')
ax.legend()

fig.tight_layout()
plt.show()
### ANCHOR_END: plot

fig.savefig('../../assets/figures/01-regression/nonlinreg_titration.svg')


def titration_sasb_model_ch(
    v_b: np.ndarray,
    c0_b: float,
    c0_a: float,
    v0: float,
) -> np.ndarray:
    k_w = 1e-14
    c_a = c0_a * v0 / (v0 + v_b)
    c_b = c0_b * v_b / (v0 + v_b)
    delta = c_a - c_b
    c_h = 0.5 * (delta + np.sqrt(delta**2 + 4 * k_w))
    return c_h


def objective_function_ch(beta, *args):
    c0_a, v0 = beta
    c0_b, v_b, ph = args
    c_h = np.power(10, -ph)
    c_h_fit = titration_sasb_model_ch(v_b, c0_b, c0_a, v0)
    return np.sum((c_h - c_h_fit)**2)


res = minimize(
    objective_function_ch, 
    beta_guess, 
    args=(C0_B, v_b, ph),
    method='Nelder-Mead',
)
ax.plot(
    v_b_interp, 
    -np.log10(titration_sasb_model_ch(v_b_interp, C0_B, res.x[0], res.x[1])),
    label='[H+] fit',
)
ax.legend()
fig.savefig('../../assets/figures/01-regression/nonlinreg_titration_ch.svg')

