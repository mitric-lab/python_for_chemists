#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from typing import Callable


### ANCHOR: dydx
def dydx(x: float, y: np.ndarray) -> np.ndarray:
    # concentrations adapted from 
    # R. J. Field, H.-D. Försterling, J. Phys. Chem. 1986, 90, 5400–5407.
    k1 = 1.3    # M^-1 s^-1
    k2 = 2.4e6  # M^-1 s^-1
    k3 = 34.0   # M^-1 s^-1
    k4 = 3.0e3  # M^-1 s^-1
    k5 = 1.0    # M^-1 s^-1
    c_a = 0.1   # M
    c_b = 0.4   # M
    
    c_x, c_y, c_z = y

    dcxdt = k1 * c_a * c_y - k2 * c_x * c_y + k3 * c_a * c_x - 2.0 * k4 * c_x**2
    dcydt = -k1 * c_a * c_y - k2 * c_x * c_y + k5 * c_b * c_z
    dczdt = k3 * c_a * c_x - k5 * c_b * c_z

    return np.array([dcxdt, dcydt, dczdt])
### ANCHOR_END: dydx


### ANCHOR: euler_step
def euler_step(
    x_n: float,
    y_n: np.ndarray,
    h: float,
    dydx: Callable[[float, np.ndarray], np.ndarray],
) -> np.ndarray:
    return y_n + h * dydx(x_n, y_n)
### ANCHOR_END: euler_step


### ANCHOR: euler_method
def euler_method(
    x0: float, 
    y0: np.ndarray,
    h: float, 
    dydx: Callable[[float, np.ndarray], np.ndarray],
    n: int,
) -> np.ndarray:
    ndim = len(y0)

    x = np.arange(0, n + 1) * h
    y = np.zeros((ndim, n + 1))
    y[:, 0] = y0

    for i in range(0, n):
        y[:, i + 1] = euler_step(x[i], y[:, i], h, dydx)

    return x, y
### ANCHOR_END: euler_method


### ANCHOR: solve_ode_bad
CX_0 = 0.0   # M
CY_0 = 0.03  # M
CZ_0 = 0.0   # M
C0 = np.array([CX_0, CY_0, CZ_0])

T0 = 0.0
STEP = 1.0
TMAX = 200.0

nsteps = int(TMAX / STEP)
x, y = euler_method(0, C0, STEP, dydx, nsteps)
### ANCHOR_END: solve_ode_bad

### ANCHOR: solve_ode
CX_0 = 0.0   # M
CY_0 = 0.03  # M
CZ_0 = 0.0   # M
C0 = np.array([CX_0, CY_0, CZ_0])

T0 = 0.0
STEP = 0.00002
TMAX = 200.0

nsteps = int(TMAX / STEP)
x, y = euler_method(0, C0, STEP, dydx, nsteps)
### ANCHOR_END: solve_ode

### ANCHOR: plot
c_x, c_y, c_z = y

fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(x[::100], c_x[::100], label='[HBrO2]')
ax.plot(x[::100], c_y[::100], label='[Br-]')
ax.plot(x[::100], c_z[::100], label='2 [Ce4+]')

ax.set_xlabel('time / s')
ax.set_ylabel('concentration / M')
ax.set_xlim(0, 200)
ax.set_ylim(-0.0001, 0.002)

fig.tight_layout()
ax.legend(loc='upper right')

plt.show()
### ANCHOR_END: plot

fig.savefig('../../assets/figures/02-differential_equations/euler_bz.svg')

