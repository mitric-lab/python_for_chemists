#!/usr/bin/env python

### ANCHOR: imports
import numpy as np
import matplotlib.pyplot as plt
from typing import Callable
### ANCHOR_END: imports


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


### ANCHOR: rk4_step
def rk4_step(
    x_n: float,
    y_n: np.ndarray,
    h: float,
    dydx: Callable[[float, np.ndarray], np.ndarray],
) -> np.ndarray:
    a21 = 1.0 / 3.0
    a31 = -1.0 / 3.0
    a32 = 1.0
    a41 = 1.0
    a42 = -1.0
    a43 = 1.0
    b1 = 1.0 / 8.0
    b2 = 3.0 / 8.0
    b3 = 3.0 / 8.0
    b4 = 1.0 / 8.0
    c2 = 1.0 / 3.0
    c3 = 2.0 / 3.0
    c4 = 1.0

    k1 = dydx(x_n, y_n)
    k2 = dydx(x_n + h * c2, y_n + h * a21 * k1)
    k3 = dydx(x_n + h * c3, y_n + h * (a31 * k1 + a32 * k2))
    k4 = dydx(x_n + h * c4, y_n + h * (a41 * k1 + a42 * k2 + a43 * k3))

    return y_n + h * (b1 * k1 + b2 * k2 + b3 * k3 + b4 * k4)
### ANCHOR_END: rk4_step


### ANCHOR: rk4_method
def rk4_method(
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
        y[:, i + 1] = rk4_step(x[i], y[:, i], h, dydx)

    return x, y
### ANCHOR_END: rk4_method


### ANCHOR: solve_ode
CX_0 = 0.0    # M
CY_0 = 0.001  # M
CZ_0 = 0.0    # M
C0 = np.array([CX_0, CY_0, CZ_0])

T0 = 0.0
STEP = 0.001
TMAX = 200.0

nsteps = int(TMAX / STEP)
x, y = rk4_method(0, C0, STEP, dydx, nsteps)
### ANCHOR_END: solve_ode

### ANCHOR: plot
c_x, c_y, c_z = y * 1000.0  # convert to mM

fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(x, c_x, label='[HBrO2]')
ax.plot(x, c_y, label='[Br-]')
ax.plot(x, c_z, label='2 [Ce4+]')

ax.set_xlabel('time / s')
ax.set_ylabel('concentration / mM')
ax.set_xlim(0, 200)
ax.set_ylim(-0.1, 2.0)

fig.tight_layout()
ax.legend(loc='upper right')

plt.show()
### ANCHOR_END: plot

fig.savefig('../../assets/figures/02-differential_equations/rk4_bz.svg')

