#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from typing import Callable


### ANCHOR: dydx
def dydx(x: float, y: float) -> float:
    k = 0.0039022970
    return -k * y
### ANCHOR_END: dydx


### ANCHOR: euler_step
def euler_step(
    x_n: float,
    y_n: float,
    h: float,
    dydx: Callable[[float, float], float],
) -> float:
    return y_n + h * dydx(x_n, y_n)
### ANCHOR_END: euler_step


### ANCHOR: euler_method
def euler_method(
    x0: float, 
    y0: float, 
    h: float, 
    dydx: Callable[[float, float], float], 
    n: int,
) -> np.ndarray:
    x = np.arange(0, n + 1) * h
    y = np.zeros(n + 1)
    y[0] = y0

    for i in range(0, n):
        y[i + 1] = euler_step(x[i], y[i], h, dydx)

    return x, y
### ANCHOR_END: euler_method

### ANCHOR: solve_ode
C0 = 1.0  # M
T0 = 0.0
STEP = 1.0  # s
MAXTIME = 900.0  # s

nsteps = int(MAXTIME / STEP)
x, y = euler_method(0, C0, STEP, dydx, nsteps)
### ANCHOR_END: solve_ode


### ANCHOR: plot
fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(x, C0 * np.exp(-0.0039022970 * x), label='analytical solution')
ax.plot(x, y, label='numerical solution')

ax.set_xlabel('time / s')
ax.set_ylabel('concentration / M')
ax.legend(loc='upper right')

fig.tight_layout()
plt.show()
### ANCHOR_END: plot

fig.savefig('../../assets/figures/02-differential_equations/euler_mn.svg')

