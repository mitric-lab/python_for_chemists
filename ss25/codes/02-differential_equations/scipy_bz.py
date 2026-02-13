#!/usr/bin/env python

### ANCHOR: imports
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
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


### ANCHOR: solve_setup
CX_0 = 0.0    # M
CY_0 = 0.001  # M
CZ_0 = 0.0    # M
C0 = np.array([CX_0, CY_0, CZ_0])

T0 = 0.0
TMAX = 200.0
MAXSTEP = 0.1
### ANCHOR_END: solve_setup
'''
print('RK45:')
### ANCHOR: solve_rk45
res = solve_ivp(
    dydx, 
    (T0, TMAX), 
    C0, 
    method='RK45', 
    max_step=MAXSTEP,
)

x, y = res.t, res.y
nsteps = len(x) - 1
minstep = np.min(np.diff(x))

print(nsteps)
print(minstep)
### ANCHOR_END: solve_rk45

### ANCHOR: verification_rk45
assert nsteps == 27611
assert np.isclose(minstep, 0.00121700)
### ANCHOR_END: verification_rk45
print()

print('DOP853:')
MAXSTEP = 0.02
res = solve_ivp(
    dydx, 
    (T0, TMAX), 
    C0, 
    method='DOP853',
    max_step=MAXSTEP,
)
print(len(res.t) - 1)
print(np.min(np.diff(res.t)))
print()

print('Radau:')
MAXSTEP = 0.1
### ANCHOR: solve_radau
res = solve_ivp(
    dydx, 
    (T0, TMAX), 
    C0, 
    method='Radau',
    max_step=MAXSTEP,
)

x, y = res.t, res.y
nsteps = len(x) - 1
minstep = np.min(np.diff(x))

print(nsteps)
print(minstep)
### ANCHOR_END: solve_radau

### ANCHOR: verification_radau
assert nsteps == 2014
assert np.isclose(minstep, 0.01494846)
### ANCHOR_END: verification_radau
'''

### ANCHOR: solve_dense
res = solve_ivp(
    dydx,
    (T0, TMAX),
    C0,
    method='Radau',
    max_step=MAXSTEP,
    dense_output=True,
)
x_plot = np.linspace(T0, TMAX, 5000)
y_plot = res.sol(x_plot)
### ANCHOR_END: solve_dense

### ANCHOR: configuration_space_plot
c_x, c_y, c_z = y_plot * 1000.0  # convert to mM

fig1, ax1 = plt.subplots(figsize=(6, 6), subplot_kw={'projection': '3d'})
ax1.scatter(c_x, c_y, c_z, s=10, alpha=0.1)
ax1.set_xlabel('[HBrO2] / mM')
ax1.set_ylabel('[Br-] / mM')
ax1.set_zlabel('2 [Ce4+] / mM')

ax1.set_xlim(0.0, 0.4)
ax1.set_ylim(0.0, 1.0)
ax1.set_zlim(0.0, 1.8)

fig1.tight_layout(rect=[0, 0, 0.95, 1.00])

plt.show()
### ANCHOR_END: configuration_space_plot

fig1.savefig('../../assets/figures/02-differential_equations/bz_configuration_space.svg')

### ANCHOR: phase_space_plot
dzdt = dydx(x_plot, y_plot)[2] * 1000.0  # convert to mM/s

fig2, ax2 = plt.subplots(figsize=(6, 6))
ax2.scatter(c_z, dzdt, c='tab:green', alpha=0.1)
ax2.set_xlabel('2 [Ce4+] / mM')
ax2.set_ylabel('(2 d[Ce4+] / dt) / (mM / s)')

fig2.tight_layout()

plt.show()
### ANCHOR_END: phase_space_plot

fig2.savefig('../../assets/figures/02-differential_equations/bz_phase_space.svg')

