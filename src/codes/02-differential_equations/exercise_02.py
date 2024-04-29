#!/usr/bin/env python

### ANCHOR: exercise_01_b
import numpy as np
import matplotlib.pyplot as plt

# Implement Euler's method
def dfdt(t, f, omega):
    x, v = f
    return np.array([v, -omega**2 * x])

def euler_step(t_n, f_n, dt, dfdt, omega):
    return f_n + dt * dfdt(t_n, f_n, omega)

def euler_method(t0, f0, dt, dfdt, omega, n):
    ndim = len(f0)

    t = t0 + np.arange(0, n + 1) * dt
    f = np.zeros((ndim, n + 1))
    f[:, 0] = f0

    for i in range(0, n):
        f[:, i + 1] = euler_step(t[i], f[:, i], dt, dfdt, omega)

    return t, f
### ANCHOR_END: exercise_01_b

### ANCHOR: exercise_01_c
# Implement analytical solution
def harm_osc(t, f_0, omega):
    x0, v0 = f_0
    x = v0 / omega * np.sin(omega * t) + x0 * np.cos(omega * t)
    v = v0 * np.cos(omega * t) - x0 * omega * np.sin(omega * t)
    return np.array([x, v])

# Set the initial conditions
t0 = 0.0
x0 = 1.0
v0 = 0.0
f0 = np.array([x0, v0])
omega = 1.0

# Set the time step and the number of steps
dt = 0.1
t_max = 10.0
n = int(t_max / dt)

# Solve the differential equation using Euler's method
t, f_euler = euler_method(t0, f0, dt, dfdt, omega, n)

# Calculate the analytical solution
f_exact = harm_osc(t, f0, omega)

# Plot the results
fig, [ax1, ax2] = plt.subplots(2, 1, figsize=(12, 6))

ax1.plot(t, f_exact[0,:], label='Exact')
ax1.plot(t, f_euler[0,:], label='Euler')

ax1.set_ylabel('position')
ax1.set_xticks([])

ax1.legend()

ax2.plot(t, f_exact[1,:], label='Exact')
ax2.plot(t, f_euler[1,:], label='Euler')

ax2.set_xlabel('time / s')
ax2.set_ylabel('velocity')

plt.subplots_adjust(hspace=0.0)
plt.show()
### ANCHOR_END: exercise_01_c

fig.savefig('../../assets/figures/02-differential_equations/euler_harm_osc.svg')

### ANCHOR: exercise_02
# Implement the Runge-Kutta method
def rk4_step(t_n, f_n, dt, dfdt, omega):

    a21 = 1.0 / 2.0
    a31 = 0.0
    a32 = 1.0 / 2.0
    a41 = 0.0
    a42 = 0.0
    a43 = 1.0
    b1 = 1.0 / 6.0
    b2 = 1.0 / 3.0
    b3 = 1.0 / 3.0
    b4 = 1.0 / 6.0
    c2 = 1.0 / 2.0
    c3 = 1.0 / 2.0
    c4 = 1.0

    k1 = dfdt(t_n, f_n, omega)
    k2 = dfdt(t_n + dt * c2, f_n + dt * a21 * k1, omega)
    k3 = dfdt(t_n + dt * c3, f_n + dt * (a31 * k1 + a32 * k2), omega)
    k4 = dfdt(t_n + dt * c4, f_n + dt * (a41 * k1 + a42 * k2 + a43 * k3), omega)

    return f_n + dt * (b1 * k1 + b2 * k2 + b3 * k3 + b4 * k4)

def rk4_method(t0, f0, dt, dfdt, omega, n):
    ndim = len(f0)

    t = t0 + np.arange(0, n + 1) * dt
    f = np.zeros((ndim, n + 1))
    f[:, 0] = f0

    for i in range(0, n):
        f[:, i + 1] = rk4_step(t[i], f[:, i], dt, dfdt, omega)

    return t, f

# Solve the differential equation using the Runge-Kutta method
t, f_rk4 = rk4_method(t0, f0, dt, dfdt, omega, n)

# Plot the results
fig, [ax1, ax2] = plt.subplots(2, 1, figsize=(12, 6))

ax1.plot(t, f_exact[0,:], label='Exact')
ax1.plot(t, f_euler[0,:], label='Euler')
ax1.plot(t, f_rk4[0,:], label='RK4')

ax1.set_ylabel('position')
ax1.set_xticks([])

ax1.legend()

ax2.plot(t, f_exact[1,:], label='Exact')
ax2.plot(t, f_euler[1,:], label='Euler')
ax2.plot(t, f_rk4[1,:], label='RK4')

ax2.set_xlabel('time / s')
ax2.set_ylabel('velocity')

plt.subplots_adjust(hspace=0.0)
plt.show()
### ANCHOR_END: exercise_02

fig.savefig('../../assets/figures/02-differential_equations/euler_rk4_harm_osc.svg')

