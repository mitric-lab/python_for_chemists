#!/usr/bin/env python

### ANCHOR: exercise_01_b
import numpy as np
import matplotlib.pyplot as plt

# Implement Euler's method
def dfdt(x, v, omega):
    return v, -omega**2 * x

def euler_step(x_n, v_n, dt, dfdt, omega):
    dx, dv = dfdt(x_n, v_n, omega)
    x_n += dx * dt
    v_n += dv * dt
    return x_n, v_n

def euler_method(x0, v0, dt, dfdt, omega, n):
    t = np.arange(0, n + 1) * dt
    x = np.zeros(n + 1)
    v = np.zeros(n + 1)
    x[0] = x0
    v[0] = v0

    for i in range(0, n):
        x[i + 1], v[i + 1] = euler_step(x[i], v[i], dt, dfdt, omega)

    return t, x, v
### ANCHOR_END: exercise_01_b

### ANCHOR: exercise_01_c
# Implement analytical solution
def harmonic_oscillator(t, x0, v0, omega):
    x = v0 / omega * np.sin(omega * t) + x0 * np.cos(omega * t)
    v = v0 * np.cos(omega * t) - x0 * omega * np.sin(omega * t)
    return x, v

# Set the initial conditions
x0 = 1
v0 = 0
omega = 1

# Set the time step and the number of steps
dt = 0.1
t_max = 10.0
n = int(t_max / dt)

# Solve the differential equation using Euler's method
t, x_euler, v_euler = euler_method(x0, v0, dt, dfdt, omega, n)

# Calculate the analytical solution
x_exact, v_exact = harmonic_oscillator(t, x0, v0, omega)

# Plot the results
fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(16, 6))

ax1.plot(t, x_euler, label='Euler')
ax1.plot(t, x_exact, label='Exact')

ax1.set_xlabel('time')
ax1.set_ylabel('position')

ax1.legend()

ax2.plot(t, v_euler, label='Euler')
ax2.plot(t, v_exact, label='Exact')

ax2.set_xlabel('time')
ax2.set_ylabel('velocity')

ax2.legend()
plt.show()
### ANCHOR_END: exercise_01_c

fig.savefig('../../assets/figures/02-differential_equations/euler_harm_osc.svg')