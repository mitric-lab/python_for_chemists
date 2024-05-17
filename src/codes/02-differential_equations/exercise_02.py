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

### ANCHOR: exercise_03_a
from scipy.integrate import solve_ivp

# Define differential equation for particle in a box
def dfdx(x, f, E):
    psi, phi = f
    dpsi_dx = phi
    dphi_dx = -2 * E * psi
    return [dpsi_dx, dphi_dx]

# Solve Schrodinger equation with shooting method
def solve_for_energy(num_states, L, step_size=0.01, tolerance=0.01):
    wavefunctions, energies = [], []
    E = 0
    while len(energies) < num_states:

        # Solve Schrodinger equation
        sol = solve_ivp(dfdx, [0, L], [0, 1.0], args=(E,), dense_output=True, t_eval=np.linspace(0, L, 1000))
        psi = sol.y[0]
        error = psi[-1]

        # Check if wavefunction is zero at boundary
        if np.abs(error) < tolerance:
            energies.append(E)
            wavefunctions.append(psi)
            E += 0.5
        else:
            E += step_size

    return np.array(energies), np.array(wavefunctions)

# Parameters for the problem
L = 4.0
num_states = 5
step_size = 0.0005
tolerance = 0.002

# Solve Schroedinger equation with shooting method
E, Psi = solve_for_energy(num_states, L, step_size, tolerance)

# Normalize wavefunctions 
dx = L / len(Psi[0])
for i in range(len(Psi)):
    Psi[i] /= np.sqrt(np.sum(Psi[i]**2) * dx)

# Plot energy levels and wavefunctions
fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(8, 4))

ax1.plot(np.arange(num_states), E, 'o', label='numerical eigenenergies')
ax1.plot(np.arange(num_states), np.arange(1, num_states + 1)**2 * np.pi**2 / (2 * L**2), label='analytical eigenenergies')
ax1.set_xlabel('state')
ax1.set_ylabel('energy')
ax1.legend()

for i in range(num_states):
    ax2.plot(np.linspace(0, L, len(Psi[i])), Psi[i] + E[i], label=f'state {i}')
ax2.set_xlabel('x')
ax2.set_ylabel('energy')
ax2.legend()

plt.tight_layout()
plt.show()
### ANCHOR_END: exercise_03_a

fig.savefig('../../assets/figures/02-differential_equations/shooting_method_pib.svg')

### ANCHOR: exercise_03_b
# Define differential equation for particle in a box with step potential
def dfdx(x, f, E, V0, a):
    psi, phi = f
    V = V0 if x > a else 0
    dpsi_dx = phi
    dphi_dx = 2 * (V - E) * psi
    return [dpsi_dx, dphi_dx]

# Solve Schrodinger equation with shooting method
def solve_for_energy(num_states, L, v0, a, step_size=0.01, tolerance=0.01):
    wavefunctions, energies = [], []
    E = 0
    while len(energies) < num_states:
        # Solve Schrodinger equation
        sol = solve_ivp(dfdx, [0, L], [0, 1.0], args=(E, V0, a), dense_output=True, t_eval=np.linspace(0, L, 1000), rtol=1e-6, atol=1e-8)
        psi = sol.y[0]
        error = psi[-1]

        # Check if wavefunction is zero at boundary
        if np.abs(error) < tolerance:
            energies.append(E)
            wavefunctions.append(psi)
            E += 0.5
        else:
            E += step_size

    return np.array(energies), np.array(wavefunctions)

# Parameters for the problem
L = 4.0
V0 = 10.0
a = 2.0
num_states = 5
step_size = 0.0005 # works with 0.0002
tolerance = 0.002 # works with 0.05

# Solve Schroedinger equation with shooting method
E, Psi = solve_for_energy(num_states, L, V0, a, step_size, tolerance)

# Normalize wavefunctions 
dx = L / len(Psi[0])
for i in range(len(Psi)):
    Psi[i] /= np.sqrt(np.sum(Psi[i]**2) * dx)

# Plot energy levels and wavefunctions
fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(8, 4))

ax1.plot(np.arange(num_states), E, 'o', label='numerical eigenenergies')
ax1.set_xlabel('state')
ax1.set_ylabel('energy')
ax1.legend()

for i in range(num_states):
    ax2.plot(np.linspace(0, L, len(Psi[i])), Psi[i] + E[i], label=f'state {i}')
ax2.plot([0, a, a, L], [0, 0, V0, V0], 'k', label='potential')
ax2.set_xlabel('x')
ax2.set_ylabel('energy')
ax2.legend()

plt.tight_layout()
plt.show()
### ANCHOR_END: exercise_03_b

fig.savefig('../../assets/figures/02-differential_equations/shooting_method_pib_step.svg')

