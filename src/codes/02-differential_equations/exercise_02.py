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

def generate_d2(n, h=1.0):
    d2 = np.zeros((n, n))
    rows, cols = np.diag_indices(n)
    d2[rows, cols] = -2
    d2[rows[:-1], cols[1:]] = 1
    d2[rows[1:], cols[:-1]] = 1
    return d2 / h**2

### ANCHOR: exercise_03_a
# Get the potential for a particle in a box with arbitrary potential
def get_potential(n, x, L, a=0, v_a=0, b=None, v_b=None, v_inf=1e4):
    # Initialize potential array
    v = np.zeros(n)
    
    # Add infinite walls
    v[x <= 0] = v_inf
    v[x >= L] = v_inf
    
    # If b and v_b are provided, it has two steps
    if b is not None and v_b is not None:
        v[(x > a) & (x < b)] = v_a
        v[(x > b) & (x < L)] = v_b
    # If only v_a is non-zero, it has one step
    elif v_a != 0:
        v[(x > a) & (x < L)] = v_a
    # Otherwise it's just a simple box with no potential
    
    return v

# Build the Hamiltonian matrix for a particle in a box with arbitrary potential
def build_hamiltonian(n, x, L, a=0, v_a=0, b=None, v_b=None, v_inf=1e4):
    h = x[1] - x[0]
    d2 = generate_d2(n, h)
    v = get_potential(n, x, L, a, v_a, b, v_b, v_inf)
    return -0.5 * d2 + np.diag(v)

# Define the parameters
L = 4.0
NX = 512
X_ARRAY = np.linspace(-0.5, 4.5, NX)

# Build the Hamiltonian matrix and solve the eigenvalue problem
hamiltonian = build_hamiltonian(NX, X_ARRAY, L)
assert np.allclose(hamiltonian, hamiltonian.T)
e, v = np.linalg.eigh(hamiltonian)

# Define the number of states to plot and the potential
NSTATES = 20
eigenenergies = e[:NSTATES]
eigenfunctions = v[:, :NSTATES] / np.sqrt(X_ARRAY[1] - X_ARRAY[0])
potential = get_potential(NX, X_ARRAY, L)

# Plot the results
fig, axs = plt.subplots(1, 2, figsize=(8, 4))
axs[0].plot(np.arange(NSTATES), eigenenergies, 'o', 
            label='numerical eigenenergies')
axs[0].plot(np.arange(NSTATES), (np.pi**2 * np.arange(NSTATES)**2) / (2 * L**2),
            label='analytical eigenenergies')
axs[0].set_xlabel('state')
axs[0].set_ylabel('energy')
axs[0].legend()
axs[1].plot(X_ARRAY, potential, color='k', lw=2, label='potential')
for i in range(5):
    axs[1].plot(X_ARRAY, eigenfunctions[:, i] + eigenenergies[i], 
                label=f'state {i}')
axs[1].set_xlabel('x')
axs[1].set_ylabel('energy')
axs[1].set_xlim(X_ARRAY[0], X_ARRAY[-1])
axs[1].set_ylim(-0.5, 9.5)
axs[1].legend()

fig.tight_layout()
plt.show()
### ANCHOR_END: exercise_03_a

fig.savefig('../../assets/figures/02-differential_equations/fdm_particle_in_box.svg')

### ANCHOR: exercise_03_b
# Define the parameters
L = 4.0
A = 2.0
V_A = 10.0
NX = 512
X_ARRAY = np.linspace(-0.5, 4.5, NX)

# Build the Hamiltonian matrix and solve the eigenvalue problem
hamiltonian = build_hamiltonian(NX, X_ARRAY, L, A, V_A)
assert np.allclose(hamiltonian, hamiltonian.T)
e, v = np.linalg.eigh(hamiltonian)

# Define the number of states to plot and the potential
NSTATES = 20
eigenenergies = e[:NSTATES]
eigenfunctions = v[:, :NSTATES] / np.sqrt(X_ARRAY[1] - X_ARRAY[0])
potential = get_potential(NX, X_ARRAY, L, A, V_A)

# Plot the results
fig, ax = plt.subplots(figsize=(4, 4))
ax.plot(X_ARRAY, potential, color='k', lw=2, label='potential')
for i in range(5):
    ax.plot(X_ARRAY, eigenfunctions[:, i] + eigenenergies[i], 
                label=f'state {i}')
ax.set_xlabel('x')
ax.set_ylabel('energy')
ax.set_xlim(X_ARRAY[0], X_ARRAY[-1])
ax.set_ylim(-0.5, 17.5)
ax.legend()

fig.tight_layout()
plt.show()
### ANCHOR_END: exercise_03_b

fig.savefig('../../assets/figures/02-differential_equations/fdm_particle_in_box_step.svg')

### ANCHOR: exercise_03_c
# Define the parameters
L = 4.0
A = 1.5
B = 2.5
V_A = 10.0
V_B = 2.0
NX = 512
X_ARRAY = np.linspace(-0.5, 4.5, NX)

# Build the Hamiltonian matrix and solve the eigenvalue problem
hamiltonian = build_hamiltonian(NX, X_ARRAY, L, A, V_A, B, V_B)
assert np.allclose(hamiltonian, hamiltonian.T)
e, v = np.linalg.eigh(hamiltonian)

# Define the number of states to plot and the potential
NSTATES = 20
eigenenergies = e[:NSTATES]
eigenfunctions = v[:, :NSTATES] / np.sqrt(X_ARRAY[1] - X_ARRAY[0])
potential = get_potential(NX, X_ARRAY, L, A, V_A, B, V_B)

# Plot the results
fig, ax = plt.subplots(figsize=(4, 4))
ax.plot(X_ARRAY, potential, color='k', lw=2, label='potential')
for i in range(5):
    ax.plot(X_ARRAY, eigenfunctions[:, i] + eigenenergies[i], 
                label=f'state {i}')
ax.set_xlabel('x')
ax.set_ylabel('energy')
ax.set_xlim(X_ARRAY[0], X_ARRAY[-1])
ax.set_ylim(-0.5, 15.0)
ax.legend()

fig.tight_layout()
plt.show()
### ANCHOR_END: exercise_03_c

fig.savefig('../../assets/figures/02-differential_equations/fdm_particle_in_box_double_well.svg')
