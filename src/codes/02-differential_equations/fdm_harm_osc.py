#!/usr/bin/env python

### ANCHOR: imports
import numpy as np
import matplotlib.pyplot as plt
### ANCHOR_END: imports


### ANCHOR: generate_d2_naive
def generate_d2_naive(n: int, h: float = 1.0) -> np.ndarray:
    d2 = np.zeros((n, n))
    for i in range(0, n):
        for j in range(0, n):
            if i == j:
                d2[i, j] = -2
            elif i == j - 1 or i == j + 1:
                d2[i, j] = 1
    return d2 / h**2
### ANCHOR_END: generate_d2_naive


### ANCHOR: generate_d2
def generate_d2(n: int, h: float = 1.0) -> np.ndarray:
    d2 = np.zeros((n, n))
    rows, cols = np.diag_indices(n)
    d2[rows, cols] = -2
    d2[rows[:-1], cols[1:]] = 1
    d2[rows[1:], cols[:-1]] = 1
    return d2 / h**2
### ANCHOR_END: generate_d2


import time

start = time.perf_counter()
for i in range(10):
    generate_d2_naive(1000)
end = time.perf_counter()
print("Naive: ", end - start)

start = time.perf_counter()
for i in range(10):
    generate_d2(1000)
end = time.perf_counter()
print("Optimized: ", end - start)

### ANCHOR: build_hamiltonian
def build_hamiltonian(n: int, x: np.ndarray, k: float = 1.0) -> np.ndarray:
    h = x[1] - x[0]
    d2 = generate_d2(n, h)
    vx = np.diag(0.5 * k * x**2)
    l_mat = -0.5 * d2 + vx

    return l_mat
### ANCHOR_END: build_hamiltonian


### ANCHOR: define_parameters
K = 1.0
NX = 512
X_ARRAY = np.linspace(-5, 5, NX)
### ANCHOR_END: define_parameters

### ANCHOR: solve_matrix_equation
hamiltonian = build_hamiltonian(NX, X_ARRAY, K)
assert np.allclose(hamiltonian, hamiltonian.T)
e, v = np.linalg.eigh(hamiltonian)
### ANCHOR_END: solve_matrix_equation

### ANCHOR: solve_harmonic_oscillator
NSTATES = 20
eigenenergies = e[:NSTATES]
eigenfunctions = v[:, :NSTATES] / np.sqrt(X_ARRAY[1] - X_ARRAY[0])
### ANCHOR_END: solve_harmonic_oscillator

### ANCHOR: plot_results_subplots
fig, axs = plt.subplots(1, 2, figsize=(8, 4))
### ANCHOR_END: plot_results_subplots

### ANCHOR: plot_eigenenergies
axs[0].plot(np.arange(NSTATES), eigenenergies, 'o', 
            label='numerical eigenenergies')
axs[0].plot(np.arange(NSTATES), np.sqrt(K) * (np.arange(NSTATES) + 0.5),
            label='analytical eigenenergies')
axs[0].set_xlabel('state')
axs[0].set_ylabel('energy')
axs[0].legend()
### ANCHOR_END: plot_eigenenergies

### ANCHOR: plot_eigenfunctions
axs[1].plot(X_ARRAY, 0.5 * K * X_ARRAY**2, color='k', lw=2, label='potential')
for i in range(5):
    axs[1].plot(X_ARRAY, eigenfunctions[:, i] + eigenenergies[i], 
                label=f'state {i}')
axs[1].set_xlabel('x')
axs[1].set_ylabel('energy')
axs[1].set_xlim(X_ARRAY[0], X_ARRAY[-1])
axs[1].set_ylim(-0.5, 9.5)
axs[1].legend()
### ANCHOR_END: plot_eigenfunctions

### ANCHOR: show_plot
fig.tight_layout()
plt.show()
### ANCHOR_END: show_plot

fig.savefig('../../assets/figures/02-differential_equations/fdm_harm_osc.svg')

