#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

AC2O_PEAKS = [(1827, 10, 0.95), (1766, 10, 0.92)]
ACOH_PEAKS = [(1715, 10, 0.96), (3000, 100, 0.30)]
K = 0.002631  # s-1
NU_MIN, NU_MAX = 1500, 3500  # cm-1
NU_POINTS = 501
T_MIN, T_MAX = 0, 1800
T_POINTS = 361


def gaussian(x, mu, sigma):
    return np.exp(-(x - mu) ** 2 / (2 * sigma ** 2))


t_grid = np.linspace(T_MIN, T_MAX, T_POINTS)
nu_grid = np.linspace(NU_MIN, NU_MAX, NU_POINTS)
trir = np.zeros((T_POINTS, NU_POINTS))

for i, t in enumerate(t_grid):
    ca = np.exp(-K * t)  # [Ac2O]
    cp = 1 - ca  # [AcOH]
    for nu0, sigma, eps in AC2O_PEAKS:
        trir[i, :] += ca * eps * gaussian(nu_grid, nu0, sigma)
    for nu0, sigma, eps in ACOH_PEAKS:
        trir[i, :] += cp * eps * gaussian(nu_grid, nu0, sigma)

rng = np.random.default_rng(seed=42)

# 1) White (detector) noise
sigma_white = 5.0e-3
trir += rng.normal(0, sigma_white, size=trir.shape)

# 2) Slow baseline drift
t_norm = np.linspace(-1, 1, trir.shape[0])
baseline = 3e-3 * t_norm[:, np.newaxis] + 1e-3 * (t_norm**2)[:, np.newaxis]
trir += baseline

# 3) Occasional impulsive “glitches”
n_glitch = int(0.001 * trir.size)
rows = rng.integers(0, t_grid.size, n_glitch)
cols = rng.integers(0, nu_grid.size, n_glitch)
trir[rows, cols] += rng.normal(0, 5e-2, n_glitch)


np.save('trir.npy', trir)


u, s, vh = np.linalg.svd(trir, full_matrices=False)
print(s[:100])  # Print the first 10 singular values

fig1, ax1 = plt.subplots(1, 1, figsize=(6, 6))
ax1.imshow(trir, aspect='auto', extent=(NU_MIN, NU_MAX, T_MAX, T_MIN))
ax1.set_xlabel('Wavenumber / cm$^{-1}$')
ax1.set_ylabel('Time / s')
cbar = plt.colorbar(ax1.images[0], ax=ax1)


fig2, axs2 = plt.subplots(2, 1, figsize=(6, 3))
axs2[0].plot(nu_grid, trir[0, :], label='t=0')
axs2[0].plot(nu_grid, trir[-1, :], label='t=1800')

axs2[1].plot(nu_grid, vh[0], label='1st SVD mode')
axs2[1].plot(nu_grid, vh[1], label='2nd SVD mode')
axs2[1].set_xlabel('Wavenumber / cm$^{-1}$')
axs2[1].set_ylabel('SVD mode amplitude')


plt.show()
