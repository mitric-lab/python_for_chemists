#!/usr/bin/env python

### ANCHOR: imports
import numpy as np
import matplotlib.pyplot as plt

XSTEP = 0.3165e-6  # m
### ANCHOR_END: imports

### ANCHOR: import_data
bg_dx, bg_int = np.loadtxt('ir_bg.txt', unpack=True)
spl_dx, spl_int = np.loadtxt('ir_spl.txt', unpack=True)

assert np.allclose(bg_dx, spl_dx)
### ANCHOR_END: import_data

### ANCHOR: plot_interferograms
fig1, ax1 = plt.subplots(figsize=(8, 4))

ax1.plot(bg_dx, bg_int, label='background')
ax1.plot(spl_dx, spl_int, label='sample')

ax1.set_xlabel('relative shift / step')
ax1.set_ylabel('intensity / arb. u.')
ax1.set_xlim(-250, 250)

ax1.legend()
fig1.tight_layout()

plt.show()
### ANCHOR_END: plot_interferograms

#fig1.savefig('../../assets/figures/03-fourier_analysis/ir_interferograms.svg')

### ANCHOR: dft_signal
int_x = spl_int - bg_int

nx = len(spl_dx)
int_nu = np.zeros(nx, dtype=complex)
n_array = np.arange(nx)
for k in range(0, nx):
    int_nu[k] = np.sum(int_x * np.exp(-1j * 2 * np.pi * k / nx * n_array))
### ANCHOR_END: dft_signal

### ANCHOR: dft_freq
x_grid = bg_dx * XSTEP
dx = x_grid[1] - x_grid[0]
if nx % 2 == 0:
    nu_pos = (1.0 / (nx * dx)) * np.arange(0, nx // 2)
    nu_neg = (1.0 / (nx * dx)) * np.arange(-nx // 2, 0)
else:
    nu_pos = (1.0 / (nx * dx)) * np.arange(0, (nx - 1) // 2 + 1)
    nu_neg = (1.0 / (nx * dx)) * np.arange(-(nx - 1) // 2, 0)
nu = np.concatenate((nu_pos, nu_neg))
### ANCHOR_END: dft_freq

### ANCHOR: dft_shift
sort_idx = np.argsort(nu)
nu = nu[sort_idx]
int_nu = int_nu[sort_idx]
### ANCHOR_END: dft_shift

### ANCHOR: signal_shift
int_nu *= np.exp(-1j * 2 * np.pi * nu * x_grid[0])
assert np.max(np.abs(int_nu.imag)) / np.max(np.abs(int_nu)) < 0.05
### ANCHOR_END: signal_shift

### ANCHOR: nu_conversion
nu /= 100  # cm^-1
### ANCHOR_END: nu_conversion

### ANCHOR: plot_spectrum
fig2, ax2 = plt.subplots(figsize=(8, 4))

ax2.plot(nu, np.abs(int_nu.real))

ax2.set_xlabel('wavenumber / cm⁻¹')
ax2.set_ylabel('intensity / arb. u.')
ax2.set_xlim(0, 4000)

fig2.tight_layout()

plt.show()
### ANCHOR_END: plot_spectrum

#np.savetxt('ir_spec.txt', np.column_stack((nu, np.abs(int_nu.real))))
#fig2.savefig('../../assets/figures/03-fourier_analysis/ir_spectrum.svg')

### ANCHOR: numpy_fft
int_nu_np = np.fft.fft(int_x)
nu_np = np.fft.fftfreq(nx, d=dx)
### ANCHOR_END: numpy_fft

### ANCHOR: numpy_fftshift
int_nu_np = np.fft.fftshift(int_nu_np)
nu_np = np.fft.fftshift(nu_np)
### ANCHOR_END: numpy_fftshift

### ANCHOR: numpy_verification
int_nu_np *= np.exp(-1j * 2 * np.pi * nu_np * x_grid[0])
nu_np /= 100  # cm^-1

assert np.allclose(int_nu, int_nu_np)
assert np.allclose(nu, nu_np)
### ANCHOR_END: numpy_verification

