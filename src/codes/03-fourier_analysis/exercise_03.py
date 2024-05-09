#!/usr/bin/env python

### ANCHOR: exercise_a
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import c

# Import dipole data
data = np.loadtxt("dipoles.txt")
time = data[:, 0]
dipoles = data[:, 1:]

# Calculate the frequency
dt = (time[1] - time[0]) * 1e3  # fs
freq = np.fft.fftfreq(len(time), dt)
freq = np.fft.fftshift(freq)
wn = (freq * 1e15 / c) / 100.0  # cm^{-1}

# Perform Fourier transform
dipole_ft = np.array([np.fft.fftshift(np.fft.fft(d)) for d in dipoles.T])
ac_ift = np.linalg.norm(dipole_ft, axis=0)**2
spec = wn**2 * ac_ift

# Normalize the spectrum
spec /= np.max(spec)

# Import experimental spectrum
data_exp = np.loadtxt("ir_spec.txt")
wn_exp = data_exp[:, 0]
spec_exp = data_exp[:, 1]

# Normalize the experimental spectrum
spec_exp /= np.max(spec_exp)

# Plot the spectra
fig, ax = plt.subplots(figsize=(8, 4))

ax.plot(wn_exp, spec_exp, label='Experimental')
ax.plot(wn, spec, label='Calculated')

ax.set_xlabel('wavenumber / cm⁻¹')
ax.set_ylabel('intensity / arb. u.')
ax.set_xlim(0, 4000)
ax.legend()

fig.tight_layout()

plt.show()
### ANCHOR_END: exercise_a

#fig.savefig('../../assets/figures/03-fourier_analysis/ir_spectrum_compare.svg')

