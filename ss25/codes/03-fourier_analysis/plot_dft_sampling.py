#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import os

X_MIN = 0.0
X_MAX = 40.0
DX = 1.0
NX = 1000

NOMEGA = 50
OMEGA_MIN = 0.0
OMEGA_MAX = 2.0 * np.pi * (NOMEGA - 1) / NOMEGA

TARGET_DIR = '../../assets/figures/03-fourier_analysis/dft_sampling_figures'

omegas = np.linspace(OMEGA_MIN, OMEGA_MAX, NOMEGA)
omegas = np.append(omegas, omegas[1:-1][::-1])
ks = np.arange(NOMEGA)
ks = np.append(ks, ks[1:-1][::-1])

x_discrete = np.arange(X_MIN, X_MAX + DX, DX)
y_discrete = np.sin(omegas[0] * x_discrete)

x_plot = np.linspace(X_MIN, X_MAX, NX)
y_plot = np.sin(omegas[0] * x_plot)


fig, ax = plt.subplots(figsize=(8, 4))
l_plot, = ax.plot(x_plot, y_plot, lw=0.5, label='continuous')
l_discrete, = ax.plot(x_discrete, y_discrete, 'o', label='discretised')
ax.set_xlabel('$t$')
ax.set_ylabel('$f(t)$')
ax.set_title(f'$k = {ks[0]}$')
ax.legend(loc='upper right')

ax.set_xlim([X_MIN, X_MAX])
ax.set_ylim([-1.1, 1.1])

fig.tight_layout()


os.system(f'rm -rf {TARGET_DIR}')
os.system(f'mkdir -p {TARGET_DIR}')
for i, omega in enumerate(omegas):
    y_discrete = np.sin(omega * x_discrete)
    y_plot = np.sin(omega * x_plot)
    l_discrete.set_ydata(y_discrete)
    l_plot.set_ydata(y_plot)
    ax.set_title(f'$k = {ks[i]}$')
    fig.savefig(f'{TARGET_DIR}/{i:04d}.png', dpi=150)

os.system(f'ffmpeg -y -i {TARGET_DIR}/%04d.png '
          f'-vf palettegen {TARGET_DIR}/palette.png')
os.system(f'ffmpeg -framerate 5 -y -i {TARGET_DIR}/%04d.png -i {TARGET_DIR}/palette.png '
          f'-lavfi paletteuse {TARGET_DIR}/dft_sampling.gif')


plt.show()

