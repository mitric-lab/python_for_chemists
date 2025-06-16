#!/usr/bin/env python

### ANCHOR: imports
import numpy as np
import matplotlib.pyplot as plt
### ANCHOR_END: imports

### ANCHOR: load_dataset
trir = np.load('trir.npy')
assert trir.shape == (361, 501)
### ANCHOR_END: load_dataset

### ANCHOR: create_grids
t_grid = np.linspace(0.0, 1800.0, 361)
nu_grid = np.arange(1500, 3500 + 1, 4.0)

### ANCHOR: plot_slices
fig1, ax1 = plt.subplots(figsize=(6, 3))
ax1.plot(nu_grid, trir[0], label='t = 0 s')
ax1.plot(nu_grid, trir[-1], label='t = 1800 s')
ax1.set_xlabel('wavenumber / cm^-1')
ax1.set_ylabel('intensity / arb. units')
ax1.legend()
fig1.tight_layout()
### ANCHOR_END: plot_slices

fig1.savefig('../../assets/figures/04-evd_and_svd/trir_slices.svg')

### ANCHOR: plot_trir
fig2, ax2 = plt.subplots(figsize=(6, 4))
im = ax2.imshow(
    trir, aspect='auto', origin='lower', cmap='viridis',
    extent=(nu_grid[0], nu_grid[-1], t_grid[0], t_grid[-1]),
)

ax2.set_xlabel('wavenumber / cm^-1')
ax2.set_ylabel('time / s')
fig2.colorbar(im, ax=ax2, label='intensity / arb. units')
fig2.tight_layout()
### ANCHOR_END: plot_trir

fig2.savefig('../../assets/figures/04-evd_and_svd/trir.svg')

### ANCHOR: elbow_plot
u, s, vh = np.linalg.svd(trir, full_matrices=False)

fig3, ax3 = plt.subplots(figsize=(6, 3))
ax3.plot(s[:10], marker='o')
ax3.set_xlabel('index')
ax3.set_ylabel('singular value')
ax3.set_yscale('log')
fig3.tight_layout()
### ANCHOR_END: elbow_plot

fig3.savefig('../../assets/figures/04-evd_and_svd/trir_elbow.svg')

### ANCHOR: plot_singular_vectors
fig4, ax4 = plt.subplots(figsize=(6, 3))
ax4.plot(nu_grid, vh[0], label='v_1')
ax4.plot(nu_grid, vh[1], label='v_2')
ax4.set_xlabel('wavenumber / cm^-1')
ax4.set_ylabel('intensity / arb. units')
ax4.legend()
fig4.tight_layout()
### ANCHOR_END: plot_singular_vectors

fig4.savefig('../../assets/figures/04-evd_and_svd/trir_rsv.svg')

plt.show()
