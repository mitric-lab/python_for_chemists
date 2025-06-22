#!/usr/bin/env python

### ANCHOR: imports
import numpy as np
import matplotlib.pyplot as plt
### ANCHOR_END: imports

### ANCHOR: load_xyz
GEOM_XYZ = 'tetracene.xyz'
symbols = []
coords = []
with open(GEOM_XYZ, 'r') as f:
    for line in f:
        tmp = line.split()
        if len(tmp) > 3:
            symbols.append(tmp[0])
            coords.append([float(x) for x in tmp[1:4]])
symbols = np.array(symbols)
coords = np.array(coords)
### ANCHOR_END: load_xyz

### ANCHOR: pca
centroid = np.mean(coords, axis=0)
coords -= centroid

u, s, vh = np.linalg.svd(coords, full_matrices=False)
### ANCHOR_END: pca

### ANCHOR: expl_var
expl_var = s**2 / np.sum(s**2)

print('Explained variance:', expl_var)
assert np.allclose(expl_var, [8.30213895e-01, 1.69730510e-01, 5.55956455e-05])
### ANCHOR_END: expl_var

### ANCHOR: plot_projection
fig, ax = plt.subplots(figsize=(6, 3))
proj = u @ np.diag(s)
ax.scatter(proj[:,0], proj[:,1], s=60, marker='o')

ax.set_aspect('equal')
ax.set_xlabel(f'PC1 ({expl_var[0]*100:.1f}%)')
ax.set_ylabel(f'PC2 ({expl_var[1]*100:.1f}%)')

fig.tight_layout()
plt.show()
### ANCHOR_END: plot_projection

fig.savefig('../../assets/figures/04-evd_and_svd/tetracene_projection.svg')

# plot atoms using 3D scatter
ATOM_PROPS = {
    'H': (50.0, 'silver'),
    'C': (100.0, 'forestgreen'),
}
XLIM = (-6, 4)
YLIM = (-7.5, 2.5)
ZLIM = (-4.2, 5.8)
VIEW = {'elev': 60, 'azim': 10, 'roll': 0}
PC_ARROW_SCALE = 2.0
PC_ARROW_COLORS = ('red', 'cyan', 'purple')

marker_sizes = np.array([ATOM_PROPS[s][0] for s in symbols])
marker_colors = np.array([ATOM_PROPS[s][1] for s in symbols])
coords += centroid  # shift back to original coordinates

fig2 = plt.figure(figsize=(6, 6))
ax2 = fig2.add_subplot(111, projection='3d')
ax2.set_aspect('equal')
ax2.view_init(**VIEW)
ax2.scatter(
    coords[:,0], coords[:,1], coords[:,2], 
    s=marker_sizes, c=marker_colors, marker='o',
)
for ci in coords:
    for cj in coords:
        rij = np.linalg.norm(ci - cj)
        if rij < 1.5:
            ax2.plot(
                [ci[0], cj[0]], [ci[1], cj[1]], [ci[2], cj[2]],
                color='tab:blue', linestyle='-', linewidth=2.5,
            )
ax2.set_xlim(XLIM)
ax2.set_ylim(YLIM)
ax2.set_zlim(ZLIM)
ax2.set_xlabel('$x$')
ax2.set_ylabel('$y$')
ax2.set_zlabel('$z$')

for i, (v, c) in enumerate(zip(vh, PC_ARROW_COLORS)):
    x, y, z = centroid
    u, v, w = v
    ax2.quiver(
        x, y, z, u, v, w, color=c, 
        normalize=True, length=PC_ARROW_SCALE, linewidths=4.0,
        label=f'PC{i+1} ({expl_var[i]*100:.1f}%)',
    )

ax2.legend()
fig2.tight_layout()

plt.show()

fig2.savefig('../../assets/figures/04-evd_and_svd/tetracene_3d.svg')

