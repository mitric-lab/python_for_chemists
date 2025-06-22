#!/usr/bin/env python

### ANCHOR: imports
import numpy as np
import matplotlib.pyplot as plt
### ANCHOR_END: imports

### ANCHOR: load_data
data = np.loadtxt('wine.csv', delimiter=',', skiprows=1)
categories = data[:, 0].astype(int)
features = data[:, 3:].astype(float)
### ANCHOR_END: load_data

### AHCNOR: filter_data
mask = categories != 0
categories = categories[mask]
features = features[mask]
n_samples, n_features = features.shape
### ANCHOR_END: filter_data

### ANCHOR: labels
CATEGORY_LABELS = ['Wappenlese Weiß', 'Wappenlese Rot', 'Soave', 'Bardolino']
FEATURE_LABELS = ['V_NO2', 'V_EtOH', 'V_VOC', 'V_CO']
### ANCHOR_END: labels

### ANCHOR: standardise
features -= features.mean(axis=0)
features /= features.std(axis=0)
### ANCHOR_END: standardise

### ANCHOR: svd
u, s, vh = np.linalg.svd(features, full_matrices=False)
pcs = vh.T
proj = u @ np.diag(s)
expl_var = s**2 / np.sum(s**2)
### ANCHOR_END: svd

### ANCHOR: verify_pca
print('First principal component:')
print(pcs[:, 0])
print('Explained variance:')
print(expl_var)

assert np.allclose(
    pcs[:, 0], 
    [0.57725187, 0.59546988, 0.55488817, 0.06553647],
)
assert np.allclose(
    expl_var,
    [0.61840175, 0.27550301, 0.06891797, 0.03717727],
)
### ANCHOR_END: verify_pca

### ANCHOR: plot_variance
fig1, ax1 = plt.subplots(figsize=(8, 6))

ax1.set_xticks(range(1, n_features + 1))
ax1.set_xlabel('principal component index')
ax1.set_ylabel('explained variance')

ax1.bar(range(1, n_features + 1), expl_var, color='tab:blue')
ax1.plot(range(1, n_features + 1), np.cumsum(expl_var), 'o-', c='tab:orange')

fig1.tight_layout()
plt.show()
### ANCHOR_END: plot_variance

fig1.savefig('../../assets/figures/04-evd_and_svd/pca_wine_variance.svg')

### ANCHOR: plot_pca
fig2, ax2 = plt.subplots(figsize=(8, 6))
ax2.set_aspect('equal')
ax2.set_xlabel(f'PC 1 ({expl_var[0]:.2f})')
ax2.set_ylabel(f'PC 2 ({expl_var[1]:.2f})')
scat = ax2.scatter(proj[:, 0], proj[:, 1])

fig2.tight_layout()
plt.show()
### ANCHOR_END: plot_pca

fig2.savefig('../../assets/figures/04-evd_and_svd/pca_wine_projection.svg')

### ANCHOR: plot_pca_coloured
CATEGORY_COLORS = ['#ffffb3', '#bebada', '#8dd3c7', '#fb8072']
colors = [CATEGORY_COLORS[c - 1] for c in categories]
scat.set_color(colors)

for c, l in zip(CATEGORY_COLORS, CATEGORY_LABELS):
    ax2.scatter([], [], c=c, label=l)
ax2.legend()

plt.show()
### ANCHOR_END: plot_pca_coloured

fig2.savefig('../../assets/figures/04-evd_and_svd/pca_wine_projection_coloured.svg')

