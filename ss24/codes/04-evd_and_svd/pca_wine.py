#!/usr/bin/env python

### ANCHOR: imports
import numpy as np
import matplotlib.pyplot as plt
### ANCHOR_END: imports

### ANCHOR: load_data
data = np.loadtxt('wine.csv', delimiter=',')
categories = data[:, 0].astype(int) - 1
features = data[:, 1:].astype(float)
### ANCHOR_END: load_data

### ANCHOR: labels
CATEGORY_LABELS = ['Barolo', 'Grignolino', 'Barbera']
FEATURE_LABELS = [
    'alcohol', 'malic acid', 'ash', 'alcalinity of ash', 'magnesium',
    'total phenols', 'flavanoids', 'nonflavanoid phenols',
    'proanthocyanins', 'color intensity', 'hue',
    'OD280/OD315', 'proline',
]
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
    [-0.1443294, 0.24518758, 0.00205106, 0.23932041, -0.14199204, -0.39466085,
     -0.4229343, 0.2985331, -0.31342949, 0.0886167, -0.29671456, -0.37616741,
     -0.28675223],
)
assert np.allclose(
    expl_var,
    [0.36198848, 0.1920749, 0.11123631, 0.0706903, 0.06563294, 0.04935823,
     0.04238679, 0.02680749, 0.02222153, 0.01930019, 0.01736836, 0.01298233,
     0.00795215],
)
### ANCHOR_END: verify_pca

### ANCHOR: plot_variance
fig1, ax1 = plt.subplots(figsize=(8, 6))

ax1.set_xticks(range(1, 14, 2))
ax1.set_xlabel('principal component index')
ax1.set_ylabel('explained variance')

ax1.bar(range(1, 14), expl_var, color='tab:blue')
ax1.plot(range(1, 14), np.cumsum(expl_var), 'o-', c='tab:orange')

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
CATEGORY_COLORS = ['#66c2a5', '#fc8d62', '#8da0cb']
colors = [CATEGORY_COLORS[c] for c in categories]
scat.set_color(colors)

for c, l in zip(CATEGORY_COLORS, CATEGORY_LABELS):
    ax2.scatter([], [], c=c, label=l)
ax2.legend()

plt.show()
### ANCHOR_END: plot_pca_coloured

fig2.savefig('../../assets/figures/04-evd_and_svd/pca_wine_projection_coloured.svg')

