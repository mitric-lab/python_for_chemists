#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import distance_matrix

NPOINTS = 100
SEED = 42

np.random.seed(SEED)
points = np.random.rand(NPOINTS, 3)

distances = distance_matrix(points, points)

c_mat = np.eye(NPOINTS) - np.ones((NPOINTS, NPOINTS)) / NPOINTS
b_mat = -0.5 * np.linalg.multi_dot([c_mat, distances**2, c_mat])
e, v = np.linalg.eigh(b_mat)
embedding = np.dot(v[:, -2:], np.diag(np.sqrt(e[-2:])))

d2 = distances
d2 = d2 - np.mean(d2, axis=0)
u, s, vt = np.linalg.svd(d2)
embedding2 = np.dot(u[:, :2], np.diag(np.sqrt(s[:2])))

fig, axs = plt.subplots(1, 2, figsize=(10, 5))

for ax in axs:
    ax.set_aspect('equal')

axs[0].scatter(embedding[:, 0], embedding[:, 1])
axs[1].scatter(embedding2[:, 0], embedding2[:, 1])

plt.show()

