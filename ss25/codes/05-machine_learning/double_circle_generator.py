#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

R1 = 1.5
R2 = 6
SIGMA1 = 0.5
SIGMA2 = 0.8


def generate_circle(radius, sigma, num_points):
    thetas = np.random.uniform(0, 2 * np.pi, num_points)
    radii = np.random.normal(radius, sigma, num_points)
    xs = radii * np.cos(thetas)
    ys = radii * np.sin(thetas)
    return np.array([xs, ys]).T


np.random.seed(0)

circle1 = generate_circle(R1, SIGMA1, 100)
circle2 = generate_circle(R2, SIGMA2, 100)
labels = np.array([1] * 100 + [-1] * 100)

np.savetxt(
    'circles.csv',
    np.hstack([np.vstack([circle1, circle2]), labels[:, None]]),
    delimiter=';',
    comments='',
    header='x_1;x_2;label',
)

fig, ax = plt.subplots(figsize=(4, 4))

ax.scatter(circle1[:, 0], circle1[:, 1], c='r', s=20)
ax.scatter(circle2[:, 0], circle2[:, 1], c='b', s=20)

ax.set_xlim(-8, 8)
ax.set_ylim(-8, 8)

fig.tight_layout()

plt.show()

fig.savefig('../../assets/figures/05-machine_learning/double_circle.svg')

