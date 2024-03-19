#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

EPSILON_ARRAY = np.array(
    [1000.0, 800.0, 600.0, 400.0, 200.0, 100.0]
)
THICKNESS = 1.0
NOISE_LEVEL = 0.002
SEED = 42

concentrations = np.array(
    [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 
     1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]
) * 1e-3



def get_absorbance(concentration, epsilons, thickness):
    intensities_in = np.ones_like(epsilons)
    intensities_out = intensities_in \
        * np.exp(-epsilons * concentration * thickness)
    return np.log(np.sum(intensities_in) / np.sum(intensities_out))


def add_noise(intensity, noise_level):
    return intensity \
        + np.random.normal(scale=noise_level, size=len(intensity))


np.random.seed(SEED)
absorbances = np.array([
    get_absorbance(c, EPSILON_ARRAY, THICKNESS) for c in concentrations
])
absorbances = add_noise(absorbances, NOISE_LEVEL)

print(concentrations)
print(absorbances)

fig, ax = plt.subplots()
ax.plot(concentrations, absorbances)
plt.show()

