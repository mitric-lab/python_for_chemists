#!/usr/bin/env python

### ANCHOR: data_list
concentrations = [
    0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
    1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
]
absorbances = [
    0.0522, 0.1010, 0.1518, 0.2017, 0.2453, 
    0.2920, 0.3404, 0.3831, 0.4240, 0.4685,
    0.5080, 0.5486, 0.5899, 0.6245, 0.6630,
    0.7020, 0.7382, 0.7766, 0.8091, 0.8424,
]
### ANCHOR_END: data_list

### ANCHOR: data_array
import numpy as np
import matplotlib.pyplot as plt

concentrations = np.array(concentrations)
absorbances = np.array(absorbances)
### ANCHOR_END: data_array

### ANCHOR: exercise_01_1
# Umbennenung der Variablen
x = concentrations
y = absorbances

# Berechnung der Mittelwerte
x_mean = np.mean(x)
y_mean = np.mean(y)

# Berechnung der Parameter der linearen Regression
beta_1 = np.sum((x - x_mean) * (y - y_mean)) / np.sum((x - x_mean)**2)
beta_0 = y_mean - beta_1 * x_mean

assert np.isclose(beta_0, 0.03735158)
assert np.isclose(beta_1, 0.41501278)
### ANCHOR_END: exercise_01_1

### ANCHOR: exercise_01_2
# Aufstellen des Gleichungssystems
a_arr = np.array([
    [len(x),        np.sum(x),      np.sum(x**2)],
    [np.sum(x),     np.sum(x**2),   np.sum(x**3)],
    [np.sum(x**2),  np.sum(x**3),   np.sum(x**4)]
])
b_arr = np.array([np.sum(y), np.sum(x * y), np.sum(x**2 * y)])

# Lösen des Gleichungssystems
beta = np.linalg.solve(a_arr, b_arr)
print(beta)

# Extrahieren der Parameter
beta_0, beta_1, beta_2 = beta

# Plotten der Daten
fig, ax = plt.subplots(figsize=(8, 6))

ax.plot(x, y, 'o', label='data')
ax.plot(x, beta_0 + beta_1 * x + beta_2 * x**2, label='fit')

ax.set_xlabel('concentration / mM')
ax.set_ylabel('absorbance')

ax.legend()
plt.show()
### ANCHOR_END: exercise_01_2

fig.savefig('../../assets/figures/01-regression/quadreg_lambert_beer.svg')