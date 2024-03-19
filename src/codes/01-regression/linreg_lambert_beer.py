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

concentrations = np.array(concentrations)
absorbances = np.array(absorbances)
### ANCHOR_END: data_array

### ANCHOR: create_system_array_elements
# make sure the number of data points is the same
assert len(concentrations) == len(absorbances)

number_of_points = len(concentrations)
sum_x = np.sum(concentrations)
sum_x_sq = np.sum(concentrations**2)
sum_y = np.sum(absorbances)
sum_xy = np.sum(concentrations * absorbances)
### ANCHOR_END: create_system_array_elements

### ANCHOR: create_system_arrays
a_arr = np.array([
    [number_of_points, sum_x],
    [sum_x, sum_x_sq],
])
b_arr = np.array([sum_y, sum_xy])
### ANCHOR_END: create_system_arrays

### ANCHOR: solve_system
beta = np.linalg.solve(a_arr, b_arr)
print(beta)
### ANCHOR_END: solve_system

### ANCHOR: beta_verification
beta0 = beta[0]
beta1 = beta[1]
assert np.isclose(beta0, 0.03735158)
assert np.isclose(beta1, 0.41501278)
### ANCHOR_END: beta_verification

### ANCHOR: import_mpl
import matplotlib.pyplot as plt
### ANCHOR_END: import_mpl

### ANCHOR: make_axes
fig, ax = plt.subplots(figsize=(8, 6))
### ANCHOR_END: make_axes

### ANCHOR: plot_data
ax.plot(concentrations, absorbances, 'o', label='data')
ax.plot(concentrations, beta0 + beta1 * concentrations, label='fit')
### ANCHOR_END: plot_data

### ANCHOR: customize_ax
ax.set_xlabel('concentration / mM')
ax.set_ylabel('absorbance')

# automatically create a legend
ax.legend()
### ANCHOR_END: customize_ax

### ANCHOR: show_plot
plt.show()
### ANCHOR_END: show_plot

fig.savefig('../../assets/figures/01-regression/linreg_lambert_beer.svg')

