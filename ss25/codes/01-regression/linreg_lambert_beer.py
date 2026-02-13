#!/usr/bin/env python

### ANCHOR: data_list
concentrations = [
    2.125, 4.250, 6.375, 8.500, 10.63, 12.75, 14.88, 17.00, 19.13, 21.25,
    23.38, 25.50, 27.63, 29.75, 31.88, 34.00, 36.13, 38.25, 40.38, 42.50,
]
absorbances = [
    0.0572, 0.1391, 0.2049, 0.2754, 0.3420, 
    0.4139, 0.4956, 0.5815, 0.6806, 0.7481,
    0.8242, 0.9130, 1.0043, 1.0809, 1.1511,
    1.2483, 1.3373, 1.4027, 1.4927, 1.5853,
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
assert np.isclose(beta0, -0.04907034)
assert np.isclose(beta1, 0.03800109)
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
ax.set_xlabel('concentration / µM')
ax.set_ylabel('absorbance')

# automatically create a legend
ax.legend()
### ANCHOR_END: customize_ax

### ANCHOR: show_plot
plt.show()
### ANCHOR_END: show_plot

### ANCHOR: plot_residuals
residuals = absorbances - (beta0 + beta1 * concentrations)

fig2, ax2 = plt.subplots(figsize=(8, 6))
ax2.bar(concentrations, residuals)
ax2.set_xlabel('concentration / µM')
ax2.set_ylabel('absorbance residuals')
plt.show()
### ANCHOR_END: plot_residuals

fig.savefig('../../assets/figures/01-regression/linreg_lambert_beer.svg')
fig2.savefig('../../assets/figures/01-regression/linreg_lambert_beer_residuals.svg')

