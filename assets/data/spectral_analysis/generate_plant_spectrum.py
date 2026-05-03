#!/usr/bin/env python

import numpy as np

PIGMENTS = {
    # pigment: (column, molar mass)
    'anthocyanin':   (1, 484.84),  # cyanidin-3-glucoside equivalent
    'betacyanin':    (2, 550.47),  # betanin equivalent
    'betaxanthin':   (3, 339.30),  # vulgaxanthin I equivalent
    'betacarotene':  (4, 536.888),
    'chlorophyll a': (5, 893.49),  # combined a and b, scaled to a
    'chlorophyll b': (6, 907.47),  # combined a and b, scaled to b
}
PLANTS = {
    # plant: [pigment concentrations in g l^-1]
    'spinach': [0.00, 0.00, 0.00, 0.05, 0.25, 0.20],
    'carrot': [0.00, 0.00, 0.00, 0.25, 0.00, 0.00],
    'blueberry': [1.00, 0.00, 0.00, 0.00, 0.00, 0.00],
    'cactus pear': [0.00, 0.05, 0.30, 0.00, 0.00, 0.00],
}
LAMBDA_MIN = 300.0
LAMBDA_MAX = 800.0

SPEC_CSV_PATH = 'plant_pigments.csv'

data = np.loadtxt(SPEC_CSV_PATH, skiprows=1, delimiter=',')

wavelengths = data[:, 0]
lambda_min_mask = (wavelengths - LAMBDA_MIN) > -1e-6
lambda_max_mask = (wavelengths - LAMBDA_MAX) < 1e-6
lambda_mask = lambda_min_mask & lambda_max_mask
data = data[lambda_mask, :]

wavelengths = data[:, 0]
pigment_spectra = np.zeros((len(wavelengths), len(PIGMENTS)))
for i, (col, molar_mass) in enumerate(PIGMENTS.values()):
    # convert from epsilon (l mol^-1 cm^-1) to alpha (l g^-1 cm^-1)
    pigment_spectra[:, i] = data[:, col] / molar_mass

plant_spectra = {}
for plant, pigment_concs in PLANTS.items():
    coeffs = np.sum(pigment_spectra * pigment_concs, axis=1)
    plant_spectra[plant] = np.clip(coeffs, a_min=0.0, a_max=None)

np.savetxt(
    'plant_spectra.csv',
    np.column_stack((wavelengths, *plant_spectra.values())),
    fmt='%.8f',
    header='wavelength,' + ','.join(plant_spectra.keys()),
    delimiter=',',
    comments=''
)

import matplotlib.pyplot as plt
for plant, alpha in plant_spectra.items():
    plt.plot(wavelengths, alpha, label=plant)
plt.xlabel('wavelength / nm')
plt.ylabel('absorption coefficient / cm^-1')
plt.legend()
plt.show()
