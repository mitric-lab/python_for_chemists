#!/usr/bin/env python

### ANCHOR: imports
import numpy as np
### ANCHOR_END: imports

### ANCHOR: water_linear_system
a_mat = np.array([
    [1.0, 1.0, 1.0],  # monopole term
    [0.0000, 0.5803, 0.5803],  # dipole z-term
])
b_vec = np.array([
     0.0,  # monopole term
    -0.434331,  # dipole z-term
])
### ANCHOR_END: water_linear_system

### ANCHOR: water_solution
a_pinv = np.linalg.pinv(a_mat, rcond=1e-12)
atom_charges = a_pinv @ b_vec

print(atom_charges)
assert np.allclose(atom_charges, [0.74845942, -0.37422971, -0.37422971])
### ANCHOR_END: water_solution

### ANCHOR: load_ch3cl_xyz
GEOM_XYZ = 'ch3cl.xyz'
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
### ANCHOR_END: load_ch3cl_xyz

### ANCHOR: ch3cl_linear_system
DIPOLE = np.array([
    -0.00000000,  # x
     0.00000000,  # y
    -0.40671888,  # z
])  # e * Angstrom

natoms = len(symbols)
a_mat = np.vstack([
    np.ones(natoms),  # monopole term
    coords[:, 0],  # dipole x-term
    coords[:, 1],  # dipole y-term
    coords[:, 2],  # dipole z-term
])
b_vec = np.concatenate(([0.0], DIPOLE))

assert a_mat.shape == (4, 5)
### ANCHOR_END: ch3cl_linear_system

### ANCHOR: ch3cl_solution
a_pinv = np.linalg.pinv(a_mat, rcond=1e-12)
atom_charges = a_pinv @ b_vec

print(symbols)
print(atom_charges)
assert np.allclose(
    atom_charges, 
    [0.01869843, -0.20764081, 0.06298079, 0.06298079, 0.06298079],
)
### ANCHOR_END: ch3cl_solution

