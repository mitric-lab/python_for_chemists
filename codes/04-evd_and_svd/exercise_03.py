#!/usr/bin/env python

### ANCHOR: exercise_01_b
import numpy as np
import matplotlib.pyplot as plt
from utils import read_xyz, get_connectivity_matrix, CubeFile

# Read benzene structure, ignore hydrogens and get connectivity matrix
atoms = read_xyz("benzene.xyz", use_bohr=True)
heavy_atoms = [atom for atom in atoms if atom[0] != 1]
connectivity = get_connectivity_matrix(heavy_atoms)
print(connectivity)

# Set up Hückel matrix parameters
alpha = 0.0  # Diagonal elements
beta = -1.0  # Off-diagonal elements for bonded atoms

# Convert connectivity matrix to Hückel matrix
# H_μν = α if μ = ν (diagonal)
# H_μν = β if atoms μ and ν are bonded
# H_μν = 0 otherwise
huckel_matrix = alpha * np.eye(connectivity.shape[0])
print(huckel_matrix)
huckel_matrix[connectivity == 1] = beta
print(huckel_matrix)
### ANCHOR_END: exercise_01_b

### ANCHOR: exercise_01_c
# Calculate eigenvalues and eigenvectors
energies, coefficients = np.linalg.eigh(huckel_matrix)

# Set up 3D grid for orbital calculation
spacing = 0.2  # Bohr
max_extent = 6.0  # Maximum extent in Bohr
npoints = int(2 * max_extent / spacing) + 1

# Create grid
x = np.linspace(-max_extent, max_extent, npoints)
y = np.linspace(-max_extent, max_extent, npoints)
z = np.linspace(-max_extent, max_extent, npoints)
X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

# Function to calculate p_z orbital
def pz_orbital(x, y, z, center):
    x0, y0, z0 = center
    r = np.sqrt((x - x0)**2 + (y - y0)**2 + (z - z0)**2)
    return (z - z0) * np.exp(-r)

# Calculate molecular orbitals
orbitals = []
for orbital_idx in range(len(energies)):
    # Initialize orbital grid
    orbital = np.zeros((npoints, npoints, npoints))
    
    # Sum over atomic contributions (only heavy atoms)
    for atom_idx, (_, x0, y0, z0) in enumerate(heavy_atoms):
        # Calculate p_z orbital for this center
        pz = pz_orbital(X, Y, Z, (x0, y0, z0))
        # Add contribution with appropriate coefficient
        orbital += coefficients[atom_idx, orbital_idx] * pz
    
    orbitals.append(orbital)
### ANCHOR_END: exercise_01_c

### ANCHOR: exercise_01_d
# Save each orbital as a cube file
for i, orbital in enumerate(orbitals):
    cube = CubeFile(atoms=atoms, data=orbital, spacing=spacing)
    cube.dump(f"benzene_orbital_{i+1}.cube", 
             comment=f"Benzene molecular orbital {i+1}, E = {energies[i]:.3f}β")
### ANCHOR_END: exercise_01_d