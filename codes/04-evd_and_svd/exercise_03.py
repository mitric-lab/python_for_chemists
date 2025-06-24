#!/usr/bin/env python

### ANCHOR: exercise_01_b
import numpy as np
from utils import read_xyz, get_connectivity_matrix

# Read benzene structure and get connectivity matrix for heavy atoms
atoms = read_xyz("benzene.xyz")
heavy_atoms = [atom for atom in atoms if atom[0] != 1]
connectivity = get_connectivity_matrix(heavy_atoms)

# Set up Hückel matrix parameters
alpha = 0.0  # Diagonal elements
beta = -1.0  # Off-diagonal elements for bonded atoms

# Convert connectivity matrix to Hückel matrix
# H_μν = α if μ = ν (diagonal)
# H_μν = β if atoms μ and ν are bonded
# H_μν = 0 otherwise
huckel_matrix = alpha * np.eye(connectivity.shape[0])
huckel_matrix[connectivity == 1] = beta
### ANCHOR_END: exercise_01_b

### ANCHOR: exercise_01_c
# Set up 3D grid for orbital calculation
spacing = 0.2  # Grid spacing
max_extent = 6.0  # Maximum extent from origin
npoints = int(2 * max_extent / spacing) + 1

# Create 1D coordinate arrays
x = np.linspace(-max_extent, max_extent, npoints)
y = np.linspace(-max_extent, max_extent, npoints)
z = np.linspace(-max_extent, max_extent, npoints)

# Create 3D grid using meshgrid
X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

# Function to calculate p_z orbital
def pz_orbital(x, y, z, center):
    """Calculate p_z orbital centered at given coordinates."""
    x0, y0, z0 = center
    r = np.sqrt((x - x0)**2 + (y - y0)**2 + (z - z0)**2)
    return (z - z0) * np.exp(-r)

# Calculate atomic orbitals on the grid
atomic_orbitals = []
for _, x0, y0, z0 in heavy_atoms:
    # Calculate p_z orbital for this center
    pz = pz_orbital(X, Y, Z, (x0, y0, z0))
    atomic_orbitals.append(pz)
### ANCHOR_END: exercise_01_c

### ANCHOR: exercise_01_d
# Calculate eigenvalues and eigenvectors of Hückel matrix
energies, coefficients = np.linalg.eigh(huckel_matrix)

# Construct molecular orbitals
molecular_orbitals = []
for orbital_idx in range(len(energies)):
    # Initialize orbital grid
    mo = np.zeros((npoints, npoints, npoints))
    
    # Sum over atomic contributions according to LCAO
    for atom_idx in range(len(heavy_atoms)):
        mo += coefficients[atom_idx, orbital_idx] * atomic_orbitals[atom_idx]
    
    molecular_orbitals.append(mo)

# Optional: Save molecular orbitals as cube files
from utils import CubeFile

for i, orbital in enumerate(molecular_orbitals):
    cube = CubeFile(atoms=atoms, data=orbital, spacing=spacing)
    cube.dump(f"benzene_orbital_{i+1}.cube", 
             comment=f"Benzene molecular orbital {i+1}, E = {energies[i]:.3f}β")
### ANCHOR_END: exercise_01_d