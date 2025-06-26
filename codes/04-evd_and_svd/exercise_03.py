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
max_extent = 6.0  # Maximum extent from origin (both in Bohr)
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
### ANCHOR_END: exercise_01_d

### ANCHOR: exercise_01_e
# Optional: Save molecular orbitals as cube files
from utils import CubeFile

for i, orbital in enumerate(molecular_orbitals):
    cube = CubeFile(atoms=atoms, data=orbital, spacing=spacing)
    cube.dump(f"benzene_mo_{i+1}.cube", 
             comment=f"Benzene molecular orbital {i+1}, E = {energies[i]:.3f}β")
### ANCHOR_END: exercise_01_e

### ANCHOR: exercise_02_a
def get_Hueckel_energies_and_coefficients(path_to_xyz_file, alpha=0.0, beta=-1.0):
    """Calculate Hückel energies and MO coefficients for a given molecule."""

    # Read structure from xyz file and get connectivity matrix
    atoms = np.array(read_xyz(path_to_xyz_file))
    heavy_atoms = [atom for atom in atoms if atom[0] != 1]
    connectivity = get_connectivity_matrix(heavy_atoms)

    # Get index of nitrogen atom (atomic number 7)
    nitrogen_idx = np.where(atoms[:, 0] == 7)[0][0]
    
    # Create Hückel matrix
    huckel_matrix = alpha * np.eye(connectivity.shape[0])
    huckel_matrix[connectivity == 1] = beta
    if nitrogen_idx:
        huckel_matrix[nitrogen_idx, nitrogen_idx] = -1.5

    # Calculate eigenvalues and eigenvectors of Hückel matrix
    energies, coefficients = np.linalg.eigh(huckel_matrix)

    return energies, coefficients

def get_AOs_on_grid(path_to_xyz_file, spacing=0.2):
    """Calculate atomic orbitals on a 3D grid for all heavy atoms."""

    # Read structure from xyz file
    atoms = read_xyz(path_to_xyz_file)
    heavy_atoms = [atom for atom in atoms if atom[0] != 1]

    # Set grid size based on molecule size
    max_extent = 3.0 * max(np.max(np.abs(atom[1:4])) for atom in heavy_atoms)
    npoints = int(2 * max_extent / spacing) + 1

    # Create 1D coordinate arrays
    x = np.linspace(-max_extent, max_extent, npoints)
    y = np.linspace(-max_extent, max_extent, npoints)
    z = np.linspace(-max_extent, max_extent, npoints)

    # Create 3D grid using meshgrid
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

    # Calculate atomic orbitals on the grid
    atomic_orbitals = []
    for _, x0, y0, z0 in heavy_atoms:
        # Calculate p_z orbital for this center
        pz = pz_orbital(X, Y, Z, (x0, y0, z0))
        atomic_orbitals.append(pz)

    return atomic_orbitals

def get_Hueckel_MOs_on_grid(coefficients, atomic_orbitals):
    """Calculate Hückel molecular orbitals on a 3D grid."""

    # Get number of points in the grid
    npoints = atomic_orbitals[0].shape[0]
    molecular_orbitals = []
    for orbital_idx in range(len(coefficients)):
        # Initialize orbital grid
        mo = np.zeros((npoints, npoints, npoints))

        # Sum over atomic contributions according to LCAO
        for atom_idx in range(len(coefficients)):
            mo += coefficients[atom_idx, orbital_idx] * atomic_orbitals[atom_idx]

        molecular_orbitals.append(mo)

    return molecular_orbitals

# Set parameters
spacing = 0.5
path_to_xyz_file = "benzoindol.xyz"

# Read structure from test.xyz
atoms = read_xyz(path_to_xyz_file)

# Calculate Hückel energies and coefficients
energies, coefficients = get_Hueckel_energies_and_coefficients(path_to_xyz_file)

# Calculate atomic orbitals on grid
atomic_orbitals = get_AOs_on_grid(path_to_xyz_file, spacing=spacing)

# Calculate Hückel molecular orbitals on grid
molecular_orbitals = get_Hueckel_MOs_on_grid(coefficients, atomic_orbitals)

# Save molecular orbitals as cube files
for i, orbital in enumerate(molecular_orbitals):
    cube = CubeFile(atoms=atoms, data=orbital, spacing=spacing)
    cube.dump(f"benzoindol_mo_{i+1}.cube", 
             comment=f"Molecular orbital {i+1}, E = {energies[i]:.3f}β")
### ANCHOR_END: exercise_02_a

T = np.array([[0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0.11, 0],
              [0, 0, 0, 0.41, -0.41, 0],
              [0, 0, 0, -0.55, -0.51, 0]])

np.save("tdm.npy", T)

### ANCHOR: exercise_02_b
# Load transition density matrix
T = np.load("tdm.npy")

# Print shape and some information about T
print(f"Shape of transition density matrix: {T.shape}")
print("This means we have:")
print(f"- {T.shape[0]} occupied orbitals")
print(f"- {T.shape[1]} unoccupied orbitals")
### ANCHOR_END: exercise_02_b

### ANCHOR: exercise_02_c
# Perform SVD on transition density matrix
U, sigma, VT = np.linalg.svd(T)

# Calculate NTOs
def calculate_NTOs(molecular_orbitals, U, VT, n_occ, n_virt):
    """Calculate hole and electron NTOs from MOs and SVD matrices."""
    ntos_hole = []
    ntos_electron = []
    
    # Calculate hole NTOs (from occupied MOs)
    for j in range(U.shape[1]):  # Loop over NTO pairs
        nto = np.zeros_like(molecular_orbitals[0])
        for i in range(n_occ):  # Sum over occupied MOs
            nto += U[i,j] * molecular_orbitals[i]
        ntos_hole.append(nto)
    
    # Calculate electron NTOs (from unoccupied MOs)
    for j in range(VT.shape[0]):  # Loop over NTO pairs
        nto = np.zeros_like(molecular_orbitals[0])
        for a in range(n_virt):  # Sum over virtual MOs
            nto += VT[j,a] * molecular_orbitals[n_occ + a]
        ntos_electron.append(nto)
    
    return ntos_hole, ntos_electron

# Get number of occupied and virtual orbitals
n_occ = T.shape[0]  # Number of rows in T
n_virt = T.shape[1]  # Number of columns in T

# Calculate NTOs
ntos_hole, ntos_electron = calculate_NTOs(molecular_orbitals, U, VT, n_occ, n_virt)
### ANCHOR_END: exercise_02_c

### ANCHOR: exercise_02_d
# Print singular values and their contributions
print("\nNTO pair contributions:")
for i, s in enumerate(sigma):
    contribution = s**2 / np.sum(sigma**2) * 100
    print(f"NTO pair {i+1}: {contribution:.1f}%")

# Save NTOs as cube files
for i, (hole, electron) in enumerate(zip(ntos_hole, ntos_electron)):
    # Only save NTOs with significant contributions (e.g., > 10%)
    contribution = sigma[i]**2 / np.sum(sigma**2) * 100
    if contribution > 10.0:
        # Save hole NTO
        cube = CubeFile(atoms=atoms, data=hole, spacing=spacing)
        cube.dump(f"benzoindol_nto_hole_{i+1}.cube",
                comment=f"Hole NTO {i+1}, contribution = {contribution:.1f}%")
        
        # Save electron NTO
        cube = CubeFile(atoms=atoms, data=electron, spacing=spacing)
        cube.dump(f"benzoindol_nto_electron_{i+1}.cube",
                comment=f"Electron NTO {i+1}, contribution = {contribution:.1f}%")
### ANCHOR_END: exercise_02_d


