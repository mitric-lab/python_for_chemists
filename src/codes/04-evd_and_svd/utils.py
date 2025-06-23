#!/usr/bin/env python
"""Utility functions for working with molecular structure files and volumetric data.

This module provides functions for reading XYZ files, analyzing molecular connectivity,
and handling Gaussian cube files for volumetric data visualization.
"""

import numpy as np

# Conversion constant from Ångström to Bohr (exact value)
ANGSTROM_TO_BOHR = 1.8897259886

# Dictionary mapping atomic symbols to atomic numbers
ATOMIC_NUMBERS = {
    "H": 1,
    "C": 6,
    "N": 7,
    "O": 8,
    "F": 9,
}

def read_xyz(fname, use_bohr=True):
    """Read atomic coordinates from an XYZ file.
    
    Parameters
    ----------
    fname : str
        Path to the XYZ file
    use_bohr : bool, optional
        If True, converts coordinates from Ångström to Bohr (default: True)

    Returns
    -------
    list of tuples
        List of (atomic_number, x, y, z) tuples for each atom.
        Coordinates will be in Bohr if use_bohr=True, otherwise in Ångström.
    """
    with open(fname, 'r') as f:
        # Read number of atoms
        natoms = int(f.readline())
        # Skip comment line
        f.readline()
        # Read atomic coordinates
        atoms = []
        for _ in range(natoms):
            symbol, *coords = f.readline().split()
            coords = [float(x) * (ANGSTROM_TO_BOHR if use_bohr else 1.0) for x in coords]
            atoms.append([ATOMIC_NUMBERS[symbol]] + coords)
        return atoms

def get_connectivity_matrix(atoms, threshold=3.0):
    """Calculate the connectivity matrix for a molecular structure.
    
    Two atoms are considered bonded if their distance is less than the threshold.
    
    Parameters
    ----------
    atoms : list of tuples
        List of (atomic_number, x, y, z) tuples for each atom.
        Coordinates must be in the same units as the threshold.
    threshold : float, optional
        Distance threshold for bonding in the same units as the atomic coordinates
        (default: 1.5)

    Returns
    -------
    numpy.ndarray
        Square connectivity matrix where 1 indicates a bond between atoms
        and 0 indicates no bond
    """    
    natoms = len(atoms)
    connectivity_matrix = np.zeros((natoms, natoms))
    
    for i in range(natoms):
        for j in range(i+1, natoms):
            # Calculate distance between atoms i and j
            dist = np.linalg.norm(np.array(atoms[i][1:4]) - np.array(atoms[j][1:4]))
            if dist < threshold:
                connectivity_matrix[i, j] = connectivity_matrix[j, i] = 1
                
    return connectivity_matrix

class CubeFile:
    """Class for handling Gaussian cube files.
    
    This class provides functionality for reading, writing, and manipulating
    Gaussian cube files, which store volumetric data and molecular geometries.
    All coordinates are handled in Bohr units internally.
    """
    
    def __init__(self, atoms=None, data=None, origin=None, spacing=0.2):
        """Initialize a CubeFile object.
        
        Parameters
        ----------
        atoms : list of tuples, optional
            List of (atomic_number, x, y, z) tuples for each atom.
            Coordinates must be in Bohr.
        data : numpy.ndarray, optional
            3D array containing the volumetric data
        origin : array-like, optional
            Origin of the grid (x, y, z) in Bohr. If None, will be calculated
            to center the grid around zero.
        spacing : float, optional
            Grid spacing in Bohr (default: 0.2)
        """
        if atoms is not None and data is not None:
            # Store atomic data
            self.atoms = [[Z, 0.0, x, y, z] for Z, x, y, z in atoms]
            self.natoms = len(atoms)
            
            # Store volumetric data
            self.data = np.array(data)
            self.NX, self.NY, self.NZ = data.shape
            
            # Set up grid vectors (orthogonal for simplicity)
            self.dx = self.dy = self.dz = spacing
            self.X = np.array([spacing, 0.0, 0.0])
            self.Y = np.array([0.0, spacing, 0.0])
            self.Z = np.array([0.0, 0.0, spacing])
            
            # Calculate or set origin
            if origin is None:
                # Center grid around zero
                extent = np.array([self.NX-1, self.NY-1, self.NZ-1]) * spacing
                self.origin = -extent/2
            else:
                self.origin = np.array(origin)
            
            # Set up coordinate arrays
            self._setup_grid()

    def _setup_grid(self):
        """Set up coordinate arrays for the grid points."""
        self.rx = self.origin[0] + np.arange(self.NX) * self.dx
        self.ry = self.origin[1] + np.arange(self.NY) * self.dy
        self.rz = self.origin[2] + np.arange(self.NZ) * self.dz

    @classmethod
    def read(cls, fname):
        """Read a Gaussian cube file.
        
        Parameters
        ----------
        fname : str
            Path to the cube file
            
        Returns
        -------
        CubeFile
            New CubeFile instance with data read from file.
            All coordinates are in Bohr.
        """
        instance = cls()
        
        with open(fname, 'r') as f:
            # Skip comment lines
            f.readline()
            f.readline()

            # Read number of atoms and origin
            tkns = f.readline().split()
            instance.natoms = int(tkns[0])
            instance.origin = np.array([float(x) for x in tkns[1:4]])

            # Read grid information
            tkns = f.readline().split()
            instance.NX = abs(int(tkns[0]))
            instance.X = np.array([float(x) for x in tkns[1:4]])
            
            tkns = f.readline().split()
            instance.NY = abs(int(tkns[0]))
            instance.Y = np.array([float(x) for x in tkns[1:4]])
            
            tkns = f.readline().split()
            instance.NZ = abs(int(tkns[0]))
            instance.Z = np.array([float(x) for x in tkns[1:4]])

            # Read atomic positions
            instance.atoms = []
            for _ in range(instance.natoms):
                tkns = f.readline().split()
                instance.atoms.append([int(tkns[0])] + [float(x) for x in tkns[1:5]])

            # Read volumetric data
            instance.data = np.zeros((instance.NX, instance.NY, instance.NZ))
            i = 0
            for line in f:
                for val in line.split():
                    instance.data[
                        i//(instance.NY*instance.NZ),
                        (i//instance.NZ)%instance.NY,
                        i%instance.NZ
                    ] = float(val)
                    i += 1

            if i != instance.NX*instance.NY*instance.NZ:
                raise ValueError("Error reading cube file: incorrect number of data points")

        # Set up grid parameters
        instance.dx, instance.dy, instance.dz = [
            np.linalg.norm(x) for x in (instance.X, instance.Y, instance.Z)
        ]
        instance._setup_grid()
        
        return instance

    def dump(self, fname, comment=""):
        """Write data to a Gaussian cube file.
        
        Parameters
        ----------
        fname : str
            Path to write the cube file
        comment : str, optional
            Comment to include in the cube file header
        
        Notes
        -----
        All coordinates are written in Bohr units with positive voxel counts.
        """
        with open(fname, 'w') as f:
            # Write header
            f.write(f" {comment}\n")
            f.write(" Cube file generated by utils.py (coordinates in Bohr)\n")
            
            # Write number of atoms and origin
            f.write(f" {self.natoms:4d} {self.origin[0]:12.6f} {self.origin[1]:12.6f} {self.origin[2]:12.6f}\n")
            
            # Write grid information
            f.write(f" {self.NX:4d} {self.X[0]:12.6f} {self.X[1]:12.6f} {self.X[2]:12.6f}\n")
            f.write(f" {self.NY:4d} {self.Y[0]:12.6f} {self.Y[1]:12.6f} {self.Y[2]:12.6f}\n")
            f.write(f" {self.NZ:4d} {self.Z[0]:12.6f} {self.Z[1]:12.6f} {self.Z[2]:12.6f}\n")
            
            # Write atomic coordinates
            for atom in self.atoms:
                f.write(f" {atom[0]:4d} {atom[1]:12.6f} {atom[2]:12.6f} {atom[3]:12.6f} {atom[4]:12.6f}\n")
            
            # Write volumetric data
            for ix in range(self.NX):
                for iy in range(self.NY):
                    for iz in range(self.NZ):
                        f.write(f" {self.data[ix, iy, iz]:14.6e}")
                        if (iz % 6 == 5):
                            f.write("\n")
                    f.write("\n")


if __name__ == "__main__":
    # Example: Create a cube file for benzene with Gaussian density on carbon atoms
    
    # Read benzene structure (automatically converts to Bohr)
    benzene_atoms = read_xyz("benzene.xyz", use_bohr=True)
    heavy_atoms = [atom for atom in benzene_atoms if atom[0] != 1]
    
    # Set up grid parameters
    spacing = 0.2  # Bohr
    max_extent = 6.0  # Maximum extent in any direction in Bohr
    
    # Calculate number of points needed to cover max_extent with given spacing
    # We want the grid to go from -max_extent to +max_extent
    # Number of points = total_extent/spacing + 1 to include both endpoints
    npoints = int(2 * max_extent / spacing) + 1
    
    # Create centered grid from -max_extent to +max_extent
    x = np.linspace(-max_extent, max_extent, npoints)
    y = np.linspace(-max_extent, max_extent, npoints)
    z = np.linspace(-max_extent, max_extent, npoints)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    
    # Create Gaussian density on carbon atoms
    density = np.zeros((npoints, npoints, npoints))
    for atomic_num, x0, y0, z0 in heavy_atoms:
        r = np.sqrt((X - x0)**2 + (Y - y0)**2 + (Z - z0)**2)
        density += np.exp(-r**2 / 2)
    
    # Create and save cube file
    cube = CubeFile(atoms=benzene_atoms, data=density, spacing=spacing)
    cube.dump("benzene.cube", comment="Benzene molecule with Gaussian density on carbons")
    
    # Verify by reading back
    cube2 = CubeFile.read("benzene.cube")
    
    # Print information
    print("\nCube file information:")
    print(f"Dimensions: {density.shape}")
    print(f"Grid spacing: {cube.dx:.3f} Bohr")
    print(f"Origin: {cube.origin}")
    print(f"Data matches: {np.allclose(cube.data, cube2.data)}")
