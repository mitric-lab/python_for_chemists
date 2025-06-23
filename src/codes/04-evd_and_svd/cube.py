#!/usr/bin/env python

import numpy as np

# Conversion factors
BOHR_TO_ANGSTROM = 0.529177249
ANGSTROM_TO_BOHR = 1.0 / BOHR_TO_ANGSTROM

class CubeFile:
    def __init__(self, atoms=None, data=None, origin=None, spacing=0.2, coords_in_angstrom=False):
        """Initialize a CubeFile object.
        
        Parameters
        ----------
        atoms : list of tuples, optional
            List of (atomic_number, x, y, z) tuples for each atom.
            Coordinates should be in Bohr by default, or in Angstrom if coords_in_angstrom=True
        data : ndarray, optional
            3D numpy array containing the volumetric data
        origin : array-like, optional
            Origin of the grid (x, y, z) in same units as atoms. If None, will be calculated automatically
        spacing : float, optional
            Grid spacing in same units as atoms (default: 0.2)
        coords_in_angstrom : bool, optional
            If True, input coordinates and spacing are assumed to be in Angstrom
        """
        if atoms is not None and data is not None:
            # Store the atoms in the input units (conversion happens in dump)
            self.atoms = [[Z, 0.0, x, y, z] for Z, x, y, z in atoms]
            self.natoms = len(atoms)
            
            # Store the volumetric data and its dimensions
            self.data = np.array(data)
            self.NX, self.NY, self.NZ = data.shape
            
            # Store whether coordinates are in Angstrom
            self.coords_in_angstrom = coords_in_angstrom
            
            # Set up the grid vectors (orthogonal for simplicity)
            self.dx = self.dy = self.dz = spacing
            self.X = np.array([spacing, 0.0, 0.0])
            self.Y = np.array([0.0, spacing, 0.0])
            self.Z = np.array([0.0, 0.0, spacing])
            
            # Calculate origin if not provided
            if origin is None:
                # For a grid of N points with spacing d, the total extent is (N-1)*d
                # To center it around 0, we start at -extent/2
                extent = np.array([self.NX-1, self.NY-1, self.NZ-1]) * spacing
                self.origin = -extent/2
            else:
                self.origin = np.array(origin)
            
            # Calculate the coordinate arrays
            self._setup_grid()

    def _setup_grid(self):
        """Set up the coordinate grid arrays."""
        self.rx = self.origin[0] + np.arange(0, self.NX)*self.dx
        self.ry = self.origin[1] + np.arange(0, self.NY)*self.dy
        self.rz = self.origin[2] + np.arange(0, self.NZ)*self.dz

    @classmethod
    def read(cls, fname):
        """Read a cube file.
        
        Parameters
        ----------
        fname : str
            Path to the cube file
            
        Returns
        -------
        CubeFile
            New instance of CubeFile. The coords_in_angstrom attribute will be set
            based on the sign of the voxel counts in the file.
        """
        instance = cls()
        
        with open(fname, 'r') as f:
            # echo comment
            for i in range(2):
                f.readline()

            # Number of atoms included in the file followed by the position of
            # the origin of the volumetric data
            tkns = f.readline().split()
            instance.natoms = int(tkns[0])
            instance.origin = np.array([float(x) for x in tkns[1:4]])

            # The next three lines give the number of voxels along
            # each axis (x, y, z) followed by the axis vector.
            # If the number of voxels is negative, coordinates are in Angstrom
            tkns = f.readline().split()
            instance.NX = abs(int(tkns[0]))
            instance.X = np.array([float(x) for x in tkns[1:4]])
            tkns = f.readline().split()
            instance.NY = abs(int(tkns[0]))
            instance.Y = np.array([float(x) for x in tkns[1:4]])
            tkns = f.readline().split()
            instance.NZ = abs(int(tkns[0]))
            instance.Z = np.array([float(x) for x in tkns[1:4]])
            
            # Determine units from the sign of the first non-zero voxel count
            for n in [int(tkns[0]) for tkns in [f.readline().split() for _ in range(3)]]:
                if n != 0:
                    instance.coords_in_angstrom = (n < 0)
                    break
            else:
                # If all voxel counts are 0 (shouldn't happen), assume Bohr
                instance.coords_in_angstrom = False

            # The last section in the header is one line for each atom
            instance.atoms = []
            for i in range(instance.natoms):
                tkns = f.readline().split()
                instance.atoms.append(
                    [int(tkns[0])] + [float(x) for x in tkns[1:5]]
                )

            # Volumetric data
            instance.data = np.zeros((instance.NX, instance.NY, instance.NZ))
            i = 0
            for s in f:
                for v in s.split():
                    instance.data[
                        i//(instance.NY*instance.NZ),
                        (i//instance.NZ)%instance.NY,
                        i%instance.NZ
                    ] = float(v)
                    i += 1

        if i != instance.NX*instance.NY*instance.NZ:
            raise ValueError("Error reading cube file: incorrect number of data points")

        # Set up the coordinate arrays
        instance.dx, instance.dy, instance.dz = [
            np.linalg.norm(x) for x in (instance.X, instance.Y, instance.Z)
        ]
        instance._setup_grid()
        
        return instance

    def dump(self, fname, comment=""):
        """Write the cube file.
        
        Parameters
        ----------
        fname : str
            Path to write the cube file
        comment : str, optional
            Comment to include in the cube file
            
        Notes
        -----
        The sign of the number of voxels in the output file indicates the coordinate units:
        - Positive: coordinates are in Bohr
        - Negative: coordinates are in Angstrom
        """
        with open(fname, 'w') as f:
            f.write(" CUBE File {}\n".format(comment))
            f.write(" generated by python (coordinates in {})\n".format(
                "Angstrom" if self.coords_in_angstrom else "Bohr"))
            
            # Convert coordinates to Bohr if necessary for output
            scale = 1.0 if not self.coords_in_angstrom else ANGSTROM_TO_BOHR
            
            # Write origin
            f.write(" {:4d} {:12.6f} {:12.6f} {:12.6f}\n".format(
                self.natoms, 
                self.origin[0] * scale,
                self.origin[1] * scale,
                self.origin[2] * scale
            ))
            
            # Write grid information - use sign to indicate units
            voxel_sign = -1 if self.coords_in_angstrom else 1
            f.write(" {:4d} {:12.6f} {:12.6f} {:12.6f}\n".format(
                self.NX * voxel_sign, 
                self.X[0] * scale, self.X[1] * scale, self.X[2] * scale
            ))
            f.write(" {:4d} {:12.6f} {:12.6f} {:12.6f}\n".format(
                self.NY * voxel_sign,
                self.Y[0] * scale, self.Y[1] * scale, self.Y[2] * scale
            ))
            f.write(" {:4d} {:12.6f} {:12.6f} {:12.6f}\n".format(
                self.NZ * voxel_sign,
                self.Z[0] * scale, self.Z[1] * scale, self.Z[2] * scale
            ))
            
            # Write atomic coordinates
            for atom in self.atoms:
                f.write(" {:4d} {:12.6f} {:12.6f} {:12.6f} {:12.6f}\n".format(
                    atom[0],
                    atom[1] * scale,
                    atom[2] * scale,
                    atom[3] * scale,
                    atom[4] * scale
                ))
            
            # Write volumetric data
            for ix in range(self.NX):
                for iy in range(self.NY):
                    for iz in range(self.NZ):
                        f.write(" {:14.6e}".format(self.data[ix, iy, iz]))
                        if (iz % 6 == 5):
                            f.write("\n")
                    f.write("\n")


if __name__ == "__main__":
    # Example: Create a cube file for a water molecule with a simple Gaussian density
    
    # Define water molecule geometry (in Angstrom)
    water_atoms = [
        (8, 0.0, 0.0, 0.0),      # O at origin
        (1, 0.96, 0.0, 0.0),     # H1 at 0.96 Å from O
        (1, -0.24, 0.93, 0.0),   # H2
    ]
    
    # Set up grid parameters (in Angstrom)
    spacing = 0.1  # Angstrom
    nx = ny = nz = 40
    
    # Create grid centered at origin
    # For N points with spacing d, the grid runs from -(N-1)d/2 to +(N-1)d/2
    extent = (nx - 1) * spacing
    x = np.linspace(-extent/2, extent/2, nx)
    y = np.linspace(-extent/2, extent/2, ny)
    z = np.linspace(-extent/2, extent/2, nz)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    
    # Create a simple Gaussian density centered on each atom
    density = np.zeros((nx, ny, nz))
    for _, x0, y0, z0 in water_atoms:
        r = np.sqrt((X - x0)**2 + (Y - y0)**2 + (Z - z0)**2)
        density += (Z - z0) * np.exp(-r)
    
    # Create and save the cube file (specifying that input coords are in Angstrom)
    cube = CubeFile(atoms=water_atoms, data=density, spacing=spacing, coords_in_angstrom=True)
    cube.dump("water.cube", comment="Water molecule with Gaussian density")
    
    # Read it back
    cube2 = CubeFile.read("water.cube")
    
    print("Created and read water.cube with dimensions:", density.shape)
    print(f"Grid spacing: {cube.dx:.3f} {'Angstrom' if cube.coords_in_angstrom else 'Bohr'}")
    print(f"Coordinates are in: {'Angstrom' if cube.coords_in_angstrom else 'Bohr'}")
    print("Origin:", cube.origin)
    print("Data matches:", np.allclose(cube.data, cube2.data))
    
    # Print grid extents to verify centering
    print("\nGrid extents:")
    print(f"x: {cube.rx.min():.2f} to {cube.rx.max():.2f}")
    print(f"y: {cube.ry.min():.2f} to {cube.ry.max():.2f}")
    print(f"z: {cube.rz.min():.2f} to {cube.rz.max():.2f}")
    
    # Verify symmetry
    print("\nGrid symmetry check:")
    print(f"x: {abs(cube.rx.min() + cube.rx.max()):.2e}")
    print(f"y: {abs(cube.ry.min() + cube.ry.max()):.2e}")
    print(f"z: {abs(cube.rz.min() + cube.rz.max()):.2e}")
    print("(Values should be close to 0 for perfect symmetry)")
    
    # Print coordinates
    print("\nAtomic coordinates:")
    for i, (Z, _, x, y, z) in enumerate(cube.atoms):
        print(f"Atom {i+1} (Z={Z}): ({x:.3f}, {y:.3f}, {z:.3f})") 