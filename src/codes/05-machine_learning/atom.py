#!/usr/bin/env python

### ANCHOR: imports
import numpy as np
### ANCHOR_END: imports

### ANCHOR: atom_init
class Atom:
    def __init__(self, symbol, charge):
        self.symbol = symbol
        self.charge = charge
### ANCHOR_END: atom_init

### ANCHOR: atom_methods
    def get_atomic_number(self):
        atomic_numbers = {'H': 1, 'C': 6, 'N': 7, 'O': 8, 'F': 9}
        return atomic_numbers[self.symbol]

    def set_charge(self, charge):
        self.charge = charge

    def get_electron_config(self):
        num_electrons = self.get_atomic_number() - self.charge
        config = []
        orbitals = [("1s", 2), ("2s", 2), ("2p", 6)]

        for orbital, max_electrons in orbitals:
            if num_electrons <= 0:
                break
            electrons_in_orbital = min(num_electrons, max_electrons)
            config.append(f"{orbital}^{electrons_in_orbital}")
            num_electrons -= electrons_in_orbital

        return ' '.join(config)
### ANCHOR_END: atom_methods

### ANCHOR: atom_print
    def __str__(self):
        return f"{self.symbol} ({self.charge})"
### ANCHOR_END: atom_print

### ANCHOR: atom_example
# Create an Atom object
atom = Atom('C', 0)
print(atom) # C (0)
print(atom.get_atomic_number()) # 6
print(atom.get_electron_config()) # 1s^2 2s^2 2p^2
atom.set_charge(1)
print(atom) # C (1)
print(atom.get_electron_config()) # 1s^2 2s^2 2p^1
### ANCHOR_END: atom_example

### ANCHOR: atom_example_2
print(atom.symbol) # C
print(atom.charge) # 1
### ANCHOR_END: atom_example_2

