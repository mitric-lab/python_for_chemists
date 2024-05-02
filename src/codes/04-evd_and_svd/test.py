#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from mol import parse_mol

DATASET = 'gdb9_subset.sdf'
BOND_TYPES = {
    'CC1': 0, 'CC2': 1, 'CC3': 2, 'CC4': 3,
    'CN1': 4, 'CN2': 5, 'CN3': 6, 'CN4': 7,
    'CO1': 8, 'CO2': 9, 'CO3': 10, 'CO4': 11,
    'NN1': 12, 'NN2': 13, 'NN3': 14, 'NN4': 15,
    'NO1': 16, 'NO2': 17, 'NO3': 18, 'NO4': 19,
    'OO1': 20, 'OO2': 21, 'OO3': 22, 'OO4': 23,
}


def get_fingerprint(mol):
    fingerprint = np.zeros(len(BOND_TYPES), dtype=np.int32)
    for btype, count in mol['bonds'].items():
        if btype in BOND_TYPES:
            fingerprint[BOND_TYPES[btype]] = count
    return fingerprint


mols = []
with open(DATASET, 'r') as f:
    mol_lines = []
    for line in f:
        if line.startswith('$$$$'):
            mol = parse_mol(mol_lines)
            mols.append(mol)
            mol_lines = []
        else:
            mol_lines.append(line)

fingerprints = np.array([get_fingerprint(mol) for mol in mols])
formulas = [mol['comment'].strip() for mol in mols]

n = len(mols)
distances = np.zeros((n, n))
for i in range(len(mols)):
    for j in range(i + 1, len(mols)):
        distances[i, j] = np.sum(np.abs(fingerprints[i] - fingerprints[j]))
        distances[j, i] = distances[i, j]

c_mat = np.eye(n) - np.ones((n, n)) / n
b_mat = -0.5 * np.linalg.multi_dot([c_mat, distances**2, c_mat])
e, v = np.linalg.eigh(b_mat)

embedding = np.dot(v[:, -2:], np.diag(np.sqrt(e[-2:])))

fig, ax = plt.subplots()
ax.scatter(embedding[:, 0], embedding[:, 1])

for i in range(0, n):
    ax.annotate(formulas[i], (embedding[i, 0], embedding[i, 1] + np.random.normal(0, 0.5)))

plt.show()

