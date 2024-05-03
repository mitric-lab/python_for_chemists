#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnnotationBbox, OffsetImage
from rdkit import Chem
from rdkit.Chem import Draw
from collections import Counter
from mol import parse_mol

DATASET = 'gdb9_subset_7.sdf'
BOND_TYPES = {
    'CC1': 0, 'CC2': 1, 'CC3': 2, 'CC4': 3,
    'CN1': 4, 'CN2': 5, 'CN3': 6, 'CN4': 7,
    'CO1': 8, 'CO2': 9, 'CO3': 10, 'CO4': 11,
    'NN1': 12, 'NN2': 13, 'NN3': 14, 'NN4': 15,
    'NO1': 16, 'NO2': 17, 'NO3': 18, 'NO4': 19,
    'OO1': 20, 'OO2': 21, 'OO3': 22, 'OO4': 23,
    'CH1': 24, 'NH1': 25, 'OH1': 26,
}


def get_fingerprint(mol):
    fingerprint = np.zeros(len(BOND_TYPES), dtype=np.int32)
    for btype, count in mol['bonds'].items():
        if btype in BOND_TYPES:
            fingerprint[BOND_TYPES[btype]] = count
    return fingerprint


mols = []
mol_strings = []
with open(DATASET, 'r') as f:
    mol_lines = []
    for line in f:
        if line.startswith('$$$$'):
            mol = parse_mol(mol_lines)
            mol_string = ''.join(mol_lines)
            mols.append(mol)
            mol_strings.append(mol_string)
            mol_lines = []
        else:
            mol_lines.append(line)

fingerprints = np.array([get_fingerprint(mol) for mol in mols])
# formulas = [mol['comment'].strip() for mol in mols]

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
for i in range(0, 10):
    print(np.sum(e[::-1][:i]) / np.sum(e[e > 0]))

fig, ax = plt.subplots()
sc = ax.scatter(embedding[:, 0], embedding[:, 1])

imagebox = OffsetImage(np.zeros((100, 100, 3)), zoom=0.5)
imagebox.image.axes = ax

ab = AnnotationBbox(
    imagebox, (40, 40), xycoords='data', 
    boxcoords="offset points", arrowprops=dict(arrowstyle="->"), 
)
ax.add_artist(ab)


def update_annotation_box(idx):
    ab.xy = (embedding[idx, 0], embedding[idx, 1])
    mol = Chem.MolFromMolBlock(mol_strings[idx])
    Chem.rdDepictor.Compute2DCoords(mol)
    img = Draw.MolToImage(mol, size=(100, 100), wedgeBonds=False)
    ab.offsetbox.set_data(img)
    # ab.offsetbox.get_bbox_patch().set_facecolor('white')
    # ab.offsetbox.get_bbox_patch().set_alpha(0.5)


def hover(event):
    vis = ab.get_visible()
    if event.inaxes == ax:
        cont, ind = sc.contains(event)
        if cont:
            update_annotation_box(ind['ind'][0])
            ab.set_visible(True)
            fig.canvas.draw_idle()
        else:
            if vis:
                ab.set_visible(False)
                fig.canvas.draw_idle()


fig.canvas.mpl_connect('motion_notify_event', hover)


plt.show()

