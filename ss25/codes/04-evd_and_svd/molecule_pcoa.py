#!/usr/bin/env python

### ANCHOR: imports
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnnotationBbox, OffsetImage
from rdkit import Chem
from rdkit.Chem import Draw
### ANCHOR_END: imports

### ANCHOR: constants
DATASET = 'gdb9_subset_5.sdf'
BOND_TYPES = {
    'CC1': 0, 'CC2': 1, 'CC3': 2, 'CC4': 3,
    'CN1': 4, 'CN2': 5, 'CN3': 6, 'CN4': 7,
    'CO1': 8, 'CO2': 9, 'CO3': 10, 'CO4': 11,
    'NN1': 12, 'NN2': 13, 'NN3': 14, 'NN4': 15,
    'NO1': 16, 'NO2': 17, 'NO3': 18, 'NO4': 19,
    'OO1': 20, 'OO2': 21, 'OO3': 22, 'OO4': 23,
}
### ANCHOR_END: constants


### ANCHOR: get_fingerprint_function
def get_fingerprint(mol):
    fingerprint = np.zeros(len(BOND_TYPES), dtype=np.int32)
    for bond in mol.GetBonds():
        sym1 = bond.GetBeginAtom().GetSymbol()
        sym2 = bond.GetEndAtom().GetSymbol()
        rd_btype = bond.GetBondType()
        
        if rd_btype == Chem.rdchem.BondType.DOUBLE:
            btype_num = 2
        elif rd_btype == Chem.rdchem.BondType.TRIPLE:
            btype_num = 3
        elif rd_btype == Chem.rdchem.BondType.AROMATIC:
            btype_num = 4
        else:
            btype_num = 1
        
        btype = ''.join(sorted([sym1, sym2])) + str(btype_num)

        if btype in BOND_TYPES:
            fingerprint[BOND_TYPES[btype]] += 1
    
    return fingerprint
### ANCHOR_END: get_fingerprint_function


### ANCHOR: get_fingerprints
mols = [mol for mol in Chem.SDMolSupplier(DATASET)]
fingerprints = np.array([get_fingerprint(mol) for mol in mols])
### ANCHOR_END: get_fingerprints


### ANCHOR: calculate_distances
n = len(fingerprints)
distances = np.zeros((n, n))
for i in range(0, n):
    for j in range(i + 1, n):
        distances[i, j] = np.sum(np.abs(fingerprints[i] - fingerprints[j]))
        distances[j, i] = distances[i, j]
### ANCHOR_END: calculate_distances


### ANCHOR: pcoa
c_mat = np.eye(n) - np.ones((n, n)) / n
b_mat = -0.5 * np.linalg.multi_dot([c_mat, distances**2, c_mat])
e, v = np.linalg.eigh(b_mat)
embedding = np.dot(v[:, -2:], np.diag(np.sqrt(e[-2:])))
### ANCHOR_END: pcoa


### ANCHOR: explained_variance
eta = np.sum(e[::-1][:2]) / np.sum(e[e > 0])
assert np.isclose(eta, 0.50722839)
print(eta)
### ANCHOR_END: explained_variance

### ANCHOR: plot
fig, ax = plt.subplots(figsize=(8, 6))

ax.scatter(embedding[:, 1], embedding[:, 0])

ax.set_aspect('equal')
ax.set_xlabel('PC1')
ax.set_ylabel('PC2')
fig.tight_layout()

plt.show()
### ANCHOR_END: plot

fig.savefig('../../assets/figures/04-evd_and_svd/molecule_pcoa.svg')

### ANCHOR: interactive_plot
fig, ax = plt.subplots(figsize=(8, 6))

sc = ax.scatter(embedding[:, 1], embedding[:, 0])

ax.set_aspect('equal')
ax.set_xlabel('PC1')
ax.set_ylabel('PC2')
fig.tight_layout(rect=[0.0, 0.0, 0.9, 0.9])

imagebox = OffsetImage(np.zeros((100, 100, 3)), zoom=0.5)
imagebox.image.axes = ax

ab = AnnotationBbox(
    imagebox, (40, 40), xycoords='data', 
    boxcoords="offset points", arrowprops=dict(arrowstyle="->"), 
)
ax.add_artist(ab)


def update_annotation_box(idx):
    ab.xy = (embedding[idx, 1], embedding[idx, 0])
    mol = mols[idx]
    Chem.rdDepictor.Compute2DCoords(mol)
    img = Draw.MolToImage(mol, size=(100, 100), wedgeBonds=False)
    ab.offsetbox.set_data(img)


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
### ANCHOR_END: interactive_plot

