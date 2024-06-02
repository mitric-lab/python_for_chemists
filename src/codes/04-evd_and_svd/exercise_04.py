#!/usr/bin/env python

### ANCHOR: exercise_01_b
import numpy as np
import matplotlib.pyplot as plt

# Build the Hamiltonian matrix of Hexatriene
num_atoms = 6
alpha = 0.0
beta = -1.0
hamiltonian = np.zeros((num_atoms, num_atoms))
rows, cols = np.diag_indices(num_atoms)
hamiltonian[rows, cols] = alpha
hamiltonian[rows[:-1], cols[1:]] = beta
hamiltonian[rows[1:], cols[:-1]] = beta

# Calculate the eigenvalues and eigenvectors
energies, coefficients = np.linalg.eigh(hamiltonian)

# Plot the orbitals and energies
fig, ax = plt.subplots(figsize=(4, 8))

for i, e in enumerate(energies):
    ax.axhline(e, xmin=0.1, xmax=0.9, color='k', ls='--')
    for j in range(num_atoms):
        color = 'blue' if coefficients[j, i] > 0 else 'red'
        ms = coefficients[j, i]**2 * 100
        ax.plot(j/(num_atoms-1), e, 'o', color=color, ms=ms)

ax.set_xlim(-0.1, 1.1)
ax.set_ylim(-2.5, 2.5)

ax.set_xticks([])
ax.set_xlabel('')
ax.set_yticks(energies)
ax.set_ylabel('energy')

# Show only y-axis
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')

fig.tight_layout()

plt.show()
### ANCHOR_END: exercise_01_b

#fig.savefig('../../assets/figures/04-evd_and_svd/hmo_orbitals.svg')

### ANCHOR: exercise_01_c
# Build the Hamiltonian matrix for Benzene
hamiltonian[0, -1] = beta
hamiltonian[-1, 0] = beta

# Calculate the eigenvalues and eigenvectors
energies, coefficients = np.linalg.eigh(hamiltonian)

# Calculate total energy of Benzene
E = 2 * np.sum(energies[:num_atoms//2]) # 6a + 8b vs 6a + 6b for three ethylene

assert np.isclose(E, -8.0)
### ANCHOR_END: exercise_01_c

### ANCHOR: exercise_02_import
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

# Import Yale Faces dataset
path = 'allFaces.mat'
data = sp.io.loadmat(path)

# Extract the data
faces = data['faces']
n = int(data['n'][0][0])
m = int(data['m'][0][0])
nfaces = data['nfaces'].flatten() # number of faces per person

# Plot first face
plt.imshow(faces[:,0].reshape(m,n).T, cmap='gray')
### ANCHOR_END: exercise_02_import

### ANCHOR: exercise_02_a
# Indices of the first face of each person
idx = np.cumsum(nfaces) - nfaces

# Define helper function for plotting faces
def plot_faces(data, n, m, ncols=5):
    nimages = data.shape[1]
    fig, axes = plt.subplots(nrows=(nimages-1)//ncols + 1, ncols=ncols, figsize=(2*ncols, 2*((nimages-1)//ncols + 1)))
    axes = axes.flatten()  # flatten the axes array

    for i in range(nimages):
        # Reshape the data point to a square image
        axes[i].imshow(data[:, i].reshape(m,n).T, cmap='gray')
        axes[i].axis('off')  # hide the axes

    # Turn off remaining axes
    for j in range(i+1, len(axes)):
        axes[j].axis('off')

    plt.tight_layout()
    plt.show()

# Plot 5 random faces from the dataset
np.random.seed(0)
rand_idx = np.random.choice(faces.shape[1], 5)
plot_faces(faces[:,rand_idx], n, m, ncols=5)
### ANCHOR_END: exercise_02_a

#plt.savefig('../../assets/figures/04-evd_and_svd/eigenfaces_samples.svg')

# Plot the first face of each person
plot_faces(faces[:,idx], n, m, ncols=6)

#plt.savefig('../../assets/figures/04-evd_and_svd/eigenfaces_first_faces.svg')

### ANCHOR: exercise_02_b
# Define training set
X = faces[:,:idx[36]]
print(idx[36])

# Calculate the mean face and subtract it from the training set
mean_face = np.mean(X, axis=1)
X = X - mean_face[:, np.newaxis]

# Perform the SVD
U, S, Vt = np.linalg.svd(X, full_matrices=False)

# Plot the first 5 eigenfaces
plot_faces(U[:,:5], n, m, ncols=5)
### ANCHOR_END: exercise_02_b

#plt.savefig('../../assets/figures/04-evd_and_svd/eigenfaces_eigenfaces.svg')

### ANCHOR: exercise_02_c
# Define test face
test_face = faces[:,idx[36]] - mean_face

# Reconstruct the test face using different number of eigenfaces
r_list = [25, 50, 100, 200, 400, 800, 1600]
recon_faces = []
recon_faces.append(test_face + mean_face)

for r in r_list:
    alpha = U[:,:r].T @ test_face
    recon_face = U[:,:r] @ alpha + mean_face
    recon_faces.append(recon_face)

# Plot the reconstructed faces with the original test face
plot_faces(np.array(recon_faces).T, n, m, ncols=4)
### ANCHOR_END: exercise_02_c

#plt.savefig('../../assets/figures/04-evd_and_svd/eigenfaces_reconstruction.svg')

### ANCHOR: exercise_02_d
# Import dog image
dog_image = plt.imread('dog.jpg', format='jpeg')[:,:,0].T.flatten() - mean_face

# Reconstruct the dog image using different number of eigenfaces
r_list = [25, 50, 100, 200, 400, 800, 1600]
recon_images = []
recon_images.append(dog_image + mean_face)

for r in r_list:
    alpha = U[:,:r].T @ dog_image
    recon_image = U[:,:r] @ alpha + mean_face
    recon_images.append(recon_image)

# Plot the reconstructed dog images with the original dog image
plot_faces(np.array(recon_images).T, n, m, ncols=4)
### ANCHOR_END: exercise_02_d

#plt.savefig('../../assets/figures/04-evd_and_svd/eigenfaces_dog_reconstruction.svg')

### ANCHOR: exercise_03_a
from matplotlib.offsetbox import AnnotationBbox, OffsetImage
from rdkit import Chem
from rdkit.Chem import Draw

DATASET = 'gdb9_subset_5.sdf'
BOND_TYPES = {
    'CC1': 0, 'CC2': 1, 'CC3': 2, 'CC4': 3,
    'CN1': 4, 'CN2': 5, 'CN3': 6, 'CN4': 7,
    'CO1': 8, 'CO2': 9, 'CO3': 10, 'CO4': 11,
    'NN1': 12, 'NN2': 13, 'NN3': 14, 'NN4': 15,
    'NO1': 16, 'NO2': 17, 'NO3': 18, 'NO4': 19,
    'OO1': 20, 'OO2': 21, 'OO3': 22, 'OO4': 23,
    'CH1': 24, 'HN1': 25, 'HO1': 26,
}

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

# Get the fingerprints of the molecules with hydrogen atoms
mols = [mol for mol in Chem.SDMolSupplier(DATASET)]
mols_with_h = [Chem.AddHs(mol) for mol in mols]
fingerprints = np.array([get_fingerprint(mol) for mol in mols_with_h])

# Calculate the pairwise distances
n = len(fingerprints)
distances = np.zeros((n, n))
for i in range(0, n):
    for j in range(i + 1, n):
        distances[i, j] = np.sum(np.abs(fingerprints[i] - fingerprints[j]))
        distances[j, i] = distances[i, j]

# Perform classical MDS
c_mat = np.eye(n) - np.ones((n, n)) / n
b_mat = -0.5 * np.linalg.multi_dot([c_mat, distances**2, c_mat])
e, v = np.linalg.eigh(b_mat)
embedding = np.dot(v[:, -2:], np.diag(np.sqrt(e[-2:])))

# Plot the embedding as interactive scatter plot
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
### ANCHOR_END: exercise_03_a