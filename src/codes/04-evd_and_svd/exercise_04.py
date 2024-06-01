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

### ANCHOR: exercise_02_a
import scipy as sp

# Import Yale Faces dataset
path = 'allFaces.mat'
data = sp.io.loadmat(path)

# Extract the data
faces = data['faces']
n = int(data['n'][0][0])
m = int(data['m'][0][0])
nfaces = data['nfaces'].flatten() # number of faces per person
idx = np.cumsum(nfaces) - nfaces # indices of the first face in each person

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