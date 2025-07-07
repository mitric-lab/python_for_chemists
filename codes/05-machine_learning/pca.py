#!/usr/bin/env python

### ANCHOR: pca_init
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

class PCA:
    def __init__(self, k=2):
        self.k = k
        self.components = None
        self.explained_variance = None
### ANCHOR_END: pca_init

### ANCHOR: pca_fit    
    def fit(self, X):
        # Center the data
        X_centered = X - np.mean(X, axis=0)
        
        # Compute the covariance matrix
        cov_matrix = np.cov(X_centered, rowvar=False)
        
        # Compute the eigenvalues and eigenvectors
        eigenvalues, eigenvectors = np.linalg.eigh(cov_matrix)
        
        # Sort the eigenvalues and corresponding eigenvectors in descending order
        sorted_indices = np.argsort(eigenvalues)[::-1]
        eigenvalues = eigenvalues[sorted_indices]
        eigenvectors = eigenvectors[:, sorted_indices]
        
        # Store the top k eigenvectors (principal components)
        self.components = eigenvectors[:, :self.k]
        
        # Calculate the explained variances of all principal components
        self.explained_variance = eigenvalues / np.sum(eigenvalues)
### ANCHOR_END: pca_fit
    
### ANCHOR: pca_transform
    def transform(self, X):
        # Project the data onto the first k principal components
        X_centered = X - np.mean(X, axis=0)
        return np.dot(X_centered, self.components)
    
    def fit_transform(self, X):
        # Fit the model and return the transformed data
        self.fit(X)
        return self.transform(X)
### ANCHOR_END: pca_transform

### ANCHOR: load_data_from_csv
path_to_csv = "aptamer_fingerprints_data.csv"
df = pd.read_csv(path_to_csv)
print(df.head())    
### ANCHOR_END: load_data_from_csv

### ANCHOR: rdkit_fingerprints
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprintAsBitVect

# Create molecule from SMILES
smiles = "CCO"  # ethanol
mol = Chem.MolFromSmiles(smiles)

# Basic molecular properties
print(f"Number of atoms: {mol.GetNumAtoms()}")
print(f"Number of heavy atoms: {mol.GetNumHeavyAtoms()}")

# Generate Morgan fingerprint
morgan_fp = GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
fp_array = np.array(morgan_fp)
print(f"Morgan fingerprint shape: {fp_array.shape}")
print(f"Number of bits set: {np.sum(fp_array)}")
### ANCHOR_END: rdkit_fingerprints

### ANCHOR: process_data
target_column = "lambda_abs_class"
X = df.drop(target_column, axis=1).values # Ignore the target column
print(X.shape)
### ANCHOR_END: process_data

### ANCHOR: pca_fit_transform
pca = PCA(k=2)
X_pca = pca.fit_transform(X)
print(X_pca.shape) # (N, 2)
### ANCHOR_END: pca_fit_transform

### ANCHOR: pca_explained_variance
n = 5

fig, ax = plt.subplots(figsize=(8, 6))
ax.bar(range(1, n + 1), pca.explained_variance[:n], color='tab:blue', label='Explained Variance')
ax.plot(range(1, n + 1), np.cumsum(pca.explained_variance[:n]), 'o-', c='tab:orange', label='Cumulative Explained Variance')
ax.set_xlabel('Principal Component')
ax.set_ylabel('Explained Variance')
ax.legend()
plt.show()
### ANCHOR_END: pca_explained_variance

print(np.cumsum(pca.explained_variance[:n]))

# fig.savefig('../../assets/figures/05-machine_learning/pca_aptamers_explained_variance.svg')

### ANCHOR: plot_pca
fig, ax = plt.subplots(figsize=(6, 6))
ax.scatter(X_pca[:, 0], X_pca[:, 1], color='blue', label='Data')
ax.set_xlabel('Principal Component 1')
ax.set_ylabel('Principal Component 2')
ax.legend()
plt.show()
### ANCHOR_END: plot_pca

# fig.savefig('../../assets/figures/05-machine_learning/pca_aptamers.svg')

### ANCHOR: plot_pca_hover
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from rdkit.Chem import Draw

# Load original data to get SMILES
df_original = pd.read_csv("aptamers_not_processed.csv")
smiles_list = df_original["smiles"].values

# Create RDKit molecule objects
mols = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]

# Create matplotlib scatter plot with minimal styling to match plot_pca
fig, ax = plt.subplots(figsize=(6, 6))

sc = ax.scatter(X_pca[:, 0], X_pca[:, 1], color='blue', label='Data')

ax.set_xlabel('Principal Component 1')
ax.set_ylabel('Principal Component 2')
ax.legend()

# Create annotation box for molecular structures
imagebox = OffsetImage(np.zeros((100, 100, 3)), zoom=0.5)
imagebox.image.axes = ax

ab = AnnotationBbox(
    imagebox, (30, 30), xycoords='data', 
    boxcoords="offset points", arrowprops=dict(arrowstyle="->"), 
)
ax.add_artist(ab)
ab.set_visible(False)  # Initially hidden

def update_annotation_box(idx):
    """Update annotation box position and molecular image."""
    ab.xy = (X_pca[idx, 0], X_pca[idx, 1])
    mol = mols[idx]
    if mol is not None:
        img = Draw.MolToImage(mol, size=(100, 100), wedgeBonds=False)
        ab.offsetbox.set_data(img)
    else:
        # Show empty image for invalid molecules
        ab.offsetbox.set_data(np.zeros((100, 100, 3)))

def hover(event):
    """Handle hover events to show/hide molecular structures."""
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

# Connect hover event to the figure
fig.canvas.mpl_connect('motion_notify_event', hover)

plt.show()
### ANCHOR_END: plot_pca_hover

# fig.savefig('../../assets/figures/05-machine_learning/pca_aptamers_interactive.svg')