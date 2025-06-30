#!/usr/bin/env python

from rdkit import Chem
from rdkit.Chem import Descriptors, Draw
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
import random

# --- Constants ---
INPUT_CSV_PATH = "aptamers_not_processed.csv"
OUTPUT_CSV_PATH = "aptamer_data.csv"
TARGET_COLUMN = "fl_int"

# --- Feature Engineering ---

def get_molecular_features(mol):
    """Calculates all molecular descriptors for a given RDKit molecule."""
    if mol is None:
        return None

    # Define SMARTS patterns once
    aromatic_ring_pattern = Chem.MolFromSmarts("a1aaaaa1")
    
    features = {
        "MolWeight": Descriptors.MolWt(mol),
        "TPSA": Descriptors.TPSA(mol),
        "NumHDonors": Descriptors.NumHDonors(mol),
        "NumHAcceptors": Descriptors.NumHAcceptors(mol),
        "NumRotatableBonds": Descriptors.NumRotatableBonds(mol),
        "RingCount": Descriptors.RingCount(mol),
        "NumHeavyAtoms": Descriptors.HeavyAtomCount(mol),
        "FractionCSP3": Descriptors.FractionCSP3(mol),
        "NumAromaticRings": len(mol.GetSubstructMatches(aromatic_ring_pattern)),
    }
    return features

# --- Main Script ---

# Load the initial dataset
df_raw = pd.read_csv(INPUT_CSV_PATH)

# Drop rows where the target is missing
df_raw = df_raw.dropna(subset=[TARGET_COLUMN])

# Create RDKit molecule objects
mols = [Chem.MolFromSmiles(smi) for smi in df_raw["smiles"]]

# --- Visualize Random Molecules ---
if len(mols) >= 5:
    sampled_mols = random.sample(mols, 5)
    img = Draw.MolsToGridImage(
        sampled_mols,
        molsPerRow=5,
        subImgSize=(200, 200),
    )
    fig, ax = plt.subplots(figsize=(15, 3))
    ax.imshow(img)
    ax.axis('off')
    
    output_path = '../../assets/figures/05-machine_learning/random_molecules.svg'
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig) # Prevent the plot from displaying
    print(f"Saved random molecule visualization to {output_path}")

# Calculate features for each molecule
feature_list = [get_molecular_features(mol) for mol in mols]

# Create a DataFrame from the calculated features
features_df = pd.DataFrame(feature_list)
features_df[TARGET_COLUMN] = df_raw[TARGET_COLUMN].values

# --- Detect and Remove Correlated Features ---
# Use a copy for correlation analysis to keep the original features safe
features_to_check = features_df.drop(TARGET_COLUMN, axis=1)
corr_matrix = features_to_check.corr().abs()

# Select upper triangle of correlation matrix
upper = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(bool))

# Find features with correlation greater than 0.95
to_drop = [column for column in upper.columns if any(upper[column] > 0.95)]

print(f"Dropping highly correlated features: {to_drop}")
# features_df = features_df.drop(to_drop, axis=1)

# Plot the correlation matrix of the remaining features
plt.figure(figsize=(8, 8))
sns.heatmap(features_df.corr(), annot=True, cmap="coolwarm", fmt=".2f")
plt.title("Correlation Matrix of Processed Features")
plt.show()

# --- Scale Data ---
scaler = StandardScaler()
scaled_data = scaler.fit_transform(features_df)
features_df = pd.DataFrame(scaled_data, columns=features_df.columns)


# --- Save Processed Data ---
features_df.to_csv(OUTPUT_CSV_PATH, index=False)

print(f"\nProcessed data with {len(features_df.columns)-1} features saved to '{OUTPUT_CSV_PATH}'")
print(features_df.head())