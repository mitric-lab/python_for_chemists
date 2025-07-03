#!/usr/bin/env python

from rdkit import Chem
from rdkit.Chem import Descriptors, Draw
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprintAsBitVect
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import random

# --- Constants ---
INPUT_CSV_PATH = "aptamers_not_processed.csv"
REGRESSION_OUTPUT_CSV_PATH = "aptamer_regression_data.csv"
CLASSIFICATION_OUTPUT_CSV_PATH = "aptamer_classification_data.csv"
REGRESSION_TARGET_COLUMN = "fl_int"
CLASSIFICATION_TARGET_COLUMN = "lambda_abs"

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

def molecules_to_fingerprints(mols, radius=2, n_bits=2048):
    """Convert molecules to Morgan fingerprints."""
    fingerprints = []
    for mol in mols:
        if mol is not None:
            fp = GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
            fingerprints.append(np.array(fp))
        else:
            fingerprints.append(np.zeros(n_bits))
    
    fingerprints_array = np.array(fingerprints)
    print(f"Generated {len(fingerprints)} fingerprints with {n_bits} bits each")
    return fingerprints_array

# --- Data Processing Functions ---

def load_data(input_path, target_column):
    """Load and clean the initial dataset."""
    df_raw = pd.read_csv(input_path)
    # Drop rows where the target is missing
    df_raw = df_raw.dropna(subset=[target_column])
    print(f"Loaded {len(df_raw)} samples from {input_path}")
    return df_raw

def create_molecules(smiles_list):
    """Create RDKit molecule objects from SMILES strings."""
    mols = [Chem.MolFromSmiles(smi) for smi in smiles_list]
    valid_mols = [mol for mol in mols if mol is not None]
    print(f"Created {len(valid_mols)} valid molecules from {len(smiles_list)} SMILES")
    return mols

def visualize_random_molecules(mols, num_samples=5, output_path=None):
    """Visualize random molecules and save to file."""
    if len(mols) < num_samples:
        print(f"Not enough molecules to sample {num_samples}, using all {len(mols)}")
        num_samples = len(mols)
    
    valid_mols = [mol for mol in mols if mol is not None]
    if len(valid_mols) < num_samples:
        print(f"Not enough valid molecules to sample {num_samples}, using all {len(valid_mols)}")
        num_samples = len(valid_mols)
    
    if num_samples == 0:
        print("No valid molecules to visualize")
        return
    
    sampled_mols = random.sample(valid_mols, num_samples)
    img = Draw.MolsToGridImage(
        sampled_mols,
        molsPerRow=num_samples,
        subImgSize=(200, 200),
    )
    fig, ax = plt.subplots(figsize=(15, 3))
    ax.imshow(img)
    ax.axis('off')
    
    if output_path is None:
        output_path = '../../assets/figures/05-machine_learning/random_molecules.svg'
    
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved random molecule visualization to {output_path}")

def calculate_molecular_features(mols, target_values, target_column):
    """Calculate molecular features for all molecules."""
    feature_list = [get_molecular_features(mol) for mol in mols]
    
    # Create a DataFrame from the calculated features
    features_df = pd.DataFrame(feature_list)
    features_df[target_column] = target_values
    
    print(f"Calculated {len(features_df.columns)-1} features for {len(features_df)} molecules")
    return features_df

def remove_correlated_features(features_df, target_column, threshold=0.95):
    """Detect and remove highly correlated features."""
    features_to_check = features_df.drop(target_column, axis=1)
    corr_matrix = features_to_check.corr().abs()
    
    # Select upper triangle of correlation matrix
    upper = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(bool))
    
    # Find features with correlation greater than threshold
    to_drop = [column for column in upper.columns if any(upper[column] > threshold)]
    
    print(f"Dropping highly correlated features (>{threshold}): {to_drop}")
    
    # Uncomment the following line to actually drop the features
    # features_df = features_df.drop(to_drop, axis=1)
    
    return features_df, to_drop

def plot_correlation_matrix(features_df, title="Correlation Matrix of Processed Features"):
    """Plot the correlation matrix of features."""
    plt.figure(figsize=(8, 8))
    sns.heatmap(features_df.corr(), annot=True, cmap="coolwarm", fmt=".2f")
    plt.title(title)
    plt.show()

def scale_features(features_df):
    """Scale features using StandardScaler."""
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(features_df)
    scaled_df = pd.DataFrame(scaled_data, columns=features_df.columns)
    print("Features scaled using StandardScaler")
    return scaled_df, scaler

def save_processed_data(features_df, output_path):
    """Save processed data to CSV file."""
    features_df.to_csv(output_path, index=False)
    print(f"Processed data with {len(features_df.columns)-1} features saved to '{output_path}'")
    print(features_df.head())

# --- Classification Analysis Functions ---

def perform_pca(fingerprints, n_components=2):
    """Perform PCA on fingerprints."""
    pca = PCA(n_components=n_components)
    pca_result = pca.fit_transform(fingerprints)
    
    explained_variance = pca.explained_variance_ratio_
    print(f"PCA explained variance: {explained_variance}")
    print(f"Total explained variance: {np.sum(explained_variance):.3f}")
    
    return pca_result, pca

def create_target_classes(target_values, threshold):
    """Create binary classes from continuous target values using given threshold."""
    classes = np.where(target_values > threshold, 1, -1)
    
    print(f"Created classes with threshold: {threshold:.3f}")
    print(f"Class distribution: {np.bincount(classes + 1)}")  # Add 1 to handle negative labels
    
    return classes

def plot_pca_2d(pca_result, target_classes, title="PCA Visualization"):
    """Plot 2D PCA results colored by target classes."""
    fig, ax = plt.subplots(figsize=(10, 8))
    
    class_names = ['Low', 'High']
    colors = ['blue', 'red']
    
    # Create scatter plot for each class
    for i, (class_name, color) in enumerate(zip(class_names, colors)):
        mask = target_classes == i
        ax.scatter(pca_result[mask, 0], pca_result[mask, 1], 
                  label=f'{class_name} (n={np.sum(mask)})', 
                  alpha=0.7, s=50, c=color)
    
    ax.set_xlabel('First Principal Component')
    ax.set_ylabel('Second Principal Component')
    ax.set_title(title)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.show()

# --- Main Pipeline ---

def main_regression(input_path, output_path, target_column):
    """Main processing pipeline."""
    # Load and clean data
    df_raw = load_data(input_path, target_column)
    
    # Create RDKit molecule objects
    mols = create_molecules(df_raw["smiles"])
    
    # Visualize random molecules
    visualize_random_molecules(mols)
    
    # Calculate molecular features
    features_df = calculate_molecular_features(mols, df_raw[target_column].values, target_column)
    
    # Remove highly correlated features
    features_df, dropped_features = remove_correlated_features(features_df, target_column)
    
    # Plot correlation matrix
    plot_correlation_matrix(features_df)
    
    # Scale features
    scaled_df, scaler = scale_features(features_df)
    
    # Save processed data
    save_processed_data(scaled_df, output_path)

def main_classification(input_path, output_path, target_column):
    """Main classification analysis pipeline."""
    # Load and clean data
    df_raw = load_data(input_path, target_column)
    
    # Create RDKit molecule objects
    mols = create_molecules(df_raw["smiles"])
    
    # Visualize random molecules
    # visualize_random_molecules(mols)
    
    # Convert molecules to fingerprints
    fingerprints = molecules_to_fingerprints(mols)
    
    # Perform PCA to 2D
    pca_result, pca_model = perform_pca(fingerprints, n_components=2)
    
    # Create target classes from continuous values
    threshold = np.median(df_raw[target_column].values)
    target_classes = create_target_classes(df_raw[target_column].values, threshold)
    
    # Plot PCA results
    # plot_pca_2d(pca_result, target_classes, title="Molecular Fingerprints in 2D PCA Space")
    
    # Create and save classification dataset
    classification_df = pd.DataFrame({
        'PC1': pca_result[:, 0],
        'PC2': pca_result[:, 1],
        f'{target_column}_continuous': df_raw[target_column].values,
        f'{target_column}_class': target_classes
    })
    
    classification_df.to_csv(output_path, index=False)
    print(f"Classification data saved to '{output_path}'")
    print(f"Dataset shape: {classification_df.shape}")
    print(classification_df.head())

if __name__ == "__main__":
    main_regression(INPUT_CSV_PATH, REGRESSION_OUTPUT_CSV_PATH, REGRESSION_TARGET_COLUMN)
    main_classification(INPUT_CSV_PATH, CLASSIFICATION_OUTPUT_CSV_PATH, CLASSIFICATION_TARGET_COLUMN)