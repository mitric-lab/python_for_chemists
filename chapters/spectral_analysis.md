# Spectral Analysis

This chapter introduces spectral methods for extracting structure from
data: eigenvalue algorithms, singular value decomposition, and
principal component analysis. Here, *spectral* is meant in the
linear-algebra sense and refers to eigenvalues, singular values, and the
associated eigenvectors and singular vectors of a matrix.
These methods help us analyse chemical transition models, explore
the electronic structure of atoms and molecules, interpret
complex molecular spectra, and uncover patterns in high-dimensional
experimental data.
This chapter includes more concepts in scientific programming, including
pseudo-random numbers, dataset loading, linear algebra tools from`numpy`, 
and advanced plotting techniques.

The accompanying problem set for this chapter is available at
{ref}`sec:pset_2`.

| Section | Covered Examples | New Concepts and Tools |
| --- | --- | --- |
| Eigenvalue Algorithms | reaction fate networks for repeated reaction and recycling steps | `@` operator, transition matrices, pseudo-random numbers, `np.linalg.eig` / `np.linalg.eigh` |
| Singular Value Decomposition | interpretation of mixture spectra | matrix norm, low-rank approximation, text file loading, `np.linalg.svd` |
| Principal Component Analysis | classification of wines | standardisation and normalisation, CSV file loadings, dimensionality reduction, advanced plotting |
