# Spectral Analysis

This chapter introduces spectral methods for extracting structure from
data: eigenvalue algorithms, singular value decomposition, and
linear systems. Here, *spectral* is meant in the
linear-algebra sense and refers to eigenvalues, singular values, and the
associated eigenvectors and singular vectors of a matrix.
These methods help us analyse chemical transition models, explore
the electronic structure of atoms and molecules, interpret
complex molecular spectra, and solve linear equations that arise in
data analysis and modelling.
This chapter includes more concepts in scientific programming, including
pseudo-random numbers, absolute and relative paths, loading external data files, 
linear algebra tools from `numpy`, and advanced plotting techniques from `matplotlib`.

The accompanying problem set for this chapter is available at
{ref}`sec:pset_2`.

| Section | Covered Examples | New Concepts and Tools |
| --- | --- | --- |
| Eigenvalue Algorithms | reaction pathway networks for repeated reaction | `@` operator, transition matrices, pseudo-random numbers, `np.linalg.eig` / `np.linalg.eigh` |
| Singular Value Decomposition | protein structure analysis from circular dichroism spectra | matrix norm, low-rank approximation, absolute and relative paths, advanced plotting, `np.linalg.svd` |
| Linear Systems | smoothie composition determination from absorption spectra | matrix inverse, Moore-Penrose pseudoinverse, matrix rank, `np.linalg.lstsq`, `np.linalg.pinv` |
