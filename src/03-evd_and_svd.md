# Eigenvalue and Singular Value Decomposition

Eigenvalues and eigenvectors of a square matrix should be familiar to you.
If all eigenvalues and eigenvectors of a diagonalisable matrix
$\bm{A}$ have been found, it can be expressed as
$\bm{A} = \bm{Q} \bm{\Lambda} \bm{Q}^{-1}$, 
where $\bm{\Lambda}$ is a diagonal matrix of eigenvalues
and $\bm{Q}$ is the matrix of eigenvectors. This representation is known as
*Eigenvalue Decomposition*.

The *Singular Value Decomposition* can be considered a generalisation of the
Eigenvalue Decomposition, as it is defined for matrices of arbitrary dimensions.
The Singular Value Decomposition of a matrix
$\bm{A}\in \R{m}{n}$ is given by $\bm{A} = \bm{U} \bm{\Sigma} \bm{V}^\intercal$, 
where $\bm{\Sigma}\in \R{m}{n}$ is a diagonal matrix of singular values
and $\bm{U}\in \R{m}{m}$ and $\bm{V}\in \R{n}{n}$ are orthogonal matrices.
The exact details will be covered later in this chapter.

In the authors' opinion, the Singular Value Decomposition is **the** 
most important matrix decomposition in (numerical) linear algebra and 
has countless applications in the natural sciences, engineering, and beyond.
At the same time, it forms the basis for many modern algorithms.
In this chapter, we will focus on the Eigenvalue and Singular Value Decomposition
and discuss their significance for quantum chemistry, data analysis, 
and machine learning.

