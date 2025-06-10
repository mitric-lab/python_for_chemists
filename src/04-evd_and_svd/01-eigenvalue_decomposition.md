## Eigenvalue Decomposition

To understand the definition of the eigenvalue decomposition (EVD), 
we first need to get acquainted with the concept of a vector space.
Although this is defined much more generally, we will restrict ourselves here
to $\mathbb{C}^n$, the $n$-dimensional vector space over the complex numbers.
Vectors in $\mathbb{C}^n$ can be represented as column vectors 
in the standard basis. The notation $\vec{v} \in \mathbb{C}^n$
denotes such a vector:
$$
  \vec{v} = \begin{pmatrix} v_1 \\ v_2 \\ \vdots \\ v_n \end{pmatrix}
$$
with components $v_i \in \mathbb{C}$. A complex conjugate row vector 
is denoted as the Hermitian conjugate of a column vector:
$$
  \vec{v}^\dag = \begin{pmatrix} v_1^* & v_2^* & \cdots & v_n^* \end{pmatrix}
$$
with the complex conjugated components $v_i^*$.
The standard scalar product of two vectors $\vec{v}, \vec{w} \in \mathbb{C}^n$
can be written in this notation as
$$
  \langle \vec{v}, \vec{w} \rangle = \vec{v}^\dag \vec{w} = \sum_{i=1}^n v_i^* w_i\,,
$$
where the rule for matrix multiplication applies.

Let $\bm{A} \in \C{n}{n}$ be a square matrix. A vector
$\vec{v}_i \in \mathbb{C}^n$ is called an eigenvector of $\bm{A}$ if
$$
  \bm{A} \vec{v}_i = \lambda_i \vec{v}_i
$$
for some $\lambda \in \mathbb{C}$. The number $\lambda_i$ is referred to
as the eigenvalue of $\bm{A}$ corresponding to the eigenvector $\vec{v}_i$.

## Theoretical Foundations

We now define the diagonal matrix of eigenvalues $\bm{\Lambda}$ as
$$
  \bm{\Lambda} = \begin{pmatrix}
    \lambda_1 & 0 & \cdots & 0 \\
    0 & \lambda_2 & \cdots & 0 \\
    \vdots & \vdots & \ddots & \vdots \\
    0 & 0 & \cdots & \lambda_n
  \end{pmatrix}
$$
and the matrix of eigenvectors $\bm{Q}$ as
$$
  \bm{Q} = \begin{pmatrix}
    \vert & \vert &  & \vert \\
    \vec{v}_1 & \vec{v}_2 & \cdots & \vec{v}_n \\
    \vert & \vert &  & \vert
  \end{pmatrix}\,.
$$
Then, the eigenvalue decomposition of $\bm{A}$ can be written as
$$
  \bm{A} = \bm{Q} \bm{\Lambda} \bm{Q}^{-1}
$$

### Eigenvalue Decomposition for Normal Matrices

Although the eigenvalue decomposition exists for all diagonalisable
matrices, we will restrict ourselves here to a subset of such matrices,
the so-called 
[*normal matrices*](https://en.wikipedia.org/wiki/Normal_matrix). 
A matrix $\bm{A}\in\C{n}{n}$ is called *normal* if it commutes with 
its Hermitian conjugate, 
i.e., if $\bm{A}^\dag \bm{A} = \bm{A} \bm{A}^\dag$ holds.
Such matrices are unitarily diagonalisable, meaning that 
the matrix $\bm{Q}$ is unitary. 
Therefore, it holds
$$
  \bm{Q}^{-1} = \bm{Q}^\dag = \begin{pmatrix}
    \text{---}\ \vec{v}_1^\dag\ \text{---} \\
    \text{---}\ \vec{v}_2^\dag\ \text{---} \\
    \vdots \\
    \text{---}\ \vec{v}_n^\dag\ \text{---}
  \end{pmatrix}\,.
$$

For a normal matrix $\bm{A}$, the eigenvalue decomposition can thus be written as
$$
  \begin{align*}
    \bm{A} &= \bm{Q} \bm{\Lambda} \bm{Q}^\dag \\
    &= \begin{pmatrix}
      \vert & \vert &  & \vert \\
      \vec{v}_1 & \vec{v}_2 & \cdots & \vec{v}_n \\
      \vert & \vert &  & \vert
    \end{pmatrix}
    \begin{pmatrix}
      \lambda_1 & 0 & \cdots & 0 \\
      0 & \lambda_2 & \cdots & 0 \\
      \vdots & \vdots & \ddots & \vdots \\
      0 & 0 & \cdots & \lambda_n
    \end{pmatrix}
    \begin{pmatrix}
      \text{---}\ \vec{v}_1^\dag\ \text{---} \\
      \text{---}\ \vec{v}_2^\dag\ \text{---} \\
      \vdots \\
      \text{---}\ \vec{v}_n^\dag\ \text{---}
    \end{pmatrix} \\
    &= \lambda_1 \vec{v}_1 \vec{v}_1^\dag + \lambda_2 \vec{v}_2 \vec{v}_2^\dag + \cdots + \lambda_n \vec{v}_n \vec{v}_n^\dag \\
    &= \sum_{i=1}^n \lambda_i \vec{v}_i \vec{v}_i^\dag \,.
    {{numeq}}{eq:eigenvalue_decomposition_normal}
  \end{align*}
$$

The product of a column vector with a row vector, such as
$\vec{v}_i \vec{v}_i^\dag$ in the above equation, is known as the
[*dyadic product*](https://en.wikipedia.org/wiki/Dyadics)
and produces a (rank-1) matrix.

To perform the eigenvalue decomposition, we first need to determine the
eigenvalues and eigenvectors of the matrix. While they can still be determined
analytically for small matrices, numerical methods are necessary for larger
matrices (in general for $n \geq 5$). We will discuss two such methods below.

#### Eigenvalue Algorithms

**Power Iteration**

The probably simplest algorithm for determining the eigenvalues of a matrix
is [*power iteration*](https://en.wikipedia.org/wiki/Power_iteration).
This algorithm finds the eigenvector of the matrix $\bm{A}\in\C{n}{n}$
corresponding to the eigenvalue with the largest absolute value, 
provided that $\bm{A}$ is diagonalisable.

We assume that $\bm{A}\in\C{n}{n}$ is a matrix with (unknown) eigenvalues
$$
  |\lambda_1| > |\lambda_2| \geq \cdots \geq |\lambda_n|
$$
Then, the iteration
$$
  \vec{x}^{k+1} = \frac{\bm{A} \vec{x}^k}{\|\bm{A} \vec{x}^k\|}\quad k=0,1,2,\ldots
  {{numeq}}{eq:power_iteration}
$$
for an arbitrary starting vector $\vec{x}^0 \in \mathbb{C}^n$ produces a sequence
of vectors $(\vec{x}^k)_{k\in\mathbb{N}}$ that converges to the eigenvector
corresponding to the eigenvalue with the largest absolute value, $\lambda_1$.

The corresponding eigenvalue can then be determined by
$$
  \lambda_1 = \frac{x^{k\dag} \bm{A} \vec{x}^k}{x^{k\dag} \vec{x}^k}
  {{numeq}}{eq:power_iteration_eigenvalue}
$$

For those interested in the mathematics, a proof of the convergence 
of the power iteration can be found here.

```admonish proof title="Proof of Convergence of Power Iteration" collapsible=true
A matrix is diagonalisable if and only if its eigenvectors form a basis.
Therefore, we can write the arbitrary starting vector $\vec{x}^0$ as a linear
combination of the eigenvectors:
$$
  \vec{x}^0 = \sum_{i=1}^n a_i \vec{v}_i
$$
with coefficients $a_i\in\mathbb{C}$ and eigenvectors $\vec{v}_i$
corresponding to the eigenvalues $\lambda_i$. Then, we have
$$
  \begin{align*}
    \vec{x}^k &= \bm{A} \vec{x}^{k-1} = \bm{A}^2 \vec{x}^{k-2} = \cdots = \bm{A}^k \vec{x}^0 \\
    &= \sum_{i=1}^n a_i \bm{A}^k \vec{v}_i \\
    &= \sum_{i=1}^n a_i \lambda_i^k \vec{v}_i \\
    &= \lambda_1^k a_1 \vec{v}_1 + \sum_{i=2}^n \lambda_i^k a_i \vec{v}_i \\
    &= \lambda_1^k \left[ 
      a_1 \vec{v}_1 + \sum_{i=2}^n \left(\frac{\lambda_i}{\lambda_1}\right)^k a_i \vec{v}_i
    \right]\,.
  \end{align*}
$$

Because $|\lambda_i| < |\lambda_1|$ for $i>1$, the second term in the parenthesis
converges to zero as $k\to\infty$, so that
$$
  \lim_{k\to\infty} \vec{x}^k = \lim_{k\to\infty} \lambda_1^k a_1 \vec{v}_1 \propto \vec{v}_1\,.
$$
With a finite number of iterations $k$, $\vec{x}^k$ can thus be considered
as an approximation of the eigenvector $\vec{v}_1$, especially since eigenvectors
can be scaled arbitrarily.

For numerical stability, it is advisable to normalise the vector $\vec{x}^k$ 
in each step, which leads to the iteration
$$
  \vec{x}^{k+1} = \frac{\bm{A} \vec{x}^k}{\|\bm{A} \vec{x}^k\|}\,.
$$
```

A generalisation of the power iteration is the so-called *inverse iteration*, 
which can be used to determine an eigenvalue close to a given 
but arbitrary number $\mu$ and its eigenvector.

**QR Algorithm**

In many cases, however, we are interested in all eigenvectors and eigenvalues
of a matrix. A frequently used method that can determine all eigenpairs
of a diagonalisable matrix is the 
[*QR algorithm*](https://en.wikipedia.org/wiki/QR_algorithm).

The key step of the QR algorithm is the decomposition of the matrix $\bm{A}$
into a unitary matrix $\bm{Q}$ and an upper triangular matrix $\bm{R}$, i.e.,
$$
  \bm{A} = \bm{Q} \bm{R}\,.
$$
This decomposition is called the 
[*QR decomposition*](https://en.wikipedia.org/wiki/QR_decomposition).
We will not go into the details of the QR decomposition here, but assume
that we can compute it for any square matrix $\bm{A} \in \C{n}{n}$.

In the QR algorithm, the initial matrix is set to $\bm{A}_0 = \bm{A}$ 
and its QR decomposition $\bm{A}_0 = \bm{Q}_0 \bm{R}_0$ is computed. 
Then, the iteration
$$
  \bm{A}_{k+1} = \bm{R}_k \bm{Q}_k
  {{numeq}}{eq:qr_algorithm}
$$
is performed. In other words, the matrix in the next step $\bm{A}_{k+1}$ is computed
by multiplying the factors of the QR decomposition in reverse order. It holds
$$
  \begin{align*}
    \bm{A}_{k+1} &= \bm{R}_k \bm{Q}_k = \textcolor{blue}{\identity} \textcolor{green}{\bm{R}_k} \bm{Q}_k \\
    &= \textcolor{blue}{\bm{Q}_k^\dag \bm{Q}_k} \textcolor{green}{\bm{R}_k} \bm{Q}_k \\
    &= \bm{Q}_k^\dag \textcolor{darkcyan}{\bm{A}_k} \bm{Q}_k \\
    &= \bm{Q}_k^\dag \textcolor{darkcyan}{\bm{R}_{k-1} \bm{Q}_{k-1}} \bm{Q}_k \\
    &= \bm{Q}_k^\dag \textcolor{orange}{\identity} \textcolor{darkcyan}{\bm{R}_{k-1} \bm{Q}_{k-1}} \bm{Q}_k \\
    &= \bm{Q}_k^\dag \textcolor{orange}{\bm{Q}_{k-1}^\dag \bm{Q}_{k-1}} \textcolor{darkcyan}{\bm{R}_{k-1} \bm{Q}_{k-1}} \bm{Q}_k \\
    &= \bm{Q}_k^\dag \textcolor{orange}{\bm{Q}_{k-1}^\dag} \textcolor{brown}{\bm{A}_{k-1}} \textcolor{darkcyan}{\bm{Q}_{k-1}} \bm{Q}_k \\
    &  \phantom{=} \vdots \\
    &= \bm{Q}_k^\dag \bm{Q}_{k-1}^\dag \cdots \bm{Q}_0^\dag \textcolor{red}{\bm{A}_0} \bm{Q}_0 \bm{Q}_1 \cdots \bm{Q}_{k-1} \bm{Q}_k \\
    &= \underbrace{\bm{Q}_k^\dag \bm{Q}_{k-1}^\dag \cdots \bm{Q}_0^\dag \textcolor{red}{\bm{Q}}}_{:=\bm{P}_k} 
      \textcolor{red}{\bm{\Lambda}}
      \underbrace{\textcolor{red}{\bm{Q}^{-1}} \bm{Q}_0 \bm{Q}_1 \cdots \bm{Q}_{k-1} \bm{Q}_k}_{=\bm{P}_k^{-1}} \\
    &= \bm{P}_k \bm{\Lambda} \bm{P}_k^{-1} \,,
  \end{align*}
$$
where $\bm{A}_0 = \bm{Q} \bm{\Lambda} \bm{Q}^{-1}$ is the exact but unknown
eigenvalue decomposition of the initial matrix.
It thus holds that the matrix $\bm{A}_k$ has the same eigenvalues 
as the initial matrix $\bm{A}$.

It can be shown that the sequence of matrices $(\bm{A}_k)_{k\in\mathbb{N}_0}$
converges to an upper triangular matrix, whose eigenvalues can be 
read off from the diagonal. These eigenvalues are then
simultaneously the eigenvalues of the initial matrix $\bm{A}$.
If the initial matrix $\bm{A}$ is also normal, the sequence of matrices converges
to a diagonal matrix, since all normal upper triangular matrices are diagonal.
Since the matrix of eigenvectors $\bm{Q}$ of a diagonal matrix 
is the identity matrix $\identity$, i.e., $\bm{P}_k \approx \identity$, 
it follows from the above definition that
$$
  \bm{P}_k = \bm{Q}_k^\dag \bm{Q}_{k-1}^\dag \cdots \bm{Q}_0^\dag \bm{Q} \approx \identity\,,
$$
which implies that
$$
  \bm{Q} \approx \bm{Q}_0 \bm{Q}_1 \cdots \bm{Q}_{k-1} \bm{Q}_k = \prod_{i=0}^k \bm{Q}_i\,,
  {{numeq}}{eq:qr_algorithm_eigenvectors}
$$
which allows us to approximate the eigenvectors of the normal matrix $\bm{A}$.

The QR algorithm can be accelerated in various ways, which will not
be discussed here.

### Implementation

We will now implement the eigenvalue decomposition using the power iteration
and the QR algorithm. For testing, we will first choose random real matrices.
To ensure that the matrix is always diagonalisable, we will add it to its
transpose to obtain a symmetric matrix, which is guaranteed to be normal.

#### Power Iteration

After importing NumPy
```python
{{#include ../codes/04-evd_and_svd/power_iteration.py:import}}
```
we implement the power iteration according to Eq. {{eqref: eq:power_iteration}}
as a function:
```python
{{#include ../codes/04-evd_and_svd/power_iteration.py:power_iteration_function}}
```

This function accepts the arguments `a` (matrix $\bm{A}$),
`eps` (error tolerance for convergence), and `maxiter` (maximum number of
iterations). It returns the eigenvalue with the largest absolute value and the
corresponding eigenvector.

First, the dimension of the matrix is determined with `a.shape[0]` and stored in
the variable `n`. Then, we use the function
[`np.random.rand`](https://numpy.org/doc/stable/reference/random/generated/numpy.random.rand.html)
to generate a random starting vector, which is then divided by its norm
and thus normalised. For this, we used the function
[`np.linalg.norm`](https://numpy.org/doc/stable/reference/generated/numpy.linalg.norm.html)
to compute the (Euclidean) norm of a vector.

A possible stopping condition is when the norm of the difference between the
vectors in two consecutive iterations is smaller than a given error tolerance 
`eps`. To implement this, we define the variable `x_prev`, which stores 
the vector from the previous iteration.

Then, the iteration begins with a `for` loop that runs up to the maximum
number of iterations `maxiter`. In each iteration, the value of `x` is first
stored in `x_prev`. Then, the new vector `x` is computed by `np.dot(a, x)`
and subsequently normalised. After that, it is checked whether the norm 
of the difference between `x` and `x_prev` is smaller than `eps`.
If this is the case, the loop is terminated.

If the loop variable `i` reaches the value `maxiter - 1`, the power iteration
may not have converged, and a warning is issued.
Then, the eigenvalue is computed using Eq. {{eqref: eq:power_iteration_eigenvalue}}.
Finally, the results are returned.

```admonish tip title="Tip for Programming Style"
Sometimes it happens that the desired variable name is already a reserved
command in Python. In this case, either a different name can be chosen or, 
according to convention, an underscore can be appended to the name.
Here, for example, the variable `lambda_` was used to store the eigenvalue.
```

In the following, we generate a symmetrised random matrix `a_mat` and call
the function `power_iteration` with `eps=1e-8`:
```python
{{#include ../codes/04-evd_and_svd/power_iteration.py:test_power_iteration}}
```
To verify our results, we use the function
[`np.linalg.eigh`](https://numpy.org/doc/stable/reference/generated/numpy.linalg.eigh.html),
which computes the eigenvalues and eigenvectors of a symmetric matrix.
The eigenvalues are returned in ascending order in the 1D array `eigvals`,
and the normalised eigenvectors as columns in the 2D array `eigvecs`.
For this reason, we take the last entry of the eigenvalues and the corresponding
eigenvector to compare them with the results of the power iteration:
```python
{{#include ../codes/04-evd_and_svd/power_iteration.py:verify_results}}
```
While we can directly compare the eigenvalues with `np.isclose`,
we must consider that the eigenvectors are only uniquely defined up to 
a scaling factor. Since our eigenvectors are normalised, as well as 
those of `np.linalg.eigh`, we can compute the projection of the eigenvectors
onto each other using `np.dot`, and the absolute value of the projection 
should be close to 1.

Repeatedly running our code would test our implementation of the power iteration
on different examples.

#### QR Algorithm
We will now implement the QR algorithm according to
Eq. {{eqref: eq:qr_algorithm}} and Eq. {{eq: eq:qr_algorithm_eigenvectors}}.
After importing NumPy, we define the function `qr_algorithm`, which has the same
signature as the function `power_iteration`:
```python
{{#include ../codes/04-evd_and_svd/qr_algorithm.py:qr_algorithm_function}}
```
After determining the dimension, the input matrix `a` is copied to `a_k`.
It is important here to make an explicit copy of `a`, as the QR algorithm
would otherwise modify the matrix `a`. To be able to multiply the future
$\bm{Q}_k$ sequentially, the variable `eigvecs` is initialised as an identity matrix.

Within the `for` loop, which runs up to `maxiter`, the QR decomposition of
`a_k` is first computed. For this, we used the function
[`np.linalg.qr`](https://numpy.org/doc/stable/reference/generated/numpy.linalg.qr.html).
Then, `a_k` is overwritten with the product of `r_k` and `q_k`, and `eigvecs`
is multiplied by `q_k`.

Since `a_k` approaches an upper triangular matrix 
(or in this case a diagonal matrix), we can use the largest absolute value 
of the elements in the strict lower triangle of the matrix as 
a termination condition. For this, the function
[`np.tril`](https://numpy.org/doc/stable/reference/generated/numpy.tril.html)
is used with the argument `k=-1` to obtain the lower triangle of the matrix.
If this value is smaller than the error tolerance `eps`, the loop is terminated.
Afterwards, as with the implementation of the power iteration, it is checked
whether the algorithm has converged and a warning is issued if necessary.

The eigenvalues are extracted from the diagonal of the matrix `a_k` using `np.diag`.
At the end, the order of the indices of the eigenvalues is computed with
`np.argsort`, and the eigenvalues and eigenvectors are sorted accordingly 
before being returned.

As with the power iteration, we test the implementation with a symmetric
random matrix:
```python
{{#include ../codes/04-evd_and_svd/qr_algorithm.py:test_qr_algorithm}}
```
and compare our results with those of `np.linalg.eigh`:
```python
{{#include ../codes/04-evd_and_svd/qr_algorithm.py:verify_results}}
```
The pairwise projection of our eigenvectors onto all eigenvectors of
`np.linalg.eigh` can be elegantly performed with the matrix multiplication
`eigvecs.T @ eigvecs_ref`. In the case of a successfully 
converged QR algorithm, this product should correspond to a diagonal matrix
with diagonal elements close to 1 or -1. 
Therefore, we take the absolute value of the elements from the product matrix
and compare them with the identity matrix, which we create with
[`np.eye`](https://numpy.org/doc/stable/reference/generated/numpy.eye.html).

The eigenvalue decomposition provides us with many essential
information about the corresponding square matrix. 
It would be nice if one could find something similar for 
non-square matrices as well. 
As you will see in the next section, the answer is a clear yes.
We will learn about the singular value decomposition, a method that can be
considered the generalisation of the eigenvalue decomposition 
for general matrices.

